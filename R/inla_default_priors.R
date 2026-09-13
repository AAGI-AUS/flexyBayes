# inla_default_priors.R -- the priors the package chooses for an INLA fit.
#
# A prior the caller states reaches INLA exactly as stated. This file decides
# the rest: the prior on every quantity the caller left unstated. Each one is
# written in the units of the quantity it governs, so multiplying the
# response or a covariate by a constant multiplies the prior with it and the
# fit stays one model:
#
#   * a random-effect standard deviation with no stated prior takes the
#     bounded uniform default, uniform(0, U), with U = 5 sd(y) on an identity
#     link and a constant on a log or logit link (.default_uniform_scale());
#   * the standard deviation of an uncorrelated slope is in response units
#     per covariate unit, so its default bound is U / sd(x).
#
# A Gaussian identity-link fit runs on a standardised response
# (R/inla_standardise.R), and each of these is stated in original units and
# restated on those coordinates with the stated priors. Before this file,
# anything unstated took INLA's own default on the standardised coordinates,
# which is scale-free in the response and blind to every covariate, and two
# failures followed. Each was silent: INLA reported mode.status 0 and the
# numerical-confirm gate passed.
#
#   * An unstated random-intercept SD degenerated. On a 400-row design with
#     four blocks (true SD 2) and only the residual prior stated, 20 of 20
#     identical per-row fits returned a block SD above 1e6 (median 1e7); on
#     the aggregated route 10 of 20 did and 5 more failed in the engine, and
#     the streamed route, which attached no prior at all, gave 10 and 4.
#   * lme4::sleepstudy with `(Dx || Subject)` and Dx = Days / 1000: the
#     default bound 5 sd(y) = 282 on the slope SD sat just above the fitted
#     279 (0.279 per day, against 6.45 at Dx = Days).
#
# A prior the caller states is never rescaled by a covariate.


# --- random-effect standard deviations ------------------------------------ #

#' The standard-deviation targets the INLA default reaches
#'
#' The walker in `.fb_default_prior_targets()` names the targets the shared
#' default reaches on every engine. INLA differs from it in two ways. A
#' multi-way interaction `a:b:c` is one iid field keyed by its colon-joined
#' name, which takes a standard-deviation prior as the two-way interaction
#' does. `us(f):g` and `(x | g)` are Wishart fields with no
#' standard-deviation slot, so a filled prior would have nowhere to go.
#'
#' @param fb An `fb_terms` object.
#' @returns A list with `shared` and `vm_ped`, group names in the key the
#'   prior DSL uses.
#' @noRd
#' @keywords internal
.fb_inla_sd_targets <- function(fb) {
  shared <- character(0)
  vm_ped <- character(0)
  for (t in fb$random_terms %||% list()) {
    ttype <- t$type %||% ""
    if (ttype %in% c("us_gxe", "simple_slope_cor")) {
      next
    }
    if (identical(ttype, "combo")) {
      shared <- c(shared, paste(.fb_nested_vars(t), collapse = ":"))
      next
    }
    walked <- .fb_default_prior_targets(list(random_terms = list(t)))
    shared <- c(shared, walked$shared)
    vm_ped <- c(vm_ped, walked$vm_ped)
  }
  list(shared = unique(shared), vm_ped = unique(vm_ped))
}

#' The standard deviation of each uncorrelated slope's covariate
#'
#' @param fb An `fb_terms` object.
#' @param data The data the fit is built from, or `NULL`.
#' @returns A named numeric vector keyed `<slope_var>_<group>`, the key the
#'   prior DSL gives the slope's standard deviation, for every uncorrelated
#'   slope on a numeric covariate with a positive standard deviation.
#' @noRd
#' @keywords internal
.fb_inla_slope_covariate_sd <- function(fb, data) {
  out <- numeric(0)
  if (is.null(data)) {
    return(out)
  }
  for (t in fb$random_terms %||% list()) {
    if (!identical(t$type %||% "", "simple_slope_uncor")) {
      next
    }
    x <- data[[t$slope_var]]
    if (!is.numeric(x)) {
      next
    }
    x <- x[is.finite(x)]
    s <- if (length(x) >= 2L) stats::sd(x) else NA_real_
    if (!is.finite(s) || s <= 0) {
      next
    }
    out[[paste0(t$slope_var, "_", t$var)]] <- s
  }
  out
}

# A uniform standard-deviation prior for a quantity divided by `s`.
.fb_inla_divide_sd_spec <- function(spec, s) {
  args <- spec$spec$args
  args$upper <- args$upper / s
  if (!is.null(args$lower)) {
    args$lower <- args$lower / s
  }
  spec$spec$args <- args
  spec
}

#' Every INLA variance-component prior, with what the caller left unstated
#' filled in
#'
#' Two additions to the `priors_to_inla()` list the emit built from what was
#' stated. First, when the fit carries a prior of any kind -- an `fb_prior()`,
#' partial or whole, or the legacy `prior_vc_sd` scalar -- each target of
#' `.fb_inla_sd_targets()` with no entry takes the bounded uniform default.
#' A fit with no prior at all is left alone: every public route resolves one
#' before it reaches an emit (flexybayes() builds the auto-default and
#' flexybayes_stream() attaches it), so that case is a direct internal call.
#' Second, a slope bound the package chose, filled here or built by the
#' auto-default, is divided by the standard deviation of its covariate. A
#' bound the caller wrote is not.
#'
#' Each added entry is marked `meta$package_chosen`, which is how the formula
#' builder tells it from a stated prior on a structure that refuses one.
#'
#' @param fb An `fb_terms` object.
#' @param hyper_ctrl The `priors_to_inla()` list of what was stated.
#' @param data The data the fit is built from; `NULL` where the route holds
#'   none.
#' @param resp_sd The response standard deviation when `data` does not hold
#'   the full response (the aggregated and streamed routes read it off the
#'   sufficient statistics).
#' @returns A list: `hyper_ctrl` with the additions, and `record`, which
#'   holds `filled` and `covariate_scaled` (lists of `fb_prior` specs, named
#'   by group) and `covariate_sd` (the divisor of each scaled bound).
#' @noRd
#' @keywords internal
.fb_inla_complete_priors <- function(fb, hyper_ctrl, data, resp_sd = NULL) {
  stated <- fb$priors
  has_prior <- inherits(stated, "fb_prior") ||
    .fb_prior_scalar_supplied(fb, "vc_sd")
  from_default <- inherits(stated, "fb_prior") &&
    !is.null(attr(stated, "fb_prior_default_basis"))
  slope_sd <- .fb_inla_slope_covariate_sd(fb, data)

  filled <- list()
  if (has_prior) {
    targets <- .fb_inla_sd_targets(fb)
    keys <- setdiff(c(targets$shared, targets$vm_ped), names(hyper_ctrl))
    if (length(keys)) {
      upper <- .default_uniform_scale(
        data,
        fb$response,
        fb$family,
        fb$link,
        resp_sd = resp_sd
      )$scale
      for (k in keys) {
        spec <- list(
          target = list(type = "sd", group = k),
          spec = list(family = "uniform", args = list(lower = 0, upper = upper))
        )
        if (k %in% names(slope_sd)) {
          spec <- .fb_inla_divide_sd_spec(spec, slope_sd[[k]])
        }
        filled[[k]] <- spec
      }
    }
  }

  covariate_scaled <- list()
  if (from_default) {
    for (spec in stated$specs) {
      k <- if (identical(spec$target$type, "sd")) spec$target$group else ""
      if (k %in% names(slope_sd) && identical(spec$spec$family, "uniform")) {
        covariate_scaled[[k]] <- .fb_inla_divide_sd_spec(spec, slope_sd[[k]])
      }
    }
  }

  added <- c(covariate_scaled, filled)
  if (length(added)) {
    entries <- priors_to_inla(structure(
      list(specs = unname(added)),
      class = c("fb_prior", "list")
    ))
    for (k in names(entries)) {
      entries[[k]]$meta$package_chosen <- TRUE
      hyper_ctrl[[k]] <- entries[[k]]
    }
  }

  list(
    hyper_ctrl = hyper_ctrl,
    record = list(
      filled = filled,
      covariate_scaled = covariate_scaled,
      covariate_sd = slope_sd[names(slope_sd) %in% names(added)]
    )
  )
}
