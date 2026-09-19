# inla_default_priors.R -- the priors the package chooses for an INLA fit.
#
# A prior the caller states reaches INLA exactly as stated. This file decides
# the rest: the prior on every quantity the caller left unstated. Each one is
# written in the units of the quantity it governs, so multiplying the
# response or a covariate by a constant multiplies the prior with it and the
# fit stays one model:
#
#   * a random-effect standard deviation with no stated prior takes the
#     half-normal default, half_normal(scale = s), with s = 5 sd(y) on an
#     identity link and a constant on a log or logit link
#     (.default_half_normal_scale());
#   * the standard deviation of an uncorrelated slope is in response units
#     per covariate unit, so its default scale is s / sd(x), on INLA and on
#     brms (.fb_default_covariate_scaled());
#   * a numeric fixed-effect column with no stated prior takes
#     normal(0, K sd(y) / sd(x_j)), with sd(y) replaced by 1 and a smaller K
#     on a log or logit link.
#
# A Gaussian identity-link fit runs on a standardised response
# (R/inla_standardise.R), and each of these is stated in original units and
# restated on those coordinates with the stated priors. Before this file,
# anything unstated took INLA's own default on the standardised coordinates,
# which is scale-free in the response and blind to every covariate, and
# three failures followed. Each was silent: INLA reported mode.status 0 and
# the numerical-confirm gate passed.
#
#   * An unstated random-intercept SD degenerated. On a 400-row design with
#     four blocks (true SD 2) and only the residual prior stated, 20 of 20
#     identical per-row fits returned a block SD above 1e6 (median 1e7); on
#     the aggregated route 10 of 20 did and 5 more failed in the engine, and
#     the streamed route, which attached no prior at all, gave 10 and 4.
#   * lme4::sleepstudy with `(Dx || Subject)` and Dx = Days / 1000: the
#     default bound 5 sd(y) = 282 on the slope SD sat just above the fitted
#     279 (0.279 per day, against 6.45 at Dx = Days).
#   * The same fit carried INLA's prec = 0.001 on the Dx coefficient, on the
#     standardised coordinates, and returned 8.64 per day against lme4's
#     10.4673; widening that one prior restored 10.4673. On Poisson and
#     Bernoulli fits with x / 1000 the coefficient fell to 41 and 4.4 per
#     cent of its value at x.
#
# A prior the caller states is never rescaled by a covariate: a response-unit
# prior such as `prior_fixed_sd = 100` is exactly that on every coefficient,
# and on a covariate recorded in small units it can bind.

# The fixed-effect multiplier K, by link class (.fb_default_link_class()).
#
# Identity link: 1e5. The default is meant to leave a coefficient where the
# data put it, as REML and brms's flat default do, and the width costs
# nothing there: with a design of full column rank the Gaussian posterior of
# the fixed effects is proper under a flat prior. Measured on INLA 26.8.7,
# lme4::sleepstudy, `Reaction ~ Dx + (Dx || Subject)` with Dx = Days * c,
# the slope bound scaled, and the coefficient read as coef(Dx) * c against
# lme4's REML 10.4672859596 (three identical calls at each of c = 0.001, 1
# and 1000):
#
#   K           at c = 1         distance from lme4   range across c
#   10          10.4665088       7.8e-4               1.1e-6
#   sqrt(1000)  10.4672082       7.8e-5               1.3e-7
#   100         10.4672782       7.8e-6               1.1e-8
#   1000        10.46728588      7.8e-8               1.1e-10
#   1e4         10.4672859588    7.8e-10              1.1e-12
#   1e5         10.4672859596    7.7e-12              0
#
# The unscaled default this replaces gave 10.4672766 at c = 1 (9.4e-6 from
# lme4), with a run-to-run range of 1.5e-10 over 20 identical calls, and
# 8.64 at c = 0.001. K = 1e5 is the smallest power of ten whose distance from
# lme4 is inside that run-to-run range at all three scales, and at c = 1 it
# shrinks less than the default it replaces.
#
# Log and logit links: sqrt(1000), about 31.6. On a covariate with unit
# standard deviation this is INLA's own prec = 0.001, so a fit on a
# unit-scale covariate does not change (Poisson and Bernoulli random-intercept
# fits, 30 groups of 8, standard-normal x: 0.508836677 against 0.508836662,
# and 0.231425423 against 0.231425542). A wider value is not free on these
# links. On a perfectly separated Bernoulli design (200 rows, 20 groups) the
# coefficient was 70 under INLA's default and 76 at sqrt(1000); at 100 the
# INLA engine failed, and at 1000, 1e4 and 1e6 it returned 369, 764 and 1476
# -- every converged fit with mode.status 0 and the numerical-confirm gate
# passing.
.FB_INLA_FIXED_SD_MULTIPLIER <- c(
  identity = 1e5,
  log = sqrt(1000),
  logit = sqrt(1000)
)

# K for a model, read off its link class.
.fb_inla_fixed_sd_multiplier <- function(fb) {
  .FB_INLA_FIXED_SD_MULTIPLIER[[.fb_default_link_class(fb$family, fb$link)]]
}


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

#' The rows a covariate-scaled default is learned on
#'
#' Every row whose working response is observed. A row with a missing
#' response adds no term to the likelihood, so its covariate values carry no
#' information about the effect a default is scaled for. Two routes hand the
#' emit such rows: a fit with `na_action = "augment"`, and an `fb_select()`
#' fold refit, which carries its assessment rows with the response set to
#' `NA`. Read over every row, the fold's default moved with the covariates of
#' the rows the fold is scored on. The response-scaled defaults already read
#' the observed rows, so both kinds of default now read the same rows.
#'
#' @param fb An `fb_terms` object.
#' @param data A data frame.
#' @returns A logical vector over the rows of `data`; all `TRUE` when the
#'   response cannot be resolved row by row.
#' @noRd
#' @keywords internal
.fb_inla_default_rows <- function(fb, data) {
  y <- .fb_working_response(fb, data)
  if (is.matrix(y) && nrow(y) == nrow(data)) {
    return(!apply(is.na(y), 1L, any))
  }
  if (is.null(y) || length(y) != nrow(data)) {
    return(rep(TRUE, nrow(data)))
  }
  !is.na(y)
}

#' The standard deviation of each uncorrelated slope's covariate
#'
#' Read over the rows `.fb_inla_default_rows()` names.
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
  rows <- .fb_inla_default_rows(fb, data)
  for (t in fb$random_terms %||% list()) {
    if (!identical(t$type %||% "", "simple_slope_uncor")) {
      next
    }
    x <- data[[t$slope_var]][rows]
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

# A half-normal standard-deviation prior for a quantity divided by `s`.
.fb_inla_divide_sd_spec <- function(spec, s) {
  spec$spec$args$scale <- spec$spec$args$scale / s
  spec
}

#' The auto-default's slope priors, divided by their covariates' SDs
#'
#' The standard deviation of an uncorrelated slope is in response units per
#' covariate unit, so the auto-default's scale for it is divided by the
#' standard deviation of the covariate, read over the rows
#' `.fb_inla_default_rows()` names. Both engines call this, so the slope slot
#' carries the same scale on INLA and on brms. A prior the caller wrote is
#' never divided.
#'
#' @param fb An `fb_terms` object.
#' @param data The data the fit is built from, or `NULL`.
#' @returns A list: `covariate_scaled`, the divided specs named by group, and
#'   `covariate_sd`, the divisor of each.
#' @noRd
#' @keywords internal
.fb_default_covariate_scaled <- function(fb, data) {
  stated <- fb$priors
  scaled <- list()
  from_default <- inherits(stated, "fb_prior") &&
    !is.null(attr(stated, "fb_prior_default_basis"))
  slope_sd <- if (from_default) .fb_inla_slope_covariate_sd(fb, data) else
    numeric(0)
  for (spec in if (from_default) stated$specs else list()) {
    k <- if (identical(spec$target$type, "sd")) spec$target$group else ""
    if (k %in% names(slope_sd) && identical(spec$spec$family, "half_normal")) {
      scaled[[k]] <- .fb_inla_divide_sd_spec(spec, slope_sd[[k]])
    }
  }
  list(
    covariate_scaled = scaled,
    covariate_sd = slope_sd[names(slope_sd) %in% names(scaled)]
  )
}

#' The auto-default with its slope priors divided, as brms receives it
#'
#' @param fb An `fb_terms` object.
#' @param scaled The `covariate_scaled` list of
#'   `.fb_default_covariate_scaled()`.
#' @returns `fb$priors` with each divided spec in place of its original, the
#'   attributes kept.
#' @noRd
#' @keywords internal
.fb_default_prior_scaled <- function(fb, scaled) {
  prior <- fb$priors
  if (!length(scaled)) {
    return(prior)
  }
  for (i in seq_along(prior$specs)) {
    spec <- prior$specs[[i]]
    k <- if (identical(spec$target$type, "sd")) spec$target$group else ""
    if (k %in% names(scaled)) {
      prior$specs[[i]] <- scaled[[k]]
    }
  }
  prior
}

#' Every INLA variance-component prior, with what the caller left unstated
#' filled in
#'
#' Two additions to the `priors_to_inla()` list the emit built from what was
#' stated. First, when the fit carries a prior of any kind -- an `fb_prior()`,
#' partial or whole, or the legacy `prior_vc_sd` scalar -- each target of
#' `.fb_inla_sd_targets()` with no entry takes the half-normal default.
#' A fit with no prior at all is left alone: every public route resolves one
#' before it reaches an emit (flexybayes() builds the auto-default and
#' flexybayes_stream() attaches it), so that case is a direct internal call.
#' Second, a slope scale the package chose, filled here or built by the
#' auto-default, is divided by the standard deviation of its covariate. A
#' scale the caller wrote is not.
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
  slope_sd <- .fb_inla_slope_covariate_sd(fb, data)

  # On an offset model the scale is read from y - o, as the standardisation
  # is (R/inla_standardise.R).
  if (is.null(resp_sd) && length(fb$offsets) && is.data.frame(data)) {
    resp_sd <- .fb_response_sd(fb, data)
  }

  filled <- list()
  if (has_prior) {
    targets <- .fb_inla_sd_targets(fb)
    keys <- setdiff(c(targets$shared, targets$vm_ped), names(hyper_ctrl))
    if (length(keys)) {
      scale <- .default_half_normal_scale(
        data,
        fb$response,
        fb$family,
        fb$link,
        resp_sd = resp_sd
      )$scale
      for (k in keys) {
        spec <- list(
          target = list(type = "sd", group = k),
          spec = list(family = "half_normal", args = list(scale = scale))
        )
        if (k %in% names(slope_sd)) {
          spec <- .fb_inla_divide_sd_spec(spec, slope_sd[[k]])
        }
        filled[[k]] <- spec
      }
    }
  }

  covariate_scaled <- .fb_default_covariate_scaled(fb, data)$covariate_scaled

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


# --- fixed effects --------------------------------------------------------- #

# The fixed part of the model as a one-sided formula over its term labels,
# or NULL when it has none.
.fb_inla_fixed_design_formula <- function(fb) {
  labels <- vapply(
    fb$fixed_terms %||% list(),
    function(term) {
      term$label %||%
        term$var %||%
        paste(term$vars %||% character(0), collapse = ":")
    },
    character(1)
  )
  labels <- labels[nzchar(labels)]
  if (!length(labels)) {
    return(NULL)
  }
  stats::as.formula(paste0("~ ", paste(labels, collapse = " + ")))
}

# The response unit a fixed-effect default is written in: sd(y) where the
# linear predictor carries the response's units, 1 on a log or logit link.
# sd(y - o) on a model with an offset o.
.fb_inla_fixed_prior_unit <- function(fb, data) {
  if (!identical(.fb_default_link_class(fb$family, fb$link), "identity")) {
    return(1)
  }
  y <- suppressWarnings(as.numeric(.fb_working_response(fb, data)))
  y <- y[is.finite(y)]
  s <- if (length(y) >= 2L) stats::sd(y) else NA_real_
  if (!is.finite(s) || s <= 0) 1 else s
}

#' The covariate-scaled default prior SD of each numeric fixed effect
#'
#' `normal(0, K u / sd(x_j))` for every design column that carries a numeric
#' variable -- a continuous term, a factor-by-covariate interaction, a
#' transformed covariate -- where `u` is `.fb_inla_fixed_prior_unit()` and
#' `K` is `.fb_inla_fixed_sd_multiplier()`.
#'
#' A column built from factors alone is left to INLA's own default. Its 0/1
#' coding carries no covariate unit, and its coefficient is a difference
#' between level means, in response units: INLA's `prec = 0.001` on the
#' standardised coordinates is `normal(0, 31.6 sd(y))` on the original ones,
#' already in those units. Scaling a contrast by the standard deviation of
#' its indicator would make the prior depend on how many rows each level
#' has, which is not a unit. The intercept keeps INLA's flat prior.
#'
#' Empty when the caller supplied `prior_fixed_sd`, which is a statement
#' about every coefficient; a coefficient named by a `b()` row keeps that
#' row. Empty too when `data` is not a data frame: the aggregated routes
#' hold the design as cell contrasts, and a numeric column never reaches
#' them.
#'
#' `sd(x_j)` is read over the rows `.fb_inla_default_rows()` names. The
#' design is built on every row, so its columns are the ones the fit
#' carries.
#'
#' @param fb An `fb_terms` object.
#' @param data The data the fit is built from.
#' @returns A named numeric vector of prior standard deviations, in the
#'   original units, keyed by design column name.
#' @noRd
#' @keywords internal
.fb_inla_default_fixed_sd <- function(fb, data) {
  if (!is.data.frame(data) || .fb_prior_scalar_supplied(fb, "fixed_sd")) {
    return(numeric(0))
  }
  form <- .fb_inla_fixed_design_formula(fb)
  if (is.null(form)) {
    return(numeric(0))
  }
  mf <- stats::model.frame(form, data = data)
  tt <- attr(mf, "terms")
  mm <- stats::model.matrix(tt, mf)
  # The design rows model.frame() kept, as rows of `data`, and of those the
  # ones with an observed response.
  kept <- seq_len(nrow(data))
  omitted <- attr(mf, "na.action")
  if (!is.null(omitted)) {
    kept <- kept[-omitted]
  }
  observed <- .fb_inla_default_rows(fb, data)[kept]
  factors <- attr(tt, "factors")
  assign <- attr(mm, "assign")
  categorical <- vapply(
    mf,
    function(v) is.factor(v) || is.character(v) || is.logical(v),
    logical(1)
  )
  named <- if (inherits(fb$priors, "fb_prior")) {
    unlist(lapply(fb$priors$specs, function(s) {
      if (identical(s$target$type, "b")) s$target$name else NULL
    }))
  } else {
    character(0)
  }
  unit <- .fb_inla_fixed_prior_unit(fb, data)
  multiplier <- .fb_inla_fixed_sd_multiplier(fb)

  out <- numeric(0)
  for (j in seq_len(ncol(mm))) {
    a <- assign[[j]]
    nm <- colnames(mm)[[j]]
    if (a == 0L || nm %in% named) {
      next
    }
    vars <- rownames(factors)[factors[, a] > 0]
    if (isTRUE(all(categorical[vars]))) {
      next
    }
    s <- stats::sd(mm[observed, j])
    if (!is.finite(s) || s <= 0) {
      next
    }
    out[[nm]] <- multiplier * unit / s
  }
  out
}
