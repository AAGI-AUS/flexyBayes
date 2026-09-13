# inla_default_priors.R -- the priors the package chooses for an INLA fit.
#
# A prior the caller states reaches INLA exactly as stated. This file decides
# the rest: the prior on a quantity the caller left unstated, written in the
# units of the quantity it governs, so multiplying the response by a
# constant multiplies the prior with it and the fit stays one model:
#
#   * a random-effect standard deviation with no stated prior takes the
#     bounded uniform default, uniform(0, U), with U = 5 sd(y) on an identity
#     link and a constant on a log or logit link (.default_uniform_scale()).
#
# A Gaussian identity-link fit runs on a standardised response
# (R/inla_standardise.R), and the default is stated in original units and
# restated on those coordinates with the stated priors. Before this file, an
# unstated standard deviation took INLA's own default on the standardised
# coordinates, and it degenerated without a sign: INLA reported mode.status 0
# and the numerical-confirm gate passed. On a 400-row design with four blocks
# (true SD 2) and only the residual prior stated, 20 of 20 identical per-row
# fits returned a block SD above 1e6 (median 1e7); on the aggregated route 10
# of 20 did and 5 more failed in the engine.


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

#' Every INLA variance-component prior, with what the caller left unstated
#' filled in
#'
#' When the fit carries a prior of any kind -- an `fb_prior()`, partial or
#' whole, or the legacy `prior_vc_sd` scalar -- each target of
#' `.fb_inla_sd_targets()` with no entry in the `priors_to_inla()` list the
#' emit built from what was stated takes the bounded uniform default. A fit
#' with no prior at all is left alone: every public route resolves one
#' before it reaches an emit (flexybayes() builds the auto-default and
#' flexybayes_stream() attaches it), so that case is a direct internal call.
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
#'   holds `filled`, a list of `fb_prior` specs named by group.
#' @noRd
#' @keywords internal
.fb_inla_complete_priors <- function(fb, hyper_ctrl, data, resp_sd = NULL) {
  stated <- fb$priors
  has_prior <- inherits(stated, "fb_prior") ||
    .fb_prior_scalar_supplied(fb, "vc_sd")

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
        filled[[k]] <- list(
          target = list(type = "sd", group = k),
          spec = list(family = "uniform", args = list(lower = 0, upper = upper))
        )
      }
    }
  }

  if (length(filled)) {
    entries <- priors_to_inla(structure(
      list(specs = unname(filled)),
      class = c("fb_prior", "list")
    ))
    for (k in names(entries)) {
      entries[[k]]$meta$package_chosen <- TRUE
      hyper_ctrl[[k]] <- entries[[k]]
    }
  }

  list(hyper_ctrl = hyper_ctrl, record = list(filled = filled))
}
