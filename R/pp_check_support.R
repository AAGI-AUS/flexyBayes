# pp_check_support.R -- the posterior predictive check, or a refusal that
# names the diagnostics the fit can answer.
#
# `pp_check()` is registered for flexybayes fits on the bayesplot
# package's own generic via @exportS3Method (delayed S3 registration), so
# bayesplot stays in Suggests.
#
# A posterior predictive check simulates datasets from the fitted model
# and shows them against the observed response. That needs predictive
# draws. A brms-engine fit has them, and the call delegates to
# brms::pp_check(), whose default display overlays the densities of
# replicated datasets on the density of y. A nested Laplace approximation
# has none, and the refusal names the residual diagnostics instead of
# drawing a different display under the check's name.
#
# plot(fit, type = "pp_check") runs the same seam (see R/plot.R), so the
# two entry points cannot disagree about what a check is.

# .fb_check_bayesplot() --- the Suggests guard for the generic's home.
#
# brms::pp_check() builds its display with bayesplot, and dispatch to
# this method normally arrives through bayesplot's generic, so the guard
# is for the direct call: plot(fit, type = "pp_check") reaches the same
# code without loading bayesplot first.
#
# @noRd
# @keywords internal
.fb_check_bayesplot <- function() {
  if (.fb_bayesplot_available()) {
    return(invisible(TRUE))
  }
  stop(
    "pp_check() builds its display with the bayesplot package, which is ",
    "not installed. Install it with install.packages(\"bayesplot\"), or ",
    "use plot(fit, type = \"residuals\") for the residual diagnostics, ",
    "which need no extra package.",
    call. = FALSE
  )
}

# .fb_pp_check_impl() --- the one seam both entry points run.
#
# @noRd
# @keywords internal
.fb_pp_check_impl <- function(object, ...) {
  if (inherits(object$brms, "brmsfit")) {
    .fb_check_brms_delegation("pp_check()")
    .fb_check_bayesplot()
    return(brms::pp_check(object$brms, ...))
  }

  engine <- .fb_fit_engine(object)
  engine_clause <- if (identical(engine, "unknown")) {
    ". "
  } else {
    paste0(
      ", and its engine (",
      engine,
      ") returns an approximation to the posterior rather than simulated ",
      "datasets. "
    )
  }

  stop(.fb_refusal_condition(
    reason_code = "pp_check_requires_predictive_draws",
    message = paste0(
      "pp_check() draws replicated datasets from the posterior predictive ",
      "distribution and overlays them on the observed response, and this ",
      "fit carries no predictive draws to replicate from. The fit's class ",
      "is <",
      paste(class(object), collapse = ", "),
      ">",
      engine_clause,
      "The diagnostics it does answer are ",
      "plot(fit, type = \"residuals\"), residuals against fitted values ",
      "beside a normal quantile-quantile plot, and, on a fit carrying a ",
      "design index, plot(fit, type = \"variogram\"), the empirical ",
      "semivariance of the residuals over the array. For a posterior ",
      "predictive check on this model, re-fit it on the brms engine -- ",
      "flexybayes(..., backend = \"brms\")."
    ),
    engine = engine
  ))
}


#' Posterior predictive check for a flexyBayes fit
#'
#' A method for the \pkg{bayesplot} package's `pp_check()` generic. On a
#' fit from the **brms** engine the call delegates to
#' [brms::pp_check()]: datasets are simulated from the posterior
#' predictive distribution and shown against the observed response. The
#' default display overlays the densities of the replicated datasets on
#' the density of the response (`type = "dens_overlay"`); every other
#' \pkg{bayesplot} check type, and every argument of
#' [brms::pp_check()], is reachable through `...`.
#'
#' On a fit from the **INLA** engine the method refuses by name. A nested
#' Laplace approximation returns marginal densities, not simulated
#' datasets, so there is nothing to overlay -- and a display of observed
#' against fitted values shown under this name would be a different
#' diagnostic wearing the check's title. The refusal names the residual
#' displays the fit does answer: `plot(fit, type = "residuals")` and, on
#' a fit carrying a design index, `plot(fit, type = "variogram")`.
#'
#' `plot(fit, type = "pp_check")` runs this same code, so the two entry
#' points cannot disagree about what a check is.
#'
#' @param object A fitted `flexybayes` object.
#' @param ... Passed to [brms::pp_check()] on the brms path -- `type`,
#'   `ndraws`, `group` and the rest; ignored on the refusal path.
#' @returns On a brms-engine fit, the \pkg{ggplot2} object
#'   [brms::pp_check()] returns. Otherwise the method does not return: it
#'   raises a refusal of class
#'   `flexybayes_refusal_pp_check_requires_predictive_draws`.
#'
#' @seealso [plot.flexybayes()] for the other displays, including the
#'   residual diagnostics this refusal points at;
#'   [as_draws_df.flexybayes()] for the posterior itself;
#'   [summary.flexybayes()], whose `$converge` slot carries the engine's
#'   own diagnostics.
#'
#' @exportS3Method bayesplot::pp_check
pp_check.flexybayes <- function(object, ...) {
  .fb_pp_check_impl(object, ...)
}
