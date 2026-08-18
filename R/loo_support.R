# loo_support.R -- approximate leave-one-out cross-validation, where the
# fit carries what it needs.
#
# `loo()` is registered for flexybayes fits on the loo package's own
# generic via @exportS3Method (delayed S3 registration), so loo stays in
# Suggests. What the method can do depends on what the engine stored:
#
#   * a brms-engine fit carries the sampler's draws and Stan's pointwise
#     log-likelihood, so the call passes straight through to
#     brms::loo(), which returns loo's own object;
#   * an INLA fit is a nested Laplace approximation. It carries marginal
#     densities, not draws of the log-likelihood, so PSIS-LOO has no
#     quantity to importance-weight. The refusal names the information
#     criteria the fit did compute instead of returning something that
#     would read as a leave-one-out estimate;
#   * a fit carrying no posterior at all reaches the sibling refusal
#     fit_lacks_posterior_draws, so the two states are distinguishable
#     by the condition class rather than by reading the message.

# .fb_check_brms_delegation() --- the package a passthrough needs.
#
# A fit carrying a brmsfit normally travels with the package that built
# it. A fit saved and read back on a machine without brms does not, and
# the delegation would fail inside brms's namespace rather than here.
#
# @noRd
# @keywords internal
.fb_check_brms_delegation <- function(what) {
  if (.fb_brms_available()) {
    return(invisible(TRUE))
  }
  stop(
    what,
    " on a brms-engine fit is answered by brms itself, and the brms ",
    "package is not installed. Install it with install.packages(\"brms\"), ",
    "or work from the draws with posterior::as_draws_df(fit).",
    call. = FALSE
  )
}

# .fb_information_criteria_line() --- what an INLA fit does carry.
#
# Built from the fit rather than recited: the per-row INLA emit asks for
# DIC and WAIC at fit time, but the aggregated emits do not, so a fit
# from that route carries neither and must not be told it does.
#
# @noRd
# @keywords internal
.fb_information_criteria_line <- function(fit) {
  waic <- fit$inla$waic$waic
  dic <- fit$inla$dic$dic
  have_waic <- is.numeric(waic) && length(waic) == 1L && is.finite(waic)
  have_dic <- is.numeric(dic) && length(dic) == 1L && is.finite(dic)

  if (have_waic && have_dic) {
    return(paste0(
      "This fit does carry the two information criteria INLA computed ",
      "at fit time: WAIC ",
      format(waic, digits = 6L),
      " at `fit$inla$waic$waic`, and DIC ",
      format(dic, digits = 6L),
      " at `fit$inla$dic$dic`. Neither is a leave-one-out estimate, and ",
      "neither carries the Pareto-k diagnostics loo() reports."
    ))
  }
  if (have_waic || have_dic) {
    one <- if (have_waic) {
      paste0("WAIC ", format(waic, digits = 6L), " at `fit$inla$waic$waic`")
    } else {
      paste0("DIC ", format(dic, digits = 6L), " at `fit$inla$dic$dic`")
    }
    return(paste0(
      "This fit does carry one information criterion INLA computed at ",
      "fit time: ",
      one,
      ". It is not a leave-one-out estimate."
    ))
  }
  paste0(
    "This fit carries no information criterion either: the aggregated ",
    "route fits on cell-level sufficient statistics and asks INLA for ",
    "neither WAIC nor DIC, so there is no stored summary to read in ",
    "place of the leave-one-out estimate. Re-fit with aggregate = FALSE ",
    "for a per-row fit carrying WAIC at `fit$inla$waic$waic` and DIC at ",
    "`fit$inla$dic$dic`."
  )
}

# .fb_refuse_loo() --- the two refusals, chosen by what the fit holds.
#
# @noRd
# @keywords internal
.fb_refuse_loo <- function(x) {
  engine <- .fb_fit_engine(x)
  has_posterior <- !is.null(x$inla) ||
    !is.null(x$greta$draws) ||
    !is.null(x$draws)

  if (!has_posterior) {
    stop(.fb_refusal_condition(
      reason_code = "fit_lacks_posterior_draws",
      message = paste0(
        "loo() cannot estimate leave-one-out predictive accuracy for this ",
        "fit: it carries no posterior at all -- no draws slot and no ",
        "engine marginals. The fit's class is <",
        paste(class(x), collapse = ", "),
        ">. Check the fit completed (return_code = FALSE). A fit that ",
        "does carry a posterior, but not the pointwise log-likelihood ",
        "loo() reads, is refused as loo_requires_sampler_draws instead."
      )
    ))
  }

  detail <- if (identical(engine, "inla")) {
    paste0(
      "This fit was estimated by INLA's nested Laplace approximation, ",
      "which returns marginal densities rather than draws of the ",
      "log-likelihood. ",
      .fb_information_criteria_line(x)
    )
  } else {
    paste0(
      "This fit carries a posterior, but not the log-likelihood of each ",
      "observation at each draw, which is the quantity PSIS-LOO ",
      "importance-weights."
    )
  }

  stop(.fb_refusal_condition(
    reason_code = "loo_requires_sampler_draws",
    message = paste0(
      "loo() estimates the expected log pointwise predictive density by ",
      "importance-weighting the log-likelihood of each observation at ",
      "each posterior draw, and this fit does not store that quantity, ",
      "so there is nothing to leave one observation out of. ",
      detail,
      " For a leave-one-out estimate on this model, re-fit it on the ",
      "brms engine -- flexybayes(..., backend = \"brms\") -- whose ",
      "sampler stores the pointwise log-likelihood loo() reads. The ",
      "sibling refusal fit_lacks_posterior_draws is raised when a fit ",
      "carries no posterior at all."
    ),
    engine = engine
  ))
}


#' Approximate leave-one-out cross-validation for a flexyBayes fit
#'
#' A method for the \pkg{loo} package's `loo()` generic. On a fit from
#' the **brms** engine the call passes through to [brms::loo()], which
#' computes PSIS-LOO from the pointwise log-likelihood Stan stored, and
#' returns \pkg{loo}'s own object -- `elpd_loo`, `p_loo`, `looic` and the
#' Pareto-k diagnostics, unchanged.
#'
#' On a fit from the **INLA** engine the method refuses by name. INLA
#' returns a nested Laplace approximation to the posterior, not draws of
#' the log-likelihood, so there is no quantity for importance sampling to
#' reweight and no leave-one-out estimate to report. The refusal names
#' the information criteria the fit did compute -- WAIC at
#' `fit$inla$waic$waic` and DIC at `fit$inla$dic$dic` on a per-row fit --
#' so the alternative is in the message rather than left to be found.
#' Neither is a leave-one-out estimate and neither carries Pareto-k
#' diagnostics.
#'
#' A fit that carries no posterior at all is refused as
#' `fit_lacks_posterior_draws` instead, so the two states are told apart
#' by the condition class rather than by reading the message.
#'
#' @param x A fitted `flexybayes` object.
#' @param ... Passed to [brms::loo()] on the brms path; ignored on the
#'   refusal paths.
#' @returns On a brms-engine fit, the `loo` object [brms::loo()] returns.
#'   Otherwise the method does not return: it raises a refusal of class
#'   `flexybayes_refusal_loo_requires_sampler_draws`, or
#'   `flexybayes_refusal_fit_lacks_posterior_draws` when the fit carries
#'   no posterior.
#'
#' @seealso [summary.flexybayes()], whose `$converge` slot carries the
#'   engine's own diagnostics; [as_draws_df.flexybayes()] for the
#'   posterior itself; [triangulate()] for comparing two engines' fits of
#'   one model; [fb_refusals()] for the refusal vocabulary.
#'
#' @exportS3Method loo::loo
loo.flexybayes <- function(x, ...) {
  if (inherits(x$brms, "brmsfit")) {
    .fb_check_brms_delegation("loo()")
    return(brms::loo(x$brms, ...))
  }
  .fb_refuse_loo(x)
}
