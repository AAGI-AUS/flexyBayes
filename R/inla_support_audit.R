# inla_support_audit.R -- Post-fit check that an INLA fit respects the
# bounds of its stated priors.
#
# A uniform prior on a standard deviation reaches INLA as an expression prior
# that is flat on its support and -1e10 outside it (.inla_uniform_sd_expr()
# in R/fb_prior.R). INLA can still return a mode where that density is -1e10
# (inst/KNOWN_ISSUES.md, issue 17). The audit reads each bound from INLA's own
# record of the priors it received (`all.hyper`) and compares the internal
# mode and the reported quantiles of every bounded log precision with it, on
# the scale the engine fitted. A fit outside its prior's support is refused.

# Two log(1.001): a log precision counts as inside its bound when the
# standard deviation it implies is within 0.1 per cent of the bound.
.FB_SUPPORT_THETA_TOL <- 2 * log1p(1e-3)

#' Bounds of the uniform expression priors in an INLA prior record
#'
#' @param prior_text Character vector of prior strings as `all.hyper`
#'   reports them.
#' @returns A data frame with columns `lower` and `upper`, the standard
#'   deviation bounds of each uniform expression prior; `upper` is `NA` for
#'   every other prior.
#' @noRd
#' @keywords internal
.fb_inla_expression_bounds <- function(prior_text) {
  text <- as.character(prior_text)
  uniform <- !is.na(text) &
    startsWith(text, "expression:") &
    grepl("[ :]U=[^;]+;", text) &
    grepl("(lb|lo)=-2\\*log\\(U\\);", text)

  upper <- rep(NA_real_, length(text))
  lower <- rep(0, length(text))
  upper[uniform] <- as.numeric(
    sub(".*[ :]U=([^;]+);.*", "\\1", text[uniform])
  )
  has_lower <- uniform & grepl("[ :]L=[^;]+;", text)
  lower[has_lower] <- as.numeric(
    sub(".*[ :]L=([^;]+);.*", "\\1", text[has_lower])
  )
  data.frame(lower = lower, upper = upper)
}

#' An empty support audit
#'
#' @returns A zero-row data frame with the columns of a support audit.
#' @noRd
#' @keywords internal
.fb_inla_support_audit_frame <- function() {
  data.frame(
    parameter = character(0),
    lower_sd = numeric(0),
    upper_sd = numeric(0),
    mode_sd = numeric(0),
    interval_lower_sd = numeric(0),
    interval_upper_sd = numeric(0),
    inside = logical(0),
    stringsAsFactors = FALSE
  )
}

#' Check that bounded priors pair with INLA's hyperparameters
#'
#' INLA lists its hyperparameters in the order of its prior record, so the
#' two are paired by position. The pairing is accepted only when the counts
#' agree, a latent term's hyperparameter carries that term's name, and a
#' likelihood's hyperparameter names the observations. The expression prior
#' is a density on INLA's internal theta, whatever the hyperparameter is
#' called (a Gaussian precision, a gamma precision parameter, a beta `phi`).
#'
#' @param rows The non-fixed rows of the prior read-back.
#' @param theta_names Row names of `internal.summary.hyperpar`.
#' @param bounded Positions of the rows under a uniform expression prior.
#' @returns `TRUE`, invisibly; otherwise an internal error.
#' @noRd
#' @keywords internal
.fb_inla_support_pairing <- function(rows, theta_names, bounded) {
  paired <- nrow(rows) == length(theta_names)
  if (paired) {
    named <- vapply(
      bounded,
      function(i) {
        switch(
          rows$section[[i]],
          random = endsWith(theta_names[[i]], paste0(" for ", rows$tag[[i]])),
          family = grepl(" observations(\\[[0-9]+\\])?$", theta_names[[i]]),
          FALSE
        )
      },
      logical(1L)
    )
    paired <- all(named)
  }
  if (!paired) {
    stop(
      ".fb_inla_support_audit(): INLA's hyperparameters (",
      paste(theta_names, collapse = ", "),
      ") could not be paired with its record of the priors it received, ",
      "so the bounded priors cannot be checked. This is a flexyBayes ",
      "internal contract violation, not a user error.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The factor that maps a standard deviation back to the response scale
#'
#' @param theta_name An internal hyperparameter name.
#' @param standardisation The response standardisation, or `NULL`.
#' @param covariate The covariate standardisation, or `NULL`. A bound on an
#'   uncorrelated slope's standard deviation is stated in response units per
#'   covariate unit, so it is reported in those units whatever the engine
#'   received.
#' @returns A positive number; 1 when the fit was not standardised.
#' @noRd
#' @keywords internal
.fb_support_sd_factor <- function(
  theta_name,
  standardisation,
  covariate = NULL
) {
  if (is.null(standardisation) && is.null(covariate)) {
    return(1)
  }
  precision_name <- sub("^Log precision", "Precision", theta_name)
  factor <- 1
  if (!is.null(standardisation)) {
    factor <- factor /
      sqrt(.fb_inla_hyper_factor(precision_name, standardisation))
  }
  if (!is.null(covariate)) {
    factor <- factor /
      sqrt(.fb_inla_covariate_hyper_factor(precision_name, covariate))
  }
  factor
}

#' Audit an INLA fit against the bounds of its stated priors
#'
#' @param inla_fit The object `INLA::inla()` returned, before any
#'   back-transformation from a standardised response.
#' @param standardisation The response standardisation the fit ran on, or
#'   `NULL`. It maps only the reported standard deviations; the comparison
#'   stays on the scale INLA fitted.
#' @param covariate The covariate standardisation the fit ran on, or `NULL`.
#'   It maps the reported standard deviations of the uncorrelated slope
#'   fields, for the same reason and in the same way.
#' @returns A data frame with one row per hyperparameter under a bounded
#'   prior: `parameter`, `lower_sd`, `upper_sd`, `mode_sd`,
#'   `interval_lower_sd`, `interval_upper_sd` and `inside`. Standard
#'   deviations are on the scale of the response; the interval is the one
#'   implied by the 2.5 and 97.5 per cent quantiles of the log precision.
#'   The frame has zero rows when no prior is bounded.
#' @noRd
#' @keywords internal
.fb_inla_support_audit <- function(
  inla_fit,
  standardisation = NULL,
  covariate = NULL
) {

  # Priors INLA received, in the order of its hyperparameters -------------

  rows <- .fb_manifest_inla_readback(inla_fit)[["rows"]]
  rows <- rows[rows$section != "fixed", , drop = FALSE]
  internal <- inla_fit[["internal.summary.hyperpar"]]
  theta_names <- rownames(internal) %||% character(0)
  bounds <- .fb_inla_expression_bounds(rows$prior)
  bounded <- which(!is.na(bounds$upper))

  if (length(bounded) == 0L) {
    return(.fb_inla_support_audit_frame())
  }
  .fb_inla_support_pairing(rows, theta_names, bounded)

  # Compare on the log-precision scale ------------------------------------

  theta_min <- -2 * log(bounds$upper[bounded])
  theta_max <- rep(Inf, length(bounded))
  two_sided <- bounds$lower[bounded] > 0
  theta_max[two_sided] <- -2 * log(bounds$lower[bounded][two_sided])

  mode <- internal[bounded, "mode"]
  q_low <- internal[bounded, "0.025quant"]
  q_high <- internal[bounded, "0.975quant"]
  within <- function(v) {
    is.finite(v) &
      v >= theta_min - .FB_SUPPORT_THETA_TOL &
      v <= theta_max + .FB_SUPPORT_THETA_TOL
  }
  inside <- within(mode) & within(q_low) & within(q_high)

  # Report standard deviations on the scale of the response ---------------

  sd_factor <- vapply(
    theta_names[bounded],
    .fb_support_sd_factor,
    numeric(1L),
    standardisation = standardisation,
    covariate = covariate
  )
  audit <- data.frame(
    parameter = theta_names[bounded],
    lower_sd = bounds$lower[bounded] * sd_factor,
    upper_sd = bounds$upper[bounded] * sd_factor,
    mode_sd = exp(-mode / 2) * sd_factor,
    interval_lower_sd = exp(-q_high / 2) * sd_factor,
    interval_upper_sd = exp(-q_low / 2) * sd_factor,
    inside = inside,
    stringsAsFactors = FALSE
  )
  rownames(audit) <- NULL
  audit
}

#' Refuse an INLA fit whose audit found a value outside a bound
#'
#' @param audit The data frame `.fb_inla_support_audit()` returned.
#' @param route The emit route: `"per_row"`, `"aggregated_count"` or
#'   `"aggregated_gaussian"`.
#' @returns `audit`, invisibly, when every value lies inside its bound;
#'   otherwise a `prior_support_violated` refusal.
#' @noRd
#' @keywords internal
.fb_inla_enforce_support <- function(audit, route) {
  outside <- audit[!audit$inside, , drop = FALSE]
  if (nrow(outside) == 0L) {
    return(invisible(audit))
  }

  detail <- sprintf(
    "%s: standard deviation %s at the mode, interval %s to %s, bound %s to %s",
    sub("^Log precision for ", "", outside$parameter),
    signif(outside$mode_sd, 4L),
    signif(outside$interval_lower_sd, 4L),
    signif(outside$interval_upper_sd, 4L),
    signif(outside$lower_sd, 4L),
    signif(outside$upper_sd, 4L)
  )
  stop(.fb_refusal_condition(
    reason_code = "prior_support_violated",
    message = paste0(
      "INLA returned a hyperparameter outside the bound of its stated ",
      "prior, so the fit is refused: ",
      paste(detail, collapse = "; "),
      ". A posterior outside its prior's support is not a fit of the ",
      "model written (INLA ",
      as.character(utils::packageVersion("INLA")),
      ", ",
      route,
      " route; inst/KNOWN_ISSUES.md, issue 17). Fit the same model with ",
      "backend = \"brms\"."
    ),
    backend = "inla",
    route = route,
    audit = audit
  ))
}
