# fb_decision() --- the decision rule a selection is read off with.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 3.3 and 6.
#
# The argument this object replaces was `metric`, on `fb_select()`. Moving it
# here is not tidying. A decision has four coupled parts -- what is optimised,
# how a size is picked from the optimised curve, how the search over sizes is
# run, and whether that search is itself cross-validated -- and the four are
# only jointly meaningful. `rule = "tolerance"` with no `tolerance` is not a
# rule; `method = "posterior_threshold"` with no practical effect size is not
# a threshold. Both refuse here, at construction, where the message can name
# the missing part, rather than several minutes into a sampling run.
#
# Posterior draws under a continuous shrinkage prior contain no exact zeros,
# so "selected" is always a named decision and never a property of the draws.
# The three methods are three different named decisions, and the object
# carries which one so that no consumer has to infer it.


#' Declare how a selection decision is made
#'
#' The decision layer of [fb_select()]: which quantity is optimised, how a
#' submodel size is chosen from it, and how the search over submodels is run.
#' Validated at construction, so a rule that is missing the number it needs
#' refuses before anything is fitted.
#'
#' @details
#' The three methods answer three different questions.
#'
#' `"projection"` is the default where the reference model's structure
#' permits it. The reference posterior is projected onto each submodel and
#' the submodel's expected log predictive density is assessed on the same
#' folds the rest of the selection uses, with the standard error aggregated
#' at the resampling unit rather than at the row. Protected terms and the
#' base model's random terms are forced into every submodel.
#'
#' `"full"` returns the reference posterior with no selection, for a user who
#' wants shrinkage without a decision.
#'
#' `"posterior_threshold"` marks a candidate retained when the posterior
#' probability that the magnitude of its coefficient exceeds `effect_size` is
#' above `threshold`. It is a decision rule and nothing more: a continuous
#' shrinkage prior puts no mass at zero, so the quantity thresholded is not a
#' posterior inclusion probability and must not be read as one. `effect_size`
#' is on the scientific scale of the candidate column and has no default --
#' what counts as a practically important effect is a property of the
#' science, not of the software.
#'
#' The three rules operate on the elpd curve by submodel size. `"min"` takes
#' the best size. `"one_se"` takes the smallest size whose paired loss
#' against the best is within one standard error of that paired difference.
#' `"tolerance"` takes the smallest size losing less than `tolerance` nats
#' per unit against the reference model.
#'
#' @param method `"projection"` (default), `"full"` or
#'   `"posterior_threshold"`.
#' @param metric `"elpd"` (default), `"rmse"` or `"nlpd"`. The projection
#'   decision is on `elpd` in this release; the other two are scored by the
#'   out-of-fold path.
#' @param rule `"one_se"` (default), `"min"` or `"tolerance"`.
#' @param tolerance Nats per unit. Required by, and only used by,
#'   `rule = "tolerance"`.
#' @param search `"forward"` (default) or `"L1"`.
#' @param nterms_max Largest submodel the search considers. `NULL` (the
#'   default) uses every eligible term.
#' @param validate_search Whether the search itself is cross-validated.
#'   `TRUE` (the default) is the assessment that is not optimistic;
#'   `FALSE` is a fast exploratory mode and is labelled `optimistic` on the
#'   result, in its print method and on the certificate.
#' @param threshold Posterior probability cut, in (0, 1). Used by, and only
#'   by, `method = "posterior_threshold"`.
#' @param effect_size The practical effect size on the scientific scale, at
#'   or below which a coefficient is not of interest. Required by, and only
#'   used by, `method = "posterior_threshold"`.
#'
#' @returns An object of class `fb_decision`: a validated list carrying
#'   `method`, `metric`, `rule`, `tolerance`, `search`, `nterms_max`,
#'   `validate_search`, `threshold` and `effect_size`.
#'
#' @seealso [fb_select()], which takes this object; [fb_resamples()] for the
#'   folds the decision is assessed on.
#'
#' @examples
#' fb_decision()
#' fb_decision(rule = "tolerance", tolerance = 0.1)
#' fb_decision(method = "posterior_threshold", effect_size = 0.2)
#' @export
fb_decision <- function(
  method = c("projection", "full", "posterior_threshold"),
  metric = c("elpd", "rmse", "nlpd"),
  rule = c("one_se", "min", "tolerance"),
  tolerance = NULL,
  search = c("forward", "L1"),
  nterms_max = NULL,
  validate_search = TRUE,
  threshold = 0.5,
  effect_size = NULL
) {
  method <- match.arg(method)
  metric <- match.arg(metric)
  rule <- match.arg(rule)
  search <- match.arg(search)

  if (!is.logical(validate_search) || length(validate_search) != 1L ||
      is.na(validate_search)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_decision_invalid",
      message = paste0(
        "fb_decision(): `validate_search` must be TRUE or FALSE. FALSE is ",
        "the fast exploratory mode and is labelled optimistic wherever the ",
        "result is reported."
      )
    ))
  }

  if (identical(rule, "tolerance")) {
    if (!.fb_dec_is_number(tolerance) || tolerance < 0) {
      stop(.fb_refusal_condition(
        reason_code = "selection_decision_invalid",
        message = paste0(
          "fb_decision(): rule = \"tolerance\" needs a single non-negative ",
          "`tolerance`, in nats per resampling unit. There is no default: ",
          "how much predictive accuracy a smaller model may cost is the ",
          "question the rule exists to answer, not something to assume."
        )
      ))
    }
  } else if (!is.null(tolerance)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_decision_invalid",
      message = paste0(
        "fb_decision(): `tolerance` was given with rule = \"",
        rule,
        "\", which does not read it. Pass rule = \"tolerance\" as well, or ",
        "drop the argument."
      )
    ))
  }

  if (!.fb_dec_is_number(threshold) || threshold <= 0 || threshold >= 1) {
    stop(.fb_refusal_condition(
      reason_code = "selection_decision_invalid",
      message = paste0(
        "fb_decision(): `threshold` is a posterior probability and must lie ",
        "strictly inside (0, 1). Got: ",
        .fb_dec_show(threshold),
        "."
      )
    ))
  }

  if (identical(method, "posterior_threshold")) {
    if (!.fb_dec_is_number(effect_size) || effect_size < 0) {
      stop(.fb_refusal_condition(
        reason_code = "selection_decision_invalid",
        message = paste0(
          "fb_decision(): method = \"posterior_threshold\" needs a single ",
          "non-negative `effect_size` -- the magnitude on the scientific ",
          "scale at or below which a coefficient is not of interest. No ",
          "default is invented: a continuous shrinkage prior puts no mass ",
          "at zero, so P(|beta| > 0) is one whatever the data say, and a ",
          "rule built on it would report every candidate as retained."
        )
      ))
    }
  } else if (!is.null(effect_size)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_decision_invalid",
      message = paste0(
        "fb_decision(): `effect_size` was given with method = \"",
        method,
        "\", which does not read it. It belongs to ",
        "method = \"posterior_threshold\"."
      )
    ))
  }

  if (!is.null(nterms_max)) {
    if (!.fb_dec_is_number(nterms_max) || nterms_max < 1 ||
        nterms_max != round(nterms_max)) {
      stop(.fb_refusal_condition(
        reason_code = "selection_decision_invalid",
        message = paste0(
          "fb_decision(): `nterms_max` must be a single whole number of ",
          "terms, at least 1. Got: ",
          .fb_dec_show(nterms_max),
          "."
        )
      ))
    }
    nterms_max <- as.integer(nterms_max)
  }

  structure(
    list(
      method = method,
      metric = metric,
      rule = rule,
      tolerance = tolerance,
      search = search,
      nterms_max = nterms_max,
      validate_search = validate_search,
      threshold = threshold,
      effect_size = effect_size
    ),
    class = c("fb_decision", "list")
  )
}


#' @export
print.fb_decision <- function(x, ...) {
  cat("<fb_decision> method = \"", x$method, "\"\n", sep = "")
  if (identical(x$method, "projection")) {
    cat(
      "  metric  : ", x$metric,
      ", rule ", x$rule,
      if (identical(x$rule, "tolerance")) {
        paste0(" (", format(x$tolerance), " nats per unit)")
      } else {
        ""
      },
      "\n",
      sep = ""
    )
    cat(
      "  search  : ", x$search,
      ", nterms_max ",
      if (is.null(x$nterms_max)) "every eligible term" else x$nterms_max,
      "\n",
      sep = ""
    )
    cat(
      "  assess  : ",
      if (isTRUE(x$validate_search)) {
        "the search is cross-validated"
      } else {
        "optimistic -- the search is NOT cross-validated"
      },
      "\n",
      sep = ""
    )
  } else if (identical(x$method, "posterior_threshold")) {
    cat(
      "  retains a candidate when P(|beta| > ",
      format(x$effect_size),
      ") > ",
      format(x$threshold),
      "\n",
      sep = ""
    )
    cat(
      "  a decision rule, not an inclusion probability: a continuous\n",
      "  shrinkage prior puts no posterior mass at zero.\n",
      sep = ""
    )
  } else {
    cat("  the reference posterior, with no selection.\n", sep = "")
  }
  invisible(x)
}


# --- small validators ---------------------------------------------------- #

.fb_dec_is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.fb_dec_show <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  paste(format(utils::head(x, 3L)), collapse = ", ")
}
