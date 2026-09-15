# fb_resamples() --- the resampling scheme for a selection.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 3.3, made
# public unchanged by SPEC_bayesian_selection_v5_2026-09-12.md section 3.3.
#
# Why the order of the resampling contract is not an implementation detail.
# An eigen or Cholesky rotation mixes observations within a covariance block,
# so a row of the rotated design is a contrast among original rows rather
# than an observation. Splitting rotated rows therefore leaves cluster mates
# on both sides of the split under a different name, and the resulting score
# answers an easier question than the one asked. Folds are assigned first, on
# the original units, and every rotation happens inside a fold. The reader of
# this object is `.fb_sel_folds()` in R/selection_resamples.R.


# fb_resamples() --- the resampling scheme, as an object rather than as a
# scattering of arguments.
#
# Normative (contract section 3.3): a `unit` alone is not a specification. It
# says what may not be split, not how the splitting is done, so naming one
# without also naming a scheme or an explicit `fold_id` refuses. The defaults
# of `v` and `repeats` are the unit-free defaults and do not silently supply a
# scheme to a grouped design.
#
#' Declare the resampling scheme for a selection
#'
#' The folds a selection is scored on, as a validated object rather than as a
#' scattering of arguments to [fb_select()]. Folds fall on the original
#' observational units before any transformation, factor levels are frozen
#' inside the training fold, and an assessment fold carrying a level the
#' training fold never saw refuses by name rather than being dropped.
#'
#' A `unit` alone is not a specification. It says what may not be split
#' across a fold boundary, not how the splitting is done, so naming one
#' without also naming `v`, `repeats` or an explicit `fold_id` refuses with
#' `selection_resamples_underspecified`. The defaults of `v` and `repeats`
#' are the unit-free defaults and do not silently supply a scheme to a
#' grouped design.
#'
#' @param unit Columns whose rows must stay together in a fold, as a
#'   one-sided formula (`~ site`) or a character vector. `NULL` (the
#'   default) splits rows.
#' @param v Number of folds, at least 2.
#' @param repeats Number of repetitions of the whole partition, at least 1.
#'   Repeats are scored as repeated measurements of the same units, never as
#'   independent evidence.
#' @param fold_id An explicit partition, one entry per row of the data. Given
#'   together with `v`, or with `repeats > 1`, it refuses: an explicit
#'   partition and a scheme that would build one are two answers to the same
#'   question.
#' @param strata Columns to stratify the fold assignment on, as a one-sided
#'   formula or a character vector. A stratifier that varies inside a unit
#'   refuses.
#'
#' @returns An object of class `fb_resamples`: a list carrying `unit`, `v`,
#'   `repeats`, `fold_id`, `strata` and the resolved `scheme`
#'   (`"v_fold"` or `"explicit"`).
#'
#' @seealso [fb_select()], which takes this object; [fb_decision()] for the
#'   decision rule scored on these folds.
#'
#' @examples
#' fb_resamples()
#' fb_resamples(unit = ~ site, v = 5L)
#' @export
fb_resamples <- function(
  unit = NULL,
  v = 10L,
  repeats = 1L,
  fold_id = NULL,
  strata = NULL
) {
  # Captured before anything is assigned to `v`: a later missing(v) reads
  # FALSE once the argument has been overwritten, which would turn the
  # fold_id conflict check below into a check that always fires.
  v_named <- !missing(v)
  scheme_named <- v_named || !missing(repeats) || !is.null(fold_id)

  unit_vars <- .fb_sel_variables(unit, "unit")
  strata_vars <- .fb_sel_variables(strata, "strata")

  if (length(unit_vars) > 0L && !scheme_named) {
    stop(.fb_refusal_condition(
      reason_code = "selection_resamples_underspecified",
      message = paste0(
        "fb_resamples(): `unit` names ",
        paste(unit_vars, collapse = ", "),
        " but no resampling scheme was given. A grouping unit says what may ",
        "not be split across a fold boundary, not how the splitting is ",
        "done. Supply `v`, `repeats` or an explicit `fold_id`."
      )
    ))
  }

  v <- .fb_sel_count(v, "v", minimum = 2L)
  repeats <- .fb_sel_count(repeats, "repeats", minimum = 1L)

  if (!is.null(fold_id)) {
    if (!is.atomic(fold_id) || length(fold_id) < 2L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_fold_id_invalid",
        message = paste0(
          "fb_resamples(): `fold_id` must be an atomic vector with one entry ",
          "per row of the data."
        )
      ))
    }
    if (v_named) {
      stop(.fb_refusal_condition(
        reason_code = "selection_fold_id_invalid",
        message = paste0(
          "fb_resamples(): `fold_id` and `v` were both given. An explicit ",
          "partition and a scheme that would build one are two answers to ",
          "the same question."
        )
      ))
    }
    if (repeats > 1L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_fold_id_invalid",
        message = paste0(
          "fb_resamples(): `fold_id` is a single explicit partition, so ",
          "`repeats = ",
          repeats,
          "` would repeat the same folds and report the repetition as ",
          "independent evidence."
        )
      ))
    }
  }

  structure(
    list(
      unit = unit_vars,
      v = if (is.null(fold_id)) v else NA_integer_,
      repeats = repeats,
      fold_id = fold_id,
      strata = strata_vars,
      scheme = if (is.null(fold_id)) "v_fold" else "explicit"
    ),
    class = "fb_resamples"
  )
}


#' @export
print.fb_resamples <- function(x, ...) {
  cat("<fb_resamples> scheme = \"", x$scheme, "\"\n", sep = "")
  cat(
    "  unit    : ",
    if (length(x$unit) > 0L) paste(x$unit, collapse = ", ") else "row",
    "\n",
    sep = ""
  )
  if (identical(x$scheme, "explicit")) {
    cat(
      "  folds   : ",
      length(unique(x$fold_id)),
      " explicit partition(s) over ",
      length(x$fold_id),
      " rows\n",
      sep = ""
    )
  } else {
    cat("  folds   : ", x$v, " x ", x$repeats, " repeat(s)\n", sep = "")
  }
  if (length(x$strata) > 0L) {
    cat("  strata  : ", paste(x$strata, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
