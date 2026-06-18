# .fb_aggregation_plan() -- model-level aggregation planning.
# Decides whether a model is a candidate for
# the aggregated-likelihood path WITHOUT materialising any design
# matrix. Reads only the IR + dataset wrapper metadata.
#
# Returns a `<fb_aggregation_plan>` S3 object carrying:
#   $eligible              logical: TRUE if the model is in-scope AND
#                          aggregation is likely productive.
#   $reason_codes          character vector of refusal reasons (empty
#                          on eligible = TRUE).
#   $cell_key_terms        list of cell-key contributions. Each entry:
#                          list(label = <human-readable>, vars =
#                          <data-column names>, K = <integer level
#                          count>; K = NA_integer_ when the cell-key
#                          count is data-dependent).
#   $requires_materialisation logical: TRUE when K cannot be sized from
#                          metadata alone (a numeric fixed term
#                          contributes ~N unique values).
#   $K_est                 integer: estimated cell count from
#                          metadata; NA_integer_ when
#                          requires_materialisation = TRUE.
#   $N                     integer: total row count.
#   $compression_est       numeric: K_est / N; NA_real_ when
#                          requires_materialisation = TRUE.
#
# The plan is informational at v0.3.0.9000: the aggregated emit
# dispatch (landed in v0.3.2) consumes the plan to
# decide between the per-row and aggregated paths. `.fb_preflight()`
# attaches the plan to the result for diagnostic visibility.
#
# Internal-only.

# Compression threshold: the model is considered productive for
# aggregation when K / N <= this value (i.e. at least 2:1 compression).
# 0.5 chosen so the dispatcher does not pay the aggregation overhead
# (cell-key materialisation + sum-of-squares accumulation) for models
# that would compress trivially. Future tuning of this threshold can
# happen via a `flexyBayes.aggregation_compression_threshold` option.
.FB_AGGREGATION_PRODUCTIVITY_THRESHOLD <- 0.5


.fb_aggregation_plan <- function(fb_ir, fb_dataset) {
  if (!inherits(fb_ir, "fb_terms")) {
    stop(
      ".fb_aggregation_plan() requires an `<fb_terms>` IR; got: ",
      paste(class(fb_ir), collapse = "/"),
      call. = FALSE
    )
  }
  if (!inherits(fb_dataset, "fb_dataset")) {
    stop(
      ".fb_aggregation_plan() requires an `<fb_dataset>` ",
      "wrapper; got: ",
      paste(class(fb_dataset), collapse = "/"),
      call. = FALSE
    )
  }

  N <- as.integer(fb_dataset$n_rows)
  reason_codes <- character(0L)

  # ---- Family / link scope check ---- #
  # Aggregation is exact for gaussian-identity, binomial-logit, and
  # poisson-log. Other families are refused; a non-canonical link is
  # refused because the aggregated emit fits the canonical link. Count
  # families carrying observation weights are refused -- the per-row
  # weights are not recoverable from the cell sums.
  canon <- c(gaussian = "identity", binomial = "logit", poisson = "log")
  fam <- fb_ir$family
  if (!fam %in% names(canon)) {
    reason_codes <- c(reason_codes, "non_aggregatable_family")
  } else if (!(is.null(fb_ir$link) || identical(fb_ir$link, canon[[fam]]))) {
    reason_codes <- c(
      reason_codes,
      if (identical(fam, "gaussian")) {
        "non_identity_link"
      } else {
        "non_canonical_link"
      }
    )
  }
  if (fam %in% c("binomial", "poisson") && length(fb_ir$addition_terms) > 0L) {
    reason_codes <- c(reason_codes, "count_weights_not_aggregatable")
  }

  # ---- Fixed-effect term scope ---- #
  any_data_dependent_cell <- FALSE
  cell_key_contribs <- list()

  for (t in fb_ir$fixed_terms) {
    ttype <- if (!is.null(t$type)) t$type else "expression"

    if (ttype %in% c("smooth", "s", "t2", "smooth_mgcv", "spline")) {
      # Smooth basis is observation-row-specific; breaks cell-constant
      # linear predictor.
      reason_codes <- c(reason_codes, "smooth_term_not_aggregatable")
    } else if (
      ttype %in%
        c(
          "numeric",
          "continuous",
          "I",
          "expression",
          "interaction",
          "polynomial"
        )
    ) {
      # Continuous fixed effect: each row carries an arbitrary value
      # so the cell key would be one cell per unique value
      # combination -- approximately N cells. Aggregation cannot
      # compress without binning the continuous variable. Flag the
      # plan as data-dependent.
      any_data_dependent_cell <- TRUE
    } else if (ttype %in% c("factor", "categorical")) {
      lbl <- .agg_plan_label(t, ttype, "fixed")
      L <- .agg_plan_factor_level_count(t, fb_dataset)
      cell_key_contribs[[length(cell_key_contribs) + 1L]] <-
        list(label = lbl, vars = as.character(t$var), K = L)
    } else if (identical(ttype, "factor_interaction")) {
      lbl <- .agg_plan_label(t, ttype, "fixed")
      Ls <- vapply(
        t$vars,
        function(v) {
          .fb_dataset_levels(fb_dataset, as.character(v))
        },
        numeric(1L)
      )
      L <- if (anyNA(Ls)) NA_integer_ else as.integer(prod(Ls))
      cell_key_contribs[[length(cell_key_contribs) + 1L]] <-
        list(label = lbl, vars = as.character(t$vars), K = L)
    } else if (identical(ttype, "factor_numeric_interaction")) {
      # Mixed factor:numeric -- the numeric component breaks the
      # cell-constant linear predictor (the indexed slope still
      # carries one slope per factor level applied to the per-row
      # numeric value).
      any_data_dependent_cell <- TRUE
    }
    # other types (structured-cov fallthrough; the random-side scope
    # check below will already have flagged the relevant random-term
    # types).
  }

  # ---- Random-effect term scope ---- #
  # The previous catch-all
  # `random_slope_in_scope` conflated random-slope and structured-
  # covariance shapes; split into a taxonomy that names what
  # actually breaks aggregation. Aligned with
  # .assert_aggregate_in_scope() so the same model produces the
  # same reason_code across both eligibility paths.
  for (t in fb_ir$random_terms) {
    rtype <- if (!is.null(t$type)) t$type else "simple"
    if (identical(rtype, "simple")) {
      lbl <- .agg_plan_label(t, rtype, "random")
      L <- .agg_plan_factor_level_count(t, fb_dataset)
      cell_key_contribs[[length(cell_key_contribs) + 1L]] <-
        list(label = lbl, vars = as.character(t$var), K = L)
    } else if (rtype %in% c("smooth_mgcv", "smooth", "s", "t2", "spline")) {
      reason_codes <- c(reason_codes, "smooth_term_not_aggregatable")
    } else if (rtype %in% c("simple_slope_uncor", "slope", "random_slope")) {
      reason_codes <- c(reason_codes, "random_slope_not_aggregatable")
    } else {
      # Structured-covariance terms (vm, ped, at, us, fa, ar1, ...):
      # the latent block is non-cell-constant; aggregation closure
      # does not hold.
      reason_codes <- c(reason_codes, "structured_random_not_aggregatable")
    }
  }

  # ---- Decide eligibility ---- #
  if (length(reason_codes) > 0L) {
    return(.new_fb_aggregation_plan(
      eligible = FALSE,
      reason_codes = unique(reason_codes),
      cell_key_terms = cell_key_contribs,
      requires_materialisation = FALSE,
      K_est = NA_integer_,
      N = N,
      compression_est = NA_real_
    ))
  }

  if (any_data_dependent_cell) {
    return(.new_fb_aggregation_plan(
      eligible = FALSE,
      reason_codes = "continuous_cell_key_data_dependent",
      cell_key_terms = cell_key_contribs,
      requires_materialisation = TRUE,
      K_est = NA_integer_,
      N = N,
      compression_est = NA_real_
    ))
  }

  # All cell-key contributions are factor-shaped. K_est = product of
  # level counts. If any level count is unresolvable from the dataset
  # wrapper's dictionaries (e.g. a metadata-only dataset that did not
  # carry dictionaries for all factors), flag the plan as
  # data-dependent.
  Ks <- vapply(cell_key_contribs, function(c) c$K, numeric(1L))
  if (anyNA(Ks)) {
    return(.new_fb_aggregation_plan(
      eligible = FALSE,
      reason_codes = "compression_level_count_unresolvable",
      cell_key_terms = cell_key_contribs,
      requires_materialisation = TRUE,
      K_est = NA_integer_,
      N = N,
      compression_est = NA_real_
    ))
  }

  K_est <- as.integer(prod(Ks))
  compression_est <- as.numeric(K_est) / as.numeric(N)

  productive <- compression_est <= .FB_AGGREGATION_PRODUCTIVITY_THRESHOLD

  .new_fb_aggregation_plan(
    eligible = productive,
    reason_codes = if (productive) {
      character(0L)
    } else {
      "compression_unproductive"
    },
    cell_key_terms = cell_key_contribs,
    requires_materialisation = FALSE,
    K_est = K_est,
    N = N,
    compression_est = compression_est
  )
}


# ---------------------------------------------------------------- #
# Helpers                                                           #
# ---------------------------------------------------------------- #

# Build a human-readable label for an IR term suitable for the plan's
# cell_key_terms list. Mirrors the .preflight_term_label() conventions.
.agg_plan_label <- function(term, ttype, kind) {
  if (!is.null(term$label) && nzchar(term$label)) {
    return(term$label)
  }
  if (identical(ttype, "factor_interaction") && !is.null(term$vars)) {
    return(paste(term$vars, collapse = ":"))
  }
  if (!is.null(term$var)) {
    return(as.character(term$var))
  }
  if (identical(kind, "random")) {
    return("<unnamed_group>")
  }
  "<unnamed>"
}

# Resolve a factor's level count from either the term's enriched
# slots (n_levels / var_n) or the dataset wrapper's dictionary.
# Returns NA_integer_ when none of the sources carry a level count.
.agg_plan_factor_level_count <- function(term, fb_dataset) {
  for (slot in c("var_n", "n_levels", "K")) {
    v <- term[[slot]]
    if (!is.null(v) && !is.na(v)) return(as.integer(v))
  }
  if (!is.null(term$var)) {
    L <- .fb_dataset_levels(fb_dataset, as.character(term$var))
    if (!is.na(L)) return(as.integer(L))
  }
  NA_integer_
}

# Constructor -- centralises the slot list shape + class attachment.
.new_fb_aggregation_plan <- function(
  eligible,
  reason_codes,
  cell_key_terms,
  requires_materialisation,
  K_est,
  N,
  compression_est
) {
  structure(
    list(
      eligible = isTRUE(eligible),
      reason_codes = reason_codes,
      cell_key_terms = cell_key_terms,
      requires_materialisation = isTRUE(requires_materialisation),
      K_est = K_est,
      N = as.integer(N),
      compression_est = compression_est
    ),
    class = c("fb_aggregation_plan", "list")
  )
}


# ---------------------------------------------------------------- #
# S3 print method                                                   #
# ---------------------------------------------------------------- #

#' Print method for an internal `<fb_aggregation_plan>` summary
#'
#' Diagnostic print of the model-level aggregation plan:
#' eligibility, reason codes, cell-key contributions, estimated cell
#' count, and the compression ratio.
#'
#' @param x   an `<fb_aggregation_plan>` object.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @keywords internal
#' @export
print.fb_aggregation_plan <- function(x, ...) {
  cat(sprintf("<fb_aggregation_plan>: eligible = %s\n", x$eligible))

  if (length(x$reason_codes)) {
    cat(
      "  reason_codes:             ",
      paste(x$reason_codes, collapse = ", "),
      "\n",
      sep = ""
    )
  }

  if (length(x$cell_key_terms)) {
    cat("  cell_key_terms:\n")
    for (c in x$cell_key_terms) {
      cat(sprintf(
        "    %-24s  L = %s\n",
        c$label,
        if (is.na(c$K)) {
          "NA"
        } else {
          format(c$K, big.mark = " ", scientific = FALSE)
        }
      ))
    }
  }

  cat(sprintf("  requires_materialisation: %s\n", x$requires_materialisation))
  cat(sprintf(
    "  N = %s; K_est = %s; compression_est = %s\n",
    format(x$N, big.mark = " ", scientific = FALSE),
    if (is.na(x$K_est)) {
      "NA"
    } else {
      format(x$K_est, big.mark = " ", scientific = FALSE)
    },
    if (is.na(x$compression_est)) {
      "NA"
    } else {
      sprintf("%.3f", x$compression_est)
    }
  ))
  invisible(x)
}
