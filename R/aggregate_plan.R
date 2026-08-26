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

# .fb_aggregation_verdict_line() --- one line stating N, K (when known),
# rows per cell, and a plain verdict on whether aggregation pays, with
# the threshold stated in the same units the line reports (S6, scale
# strategy 2026-08-22: "a gen:loc:yearf model ... measured 1.06 rows per
# cell against the shipped demo's 1,200 cells" -- the printed line is
# what would have shown a user that contrast directly, on either
# model, without reading the strategy document).
#
# Rows per cell (N/K) is the quantity a breeder reads directly off a
# trial design; the productivity test is expressed internally as
# K/N <= .FB_AGGREGATION_PRODUCTIVITY_THRESHOLD (>= 2:1 compression), so
# the equivalent rows-per-cell threshold is its reciprocal. Shared by
# print.fb_aggregation_plan() and print.fb_plan()'s Aggregation: line so
# the two plan surfaces cannot drift on wording -- this is the S6 half
# of "the plan surface tells the truth" (the FS-22 half is the grammar
# fix in R/fb_plan.R).
#
# `rows_per_cell` is generally N/K precomputed by the caller (both
# callers already carry it under different field names); passing it
# explicitly rather than recomputing from N and K keeps this a pure
# formatting function and avoids a second definition of "rows per
# cell" that could drift from the one the caller actually used.
.fb_aggregation_verdict_line <- function(N, K, rows_per_cell) {
  thresh <- 1 / .FB_AGGREGATION_PRODUCTIVITY_THRESHOLD
  n_txt <- if (is.na(N)) {
    "NA"
  } else {
    format(N, big.mark = " ", scientific = FALSE)
  }
  k_txt <- if (is.na(K)) {
    "NA"
  } else {
    format(K, big.mark = " ", scientific = FALSE)
  }
  rpc_txt <- if (is.na(rows_per_cell)) {
    "NA"
  } else {
    sprintf("%.2f", rows_per_cell)
  }
  verdict <- if (is.na(rows_per_cell)) {
    "aggregation applicability unknown (cell count not resolved)"
  } else if (rows_per_cell >= thresh) {
    sprintf("aggregation will pay (threshold: rows/cell >= %.1f)", thresh)
  } else {
    sprintf(
      "aggregation will NOT pay (threshold: rows/cell >= %.1f)",
      thresh
    )
  }
  sprintf("N = %s; K = %s; rows/cell = %s; %s", n_txt, k_txt, rpc_txt, verdict)
}

# .fb_checked_cell_count() --- a cell-count estimate carried safely
# past R's integer limit, or a refusal.
#
# K_est = prod(Ks) is the product of the aggregation cell key's level
# counts. Growing a design by genotype the way a breeding programme
# does (FS-24) reaches this product well before the row count reaches
# its own limit: at 911,808 rows / 54,208 genotypes the cell-key
# product is already about 2.4e9. `as.integer(prod(Ks))` past
# `2^31 - 1` returned `NA` with only a base-R coercion warning, and
# `compression_est` then divided by that `NA` silently -- the plan
# lost both its cell count and its compression estimate rather than
# refusing. `Ks` already arrives as a double vector (`vapply(...,
# numeric(1L))`), so `prod(Ks)` is computed in double space; only the
# final cast back to `integer` for the common, well within-range case
# could overflow, and this is where the guard sits.
#
# This is the cell-count member of the same class the row-count guard
# (`.fb_checked_row_count()`, R/stream_aggregate.R) repairs for `N`.
.fb_checked_cell_count <- function(Ks, cell_key_contribs) {
  K_est_dbl <- prod(Ks)
  if (is.finite(K_est_dbl) && K_est_dbl <= .Machine$integer.max) {
    return(as.integer(round(K_est_dbl)))
  }
  labels <- vapply(
    cell_key_contribs,
    function(c) if (!is.null(c$label)) c$label else "<unlabelled>",
    character(1L)
  )
  stop(.fb_refusal_condition(
    reason_code = "cell_count_exceeds_integer",
    message = paste0(
      "This model's aggregation cell key (",
      paste(labels, collapse = " x "),
      ") has an estimated cell count of ",
      format(K_est_dbl, scientific = FALSE, big.mark = ","),
      ", past R's integer limit of ",
      format(.Machine$integer.max, scientific = FALSE),
      ". A count that large silently becomes NA when recorded, so the ",
      "aggregation plan would carry no cell count and no compression ",
      "estimate at all. Reduce the cell-key cardinality (drop or coarsen ",
      "an interacting factor) or fit on the per-row route."
    ),
    K_est = K_est_dbl
  ))
}


.fb_aggregation_plan <- function(fb_ir, fb_dataset) {
  .check_fb_terms(
    fb_ir,
    ".fb_aggregation_plan() requires an `<fb_terms>` IR; got: ",
    paste(class(fb_ir), collapse = "/")
  )
  .check_fb_dataset(
    fb_dataset,
    ".fb_aggregation_plan() requires an `<fb_dataset>` ",
    "wrapper; got: ",
    paste(class(fb_dataset), collapse = "/")
  )

  N <- .fb_checked_row_count(fb_dataset$n_rows, "This dataset")
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
  #
  # The cell key is the set of VARIABLES the cell-constant predictors
  # need, not the list of terms that need them. `y ~ a * b` produces
  # three terms -- `a`, `b` and `a:b` -- over two variables, and the
  # runtime aggregator (.fb_stream_key_cols()) has always keyed on the
  # union. The planner multiplied the term-level counts instead, so the
  # same replicated 4-by-4 factorial that compresses 320 rows to 16
  # cells was estimated at 4 * 4 * 16 = 256 and refused as
  # `compression_unproductive` -- the plan disagreed with the emitter
  # about the very quantity the decision turns on.
  #
  # `cell_key_contribs` therefore accumulates one entry per distinct
  # variable, and K_est is the product over that set.
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
      cell_key_contribs <- .agg_plan_add_var(
        cell_key_contribs,
        as.character(t$var),
        .agg_plan_factor_level_count(t, fb_dataset)
      )
    } else if (identical(ttype, "factor_interaction")) {
      for (v in as.character(t$vars)) {
        cell_key_contribs <- .agg_plan_add_var(
          cell_key_contribs,
          v,
          .fb_dataset_levels(fb_dataset, v)
        )
      }
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
      grp <- if (is.null(t$var)) {
        .agg_plan_label(t, rtype, "random")
      } else {
        as.character(t$var)
      }
      cell_key_contribs <- .agg_plan_add_var(
        cell_key_contribs,
        grp,
        .agg_plan_factor_level_count(t, fb_dataset)
      )
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

  # ---- Residual-term scope ---- #
  # Only the plain `units` residual is cell-constant, so only `units` is
  # aggregatable. An AR1 / separable AR1xAR1 residual is reinterpreted on
  # INLA as a per-observation latent field, and a sectioned
  # `dsum(~ units | f)` residual carries one variance per level of `f`;
  # neither is constant within a predictor cell, and no aggregated emitter
  # represents either. The check names what IS representable rather than
  # enumerating what is not, so a residual class added later is refused
  # here until an aggregated emitter for it exists -- the earlier
  # enumeration silently let `at_units` through, and the plan surface then
  # advertised an aggregated INLA route for a model the INLA gate refuses.
  for (t in fb_ir$residual_terms) {
    rtype <- if (!is.null(t$type)) t$type else "units"
    if (!identical(rtype, "units")) {
      reason_codes <- c(reason_codes, "structured_residual_not_aggregatable")
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

  # All cell-key contributions are factor-shaped. K_est = product of the
  # level counts over the DISTINCT key variables. If any level count is
  # unresolvable from the dataset wrapper's dictionaries (e.g. a
  # metadata-only dataset that did not carry dictionaries for all
  # factors), flag the plan as data-dependent.
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

  K_est <- .fb_checked_cell_count(Ks, cell_key_contribs)
  compression_est <- as.numeric(K_est) / as.numeric(N)

  # ---- Engine limit: interaction design columns ---- #
  # Both aggregated emitters pre-expand the fixed effects to a
  # model matrix and then name its columns in the INLA formula. An
  # interaction column is named `a2:b2`, and INLA does not resolve a
  # column whose name contains a colon even when the formula backticks
  # it -- a live probe on a replicated 4-by-4 binomial factorial fails
  # inside INLA with `object 'aa2:bb2' not found` on both the binomial
  # and the Poisson emit. Renaming the columns would put a synthesised
  # token into `summary.fixed`, `marginals.fixed` and the latent field
  # that `.inla_fixef_draws()` matches on, so the aggregated route
  # refuses and the per-row route -- which uses INLA's native `a:b`
  # notation on the raw columns and fits -- takes the model.
  #
  # K_est is reported anyway: it is the compression the model would
  # have had, and it is what makes the refusal legible.
  if (.agg_plan_has_factor_interaction(fb_ir)) {
    return(.new_fb_aggregation_plan(
      eligible = FALSE,
      reason_codes = "aggregated_interaction_column_not_representable",
      cell_key_terms = cell_key_contribs,
      requires_materialisation = FALSE,
      K_est = K_est,
      N = N,
      compression_est = compression_est
    ))
  }

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

# One sentence of context for the reason codes whose name alone does not
# tell a user what to do about them. Returns "" for the rest, so the
# refusal message reads the same as before where nothing is to be added.
.agg_refusal_note <- function(reason_codes) {
  if (
    "aggregated_interaction_column_not_representable" %in% reason_codes
  ) {
    return(paste0(
      " The aggregated emit expands the fixed effects to a model matrix",
      " and names its columns in the INLA formula, and INLA cannot",
      " reference a column whose name contains a colon, which is how an",
      " interaction column is named. The per-row route uses INLA's own",
      " `a:b` notation and fits this model."
    ))
  }
  ""
}

# TRUE when the fixed effects contain a factor-by-factor interaction, so
# the aggregated design matrix would carry a column whose name holds a
# colon. See the call site for why that is fatal on the INLA emit.
.agg_plan_has_factor_interaction <- function(fb_ir) {
  any(vapply(
    fb_ir$fixed_terms,
    function(t) identical(t$type %||% "expression", "factor_interaction"),
    logical(1L)
  ))
}

# Record one variable's contribution to the joint cell key, keeping the
# list one-entry-per-variable.
#
# A variable named by more than one term (`a` as a main effect and again
# inside `a:b`, or a group factor that is also a fixed effect) enters the
# key once. Where the two sightings disagree on the level count -- a
# term slot carrying `var_n` against a dataset dictionary, say -- the
# first is kept and the second is a no-op, except that a resolvable count
# replaces an unresolvable one so an NA from a partial dictionary does
# not poison a count another term already knew.
.agg_plan_add_var <- function(contribs, var, count) {
  var <- as.character(var)
  count <- if (is.null(count) || is.na(count)) {
    NA_integer_
  } else {
    as.integer(count)
  }
  for (i in seq_along(contribs)) {
    if (identical(contribs[[i]]$label, var)) {
      if (is.na(contribs[[i]]$K) && !is.na(count)) {
        contribs[[i]]$K <- count
      }
      return(contribs)
    }
  }
  contribs[[length(contribs) + 1L]] <- list(
    label = var,
    vars = var,
    K = count
  )
  contribs
}

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
#' @param x   An `<fb_aggregation_plan>` object, the model-level plan the
#'   aggregation layer builds.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the eligibility and
#'   cell-count lines it prints.
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
  rows_per_cell <- if (is.na(x$K_est) || identical(x$K_est, 0L)) {
    NA_real_
  } else {
    x$N / x$K_est
  }
  cat(
    "  ",
    .fb_aggregation_verdict_line(
      N = x$N, K = x$K_est, rows_per_cell = rows_per_cell
    ),
    "\n",
    sep = ""
  )
  invisible(x)
}
