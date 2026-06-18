# .fb_aggregate_gaussian() -- Gaussian exact aggregation by
# sufficient statistics (target v0.3.2).
#
# Cell construction. A cell k is the maximal subset of observations
# that share the same linear-predictor value mu_k = x_k^T beta + sum_j
# u_{g_j(k)}, where the right-hand side is constant within the cell by
# construction. The cell key is therefore the unique combination of
#   (i)  every column of the fixed-effects model.matrix();
#   (ii) every random-intercept grouping-factor level.
#
# Three sufficient statistics per cell:
#   n_k    -- cell size              (integer)
#   S1_k   -- sum_{i in k} y_i       (double)
#   S2_k   -- sum_{i in k} y_i^2     (double)
#
# Aggregated gaussian-identity log-likelihood:
#   ell(beta, sigma^2) = -N/2 * log(2 pi sigma^2)
#                        - (1 / (2 sigma^2)) * sum_k [ S2_k
#                                                       - 2 mu_k S1_k
#                                                       + n_k mu_k^2 ]
# which is algebraically identical (not approximate) to the per-row
# gaussian log-likelihood. Bit-exact equivalence is the property-based
# test gate; the cell-mean shortcut (replacing the cell by its mean
# and modelling N(mu_k, sigma^2/n_k)) is explicitly forbidden because
# it drops the sigma-dependent within-cell sum-of-squares term and
# biases the residual-variance posterior.
#
# Internal-only at v0.3.2 -- consumed by emit_gaussian_aggregated()
# (backend wiring; not in this commit). No @export.

# Top-level aggregator. Returns an <fb_aggregated> S3 list:
#   $cell_design     K x p cell-level fixed-effects design matrix.
#                    One row per cell carrying the unique fixed-effect
#                    column values for that cell. The v0.3.0.9000
#                    foundation stored the full N x p design matrix;
#                    this defeated compression at large N and did not
#                    scale to the 10M-100M-row roadmap. The cell-level
#                    form contains
#                    the same information for the aggregated emit path
#                    (the per-cell linear predictor) at a fraction of
#                    the memory cost.
#   $sufficient_stats data.table with cols: cell_id (integer), n_k,
#                     S1_k, S2_k, plus the cell-key columns
#   $cell_key_cols   character vector: the columns that define a cell
#   $fixed_cols      character vector: just the fixed-effects columns
#   $random_cols     character vector: just the random-intercept group
#                     factor columns
#   $K              integer count of distinct cells
#   $N              integer total observations (sum of n_k)
#   $compression    K / N ratio (smaller is better; 1.0 = no compression)
#   $response       response variable name (string)
#
# Per-row equivalence verification (the bit-exact log-lik gate)
# rebuilds the per-row design matrix from the underlying data via
# stats::model.matrix(); the aggregator no longer carries it.
#
# Refuses (typed condition `flexybayes_aggregate_out_of_scope`) when:
#   - family != gaussian (link identity)
#   - any smooth fixed-effect term
#   - any random slope (correlated or uncorrelated)
.fb_aggregate_gaussian <- function(fb_ir, fb_dataset) {
  if (!inherits(fb_ir, "fb_terms")) {
    stop(
      ".fb_aggregate_gaussian() requires an `<fb_terms>` IR; got: ",
      paste(class(fb_ir), collapse = "/"),
      call. = FALSE
    )
  }
  if (!inherits(fb_dataset, "fb_dataset")) {
    stop(
      ".fb_aggregate_gaussian() requires an `<fb_dataset>` ",
      "wrapper; got: ",
      paste(class(fb_dataset), collapse = "/"),
      call. = FALSE
    )
  }
  if (.fb_dataset_is_metadata(fb_dataset)) {
    stop(
      ".fb_aggregate_gaussian() requires a data-backed dataset; ",
      "metadata-only descriptors cannot be aggregated (no `y` to ",
      "sum).",
      call. = FALSE
    )
  }

  # Scope gate -- the cell-constant mu_k property only holds for the
  # aggregation envelope: gaussian-identity + fixed (numeric /
  # factor) + random intercept.
  .assert_aggregate_in_scope(fb_ir)

  data <- fb_dataset$data # data.table; always a wrapper-owned copy
  # (see .fb_dataset() header).
  response <- fb_ir$response
  if (!response %in% names(data)) {
    stop(
      "Response variable '",
      response,
      "' not found in dataset.",
      call. = FALSE
    )
  }

  y <- as.numeric(data[[response]])
  N <- length(y)

  # ----------------------- Fixed-effects design --------------------- #
  # Reconstruct a formula-like RHS for model.matrix. We build it from
  # the IR's fixed_terms labels so that the byte-identity contract
  # with the per-row emit is preserved: emit_greta() and emit_inla()
  # both feed model.matrix(<same formula>, data) into the backend.
  fixed_form <- .fb_aggregate_fixed_formula(fb_ir)
  if (is.null(fixed_form)) {
    # Intercept-only -- model.matrix(~ 1, data) is N x 1 of ones.
    X <- matrix(1, nrow = N, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  } else {
    X <- stats::model.matrix(fixed_form, data = data)
  }

  fixed_cols <- colnames(X)

  # ----------------------- Random-intercept keys -------------------- #
  random_intercept_groups <- vapply(
    fb_ir$random_terms,
    function(t) {
      if (identical(t$type, "simple")) {
        as.character(t$var)
      } else {
        NA_character_
      }
    },
    character(1L)
  )
  random_cols <- random_intercept_groups[
    !is.na(random_intercept_groups)
  ]

  # ----------------------- Build the cell-key data.table ------------ #
  # Cell key = every fixed-effect column + every random-intercept
  # grouping-factor column. We materialise the columns as integer
  # codes (factors -> .levels lookup) so data.table's by= can hash
  # them efficiently. The dictionary-frozen <fb_dataset> guarantees
  # the level set is stable across re-aggregation.
  key_dt <- data.table::as.data.table(X)
  data.table::setnames(key_dt, fixed_cols)
  for (g in random_cols) {
    data.table::set(key_dt, j = g, value = data[[g]])
  }
  data.table::set(key_dt, j = ".y", value = y)

  key_cols <- c(fixed_cols, random_cols)

  # ----------------------- Aggregate ------------------------------- #
  # data.table .N / sum / sum-of-squares by the full cell key. Each
  # row of the aggregated table is one cell. Using the env-arg form
  # (compatible with `cedta()` check in non-data.table-imported
  # namespaces by routing through evalq() on a data.table-aware env).
  agg_expr <- substitute(
    key_dt[, list(n_k = .N, S1_k = sum(.y), S2_k = sum(.y * .y)), by = key_cols]
  )
  agg <- eval(agg_expr, envir = asNamespace("data.table"))

  K <- nrow(agg)
  data.table::set(agg, j = "cell_id", value = seq_len(K))
  data.table::setcolorder(
    agg,
    c("cell_id", key_cols, "n_k", "S1_k", "S2_k")
  )

  # Extract the cell-level design matrix (K x p) from the aggregated
  # data.table. The full N x p `X` matrix is no longer carried on the
  # returned object -- it goes out of scope on function return.
  cell_design <- as.matrix(agg[, fixed_cols, with = FALSE])
  colnames(cell_design) <- fixed_cols

  structure(
    list(
      cell_design = cell_design,
      sufficient_stats = agg,
      cell_key_cols = key_cols,
      fixed_cols = fixed_cols,
      random_cols = random_cols,
      K = K,
      N = as.integer(N),
      compression = as.numeric(K) / as.numeric(N),
      response = response
    ),
    class = c("fb_aggregated", "list")
  )
}


# ---------------------------------------------------------------- #
# Scope gate                                                        #
# ---------------------------------------------------------------- #

# Refuses IRs outside the aggregation scope envelope. Raises a typed
# `flexybayes_aggregate_out_of_scope` condition so downstream tooling
# (the eventual emit_gaussian_aggregated() dispatch branch) can
# pattern-match without parsing free text.
.assert_aggregate_in_scope <- function(fb_ir) {
  if (!identical(fb_ir$family, "gaussian")) {
    .stop_aggregate_out_of_scope(
      reason_code = "non_gaussian_family",
      detail = paste0(
        "family = '",
        fb_ir$family,
        "'; aggregated emit requires the gaussian family."
      )
    )
  }
  if (!(is.null(fb_ir$link) || identical(fb_ir$link, "identity"))) {
    .stop_aggregate_out_of_scope(
      reason_code = "non_identity_link",
      detail = paste0(
        "link = '",
        fb_ir$link,
        "'; aggregated emit requires the identity link."
      )
    )
  }

  # Fixed-effect term scope: numeric + factor + interactions thereof.
  # Smooths break cell-constant mu (basis is observation-row-specific).
  for (t in fb_ir$fixed_terms) {
    ttype <- if (!is.null(t$type)) t$type else "expression"
    if (ttype %in% c("smooth", "s", "t2", "smooth_mgcv", "spline")) {
      .stop_aggregate_out_of_scope(
        reason_code = "smooth_term_not_aggregatable",
        detail = paste0(
          "smooth fixed-effect term '",
          if (!is.null(t$var)) t$var else "?",
          "' breaks the cell-constant mu property."
        )
      )
    }
  }

  # Random-effect term scope: random intercept only. The previous
  # catch-all `random_slope_in_scope` conflated three distinct
  # out-of-scope shapes; split into a
  # taxonomy that names what actually breaks aggregation:
  #   - smooth_term_not_aggregatable    (smooth / spline on random side)
  #   - random_slope_not_aggregatable   (correlated / uncorrelated slope)
  #   - structured_random_not_aggregatable (vm, ped, at, us, fa, ar1, ...)
  for (t in fb_ir$random_terms) {
    rtype <- if (!is.null(t$type)) t$type else "simple"
    if (identical(rtype, "simple")) {
      next
    }

    reason <- if (rtype %in% c("smooth", "s", "t2", "smooth_mgcv", "spline")) {
      "smooth_term_not_aggregatable"
    } else if (rtype %in% c("simple_slope_uncor", "slope", "random_slope")) {
      "random_slope_not_aggregatable"
    } else {
      "structured_random_not_aggregatable"
    }

    .stop_aggregate_out_of_scope(
      reason_code = reason,
      detail = paste0(
        "random-effect term type = '",
        rtype,
        "' breaks the cell-constant mu property."
      )
    )
  }

  invisible(TRUE)
}

.stop_aggregate_out_of_scope <- function(reason_code, detail) {
  msg <- paste0(
    "Aggregated Gaussian emit refused -- the model is out of scope.\n",
    "Reason: ",
    reason_code,
    "\n",
    detail,
    "\n\n",
    "Out-of-scope IRs route to the standard per-observation emit ",
    "(no compression); this condition is informational, not an error ",
    "in the per-row path."
  )
  cond <- structure(
    class = c("flexybayes_aggregate_out_of_scope", "error", "condition"),
    list(
      message = msg,
      call = NULL,
      reason_code = reason_code,
      detail = detail
    )
  )
  stop(cond)
}


# ---------------------------------------------------------------- #
# Likelihood evaluators (pure-R; for the bit-exact equivalence test) #
# ---------------------------------------------------------------- #

# Raw per-row gaussian log-likelihood. mu vector of length N; sigma
# positive scalar.
.gaussian_loglik_raw <- function(y, mu, sigma) {
  N <- length(y)
  -0.5 * N * log(2 * pi * sigma^2) - 0.5 / sigma^2 * sum((y - mu)^2)
}

# Aggregated per-cell gaussian log-likelihood using the three
# sufficient statistics (n_k, S1_k, S2_k). mu_cell is the
# linear-predictor value at each cell (length K); sigma positive
# scalar. Algebraically identical to the raw form.
.gaussian_loglik_aggregated <- function(n_k, S1_k, S2_k, mu_cell, sigma) {
  N <- sum(n_k)
  -0.5 *
    N *
    log(2 * pi * sigma^2) -
    0.5 / sigma^2 * sum(S2_k - 2 * mu_cell * S1_k + n_k * mu_cell^2)
}

# Cell-mean shortcut. Drops the within-cell
# sum-of-squares term `S2_k - S1_k^2 / n_k` (sigma-dependent). This is
# the BIASED form -- exposed as a helper so the property-based test
# can confirm it produces a different log-likelihood from the correct
# aggregated form.
.gaussian_loglik_cellmean <- function(n_k, S1_k, mu_cell, sigma) {
  ybar_k <- S1_k / n_k
  -0.5 *
    sum(n_k * log(2 * pi * sigma^2 / n_k)) -
    0.5 / sigma^2 * sum(n_k * (ybar_k - mu_cell)^2)
}


# ---------------------------------------------------------------- #
# Fixed-effects formula reconstruction                              #
# ---------------------------------------------------------------- #

# Build a one-sided formula from the IR's fixed_terms list. Returns
# NULL if the only fixed term is an implicit intercept (caller
# handles the matrix-of-ones case).
.fb_aggregate_fixed_formula <- function(fb_ir) {
  labels <- vapply(
    fb_ir$fixed_terms,
    function(t) {
      if (!is.null(t$label) && nzchar(t$label)) {
        return(t$label)
      }
      if (!is.null(t$var) && nzchar(t$var)) {
        return(as.character(t$var))
      }
      if (!is.null(t$deparse) && nzchar(t$deparse)) {
        return(t$deparse)
      }
      NA_character_
    },
    character(1L)
  )
  labels <- labels[!is.na(labels)]
  if (!length(labels)) {
    return(NULL)
  }
  rhs <- paste(labels, collapse = " + ")
  if (!isTRUE(fb_ir$intercept)) {
    rhs <- paste0(rhs, " - 1")
  }
  stats::as.formula(paste0("~ ", rhs))
}


# ---------------------------------------------------------------- #
# S3 print method                                                   #
# ---------------------------------------------------------------- #

#' Print method for an internal `<fb_aggregated>` summary
#'
#' Diagnostic print of the sufficient-statistics aggregation:
#' compression ratio, cell count, total N, and the per-cell key
#' columns.
#'
#' @param x   an `<fb_aggregated>` object.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @keywords internal
#' @export
print.fb_aggregated <- function(x, ...) {
  ratio <- if (x$compression > 0) {
    sprintf("%.1f:1", 1 / x$compression)
  } else {
    "n/a"
  }
  cat(sprintf(
    "<fb_aggregated>: N = %s rows -> K = %s cells (compression %s)\n",
    format(x$N, big.mark = " ", scientific = FALSE),
    format(x$K, big.mark = " ", scientific = FALSE),
    ratio
  ))
  if (length(x$fixed_cols)) {
    cat("  fixed cols:  ", paste(x$fixed_cols, collapse = ", "), "\n", sep = "")
  }
  if (length(x$random_cols)) {
    cat(
      "  random grps: ",
      paste(x$random_cols, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  cat("  response:    ", x$response, "\n", sep = "")
  invisible(x)
}
