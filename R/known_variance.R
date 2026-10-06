# known_variance.R -- a known variance per row, added to the residual.
#
# `flexybayes(known_variance = v)` fits
#
#   y_i = x_i' beta + z_i' u + e_i,   e_i ~ N(0, sigma^2 + v_i),
#
# with v_i known and sigma^2 estimated. It is the second stage of a two-stage
# analysis, where y_i is a first-stage estimate (a variety mean in one trial)
# and v_i is the variance the first stage attached to it. `weights = w` is a
# different model, Var(e_i) = sigma^2 / w_i, which scales the estimated
# variance instead of adding to it, so the two do not share an argument.
#
# The lowering, one per engine.
#
#   brms   y | se(s, sigma = TRUE) with s_i = sqrt(v_i). brms writes the
#          likelihood as normal_lpdf(Y | mu, sqrt(square(sigma) + se2)), and
#          with a sigma sub-formula (a heterogeneous residual, the weights
#          offset) sigma is that row's value, so the two add as they do on
#          INLA.
#   INLA   a row effect f(i, w, model = "iid") whose precision is fixed at 1,
#          with w_i = sqrt(v_i) / s on the standardised response (s is its
#          scale). Row i's variance on the original scale gains
#          s^2 w_i^2 = v_i. The Gaussian likelihood's own precision offset
#          cannot carry this: it is one value for every row, and INLA fixes
#          it with no prior.
#
# The INLA row effect sits in the linear predictor, so INLA's fitted values
# include its posterior mean, where brms's exclude the known variance, which
# is part of the likelihood there. fitted() on an INLA fit subtracts it
# (.fb_inla_known_variance_effect()), and the effect is left out of every
# random-effect table.
#
# Where the values live. The vector is checked against the rows it was given
# with, then carried as a column of the analysis data, so it follows every
# row operation the missing-data layer performs: a row dropped under `omit`
# drops its variance, and a cell the design grid adds, which has no response
# and so no likelihood, gets zero. The IR records the column and the vector
# as supplied, which update() re-issues.
#
# A vector of zeros is the model without a known variance under another
# spelling, so nothing is attached and the engines receive the call they
# would have received without it.
#
# A known covariance matrix. `known_variance = K`, K symmetric and positive
# definite with one row per data row, fits
#
#   y = X beta + Z u + e + k,   e ~ N(0, sigma^2 I),   k ~ N(0, K),
#
# the second stage of a two-stage analysis whose first stage returned the
# covariance of its estimates, not only their variances. A row with a zero
# diagonal entry, and so a zero row and column, carries no known part. INLA
# fits k as a generic0 effect over the rows with a positive variance, with
# Cmatrix the inverse of that block scaled to the standardised response and
# its precision fixed at 1. sigma^2 stays estimated, as it is for a vector:
# with one row per trial it is the trial term. A diagonal K is the vector
# model under a second lowering. brms refuses a matrix: its fcor() term fits
# sigma^2 K and does not combine with se() (brms 2.23.0), so sigma^2 I + K
# has no spelling there.
#
# The matrix is stored once on the IR, and each data row carries its row
# index into it in a second column, so the matrix follows the same row
# operations the vector does. A cell the design grid adds has no index.

.FB_KNOWN_VARIANCE_COL <- "flexybayes_known_variance"
.FB_BRMS_KNOWN_SE_COL <- "flexybayes_known_se"
.FB_INLA_KNOWN_INDEX <- "flexybayes_known_index"
.FB_INLA_KNOWN_SCALE <- "flexybayes_known_scale"
.FB_KNOWN_ROW_COL <- "flexybayes_known_row"
.FB_INLA_KNOWN_PRECISION <- "flexybayes_known_precision"

.fb_known_variance_columns <- function() {
  c(
    .FB_KNOWN_VARIANCE_COL,
    .FB_BRMS_KNOWN_SE_COL,
    .FB_INLA_KNOWN_INDEX,
    .FB_INLA_KNOWN_SCALE,
    .FB_KNOWN_ROW_COL,
    .FB_INLA_KNOWN_PRECISION
  )
}

# .fb_known_variance_check() --- refuse a vector that is not one finite,
# non-negative variance per row of `data`, and a matrix that is not a
# covariance over those rows.
.fb_known_variance_check <- function(known_variance, data) {
  if (is.matrix(known_variance) && ncol(known_variance) > 1L) {
    return(.fb_known_covariance_check(known_variance, data))
  }
  n <- nrow(data)
  values <- function(k) paste0(k, if (k == 1L) " value is" else " values are")
  problem <- if (!is.numeric(known_variance) || is.matrix(known_variance)) {
    paste0("a ", class(known_variance)[[1L]], " was given")
  } else if (length(known_variance) != n) {
    paste0(
      "it has ", length(known_variance), " values and `data` has ", n,
      " rows"
    )
  } else if (any(!is.finite(known_variance))) {
    paste0(values(sum(!is.finite(known_variance))), " missing or not finite")
  } else if (any(known_variance < 0)) {
    paste0(values(sum(known_variance < 0)), " negative")
  }
  if (is.null(problem)) {
    return(invisible(NULL))
  }
  .fb_known_variance_invalid(problem)
}

# .fb_known_variance_invalid(): the refusal of a `known_variance` that is
# neither form, with what was wrong with it.
.fb_known_variance_invalid <- function(problem) {
  stop(.fb_refusal_condition(
    reason_code = "known_variance_invalid",
    message = paste0(
      "`known_variance` must be a numeric vector with one finite, ",
      "non-negative variance for each row of `data`, or a symmetric ",
      "positive-definite covariance matrix with one row and column for ",
      "each row of `data`, and ", problem, ". ",
      "A row whose first-stage estimate carries no error takes 0."
    )
  ))
}

# Relative to the largest entry; a covariance read from a first-stage fit is
# symmetric to rounding.
.FB_KNOWN_SYMMETRY_TOL <- 1e-8

# .fb_known_covariance_check(): refuse a matrix that is not a covariance
# over the rows of `data`. A row with a zero variance must be zero
# throughout, and the other rows must form a positive-definite block.
.fb_known_covariance_check <- function(known_variance, data) {
  n <- nrow(data)
  k <- known_variance
  problem <- if (!is.numeric(k)) {
    paste0("a ", typeof(k), " matrix was given")
  } else if (nrow(k) != ncol(k)) {
    paste0("it is ", nrow(k), " x ", ncol(k))
  } else if (nrow(k) != n) {
    paste0("it has ", nrow(k), " rows and `data` has ", n)
  } else if (any(!is.finite(k))) {
    paste0(sum(!is.finite(k)), " entries are missing or not finite")
  } else if (max(abs(k - t(k))) > .FB_KNOWN_SYMMETRY_TOL * max(abs(k))) {
    paste0(
      "it is not symmetric: an entry and its transpose differ by up to ",
      signif(max(abs(k - t(k))), 3L)
    )
  } else if (any(diag(k) < 0)) {
    paste0(sum(diag(k) < 0), " diagonal entries are negative")
  }
  if (is.null(problem)) {
    zero <- which(diag(k) == 0)
    positive <- which(diag(k) > 0)
    linked <- zero[rowSums(k[zero, , drop = FALSE] != 0) > 0L]
    if (length(linked) > 0L) {
      problem <- paste0(
        "row ", linked[[1L]], " has a zero variance and a nonzero covariance"
      )
    } else if (
      length(positive) > 0L &&
        is.null(tryCatch(
          chol(k[positive, positive, drop = FALSE]),
          error = function(e) NULL
        ))
    ) {
      problem <- "it is not positive definite"
    }
  }
  if (is.null(problem)) {
    return(invisible(NULL))
  }
  .fb_known_variance_invalid(problem)
}

# .fb_known_variance_refuse() --- the refusal for a request that cannot carry
# a known variance, with the reason and the remedy in `detail`.
.fb_known_variance_refuse <- function(detail) {
  stop(.fb_refusal_condition(
    reason_code = "known_variance_not_supported",
    message = detail
  ))
}

# .fb_known_variance_attach() --- the vector as a column of the analysis data,
# and its record on the IR. A vector of zeros attaches nothing.
.fb_known_variance_attach <- function(fb, data, known_variance, aggregate) {
  covariance <- NULL
  if (is.matrix(known_variance) && ncol(known_variance) > 1L) {
    covariance <- (known_variance + t(known_variance)) / 2
    dimnames(covariance) <- NULL
    known_variance <- diag(covariance)
  }
  if (!any(known_variance > 0)) {
    return(list(fb = fb, data = data))
  }
  fam <- fb$family %||% ""
  link_ok <- is.null(fb$link) || identical(fb$link, "identity")
  if (!identical(fam, "gaussian") || !link_ok) {
    .fb_known_variance_refuse(paste0(
      "`known_variance` adds a known variance to the residual variance of a ",
      "Gaussian response on the identity link, Var(y_i) = sigma^2 + v_i. ",
      "This model is family = \"", fam, "\"",
      if (!link_ok) paste0(", link = \"", fb$link, "\"") else "",
      ", which has no additive residual variance for it to join."
    ))
  }
  if (isTRUE(aggregate)) {
    .fb_known_variance_refuse(paste0(
      "`aggregate = TRUE` refused: the aggregated route compresses rows into ",
      "cell statistics, and a variance that differs from row to row does ",
      "not survive the compression. Leave `aggregate` at \"auto\", which ",
      "fits a Gaussian model row by row."
    ))
  }
  taken <- intersect(.fb_known_variance_columns(), names(data))
  if (length(taken) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "nested_index_name_collision",
      message = paste0(
        "`known_variance` is carried in columns named ",
        paste0("`", .fb_known_variance_columns(), "`", collapse = ", "),
        ", and `data` already has ",
        paste0("`", taken, "`", collapse = ", "),
        ". Rename ", if (length(taken) > 1L) "them" else "it", "."
      )
    ))
  }
  data[[.FB_KNOWN_VARIANCE_COL]] <- as.numeric(known_variance)
  term <- list(
    type = "known_variance",
    column = .FB_KNOWN_VARIANCE_COL,
    values = as.numeric(known_variance)
  )
  if (!is.null(covariance)) {
    data[[.FB_KNOWN_ROW_COL]] <- seq_len(nrow(data))
    term$matrix <- covariance
    term$row_column <- .FB_KNOWN_ROW_COL
  }
  fb$addition_terms <- c(fb$addition_terms %||% list(), list(term))
  list(fb = fb, data = data)
}

# .fb_known_variance_term() --- the IR's known-variance record, or NULL.
.fb_known_variance_term <- function(fb) {
  for (t in fb$addition_terms %||% list()) {
    if (identical(t$type %||% "", "known_variance")) {
      return(t)
    }
  }
  NULL
}

# .fb_known_variance_fill() --- zero on the `n_added` cells the design grid
# appended to the data (R/fb_complete_grid.R), which carry row 1's value in
# every column the model does not name. Such a cell has no response, and no
# first stage estimated it.
.fb_known_variance_fill <- function(fb, data, n_added) {
  if (is.null(.fb_known_variance_term(fb)) || n_added < 1L) {
    return(data)
  }
  rows <- nrow(data) - seq_len(n_added) + 1L
  data[[.FB_KNOWN_VARIANCE_COL]][rows] <- 0
  if (.FB_KNOWN_ROW_COL %in% names(data)) {
    data[[.FB_KNOWN_ROW_COL]][rows] <- NA_integer_
  }
  data
}

# .fb_known_variance_values() --- the known variance of each row of `data`,
# or NULL when the model carries none.
.fb_known_variance_values <- function(fb, data) {
  if (is.null(.fb_known_variance_term(fb))) {
    return(NULL)
  }
  as.numeric(data[[.FB_KNOWN_VARIANCE_COL]])
}

# .fb_known_covariance_values(): the known covariance over the rows of
# `data`, or NULL when the model carries a vector or nothing. A row without
# an index, a cell the design grid added, has a zero row and column.
.fb_known_covariance_values <- function(fb, data) {
  term <- .fb_known_variance_term(fb)
  if (is.null(term$matrix)) {
    return(NULL)
  }
  index <- data[[.FB_KNOWN_ROW_COL]]
  if (is.null(index) || length(index) != nrow(data)) {
    .fb_known_variance_refuse(paste0(
      "The model carries a known covariance matrix, and its data has no `",
      .FB_KNOWN_ROW_COL, "` column holding each row's index into it."
    ))
  }
  out <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  keep <- which(!is.na(index))
  out[keep, keep] <- term$matrix[index[keep], index[keep]]
  out
}

# .fb_known_variance_reissue(): the `known_variance` argument that refits
# the rows of `data`: the covariance over them, or the vector.
.fb_known_variance_reissue <- function(fb, data) {
  .fb_known_covariance_values(fb, data) %||% data[[.FB_KNOWN_VARIANCE_COL]]
}

# .fb_inla_known_index(): each row's generic0 index: its position among
# the rows with a positive known variance, NA elsewhere.
.fb_inla_known_index <- function(covariance) {
  positive <- diag(covariance) > 0
  index <- rep(NA_integer_, nrow(covariance))
  index[positive] <- seq_len(sum(positive))
  index
}

# .fb_inla_known_variance_formula_term(): the INLA row effect: a weighted
# iid effect for a vector, generic0 over the known precision for a matrix.
.fb_inla_known_variance_formula_term <- function(fb = NULL) {
  if (!is.null(.fb_known_variance_term(fb)$matrix)) {
    return(paste0(
      "f(", .FB_INLA_KNOWN_INDEX, ", model = \"generic0\", Cmatrix = ",
      .FB_INLA_KNOWN_PRECISION, ", ",
      "hyper = list(prec = list(initial = 0, fixed = TRUE)))"
    ))
  }
  paste0(
    "f(", .FB_INLA_KNOWN_INDEX, ", ", .FB_INLA_KNOWN_SCALE,
    ", model = \"iid\", ",
    "hyper = list(prec = list(initial = 0, fixed = TRUE)))"
  )
}

# .fb_inla_known_variance_effect() --- the posterior mean of the row effect on
# each row, on the response scale, or NULL. INLA's linear predictor carries
# it and flexyBayes's fitted values do not.
.fb_inla_known_variance_effect <- function(object) {
  fb <- object$extras$fb_terms %||% object$fb
  v <- .fb_known_variance_values(fb, object$data)
  effect <- object$inla$summary.random[[.FB_INLA_KNOWN_INDEX]]
  if (is.null(v) || is.null(effect)) {
    return(NULL)
  }
  covariance <- .fb_known_covariance_values(fb, object$data)
  if (!is.null(covariance)) {
    # A generic0 effect carries no weight: its mapped mean is the row's part.
    out <- effect$mean[match(.fb_inla_known_index(covariance), effect$ID)]
    out[is.na(out)] <- 0
    return(out)
  }
  scale <- object$extras$response_standardisation$scale %||% 1
  # summary.random has been mapped back to the response scale, a factor of
  # `scale`, so each row adds w_i times its mean, w_i = sqrt(v_i) / scale.
  sqrt(v) / scale * effect$mean[match(seq_along(v), effect$ID)]
}


# --- fb_select() ------------------------------------------------------- #
#
# The fold covariance of the benchmark and optimiser workflows under
# `covariance = "base_model"` carries each row's known variance on V's
# diagonal, sigma^2 + v_i, and a known covariance matrix K in full,
# sigma^2 I + K, with sigma^2 read off the fold refit, which is
# itself fitted with the known variances of the rows it is handed
# (R/selection_cov_model.R). The INLA refit of each fold's chosen support
# (`refit = list(method = "bayes")`) is fitted with it, and each held-out
# row's predictive variance gains v_i (R/selection_value.R); on brms, or
# where the held-out density is read off the engine's draws, that refit
# refuses. Every other route leaves it out, and refuses: the REML plug-in
# estimates one residual variance common to every row, and the Bayesian
# workflow and the projection's reference fit score or project posterior
# draws whose predictive variance omits it.

# .fb_sel_known_variance() --- each row's known variance, aligned to the
# base's data, or NULL when the base carries none.
.fb_sel_known_variance <- function(base) {
  fb <- base$fb %||% base$extras$fb_terms
  if (is.null(.fb_known_variance_term(fb))) {
    return(NULL)
  }
  v <- base$data[[.FB_KNOWN_VARIANCE_COL]]
  if (is.null(v) || length(v) != nrow(base$data)) {
    .fb_known_variance_refuse(paste0(
      "fb_select(): the base carries a known variance per row, and its data ",
      "has no `", .FB_KNOWN_VARIANCE_COL, "` column holding one value per ",
      "row. Refit the base with flexybayes(known_variance = ) on the data ",
      "fb_select() is to read."
    ))
  }
  as.numeric(v)
}

# .fb_sel_known_variance_refuse_route() --- a route that would leave the
# known variance out.
.fb_sel_known_variance_refuse_route <- function(route) {
  .fb_known_variance_refuse(paste0(
    "fb_select(): the base carries a known variance per row, which ", route,
    " would leave out. The penalised optimiser and the unpenalised ",
    "benchmark carry it into the fold covariance with ",
    "covariance = \"base_model\"."
  ))
}

# .fb_sel_known_variance_refuse_reml() --- the REML plug-in's refusal.
.fb_sel_known_variance_refuse_reml <- function() {
  .fb_sel_known_variance_refuse_route(paste0(
    "the REML plug-in (covariance = \"reml_plugin\"), which estimates one ",
    "residual variance common to every row,"
  ))
}

# .fb_select_require_known_variance_route() --- refuse, before any fit, a
# route that would leave a base's known variance out.
.fb_select_require_known_variance_route <- function(base,
                                                    workflow,
                                                    covariance,
                                                    value_request,
                                                    block_scope = FALSE) {
  if (is.null(.fb_sel_known_variance(base))) {
    return(invisible(NULL))
  }
  # The reference fit of a block-scoped prior carries the known variance in
  # its likelihood, and with no out-of-fold score nothing else is fitted.
  if (identical(workflow, "bayes") && isTRUE(block_scope)) {
    return(invisible(NULL))
  }
  if (identical(workflow, "bayes")) {
    .fb_sel_known_variance_refuse_route(paste0(
      "the Bayesian workflow, which scores each fold refit's posterior ",
      "predictive draws,"
    ))
  }
  if (identical(covariance, "reml_plugin")) {
    .fb_sel_known_variance_refuse_reml()
  }
  if (!is.null(value_request) &&
        !identical(.fb_sel_bayes_provider(base), "inla")) {
    .fb_sel_known_variance_refuse_route(paste0(
      "the Bayesian refit of each fold's chosen support (`refit`) on brms, ",
      "whose held-out log predictive density reads posterior draws,"
    ))
  }
  invisible(NULL)
}

# .fb_brms_known_se_newdata() --- brms reads the se() column of `newdata`
# even for a posterior mean, which does not depend on it. A new row's known
# variance is not part of its mean, so an absent column is supplied as zero.
.fb_brms_known_se_newdata <- function(object, newdata) {
  if (
    is.null(.fb_known_variance_term(.fb_fit_ir(object))) ||
      .FB_BRMS_KNOWN_SE_COL %in% names(newdata)
  ) {
    return(newdata)
  }
  newdata[[.FB_BRMS_KNOWN_SE_COL]] <- 0
  newdata
}
