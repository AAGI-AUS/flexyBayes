# Genomic-prediction cross-validation (G5) --- the payoff layer.
#
# Genomic selection is judged by *prediction accuracy*: how well a model
# trained on phenotyped genotypes predicts the field performance of
# genotypes whose phenotypes were held out, using only the genomic
# relationship matrix to borrow strength from relatives. `fb_gblup_cv()`
# runs that cross-validation analytically -- it is the GBLUP best-linear-
# unbiased-prediction equation evaluated per fold, leaning on the same
# spectral primitive (G0a) and REML estimator (`.fb_reml_vc()`, G3) the
# rest of the genomics track uses, so no MCMC backend is needed. The
# expensive part (the eigendecomposition of the training relationship
# matrix) is the shared primitive; everything else is matrix algebra.
#
# This is the analytic complement to a full Bayesian CV (which would re-fit
# the actual posterior per fold and propagate prediction uncertainty); the
# analytic accuracy is the field-standard metric and is fast enough to
# repeat many times.

# --- the cross-validation ----------------------------------------- #

#' Genomic-prediction accuracy by cross-validation
#'
#' Estimate how accurately a genomic BLUP predicts the performance of
#' genotypes whose phenotypes are held out. For each fold the variance
#' components are estimated by REML on the training genotypes, the held-out
#' breeding values are predicted from the genomic relationship matrix
#' (\eqn{\hat{u}_V = K_{VT}(K_{TT} + \delta I)^{-1}(y_T - X_T\hat\beta)}),
#' and the predictions are compared with the observed phenotypes. The
#' headline metric is *prediction accuracy*, the correlation between
#' predicted and observed across the held-out genotypes; *bias* is the
#' slope of observed on predicted (1 is unbiased).
#'
#' @param y Numeric vector of phenotypes, one per genotype, aligned to the
#'   rows of `K`.
#' @param K An `n x n` genomic (or pedigree) relationship matrix.
#' @param X Optional `n x p` fixed-effect design matrix (an intercept is
#'   used when `NULL`).
#' @param folds Number of cross-validation folds (default 5).
#' @param repeats Number of independent fold partitions to average over
#'   (default 1); repeated CV stabilises the estimate.
#' @param seed Optional integer seed for reproducible fold assignment.
#' @param tol PSD / symmetry tolerance for the spectral decomposition.
#'
#' @return An `fb_gblup_cv` object: `accuracy` (pooled predicted-observed
#'   correlation), `accuracy_per_fold`, `bias` (pooled observed-on-
#'   predicted slope), `predicted` (out-of-fold predictions in input
#'   order, averaged over repeats), `h2` (REML heritability on the full
#'   data), and metadata.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 80L; m <- 400L
#' Z <- matrix(rbinom(n * m, 2L, 0.3), n, m)
#' Zc <- scale(Z, scale = FALSE)
#' G <- tcrossprod(Zc) / m + diag(n) * 1e-3
#' u <- as.vector(t(chol(G)) %*% rnorm(n))
#' y <- u + rnorm(n)
#' cv <- fb_gblup_cv(y, G, folds = 5L, seed = 1L)
#' cv$accuracy
#' }
#' @export
fb_gblup_cv <- function(
  y,
  K,
  X = NULL,
  folds = 5L,
  repeats = 1L,
  seed = NULL,
  tol = 1e-8
) {
  K <- as.matrix(K)
  n <- length(y)
  if (nrow(K) != n || ncol(K) != n) {
    stop(
      "`K` must be an ", n, " x ", n, " relationship matrix aligned to `y`; ",
      "got ", nrow(K), " x ", ncol(K), ".",
      call. = FALSE
    )
  }
  if (is.null(X)) {
    X <- matrix(1, n, 1L)
  } else {
    X <- as.matrix(X)
    if (nrow(X) != n) {
      stop("`X` must have ", n, " rows.", call. = FALSE)
    }
  }
  folds <- as.integer(folds)
  if (folds < 2L || folds > n) {
    stop("`folds` must be between 2 and ", n, ".", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  pred_acc <- matrix(NA_real_, repeats, folds)
  pred_pooled <- matrix(NA_real_, n, repeats)

  for (rep_i in seq_len(repeats)) {
    fold_id <- .fb_cv_folds(n, folds)
    preds <- rep(NA_real_, n)
    for (f in seq_len(folds)) {
      test <- which(fold_id == f)
      train <- which(fold_id != f)
      preds[test] <- .fb_gblup_predict(
        y, X, K, train, test, tol = tol
      )
      pred_acc[rep_i, f] <- suppressWarnings(
        stats::cor(preds[test], y[test])
      )
    }
    pred_pooled[, rep_i] <- preds
  }

  predicted <- rowMeans(pred_pooled, na.rm = TRUE)
  accuracy <- suppressWarnings(stats::cor(predicted, y))
  bias <- tryCatch(
    unname(stats::coef(stats::lm(y ~ predicted))[2L]),
    error = function(e) NA_real_
  )

  # REML heritability on the full data for context.
  spec_full <- .fb_spectral(K, tol = tol)
  h2_full <- .fb_reml_vc(spec_full, y, X)$h2

  structure(
    list(
      accuracy = accuracy,
      accuracy_per_fold = as.vector(t(pred_acc)),
      accuracy_mean = mean(pred_acc, na.rm = TRUE),
      accuracy_se = stats::sd(as.vector(pred_acc), na.rm = TRUE) /
        sqrt(sum(!is.na(pred_acc))),
      bias = bias,
      predicted = predicted,
      h2 = h2_full,
      n = n,
      folds = folds,
      repeats = repeats
    ),
    class = c("fb_gblup_cv", "list")
  )
}

# --- display ------------------------------------------------------ #

#' @exportS3Method print fb_gblup_cv
print.fb_gblup_cv <- function(x, ...) {
  cat("<fb_gblup_cv>  genomic-prediction cross-validation\n")
  cat(
    "  ", x$n, " genotypes, ", x$folds, "-fold",
    if (x$repeats > 1L) paste0(" x ", x$repeats, " repeats") else "", "\n",
    sep = ""
  )
  cat("  REML heritability (full data): ", format(round(x$h2, 3L)), "\n",
    sep = ""
  )
  cat(
    "  prediction accuracy: ", format(round(x$accuracy, 3L)),
    "  (per-fold mean ", format(round(x$accuracy_mean, 3L)),
    " +/- ", format(round(x$accuracy_se, 3L)), ")\n",
    sep = ""
  )
  cat("  bias (obs ~ pred slope): ", format(round(x$bias, 3L)), "\n", sep = "")
  invisible(x)
}

# --- internal helpers --------------------------------------------- #

# GBLUP prediction of held-out genotypes from the training fit. Estimates
# the variance components by REML on the training block, then predicts the
# test breeding values via the relationship matrix:
#   u_hat_V = K_VT (K_TT + delta I)^{-1} (y_T - X_T beta_hat),
# implemented through the spectral decomposition of K_TT so the inverse is
# diag(1 / (lambda + delta)) in the eigenbasis. The predicted phenotype
# adds back the fixed-effect part X_V beta_hat.
.fb_gblup_predict <- function(y, X, K, train, test, tol = 1e-8) {
  K_tt <- K[train, train, drop = FALSE]
  K_vt <- K[test, train, drop = FALSE]
  y_t <- y[train]
  X_t <- X[train, , drop = FALSE]
  X_v <- X[test, , drop = FALSE]

  spec <- .fb_spectral(K_tt, tol = tol, name = "K_train")
  reml <- .fb_reml_vc(spec, y_t, X_t)
  beta <- reml$beta
  resid_t <- y_t - as.vector(X_t %*% beta)

  # alpha = (K_TT + delta I)^{-1} resid_t via the spectrum.
  ut_r <- as.vector(crossprod(spec$vectors, resid_t))
  alpha <- spec$vectors %*% (ut_r / (spec$values + reml$delta))
  u_hat_v <- as.vector(K_vt %*% alpha)
  as.vector(X_v %*% beta) + u_hat_v
}

# Balanced random fold assignment: a permutation of recycled fold labels.
.fb_cv_folds <- function(n, folds) {
  sample(rep(seq_len(folds), length.out = n))
}
