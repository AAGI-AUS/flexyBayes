# Genomic-prediction cross-validation (G5). The per-fold prediction is the
# exact GBLUP best-linear-unbiased-prediction equation, validated against a
# brute-force full-matrix solve; the behavioural checks confirm accuracy is
# positive, rises with heritability, and collapses to zero under a null.

# Brute-force GBLUP prediction for one split (no spectral shortcut) -- the
# exact reference .fb_gblup_predict() must reproduce, at the same REML
# variance ratio.
.brute_gblup_predict <- function(y, X, K, train, test, delta, beta) {
  K_tt <- K[train, train, drop = FALSE]
  K_vt <- K[test, train, drop = FALSE]
  resid_t <- y[train] - as.vector(X[train, , drop = FALSE] %*% beta)
  alpha <- solve(K_tt + delta * diag(length(train)), resid_t)
  as.vector(X[test, , drop = FALSE] %*% beta) + as.vector(K_vt %*% alpha)
}

test_that(".fb_gblup_predict() reproduces the exact full-matrix GBLUP prediction", {
  G <- sim_kinship(n_geno = 60L, n_markers = 400L, seed = 6L)
  sim <- sim_gblup_pheno(G, var_g = 2, var_e = 1, n_rep = 1L, seed = 7L)
  y <- sim$data$y
  X <- matrix(1, length(y), 1L)
  train <- 1:45
  test <- 46:60

  spec <- flexyBayes:::.fb_spectral(G[train, train])
  reml <- flexyBayes:::.fb_reml_vc(spec, y[train], X[train, , drop = FALSE])
  spectral_pred <- flexyBayes:::.fb_gblup_predict(y, X, G, train, test)
  brute_pred <- .brute_gblup_predict(y, X, G, train, test, reml$delta, reml$beta)

  expect_equal(spectral_pred, brute_pred, tolerance = 1e-7)
})

test_that("fb_gblup_cv() gives strong accuracy when genotypes are related", {
  # Related genotypes (founder structure) are the realistic genomic-
  # selection setting: a held-out line is predicted from its relatives via
  # the relationship matrix. Accuracy is then substantial.
  Z <- sim_related_markers(n_geno = 150L, n_markers = 600L, n_founders = 10L,
    seed = 8L)
  G <- kinship_from_markers(Z)
  sim <- sim_gblup_pheno(G, var_g = 3, var_e = 1, n_rep = 1L, seed = 9L)
  cv <- fb_gblup_cv(sim$data$y, G, folds = 5L, seed = 11L)

  expect_s3_class(cv, "fb_gblup_cv")
  expect_gt(cv$accuracy, 0.4)
  expect_length(cv$accuracy_per_fold, 5L)
  # Predictions are roughly unbiased (slope near 1, generous tolerance).
  expect_gt(cv$bias, 0.4)
  expect_lt(cv$bias, 1.8)
})

test_that("fb_gblup_cv() accuracy rises with heritability", {
  G <- sim_kinship(n_geno = 160L, n_markers = 700L, seed = 12L)
  lo <- sim_gblup_pheno(G, var_g = 0.2, var_e = 1, n_rep = 1L, seed = 13L)
  hi <- sim_gblup_pheno(G, var_g = 6, var_e = 1, n_rep = 1L, seed = 13L)
  cv_lo <- fb_gblup_cv(lo$data$y, G, folds = 5L, seed = 14L)
  cv_hi <- fb_gblup_cv(hi$data$y, G, folds = 5L, seed = 14L)
  expect_gt(cv_hi$accuracy, cv_lo$accuracy)
  expect_gt(cv_hi$h2, cv_lo$h2)
})

test_that("fb_gblup_cv() accuracy is near zero when the trait is unrelated to K", {
  G <- sim_kinship(n_geno = 150L, n_markers = 600L, seed = 15L)
  set.seed(16L)
  y_null <- stats::rnorm(150L) # pure noise, no genetic signal
  cv <- fb_gblup_cv(y_null, G, folds = 5L, repeats = 2L, seed = 17L)
  expect_lt(abs(cv$accuracy), 0.2)
})

test_that("fb_gblup_cv() repeated CV runs and averages out-of-fold predictions", {
  G <- sim_kinship(n_geno = 80L, n_markers = 300L, seed = 18L)
  sim <- sim_gblup_pheno(G, var_g = 2, var_e = 1, n_rep = 1L, seed = 19L)
  cv <- fb_gblup_cv(sim$data$y, G, folds = 4L, repeats = 3L, seed = 20L)
  expect_length(cv$accuracy_per_fold, 12L)
  expect_equal(length(cv$predicted), 80L)
  expect_true(is.finite(cv$accuracy_se))
})

test_that("fb_gblup_cv() validates its arguments", {
  G <- sim_kinship(n_geno = 20L, n_markers = 100L, seed = 21L)
  y <- stats::rnorm(20L)
  expect_error(fb_gblup_cv(y, G[1:10, 1:10]), "aligned to")
  expect_error(fb_gblup_cv(y, G, folds = 1L), "between 2")
  expect_error(fb_gblup_cv(y, G, X = matrix(1, 10L, 1L)), "rows")
})

test_that("print.fb_gblup_cv() renders accuracy and heritability", {
  G <- sim_kinship(n_geno = 60L, n_markers = 300L, seed = 22L)
  sim <- sim_gblup_pheno(G, var_g = 2, var_e = 1, n_rep = 1L, seed = 23L)
  cv <- fb_gblup_cv(sim$data$y, G, folds = 5L, seed = 24L)
  out <- utils::capture.output(print(cv))
  expect_true(any(grepl("<fb_gblup_cv>", out, fixed = TRUE)))
  expect_true(any(grepl("prediction accuracy", out)))
  expect_true(any(grepl("heritability", out)))
})
