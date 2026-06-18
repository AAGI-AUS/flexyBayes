# Known-truth simulators (G0c) --- the data-generating processes the
# genomics recovery cells ground against. These tests assert the
# simulators produce data matching the variance components / QTL effects
# they author (ground the config, not just the data), so a downstream
# recovery check compares a fit against trustworthy truth.

test_that("sim_kinship() returns a PSD relationship matrix scaled to unit diagonal", {
  G <- sim_kinship(n_geno = 30L, n_markers = 250L, seed = 2L)
  expect_equal(dim(G), c(30L, 30L))
  expect_equal(mean(diag(G)), 1, tolerance = 1e-8)
  expect_true(isSymmetric(unname(G), tol = 1e-8))
  expect_gt(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  expect_equal(rownames(G), paste0("g", seq_len(30L)))
})

test_that("sim_gblup_pheno() realises the authored breeding-value variance", {
  G <- sim_kinship(n_geno = 60L, n_markers = 400L, seed = 3L)
  sim <- sim_gblup_pheno(G, var_g = 2, var_e = 1, n_rep = 1L, seed = 4L)
  # The realised genetic variance of u (with K scaled to unit diagonal)
  # should match var_g within Monte-Carlo error at this n.
  expect_equal(stats::var(sim$u_true), 2, tolerance = 0.6)
  expect_equal(sim$h2, 2 / 3, tolerance = 1e-12)
  expect_equal(nrow(sim$data), 60L)
  expect_setequal(levels(sim$data$geno), rownames(G))
})

test_that("sim_gblup_pheno() truth tracks the heritability target", {
  G <- sim_kinship(n_geno = 80L, n_markers = 500L, seed = 5L)
  lo <- sim_gblup_pheno(G, var_g = 0.25, var_e = 1, seed = 6L)
  hi <- sim_gblup_pheno(G, var_g = 4, var_e = 1, seed = 6L)
  expect_lt(lo$h2, hi$h2)
  # Empirical genotype-level signal-to-noise rises with the target h2:
  # the spread of u_true relative to the residual SD is larger at hi.
  expect_lt(stats::sd(lo$u_true), stats::sd(hi$u_true))
})

test_that("sim_gwas_pheno() QTL markers carry stronger marginal signal than nulls", {
  sim <- sim_gwas_pheno(
    n_geno = 200L,
    n_markers = 300L,
    qtl_idx = c(50L, 150L, 250L),
    qtl_effect = c(2, -1.8, 1.5),
    var_e = 1,
    seed = 7L
  )
  Zstd <- scale(sim$markers)
  Zstd[!is.finite(Zstd)] <- 0
  # Marginal single-marker association: |t| from y ~ marker.
  abs_t <- vapply(seq_len(ncol(Zstd)), function(j) {
    f <- stats::lm(sim$y ~ Zstd[, j])
    abs(summary(f)$coefficients[2L, 3L])
  }, numeric(1))
  qtl_t <- abs_t[sim$qtl_idx]
  null_t <- abs_t[-sim$qtl_idx]
  # Every authored QTL beats the 95th percentile of the null markers.
  expect_true(all(qtl_t > stats::quantile(null_t, 0.95)))
})

test_that("sim_gwas_pheno() polygenic background adds variance without inventing QTL", {
  sim0 <- sim_gwas_pheno(var_poly = 0, seed = 8L)
  simp <- sim_gwas_pheno(var_poly = 1, seed = 8L)
  expect_length(sim0$y, length(simp$y))
  expect_equal(simp$var_poly, 1)
  expect_equal(simp$qtl_idx, sim0$qtl_idx)
})

test_that("sim_gwas_pheno() refuses mismatched QTL specification", {
  expect_error(
    sim_gwas_pheno(qtl_idx = c(1L, 2L), qtl_effect = 1),
    "same length"
  )
})
