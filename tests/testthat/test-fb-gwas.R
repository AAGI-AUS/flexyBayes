# Genome-wide association scan (G3). The EMMAX / P3D scan is deterministic
# and validated against exact, dependency-free oracles -- brute-force
# full-matrix REML for the variance components and brute-force per-marker
# generalised least squares for the marker tests (both exact by the
# spectral rotation = GLS identity established for .fb_spectral()) -- and
# cross-checked against sommer's independent REML (the Independent Oracle
# Principle: a number whose oracle I did not also author).

# ---------------------------------------------------------------- #
# (a) REML variance components vs exact + independent oracles.      #
# ---------------------------------------------------------------- #

# Brute-force full-matrix REML over the variance ratio delta -- the exact
# reference the rotated .fb_reml_vc() must reproduce.
.brute_reml <- function(y, X, K) {
  n <- length(y)
  p <- ncol(X)
  nll <- function(ld) {
    delta <- exp(ld)
    V <- K + delta * diag(n)
    Vi <- solve(V)
    xtvix <- crossprod(X, Vi %*% X)
    b <- solve(xtvix, crossprod(X, Vi %*% y))
    r <- y - X %*% b
    rss <- as.numeric(crossprod(r, Vi %*% r))
    0.5 * ((n - p) * log(rss) +
      as.numeric(determinant(V, logarithm = TRUE)$modulus) +
      as.numeric(determinant(xtvix, logarithm = TRUE)$modulus))
  }
  o <- stats::optimize(nll, c(-10, 10), tol = 1e-9)
  delta <- exp(o$minimum)
  V <- K + delta * diag(n)
  Vi <- solve(V)
  b <- solve(crossprod(X, Vi %*% X), crossprod(X, Vi %*% y))
  r <- y - X %*% b
  var_g <- as.numeric(crossprod(r, Vi %*% r)) / (n - p)
  list(var_g = var_g, var_e = delta * var_g, delta = delta)
}

test_that(".fb_reml_vc() reproduces brute-force full-matrix REML exactly", {
  G <- sim_kinship(n_geno = 60L, n_markers = 400L, seed = 3L)
  sim <- sim_gblup_pheno(G, var_g = 2, var_e = 1, n_rep = 1L, seed = 4L)
  y <- sim$data$y
  X <- matrix(1, length(y), 1L)

  spec <- flexyBayes:::.fb_spectral(G)
  fb <- flexyBayes:::.fb_reml_vc(spec, y, X)
  bf <- .brute_reml(y, X, G)

  expect_equal(fb$var_g, bf$var_g, tolerance = 1e-4)
  expect_equal(fb$var_e, bf$var_e, tolerance = 1e-4)
  expect_equal(fb$delta, bf$delta, tolerance = 1e-4)
  expect_equal(fb$h2, 1 / (1 + bf$delta), tolerance = 1e-4)
})

test_that(".fb_reml_vc() agrees with sommer's independent REML (Independent Oracle)", {
  skip_if_not_installed("sommer")
  G <- sim_kinship(n_geno = 60L, n_markers = 400L, seed = 3L)
  sim <- sim_gblup_pheno(G, var_g = 2, var_e = 1, n_rep = 1L, seed = 4L)
  y <- sim$data$y
  spec <- flexyBayes:::.fb_spectral(G)
  fb <- flexyBayes:::.fb_reml_vc(spec, y, matrix(1, length(y), 1L))

  d2 <- data.frame(geno = factor(rownames(G), levels = rownames(G)), y = y)
  s <- tryCatch(
    sommer::mmer(y ~ 1,
      random = ~ sommer::vsr(geno, Gu = G),
      data = d2, verbose = FALSE, dateWarning = FALSE
    ),
    error = function(e) NULL
  )
  skip_if(is.null(s), "sommer::mmer() did not converge in this environment")
  vc <- summary(s)$varcomp[, 1L]
  # vc[1] = genetic, vc[2] = residual (units).
  expect_equal(fb$var_g, unname(vc[[1L]]), tolerance = 1e-3)
  expect_equal(fb$var_e, unname(vc[[2L]]), tolerance = 1e-3)
})

# ---------------------------------------------------------------- #
# (b) EMMAX marker test vs exact per-marker GLS.                    #
# ---------------------------------------------------------------- #

test_that("the EMMAX scan statistic equals exact per-marker GLS at the null VCs", {
  set.seed(31L)
  n <- 70L
  m <- 40L
  G <- sim_kinship(n_geno = n, n_markers = 300L, seed = 5L)
  markers <- matrix(stats::rbinom(n * m, 2L, 0.3), n, m)
  colnames(markers) <- paste0("snp", seq_len(m))
  y <- as.numeric(2 * scale(markers[, 7L]) + MASS::mvrnorm(1L, rep(0, n), 0.5 * G) +
    stats::rnorm(n))
  dat <- data.frame(y = y)

  scan <- fb_gwas(y ~ 1, data = dat, markers = markers, K = G)

  # Exact per-marker GLS at the same (fixed) null variance components.
  V <- scan$var_g * G + scan$var_e * diag(n)
  Vi <- solve(V)
  X0 <- matrix(1, n, 1L)
  brute_chisq <- vapply(seq_len(m), function(j) {
    Xj <- cbind(X0, markers[, j])
    xtvix <- crossprod(Xj, Vi %*% Xj)
    b <- solve(xtvix, crossprod(Xj, Vi %*% y))
    covb <- solve(xtvix)
    (b[2L] / sqrt(covb[2L, 2L]))^2
  }, numeric(1))

  expect_equal(scan$results$statistic, brute_chisq, tolerance = 1e-6)
})

# ---------------------------------------------------------------- #
# (c) Known-QTL recovery + genomic-control calibration.            #
# ---------------------------------------------------------------- #

test_that("fb_gwas() detects simulated QTL with a background relationship matrix", {
  # Proper GWAS practice: the relationship matrix is built from
  # *background* markers, not the markers being tested -- testing a marker
  # that is also in K causes proximal contamination (the polygenic term
  # absorbs the marker's own signal, deflating the statistic). Here K
  # comes from an independent marker panel on the same individuals.
  sim <- sim_gwas_pheno(
    n_geno = 250L,
    n_markers = 500L,
    qtl_idx = c(80L, 250L, 420L),
    qtl_effect = c(2.2, -2.0, 1.8),
    var_e = 1,
    var_poly = 0.3,
    seed = 41L
  )
  set.seed(71L)
  bg <- matrix(stats::rbinom(250L * 800L, 2L, 0.25), 250L, 800L)
  Zc <- scale(bg, center = TRUE, scale = FALSE)
  K_bg <- tcrossprod(Zc) / 800L
  K_bg <- K_bg / mean(diag(K_bg)) + diag(250L) * 1e-6

  scan <- fb_gwas(y ~ 1, data = data.frame(y = sim$y), markers = sim$markers,
    K = K_bg)

  # Every authored QTL is genome-wide (Bonferroni) significant.
  expect_true(all(scan$results$p_bonferroni[sim$qtl_idx] < 0.05))
  # The QTL dominate the top of the ranking.
  expect_true(all(sim$qtl_idx %in% order(scan$results$p_value)[1:10]))
  # With an uncontaminated K the test statistics are well-calibrated.
  expect_gt(scan$lambda_gc, 0.8)
  expect_lt(scan$lambda_gc, 1.25)
})

test_that("fb_gwas() default K-from-markers detects QTL but is conservative (proximal contamination)", {
  # The convenience default (K built from the tested markers) still finds
  # strong QTL, but the shared-marker proximal contamination deflates the
  # genome-wide statistics -- conservative, documented, and the reason an
  # independent / leave-one-chromosome-out K is preferred.
  sim <- sim_gwas_pheno(
    n_geno = 250L, n_markers = 500L,
    qtl_idx = c(80L, 250L, 420L), qtl_effect = c(2.2, -2.0, 1.8),
    var_e = 1, var_poly = 0.3, seed = 41L
  )
  scan <- fb_gwas(y ~ 1, data = data.frame(y = sim$y), markers = sim$markers)
  expect_true(all(scan$results$p_bonferroni[sim$qtl_idx] < 0.05))
  expect_lt(scan$lambda_gc, 1) # deflation, not inflation
})

test_that("fb_gwas() under a pure polygenic null gives uniform p-values (lambda_GC ~ 1)", {
  sim <- sim_gwas_pheno(
    n_geno = 250L,
    n_markers = 400L,
    qtl_idx = integer(0),
    qtl_effect = numeric(0),
    var_e = 1,
    var_poly = 0.5,
    seed = 42L
  )
  dat <- data.frame(y = sim$y)
  scan <- fb_gwas(y ~ 1, data = dat, markers = sim$markers)
  expect_gt(scan$lambda_gc, 0.7)
  expect_lt(scan$lambda_gc, 1.3)
  # No false genome-wide hits under the null.
  expect_equal(sum(scan$results$p_bonferroni < 0.05), 0L)
})

# ---------------------------------------------------------------- #
# (d) Interface: K-from-markers, map join, guards, display.        #
# ---------------------------------------------------------------- #

test_that("fb_gwas() builds K from markers when none is supplied", {
  # Fewer markers than individuals -> the auto-built relationship matrix is
  # rank-deficient; the full eigendecomposition must still drive REML.
  sim <- sim_gwas_pheno(n_geno = 120L, n_markers = 80L,
    qtl_idx = c(20L, 60L), qtl_effect = c(2, -2), seed = 43L)
  dat <- data.frame(y = sim$y)
  scan_auto <- fb_gwas(y ~ 1, data = dat, markers = sim$markers)
  expect_s3_class(scan_auto, "fb_gwas")
  expect_equal(nrow(scan_auto$results), 80L)
  expect_true(is.finite(scan_auto$lambda_gc))
  expect_true(is.finite(scan_auto$h2))
})

test_that("fb_gwas() joins a marker map and the results carry FDR + Bonferroni", {
  sim <- sim_gwas_pheno(n_geno = 100L, n_markers = 30L, qtl_idx = 10L,
    qtl_effect = 2, seed = 44L)
  colnames(sim$markers) <- paste0("snp", seq_len(30L))
  map <- data.frame(
    marker = paste0("snp", seq_len(30L)),
    chr = rep(1:3, each = 10L),
    pos = rep(seq_len(10L), 3L)
  )
  scan <- fb_gwas(y ~ 1, data = data.frame(y = sim$y), markers = sim$markers,
    marker_map = map)
  expect_true(all(c("chr", "pos") %in% names(scan$results)))
  expect_true(all(scan$results$q_value >= scan$results$p_value - 1e-9))
  expect_true(all(scan$results$p_bonferroni >= scan$results$p_value - 1e-9))
})

test_that("fb_gwas() refuses a marker matrix with the wrong row count", {
  sim <- sim_gwas_pheno(n_geno = 50L, n_markers = 20L, qtl_idx = 5L,
    qtl_effect = 2, seed = 45L)
  expect_error(
    fb_gwas(y ~ 1, data = data.frame(y = sim$y[1:40]), markers = sim$markers),
    "one row per individual"
  )
})

test_that("print.fb_gwas() renders the headline diagnostics", {
  sim <- sim_gwas_pheno(n_geno = 80L, n_markers = 50L, qtl_idx = 25L,
    qtl_effect = 2.5, seed = 46L)
  scan <- fb_gwas(y ~ 1, data = data.frame(y = sim$y), markers = sim$markers)
  out <- utils::capture.output(print(scan))
  expect_true(any(grepl("<fb_gwas>", out, fixed = TRUE)))
  expect_true(any(grepl("lambda_GC", out)))
  expect_true(any(grepl("h\\^2", out)))
})
