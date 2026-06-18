# Spectral primitive (G-I4) --- the shared efficiency machinery for the
# genomics / MET expansion. These tests pin the exact numerical
# contract base R + Matrix can verify without any backend: the
# decomposition identities, the rotation operators, and --- the
# load-bearing claim --- that the rotated weighted least squares
# reproduces full-covariance generalised least squares (the EMMAX /
# P3D fast path the GWAS scan rests on) and that the spectral
# log-determinant matches the brute-force value.

# Small SPD relationship-matrix fixture: G = Z Z' / m for random
# markers, the VanRaden shape, guaranteed PSD.
.fb_spectral_kinship <- function(n = 12L, n_markers = 60L, seed = 1L) {
  set.seed(seed)
  Z <- matrix(stats::rbinom(n * n_markers, 2L, 0.3), n, n_markers)
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  tcrossprod(Zc) / n_markers + diag(n) * 1e-6
}

# ---------------------------------------------------------------- #
# (a) Decomposition identities.                                     #
# ---------------------------------------------------------------- #

test_that(".fb_spectral() reconstructs K = U Lambda U' to machine tol", {
  K <- .fb_spectral_kinship()
  spec <- flexyBayes:::.fb_spectral(K, name = "G")
  K_hat <- spec$vectors %*% diag(spec$values) %*% t(spec$vectors)
  expect_equal(K_hat, K, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that(".fb_spectral() returns orthonormal eigenvectors and the eigen-equation", {
  K <- .fb_spectral_kinship()
  spec <- flexyBayes:::.fb_spectral(K)
  expect_equal(crossprod(spec$vectors), diag(spec$n), tolerance = 1e-9)
  # K U = U Lambda.
  expect_equal(
    K %*% spec$vectors,
    sweep(spec$vectors, 2L, spec$values, `*`),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that(".fb_spectral() eigenvalues are descending and non-negative", {
  K <- .fb_spectral_kinship()
  spec <- flexyBayes:::.fb_spectral(K)
  expect_equal(spec$values, sort(spec$values, decreasing = TRUE))
  expect_true(all(spec$values >= 0))
  expect_equal(spec$rank, spec$n)
  expect_equal(spec$rank_full, spec$n)
})

# ---------------------------------------------------------------- #
# (b) Rotation operators.                                          #
# ---------------------------------------------------------------- #

test_that(".fb_spectral_rotate() applies U' and round-trips at full rank", {
  K <- .fb_spectral_kinship()
  spec <- flexyBayes:::.fb_spectral(K)
  y <- stats::rnorm(spec$n)
  y_star <- flexyBayes:::.fb_spectral_rotate(spec, y)
  expect_equal(as.vector(y_star), as.vector(crossprod(spec$vectors, y)))
  back <- flexyBayes:::.fb_spectral_backrotate(spec, y_star)
  expect_equal(as.vector(back), y, tolerance = 1e-9)
})

test_that(".fb_spectral_rotate() handles a design matrix and rejects wrong nrow", {
  K <- .fb_spectral_kinship(n = 10L)
  spec <- flexyBayes:::.fb_spectral(K)
  X <- cbind(1, stats::rnorm(spec$n))
  X_star <- flexyBayes:::.fb_spectral_rotate(spec, X)
  expect_equal(dim(X_star), c(spec$rank, 2L))
  expect_equal(X_star, crossprod(spec$vectors, X), ignore_attr = TRUE)
  expect_error(
    flexyBayes:::.fb_spectral_rotate(spec, stats::rnorm(spec$n + 1L)),
    "rows"
  )
})

test_that(".fb_spectral_sqrt() gives B B' = K (full) and K_k (truncated)", {
  K <- .fb_spectral_kinship(n = 14L)
  spec <- flexyBayes:::.fb_spectral(K)
  B <- flexyBayes:::.fb_spectral_sqrt(spec)
  expect_equal(tcrossprod(B), K, tolerance = 1e-8, ignore_attr = TRUE)

  spec3 <- flexyBayes:::.fb_spectral(K, rank = 3L)
  B3 <- flexyBayes:::.fb_spectral_sqrt(spec3)
  K3 <- spec3$vectors %*% diag(spec3$values) %*% t(spec3$vectors)
  expect_equal(tcrossprod(B3), K3, tolerance = 1e-8, ignore_attr = TRUE)
})

# ---------------------------------------------------------------- #
# (c) The EMMAX / GLS equivalence --- the load-bearing claim.       #
# ---------------------------------------------------------------- #

test_that(".fb_spectral_dvar() equals diag(U' (var_g K + var_e I) U)", {
  K <- .fb_spectral_kinship()
  spec <- flexyBayes:::.fb_spectral(K)
  var_g <- 1.7
  var_e <- 0.8
  V <- var_g * K + var_e * diag(spec$n)
  rotated_V <- crossprod(spec$vectors, V %*% spec$vectors)
  expect_equal(
    flexyBayes:::.fb_spectral_dvar(spec, var_g, var_e),
    diag(rotated_V),
    tolerance = 1e-8
  )
  # Off-diagonals of the rotated covariance vanish (it is diagonal).
  expect_lt(max(abs(rotated_V - diag(diag(rotated_V)))), 1e-7)
})

test_that(".fb_spectral_logdet() matches the brute-force log-determinant", {
  K <- .fb_spectral_kinship()
  spec <- flexyBayes:::.fb_spectral(K)
  var_g <- 2.3
  var_e <- 0.5
  V <- var_g * K + var_e * diag(spec$n)
  brute <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
  expect_equal(flexyBayes:::.fb_spectral_logdet(spec, var_g, var_e), brute,
    tolerance = 1e-8
  )
})

test_that("rotated weighted least squares reproduces full-covariance GLS", {
  # The EMMAX claim: solving beta by weighting the rotated data with
  # 1 / (var_g lambda + var_e) equals the generalised-least-squares
  # solution with the full covariance V = var_g K + var_e I.
  K <- .fb_spectral_kinship(n = 20L, n_markers = 80L, seed = 7L)
  n <- nrow(K)
  spec <- flexyBayes:::.fb_spectral(K)
  set.seed(99L)
  X <- cbind(1, stats::rnorm(n), stats::rnorm(n))
  y <- stats::rnorm(n)
  var_g <- 1.4
  var_e <- 0.6
  V <- var_g * K + var_e * diag(n)
  Vinv <- solve(V)
  beta_gls <- solve(crossprod(X, Vinv %*% X), crossprod(X, Vinv %*% y))

  y_star <- flexyBayes:::.fb_spectral_rotate(spec, y)
  X_star <- flexyBayes:::.fb_spectral_rotate(spec, X)
  w <- 1 / flexyBayes:::.fb_spectral_dvar(spec, var_g, var_e)
  XtWX <- crossprod(X_star, X_star * w)
  XtWy <- crossprod(X_star, as.vector(y_star) * w)
  beta_wls <- solve(XtWX, XtWy)

  expect_equal(as.vector(beta_wls), as.vector(beta_gls), tolerance = 1e-7)
})

# ---------------------------------------------------------------- #
# (d) Truncation + capture metrics.                                #
# ---------------------------------------------------------------- #

test_that("truncation keeps the leading eigenpairs with correct capture", {
  K <- .fb_spectral_kinship(n = 16L)
  full <- flexyBayes:::.fb_spectral(K)
  k <- 4L
  spec <- flexyBayes:::.fb_spectral(K, rank = k)
  expect_equal(spec$rank, k)
  expect_equal(spec$values, full$values[seq_len(k)])
  expect_equal(
    spec$capture_trace,
    sum(full$values[seq_len(k)]) / sum(full$values)
  )
  expect_equal(
    spec$capture_frobenius,
    sum(full$values[seq_len(k)]^2) / sum(full$values^2)
  )
  # Frobenius capture is the complement of the relative squared
  # truncation error ||K - K_k||_F^2 / ||K||_F^2.
  K_k <- spec$vectors %*% diag(spec$values) %*% t(spec$vectors)
  rel_err <- sum((K - K_k)^2) / sum(K^2)
  expect_equal(1 - spec$capture_frobenius, rel_err, tolerance = 1e-8)
})

test_that("full-rank capture is 1 on both metrics", {
  K <- .fb_spectral_kinship()
  spec <- flexyBayes:::.fb_spectral(K)
  expect_equal(spec$capture_trace, 1, tolerance = 1e-12)
  expect_equal(spec$capture_frobenius, 1, tolerance = 1e-12)
})

# ---------------------------------------------------------------- #
# (e) PSD contract: clamp noise, refuse indefiniteness.            #
# ---------------------------------------------------------------- #

test_that(".fb_spectral() clamps numerical-noise negatives and records it", {
  # A rank-deficient PSD matrix perturbed by tiny negative noise on the
  # null directions: should clamp, not refuse, and report the repair.
  set.seed(3L)
  U <- qr.Q(qr(matrix(stats::rnorm(36L), 6L, 6L)))
  lam <- c(5, 3, 1, 0, 0, 0) + c(0, 0, 0, -1e-12, -2e-12, -5e-13)
  K <- U %*% diag(lam) %*% t(U)
  K <- (K + t(K)) / 2
  spec <- flexyBayes:::.fb_spectral(K)
  expect_true(all(spec$values >= 0))
  expect_gt(spec$negative_clamped$count, 0L)
  expect_lt(spec$negative_clamped$min, 0)
})

test_that(".fb_spectral() refuses a genuinely indefinite matrix", {
  set.seed(5L)
  U <- qr.Q(qr(matrix(stats::rnorm(36L), 6L, 6L)))
  K <- U %*% diag(c(5, 3, 1, 0.5, 0.2, -2)) %*% t(U)
  K <- (K + t(K)) / 2
  expect_error(flexyBayes:::.fb_spectral(K, name = "G"), "positive-semidefinite")
})

test_that(".fb_spectral() refuses non-square and non-symmetric input", {
  expect_error(flexyBayes:::.fb_spectral(matrix(0, 3L, 4L)), "square")
  A <- matrix(c(1, 2, 0, 1), 2L, 2L)
  expect_error(flexyBayes:::.fb_spectral(A), "symmetric")
})

# ---------------------------------------------------------------- #
# (f) Marker-matrix SVD path.                                      #
# ---------------------------------------------------------------- #

test_that(".fb_spectral_from_markers() matches eigendecomposition of G = Zc Zc' / m", {
  set.seed(11L)
  n <- 15L
  m <- 50L
  Z <- matrix(stats::rbinom(n * m, 2L, 0.25), n, m)
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  G <- tcrossprod(Zc) / m

  from_marker <- flexyBayes:::.fb_spectral_from_markers(Z, m = m)
  from_kinship <- flexyBayes:::.fb_spectral(G + diag(n) * 0)

  # Eigenvalues agree (the SVD path may carry trailing ~0 values for
  # the n - rank(Zc) null directions; compare the leading n).
  ev_marker <- sort(from_marker$values, decreasing = TRUE)
  ev_kin <- sort(from_kinship$values, decreasing = TRUE)
  expect_equal(ev_marker[seq_len(n)], ev_kin, tolerance = 1e-7)

  # Reconstruction agrees regardless of eigenvector sign convention.
  G_marker <- from_marker$vectors %*% diag(from_marker$values) %*%
    t(from_marker$vectors)
  expect_equal(G_marker, G, tolerance = 1e-7, ignore_attr = TRUE)
})

test_that(".fb_spectral_from_markers() refuses a monomorphic marker under unit scaling", {
  Z <- cbind(stats::rbinom(10L, 2L, 0.3), rep(1L, 10L))
  expect_error(
    flexyBayes:::.fb_spectral_from_markers(Z, scale = TRUE),
    "non-finite"
  )
})

# ---------------------------------------------------------------- #
# (g) Determinism, guards, display.                                #
# ---------------------------------------------------------------- #

test_that(".fb_spectral() is deterministic for a fixed input", {
  K <- .fb_spectral_kinship(seed = 42L)
  a <- flexyBayes:::.fb_spectral(K)
  b <- flexyBayes:::.fb_spectral(K)
  expect_equal(a$values, b$values)
  expect_equal(a$vectors, b$vectors)
})

test_that(".fb_spectral() validates the rank argument", {
  K <- .fb_spectral_kinship(n = 8L)
  expect_error(flexyBayes:::.fb_spectral(K, rank = 0L), "rank")
  expect_error(flexyBayes:::.fb_spectral(K, rank = 9L), "rank")
  expect_error(flexyBayes:::.fb_spectral(K, rank = 2.5), "rank")
})

test_that(".fb_spectral_logdet() refuses a truncated object", {
  K <- .fb_spectral_kinship(n = 10L)
  spec <- flexyBayes:::.fb_spectral(K, rank = 4L)
  expect_error(
    flexyBayes:::.fb_spectral_logdet(spec, 1, 1),
    "full-rank"
  )
})

test_that(".fb_spectral_dvar() rejects negative variance components", {
  K <- .fb_spectral_kinship(n = 6L)
  spec <- flexyBayes:::.fb_spectral(K)
  expect_error(flexyBayes:::.fb_spectral_dvar(spec, -1, 1), "non-negative")
  expect_error(flexyBayes:::.fb_spectral_dvar(spec, 1, NA_real_), "non-negative")
})

test_that("print.fb_spectral() renders rank and capture", {
  K <- .fb_spectral_kinship(n = 10L)
  out <- utils::capture.output(print(flexyBayes:::.fb_spectral(K, rank = 3L)))
  expect_true(any(grepl("<fb_spectral>", out, fixed = TRUE)))
  expect_true(any(grepl("truncated", out, fixed = TRUE)))
  expect_true(any(grepl("variance captured", out, fixed = TRUE)))
})

test_that("is_fb_spectral() identifies the class", {
  K <- .fb_spectral_kinship(n = 6L)
  expect_true(flexyBayes:::is_fb_spectral(flexyBayes:::.fb_spectral(K)))
  expect_false(flexyBayes:::is_fb_spectral(K))
})
