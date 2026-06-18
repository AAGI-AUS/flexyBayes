# fb_structured_cov(): identified covariance reconstruction for fa terms.
# Driven by synthetic posterior draws so the contract is exact and fast.

mk_fa_fit <- function(trueL, truePsi, n_iter = 200L, sd_noise = 0.02) {
  no <- nrow(trueL)
  k <- ncol(trueL)
  tag <- "geno_env_fa2"
  lambda_nm <- function(i, j) sprintf("Lambda_%s[%d,%d]", tag, i, j)
  psi_nm <- function(i) sprintf("psi_%s[%d,1]", tag, i)
  cols <- c(
    as.vector(vapply(
      seq_len(k), function(j) vapply(seq_len(no), lambda_nm, character(1), j = j),
      character(no)
    )),
    vapply(seq_len(no), psi_nm, character(1))
  )
  mk_chain <- function(seed) {
    set.seed(seed)
    M <- matrix(0, n_iter, length(cols), dimnames = list(NULL, cols))
    for (j in seq_len(k)) {
      for (i in seq_len(no)) {
        M[, lambda_nm(i, j)] <- trueL[i, j] + stats::rnorm(n_iter, 0, sd_noise)
      }
    }
    for (i in seq_len(no)) {
      M[, psi_nm(i)] <- truePsi[i] + abs(stats::rnorm(n_iter, 0, sd_noise))
    }
    coda::mcmc(M)
  }
  draws <- coda::as.mcmc.list(list(mk_chain(1L), mk_chain(2L)))
  fit <- list(
    greta = list(draws = draws),
    extras = list(parse_info = list(random = list(
      list(type = "fa_gxe", outer = "env", inner = "geno",
           k = k, n_outer = no, n_inner = 10L)
    )))
  )
  class(fit) <- "flexybayes"
  fit
}

test_that("fb_structured_cov reconstructs the identified fa covariance", {
  trueL <- matrix(c(0.8, 0.2, -0.5, 0.1, 0.6, 0.3), nrow = 3, ncol = 2)
  truePsi <- c(0.20, 0.30, 0.15)
  Gtrue <- trueL %*% t(trueL) + diag(truePsi)
  fit <- mk_fa_fit(trueL, truePsi)

  sc <- fb_structured_cov(fit)
  expect_named(sc, "env")
  expect_identical(dim(sc$env$cov_mean), c(3L, 3L))
  expect_equal(unname(sc$env$cov_mean), Gtrue, tolerance = 0.03)
  expect_true(isSymmetric(unname(sc$env$cov_mean)))
  # diagonal of the correlation matrix is 1
  expect_equal(unname(diag(sc$env$cor_mean)), rep(1, 3), tolerance = 1e-8)
})

test_that("fb_structured_cov gives a meaningful (low) Rhat on the identified quantity", {
  trueL <- matrix(c(0.8, 0.2, -0.5, 0.1, 0.6, 0.3), nrow = 3, ncol = 2)
  fit <- mk_fa_fit(trueL, c(0.2, 0.3, 0.15))
  sc <- fb_structured_cov(fit)
  expect_true(is.finite(sc$env$max_rhat))
  expect_lt(sc$env$max_rhat, 1.1)
})

test_that("fb_structured_cov messages and returns empty without an fa term", {
  fit <- list(
    greta = list(draws = NULL),
    extras = list(parse_info = list(random = list()))
  )
  class(fit) <- "flexybayes"
  expect_message(out <- fb_structured_cov(fit), "no factor-analytic")
  expect_length(out, 0L)
})

test_that("fb_structured_cov rejects non-flexybayes input", {
  expect_error(fb_structured_cov(list(a = 1)), "must be a flexybayes")
})
