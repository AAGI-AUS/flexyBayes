# =============================================================================
# fb_marginal_cov(): the marginal covariance of the response of a fit.
#
# Oracles, none written by the code under test.
#
#   Dense assembly. as.matrix() of the object equals V written out from the
#   fit's own variance-component table and the data, element by element:
#   E[sigma^2] = estimate^2 + std.error^2 for every variance, a correlation
#   at its estimate. Fixtures: iid intercept, nested intercepts, one
#   AR1 x AR1 field per trial, a residual sectioned by a factor, a known
#   variance per row and a known covariance matrix.
#
#   Kronecker factors. A separable block's factors, multiplied out, equal
#   that block of V.
#
#   Rows. A subset of rows gives the submatrix of the full V, and a cross
#   block (rows, cols) the matching rectangle.
#
# Tolerances, fixed before any result.
#   - V against its written-out definition: 1e-12 relative. Both sides add
#     the same few terms per entry, possibly in a different order.
#   - Kronecker factors against V: 1e-12 absolute, from subtracting the
#     other terms of V.
#   - Subsets: identical(). The same entries pass through the same
#     arithmetic.
# =============================================================================

fbmc_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE,
    flexyBayes.silence_high_missingness_warning = TRUE,
    .local_envir = envir
  )
  old <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  withr::defer(INLA::inla.setOption(num.threads = old), envir = envir)
}

fbmc_rel <- function(a, b) {
  max(abs(a - b)) / max(abs(b))
}

# E[sigma^2] for each row of a fit's variance-component table, and each
# correlation at its estimate.
fbmc_values <- function(fit) {
  vc <- summary(fit)$varcomp
  is_rho <- startsWith(vc$component, "rho_")
  list(
    var = stats::setNames(
      vc$estimate[!is_rho]^2 + vc$std.error[!is_rho]^2,
      vc$component[!is_rho]
    ),
    rho = stats::setNames(vc$estimate[is_rho], vc$component[is_rho])
  )
}

fbmc_same <- function(f) {
  outer(as.character(f), as.character(f), "==") * 1
}

fbmc_site_data <- function(seed = 3L) {
  withr::with_seed(seed, {
    d <- data.frame(
      site = factor(rep(sprintf("s%d", 1:6), each = 8L)),
      year = factor(rep(c("y1", "y2", "y3", "y4"), length.out = 48L)),
      env = factor(rep(c("a", "b"), 24L))
    )
    d$x <- stats::rnorm(48L)
    d$y <- 1 + d$x + stats::rnorm(6L)[d$site] +
      stats::rnorm(24L, sd = 0.4)[as.integer(interaction(d$site, d$year))] +
      stats::rnorm(48L, sd = ifelse(d$env == "a", 0.5, 1.5))
    d
  })
}

fbmc_field_data <- function(seed = 4L) {
  withr::with_seed(seed, {
    d <- expand.grid(row = 1:4, range = 1:3, trial = c("t1", "t2", "t3"))
    d <- d[sample.int(nrow(d)), ]
    d$trial <- factor(d$trial)
    d$range <- factor(d$range)
    d$row <- factor(d$row)
    d$x <- stats::rnorm(nrow(d))
    d$y <- 2 + d$x + stats::rnorm(3L)[d$trial] +
      0.5 * as.integer(d$range) - 0.3 * as.integer(d$row) +
      stats::rnorm(nrow(d), sd = 0.5)
    rownames(d) <- NULL
    d
  })
}

fbmc_fit <- function(...) {
  suppressWarnings(flexybayes(..., backend = "inla", verbose = FALSE))
}


# --- (a) dense assembly against V written out ------------------------------ #

test_that("an iid intercept: V = sigma^2 I + sigma_site^2 Z Z'", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  fit <- fbmc_fit(y ~ x, random = ~site, data = d)
  s <- fbmc_values(fit)$var
  want <- diag(s[["sigma"]], nrow(d)) + s[["sd_site"]] * fbmc_same(d$site)
  x <- fb_marginal_cov(fit)
  expect_s3_class(x, "fb_marginal_cov")
  expect_lt(fbmc_rel(as.matrix(x), want), 1e-12)
})

test_that("nested intercepts: V adds sigma_{site:year}^2 for the site:year cell", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  fit <- fbmc_fit(y ~ x, random = ~ site + site:year, data = d)
  s <- fbmc_values(fit)$var
  want <- diag(s[["sigma"]], nrow(d)) +
    s[["sd_site"]] * fbmc_same(d$site) +
    s[["sd_site:year"]] * fbmc_same(paste(d$site, d$year))
  x <- fb_marginal_cov(fit)
  expect_lt(fbmc_rel(as.matrix(x), want), 1e-12)
  # The two intercepts both nest within site, so V is block-diagonal by site.
  expect_true(x$groups$applicable)
  expect_identical(x$groups$column, "site")
  expect_length(x$groups$blocks, 6L)
})

test_that("one AR1 x AR1 field per trial: rho_row^|dr| rho_col^|dc| within a trial, zero across", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_field_data()
  fit <- fbmc_fit(y ~ x, random = ~ trial + id(trial):ar1(range):ar1(row),
                  data = d)
  v <- fbmc_values(fit)
  s <- v$var
  rho <- v$rho
  same <- fbmc_same(d$trial)
  want <- diag(s[["sigma"]], nrow(d)) +
    s[["sd_trial"]] * same +
    s[["sd_spatial"]] * same *
      rho[["rho_row"]]^abs(outer(as.integer(d$range), as.integer(d$range), "-")) *
      rho[["rho_col"]]^abs(outer(as.integer(d$row), as.integer(d$row), "-"))
  x <- fb_marginal_cov(fit)
  got <- as.matrix(x)
  expect_lt(fbmc_rel(got, want), 1e-12)
  expect_true(all(got[same == 0] == 0))
})

test_that("a residual sectioned by a factor: each row's residual variance is its section's", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  fit <- fbmc_fit(y ~ x, random = ~site, residual = ~ dsum(~ units | env),
                  data = d)
  s <- fbmc_values(fit)$var
  want <- diag(s[paste0("sigma_", as.character(d$env))]) +
    s[["sd_site"]] * fbmc_same(d$site)
  x <- fb_marginal_cov(fit)
  expect_identical(x$residual$by, "env")
  expect_lt(fbmc_rel(as.matrix(x), want), 1e-12)
})

test_that("a known variance per row adds v_i to the residual on the diagonal", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  known <- withr::with_seed(9L, stats::runif(nrow(d), 0.1, 0.6))
  fit <- fbmc_fit(y ~ x, random = ~site, data = d, known_variance = known)
  s <- fbmc_values(fit)$var
  want <- diag(s[["sigma"]] + known) + s[["sd_site"]] * fbmc_same(d$site)
  x <- fb_marginal_cov(fit)
  expect_equal(x$residual$known_variance, known)
  expect_lt(fbmc_rel(as.matrix(x), want), 1e-12)
})

test_that("a known covariance matrix K is added whole", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  k <- withr::with_seed(11L, {
    a <- matrix(stats::rnorm(nrow(d)^2, 0, 0.05), nrow(d))
    diag(stats::runif(nrow(d), 0.05, 0.3)) + tcrossprod(a)
  })
  fit <- fbmc_fit(y ~ x, random = ~site, data = d, known_variance = k)
  s <- fbmc_values(fit)$var
  want <- diag(s[["sigma"]], nrow(d)) + k + s[["sd_site"]] * fbmc_same(d$site)
  x <- fb_marginal_cov(fit)
  expect_lt(fbmc_rel(as.matrix(x), want), 1e-12)
  # K links rows across every site, so V has no block-diagonal grouping.
  expect_false(x$groups$applicable)
})


# --- (b) Kronecker factors -------------------------------------------------- #

test_that("a per-trial field's Kronecker factors, multiplied out, equal its block of V", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_field_data()
  fit <- fbmc_fit(y ~ x, random = ~ trial + id(trial):ar1(range):ar1(row),
                  data = d)
  x <- fb_marginal_cov(fit)
  field <- x$blocks[["spatial"]]
  f <- field$kronecker
  expect_identical(f$n_at, 3L)
  expect_identical(dim(f$row), c(3L, 3L))
  expect_identical(dim(f$col), c(4L, 4L))
  g <- kronecker(diag(f$n_at), kronecker(f$row, f$col))
  from_factors <- field$variance * g[field$index, field$index]

  # V less the residual and the trial intercept is the field's block.
  others <- diag(x$residual$variance, nrow(d)) +
    x$blocks[["trial"]]$variance * fbmc_same(d$trial)
  expect_lt(max(abs(as.matrix(x) - others - from_factors)), 1e-12)
})


# --- (c) rows --------------------------------------------------------------- #

test_that("a subset of rows gives the submatrix, and rows by cols the rectangle", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_field_data()
  fit <- fbmc_fit(y ~ x, random = ~ trial + id(trial):ar1(range):ar1(row),
                  data = d)
  full <- as.matrix(fb_marginal_cov(fit))
  rows <- c(30L, 2L, 17L, 5L, 11L)
  sub <- fb_marginal_cov(fit, rows = rows)
  expect_identical(sub$n, length(rows))
  expect_identical(sub$rows, rows)
  expect_identical(as.matrix(sub), full[rows, rows])
  x <- fb_marginal_cov(fit)
  a <- 1:10
  b <- 6:20
  expect_identical(as.matrix(x, rows = a, cols = b), full[a, b])
})


# --- (d) defaults and refusals ---------------------------------------------- #

test_that("defaults cover every row of the fit's data, and the object prints", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  fit <- fbmc_fit(y ~ x, random = ~site, data = d)
  x <- fb_marginal_cov(fit)
  expect_identical(x$n, nrow(d))
  expect_identical(x$rows, seq_len(nrow(d)))
  expect_identical(dim(as.matrix(x)), c(nrow(d), nrow(d)))
  expect_identical(names(x$blocks), "site")
  expect_output(print(x), "fb_marginal_cov")
})

test_that("variances = FALSE gives the structure alone, and as.matrix() refuses it", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  fit <- fbmc_fit(y ~ x, random = ~site, data = d)
  x <- fb_marginal_cov(fit, variances = FALSE)
  expect_true(is.na(x$blocks[["site"]]$variance))
  expect_error(as.matrix(x), class = "flexybayes_refusal_marginal_cov_not_buildable")
  # The structure re-evaluated at given variances.
  y <- fb_marginal_cov(x, variances = list(residual = 0.5, terms = c(site = 2)))
  want <- diag(0.5, nrow(d)) + 2 * fbmc_same(d$site)
  expect_identical(as.matrix(y), want)
})

test_that("a structure it does not build refuses by name", {
  skip_if_not_installed("INLA")
  fbmc_quiet()
  d <- fbmc_site_data()
  fit <- fbmc_fit(y ~ x, random = ~ site:year:env, data = d)
  expect_error(
    fb_marginal_cov(fit),
    class = "flexybayes_refusal_marginal_cov_not_buildable"
  )
})
