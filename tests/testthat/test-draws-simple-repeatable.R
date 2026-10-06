# DEF-007: fb_as_draws_simple() on an INLA fit repeats under set.seed() and
# draws from the same distribution as INLA::inla.posterior.sample().

.dsr_fit <- function() {
  set.seed(3)
  d <- data.frame(
    site = factor(rep(1:6, each = 12)),
    trt = factor(rep(c("a", "b", "c"), 24)),
    x = stats::rnorm(72)
  )
  d$y <- 1 + c(0, 0.5, -0.3)[d$trt] + 0.4 * d$x +
    stats::rnorm(6)[d$site] + stats::rnorm(72)
  suppressWarnings(suppressMessages(flexybayes(
    y ~ trt + x, random = ~site, data = d, backend = "inla", verbose = FALSE
  )))
}

.dsr_poisson_fit <- function() {
  set.seed(20260916L)
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:30), each = 8L)))
  d$x <- stats::rnorm(nrow(d), 0, 1000)
  u <- stats::rnorm(30L, 0, 0.7)
  b <- stats::rnorm(30L, 0, 2e-4)
  eta <- 0.4 + (5e-4 + b[as.integer(d$g)]) * d$x + u[as.integer(d$g)]
  d$y <- stats::rpois(nrow(d), exp(eta))
  suppressWarnings(suppressMessages(flexybayes(
    y ~ x + (x || g), data = d, family = "poisson", backend = "inla",
    verbose = FALSE
  )))
}

# Draws as a matrix with one column per named parameter.
.dsr_matrix <- function(draws, rows) {
  vapply(draws[rows], as.numeric, numeric(length(draws[[rows[1L]]])))
}

.dsr_inla_matrix <- function(fit, n, rows) {
  s <- .fb_inla_posterior_sample(fit, n, num.threads = "1:1")
  t(vapply(s, function(z) z$latent[rows, 1L], numeric(length(rows))))
}

test_that("DEF-007: two calls under one seed return the same draws", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  fit <- .dsr_fit()
  set.seed(1)
  a <- fb_as_draws_simple(fit, n_samples = 200L)
  set.seed(1)
  b <- fb_as_draws_simple(fit, n_samples = 200L)
  expect_identical(a, b)
})

test_that("DEF-007: draws keep the names and scale of INLA's sampler", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  fit <- .dsr_fit()
  set.seed(1)
  ours <- fb_as_draws_simple(fit, n_samples = 20000L)
  s <- .fb_inla_posterior_sample(fit, 20000L, num.threads = "1:1")
  expect_identical(
    names(ours),
    c(rownames(s[[1L]]$latent), names(s[[1L]]$hyperpar))
  )
  expect_true(all(lengths(ours) == 20000L))

  # Each hyperparameter draw is one of INLA's configurations, mapped alike.
  for (nm in names(s[[1L]]$hyperpar)) {
    theirs <- vapply(s, function(z) z$hyperpar[[nm]], numeric(1L))
    expect_true(all(unique(ours[[nm]]) %in% unique(theirs)), label = nm)
  }
})

test_that("DEF-007: draws agree with inla.posterior.sample() and INLA's summary", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  fit <- .dsr_fit()
  fixed <- paste0(rownames(fit$inla$summary.fixed), ":1")
  rows <- c(fixed, "site:1", "site:4", "Predictor:1", "Predictor:40")
  set.seed(1)
  ours <- .dsr_matrix(fb_as_draws_simple(fit, n_samples = 20000L), rows)
  set.seed(2)
  theirs <- .dsr_inla_matrix(fit, 20000L, rows)

  sd_ref <- sqrt(diag(stats::cov(theirs)))
  # Over 10 seed pairs: mean gap up to 0.022 SD, SD ratio within 0.027 of 1,
  # correlation gap up to 0.028; two INLA runs: 0.029, 0.027, 0.026.
  expect_lt(max(abs(colMeans(ours) - colMeans(theirs)) / sd_ref), 0.04)
  expect_lt(max(abs(sqrt(diag(stats::cov(ours))) / sd_ref - 1)), 0.04)
  expect_lt(max(abs(stats::cor(ours) - stats::cor(theirs))), 0.04)

  sf <- fit$inla$summary.fixed
  # Over 10 seeds: mean gap up to 0.016 SD; SD ratio within 0.044 of 1 (0.036
  # for INLA's own draws), as summary SDs come from INLA's marginal densities.
  expect_lt(max(abs(colMeans(ours[, fixed]) - sf$mean) / sf$sd), 0.03)
  expect_lt(max(abs(apply(ours[, fixed], 2L, stats::sd) / sf$sd - 1)), 0.06)
  sr <- fit$inla$summary.random$site[c(1L, 4L), ]
  site <- ours[, c("site:1", "site:4")]
  # Over 10 seeds: mean gap up to 0.012 SD, SD ratio within 0.039 of 1.
  expect_lt(max(abs(colMeans(site) - sr$mean) / sr$sd), 0.03)
  expect_lt(max(abs(apply(site, 2L, stats::sd) / sr$sd - 1)), 0.06)
})

test_that("DEF-007: draws on a covariate-standardised slope are on the data scale", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  fit <- .dsr_poisson_fit()
  expect_false(is.null(fit$extras$covariate_standardisation))
  rows <- c("(Intercept):1", "x:1", "g_for_slope_x:1", "g_for_slope_x:7")
  set.seed(1)
  ours <- .dsr_matrix(fb_as_draws_simple(fit, n_samples = 20000L), rows)

  ref <- rbind(
    fit$inla$summary.fixed[, c("mean", "sd")],
    fit$inla$summary.random$g_for_slope_x[c(1L, 7L), c("mean", "sd")]
  )
  # Over 10 seeds: mean gap up to 0.018 SD, SD ratio within 0.015 of 1.
  expect_lt(max(abs(colMeans(ours) - ref$mean) / ref$sd), 0.03)
  expect_lt(max(abs(apply(ours, 2L, stats::sd) / ref$sd - 1)), 0.03)
})

test_that("DEF-007: draws on a sum-to-zero random walk keep the constraint", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  set.seed(5)
  d <- data.frame(x = rep(1:25, 3))
  d$y <- sin(d$x / 4) + stats::rnorm(75, 0, 0.3)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ 1, random = ~ spl(x), data = d, backend = "inla", verbose = FALSE
  )))
  expect_false(is.null(fit$inla$misc$configs$constr))
  set.seed(1)
  draws <- fb_as_draws_simple(fit, n_samples = 20000L)
  walk <- .dsr_matrix(draws, paste0("x:", 1:25))
  # Over 10 seeds the largest per-draw sum was 9.8e-13.
  expect_lt(max(abs(rowSums(walk))), 1e-10)

  rows <- c(1L, 5L, 13L, 25L)
  ref <- fit$inla$summary.random$x[rows, ]
  # Over 10 seeds: mean gap up to 0.018 SD, SD ratio within 0.014 of 1.
  expect_lt(max(abs(colMeans(walk[, rows]) - ref$mean) / ref$sd), 0.03)
  expect_lt(max(abs(apply(walk[, rows], 2L, stats::sd) / ref$sd - 1)), 0.03)
})
