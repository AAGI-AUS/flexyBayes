# DEF-001: SD medians and 95% limits of INLA variance components agree with
# INLA's own summary and with draws from inla.hyperpar.sample().

.vq_fit <- function() {
  set.seed(6)
  d <- expand.grid(env = factor(1:6), gen = factor(1:20), blk = factor(1:4))
  d$eb <- interaction(d$env, d$blk)
  ge <- stats::rnorm(20)
  # eb has true SD 0, so its precision marginal has a long right tail.
  d$y <- 2 + ge[d$gen] + stats::rnorm(nrow(d))
  suppressWarnings(suppressMessages(flexybayes(
    y ~ env,
    random = ~ gen + eb,
    data = d,
    backend = "inla",
    verbose = FALSE
  )))
}

test_that("DEF-001: SD quantiles match INLA's summary on a near-zero component", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vq_fit()
  vc <- fit$extras$variance_comps
  hp <- fit$inla$summary.hyperpar

  ref <- data.frame(
    component = c("sigma", "sd_gen", "sd_eb"),
    q50 = 1 / sqrt(hp[, "0.5quant"]),
    q2.5 = 1 / sqrt(hp[, "0.975quant"]),
    q97.5 = 1 / sqrt(hp[, "0.025quant"])
  )
  vc <- vc[match(ref$component, vc$component), ]
  med <- attr(fit$extras$variance_comps, "posterior_median")[ref$component]

  # INLA computes its summary on its own grid. On the Barley base model it
  # differs from 20,000 draws by up to 2.8%.
  expect_lt(max(abs(med / ref$q50 - 1)), 0.03)
  expect_lt(max(abs(vc$q2.5 / ref$q2.5 - 1)), 0.03)
  expect_lt(max(abs(vc$q97.5 / ref$q97.5 - 1)), 0.03)
})

test_that("DEF-001: SD quantiles match inla.hyperpar.sample() draws", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vq_fit()
  vc <- fit$extras$variance_comps
  set.seed(1)
  draws <- 1 /
    sqrt(INLA::inla.hyperpar.sample(
      20000,
      fit$inla,
      improve.marginals = TRUE
    ))
  q <- apply(draws, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
  comp <- c("sigma", "sd_gen", "sd_eb")
  vc <- vc[match(comp, vc$component), ]
  med <- attr(fit$extras$variance_comps, "posterior_median")[comp]

  # Each ratio is checked on its own: a vector tolerance would let the large
  # components mask an error on sd_eb.
  expect_lt(max(abs(med / q[2, ] - 1)), 0.03)
  expect_lt(max(abs(vc$q2.5 / q[1, ] - 1)), 0.03)
  expect_lt(max(abs(vc$q97.5 / q[3, ] - 1)), 0.03)
})

# DEF-002: hyperparameter draws reproduce the reported medians and 95% limits
# on every component. Improved draws match the medians to about 1e-9 and the
# limits to under 1e-3 at 20,000 draws, and default draws missed by up to 1%.

.vq_draw_error <- function(fit, draws) {
  vc <- fit$extras$variance_comps
  med <- attr(vc, "posterior_median")
  comp <- intersect(names(med), names(draws))
  vc <- vc[match(comp, vc$component), ]
  q <- vapply(
    comp,
    function(k) stats::quantile(draws[[k]], c(0.025, 0.5, 0.975)),
    numeric(3L)
  )
  list(
    med = abs(q[2L, ] / med[comp] - 1),
    lim = abs(c(q[1L, ] / vc$q2.5, q[3L, ] / vc$q97.5) - 1)
  )
}

test_that("DEF-002: fb_variance_draws() reproduces the reported quantiles", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vq_fit()
  err <- .vq_draw_error(fit, fb_variance_draws(fit, 20000L, seed = 1L))
  expect_lt(max(err$med), 1e-4)
  expect_lt(max(err$lim), 2e-3)
})

test_that("DEF-002: draws on a covariate-standardised slope are on the data scale", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(20260916L)
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:30), each = 8L)))
  d$x <- stats::rnorm(nrow(d), 0, 1000)
  u <- stats::rnorm(30L, 0, 0.7)
  b <- stats::rnorm(30L, 0, 2e-4)
  eta <- 0.4 + (5e-4 + b[as.integer(d$g)]) * d$x + u[as.integer(d$g)]
  d$y <- stats::rpois(nrow(d), exp(eta))
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ x + (x || g),
    data = d,
    family = "poisson",
    backend = "inla",
    verbose = FALSE
  )))
  expect_false(is.null(fit$extras$covariate_standardisation))
  err <- .vq_draw_error(fit, fb_variance_draws(fit, 20000L, seed = 1L))
  expect_lt(max(err$med), 1e-4)
  expect_lt(max(err$lim), 2e-3)
})

test_that("DEF-002: vpredict() medians match the reported medians", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vq_fit()
  med <- attr(fit$extras$variance_comps, "posterior_median")
  set.seed(1)
  vp <- vpredict(fit, g ~ sd_gen, e ~ sd_eb, s ~ sigma, n_draws = 20000L)
  got <- vp$median[match(c("g", "e", "s"), vp$functional)]
  expect_lt(max(abs(got / med[c("sd_gen", "sd_eb", "sigma")] - 1)), 1e-4)
})
