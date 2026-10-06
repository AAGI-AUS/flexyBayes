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
    y ~ env, random = ~ gen + eb, data = d, backend = "inla", verbose = FALSE
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

  # INLA computes its summary on its own grid; on the Barley base model it
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
  draws <- 1 / sqrt(INLA::inla.hyperpar.sample(
    20000, fit$inla, improve.marginals = TRUE
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
