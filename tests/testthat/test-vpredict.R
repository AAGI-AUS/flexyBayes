# Posterior functionals of the variance components: vpredict().
#
# ASReml's vpredict() returns a heritability with a DELTA-METHOD standard
# error, which linearises at the REML estimate and therefore does not know
# where the parameter space ends. This one transforms the draws, so the
# functional's interval is a quantile interval of its own posterior and
# cannot leave the range the function maps into.
#
# The tests hold three things: the number agrees with an oracle this
# package did not write (lme4's REML fit, through the same plug-in
# formula); the boundary property holds where the delta method's does not;
# and a component name the fit does not carry refuses instead of resolving
# against whatever the caller has in their workspace.

.vp_fixture <- function(
  seed = 11L,
  n_g = 40L,
  n_rep = 8L,
  sd_g = 1.0,
  sd_resid = 0.6
) {
  set.seed(seed)
  tab <- expand.grid(rep = seq_len(n_rep), gen = factor(seq_len(n_g)))
  u <- rnorm(n_g, 0, sd_g)
  tab$y <- 10 + u[as.integer(tab$gen)] + rnorm(nrow(tab), 0, sd_resid)
  tab
}

.vp_fit <- function(tab, backend = "inla") {
  suppressMessages(flexybayes(
    y ~ 1,
    random = ~gen,
    data = tab,
    backend = backend,
    verbose = FALSE
  ))
}

test_that("a heritability agrees with lme4's REML fit of the same model", {
  # The oracle is lme4, not this package: the REML variance components put
  # through the same plug-in ratio. A posterior mean of a ratio is not the
  # ratio of REML point estimates, so the tolerance is the width of that
  # difference and not machine precision, but a mistranslated component,
  # or a functional evaluated on the precision scale, would miss by far
  # more than this.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .vp_fixture()
  fit <- .vp_fit(tab)
  res <- vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2), n_draws = 4000L)

  ref <- lme4::lmer(y ~ 1 + (1 | gen), data = tab)
  vc <- as.data.frame(lme4::VarCorr(ref))
  ref_h2 <- vc$vcov[vc$grp == "gen"] / sum(vc$vcov)

  expect_identical(res$functional, "h2")
  expect_equal(res$estimate, ref_h2, tolerance = 0.05)
  expect_lt(res$conf.low, ref_h2)
  expect_gt(res$conf.high, ref_h2)
})

test_that("the interval respects the boundary the delta method does not", {
  # A heritability lives in [0, 1] and a variance ratio in [0, Inf). Every
  # DRAW does, so the interval does. The delta method has no such
  # guarantee: it is a normal interval around a point estimate, and near a
  # boundary it crosses one.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  # A genotype variance close to zero, the case that makes a delta-method
  # heritability interval run below zero.
  tab <- .vp_fixture(seed = 5L, sd_g = 0.05, sd_resid = 1.0)
  fit <- .vp_fit(tab)
  res <- vpredict(
    fit,
    h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2),
    ratio ~ sd_gen / sigma,
    n_draws = 4000L
  )
  h2 <- res[res$functional == "h2", ]
  expect_gte(h2$conf.low, 0)
  expect_lte(h2$conf.high, 1)
  expect_gte(res$conf.low[res$functional == "ratio"], 0)

  # And the draws themselves, not only their summary.
  d <- vpredict(
    fit,
    h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2),
    n_draws = 1000L,
    draws = TRUE
  )
  expect_true(all(d[, "h2"] >= 0 & d[, "h2"] <= 1))
})

test_that("the functional's draws are continuous, not INLA's integration grid", {
  # The generic draws path takes hyperparameters from
  # inla.posterior.sample(), which resamples the stored integration
  # configurations: 1000 draws carried 15 distinct values on this fixture,
  # so a functional built on them had a posterior supported on 15 atoms
  # while its quantiles were reported as a smooth interval. The point
  # estimate was unaffected; the tails were not. vpredict() draws from the
  # joint hyperparameter posterior instead.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vp_fit(.vp_fixture())
  d <- vpredict(
    fit,
    h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2),
    n_draws = 1000L,
    draws = TRUE
  )
  expect_gt(length(unique(d[, "h2"])), 500L)
  # And the support is genuinely wider than the central interval, which a
  # 15-atom distribution cannot be: its min IS its 2.5th percentile.
  q <- stats::quantile(d[, "h2"], c(0.025, 0.975), names = FALSE)
  expect_lt(min(d[, "h2"]), q[[1L]])
  expect_gt(max(d[, "h2"]), q[[2L]])
})

test_that("draws = TRUE returns one column per functional", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vp_fit(.vp_fixture(n_g = 20L, n_rep = 4L))
  d <- vpredict(
    fit,
    h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2),
    ratio ~ sd_gen / sigma,
    n_draws = 300L,
    draws = TRUE
  )
  expect_true(is.matrix(d))
  expect_identical(colnames(d), c("h2", "ratio"))
  expect_identical(nrow(d), 300L)
})

test_that("a component the fit does not carry refuses and lists what it has", {
  # Evaluated in baseenv(), so a name that happens to exist in the caller's
  # workspace cannot supply a value. Silently returning a number computed
  # from something other than this posterior is the failure mode worth
  # refusing over.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vp_fit(.vp_fixture(n_g = 20L, n_rep = 4L))

  sd_genotype <- 999
  err <- expect_error(
    vpredict(fit, h2 ~ sd_genotype^2 / (sd_genotype^2 + sigma^2)),
    "does not carry"
  )
  expect_match(conditionMessage(err), "sd_gen", fixed = TRUE)
  expect_match(conditionMessage(err), "sigma", fixed = TRUE)
})

test_that("vpredict refuses an argument that is not a two-sided formula", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vp_fit(.vp_fixture(n_g = 20L, n_rep = 4L))
  expect_error(vpredict(fit), "at least one functional")
  expect_error(vpredict(fit, ~ sd_gen / sigma), "two-sided")
  expect_error(
    vpredict(fit, h2 ~ sd_gen / sigma, prob = 1.5),
    "between 0 and 1"
  )
})

test_that("a functional that is not one number per draw refuses", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- .vp_fit(.vp_fixture(n_g = 20L, n_rep = 4L))
  expect_error(
    vpredict(fit, h2 ~ mean(sd_gen), n_draws = 200L),
    "one number per posterior draw"
  )
})

test_that("the same functional reads the same on either engine", {
  # The point of writing this over CANONICAL draws: a heritability is one
  # question, and the answer should not depend on which engine was asked.
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  tab <- .vp_fixture(n_g = 30L, n_rep = 6L)
  h2 <- function(fit, n) {
    vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2), n_draws = n)$estimate
  }
  a <- h2(.vp_fit(tab, "inla"), 4000L)
  b <- h2(
    suppressMessages(flexybayes(
      y ~ 1,
      random = ~gen,
      data = tab,
      backend = "brms",
      chains = 2L,
      n_samples = 1000L,
      warmup = 1000L,
      seed = 1L,
      verbose = FALSE
    )),
    2000L
  )
  expect_equal(a, b, tolerance = 0.06)
})
