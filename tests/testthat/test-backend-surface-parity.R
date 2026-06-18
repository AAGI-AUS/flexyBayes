# Regression guards for the 2026-05-31 audit fixes. Each test pins a
# bug that the audit found and fixed, so a future change that
# re-introduces it fails here rather than shipping silently.
#
#   1. triangulate() must compare variance components on a common (SD)
#      scale across greta and INLA -- the INLA precision -> SD transform
#      must fire (was silently no-op: precision compared against SD).
#   2. plot() must dispatch on every backend fit class -- "effects"
#      works everywhere; MCMC-only types degrade to a message on a
#      backend that has no draws (was: plot.default crash on INLA).
#   3. Positional and named prior arguments must parse identically (was:
#      positional scale silently dropped to the default on by-name emit).
#   4. A non-binary binomial response must be refused on greta, not
#      silently fitted as Bernoulli.

silence_notes <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}


# ---------------------------------------------------------------- #
# 1. Cross-engine variance-component commensurability               #
# ---------------------------------------------------------------- #

test_that("triangulate() puts INLA variance components on the SD scale, not precision", {
  skip_on_cran()
  skip_if_greta_backend_unusable()
  skip_if_not_installed("INLA")
  silence_notes()

  # True sigma_e = 2.0 and sd_g = 1.5 are both clearly > 1, so the
  # SD scale (~2, ~1.5) is unmistakably distinct from the precision
  # scale (1/4 = 0.25, 1/2.25 = 0.44). A precision/SD mismatch would
  # collapse the INLA rows toward those small values.
  set.seed(101L)
  J <- 8L
  N <- 160L
  g <- factor(sample.int(J, N, replace = TRUE))
  x <- stats::rnorm(N)
  u <- stats::rnorm(J, 0, 1.5)
  y <- 2 + 0.8 * x + u[as.integer(g)] + stats::rnorm(N, 0, 2.0)
  d <- data.frame(y, x, g)

  fg <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    backend = "greta",
    n_samples = 400L,
    warmup = 400L,
    chains = 2L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fi <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  tri <- triangulate(fg, fi, n_samples = 1000L)

  expect_true(all(c("sigma", "sd_g") %in% tri$common))
  m <- tri$metrics

  sig <- m[m$param == "sigma", ]
  sdg <- m[m$param == "sd_g", ]

  # INLA side (mean_b) must be on the SD scale: residual SD near 2,
  # group SD near 1.5 -- well above the precision-scale values they
  # would collapse to if the transform failed.
  expect_gt(sig$mean_b, 1.0)
  expect_gt(sdg$mean_b, 0.8)

  # And the two engines must agree to a sane tolerance on the SD scale.
  expect_lt(abs(sig$mean_a - sig$mean_b) / sig$mean_a, 0.5)
  expect_lt(abs(sdg$mean_a - sdg$mean_b) / max(sdg$mean_a, 0.1), 0.8)
})


# ---------------------------------------------------------------- #
# 2. plot() dispatches across backend fit classes                   #
# ---------------------------------------------------------------- #

test_that("plot() dispatches on an INLA fit: effects + residuals render, MCMC-only types degrade to a message", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  silence_notes()

  set.seed(202L)
  d <- data.frame(
    y = stats::rnorm(60L),
    x = stats::rnorm(60L),
    g = factor(rep(seq_len(6L), length.out = 60L))
  )
  fi <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  expect_s3_class(fi, "flexybayes_inla")
  expect_false(inherits(fi, "flexybayes"))

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # The forest plot reads coef() / confint(), which INLA supports.
  expect_no_error(plot(fi, type = "effects"))

  # The residual plot reads fitted() / residuals(), which the INLA fit now
  # exposes (fitted.flexybayes_inla / residuals.flexybayes_inla, added
  # 2026-06-02 after the corner-to-corner stress run found they silently
  # returned NULL) -- so it renders rather than degrading.
  expect_no_error(plot(fi, type = "residuals"))

  # Trace / diagnostics plots genuinely need MCMC draws the INLA (Laplace)
  # fit does not have -- they must message, never crash to plot.default.
  expect_message(plot(fi, type = "diagnostics"), "not available")
})


# ---------------------------------------------------------------- #
# 3. Positional vs named prior arguments                            #
# ---------------------------------------------------------------- #

test_that("positional and named prior arguments parse to identical specs", {
  pn <- fb_prior(b("x") ~ normal(0, sd = 50))
  pp <- fb_prior(b("x") ~ normal(0, 50))
  expect_identical(pn$specs[[1]]$spec$args, pp$specs[[1]]$spec$args)
  expect_equal(pp$specs[[1]]$spec$args$sd, 50)

  qn <- fb_prior(sd(group = "g") ~ pc(upper = 2, prob = 0.05))
  qp <- fb_prior(sd(group = "g") ~ pc(2, 0.05))
  expect_identical(qn$specs[[1]]$spec$args, qp$specs[[1]]$spec$args)
  expect_equal(qp$specs[[1]]$spec$args$upper, 2)

  hn <- fb_prior(sd(group = "g") ~ half_normal(scale = 3))
  hp <- fb_prior(sd(group = "g") ~ half_normal(3))
  expect_identical(hn$specs[[1]]$spec$args, hp$specs[[1]]$spec$args)
})


# ---------------------------------------------------------------- #
# 5. predict() / emmeans / marginaleffects on an INLA mixed model   #
# ---------------------------------------------------------------- #

test_that("predict() works on an INLA mixed model (brms-grammar bar term stripped)", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  silence_notes()

  set.seed(404L)
  d <- data.frame(
    y = stats::rnorm(80L),
    x = stats::rnorm(80L),
    g = factor(rep(seq_len(8L), each = 10L))
  )

  # brms-grammar mixed model: the random-effect bar `(1 | g)` is stored
  # whole in call_info$fixed; formula() must strip it so the fixed-effect
  # design-matrix reconstruction matches the coefficient basis.
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d,
    backend = "inla",
    verbose = FALSE
  ))

  expect_identical(
    deparse(stats::formula(fit)),
    "y ~ x"
  ) # bar stripped

  p <- predict(fit)
  expect_length(p, nrow(d))
  expect_true(all(is.finite(p)))

  ps <- predict(fit, se.fit = TRUE)
  expect_length(ps$fit, nrow(d))
  expect_length(ps$se.fit, nrow(d))

  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    expect_identical(
      nrow(marginaleffects::predictions(fit)),
      nrow(d)
    )
  }
})


# ---------------------------------------------------------------- #
# 4. Non-binary binomial is refused on greta                        #
# ---------------------------------------------------------------- #

test_that("a non-binary binomial response is refused on greta, not silently Bernoulli-fit", {
  skip_on_cran()
  skip_if_greta_backend_unusable()
  silence_notes()

  set.seed(303L)
  d <- data.frame(
    y = sample(0:5, 50L, replace = TRUE), # counts, not 0/1
    x = stats::rnorm(50L)
  )

  expect_error(
    suppressMessages(fb(
      y ~ x,
      data = d,
      family = "binomial",
      backend = "greta",
      n_samples = 10L,
      warmup = 10L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    "binary"
  )
})
