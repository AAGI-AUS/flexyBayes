# =============================================================================
# Observation weights (C6): fitting side.
#
# `weights` is lowered for family = "gaussian" with the identity link only,
# in the ASReml / lme4 / glm(weights =) precision sense, Var(y_i) =
# sigma^2 / w_i:
#
#   - INLA: `scale = w`, INLA's own per-observation precision multiplier
#     (?INLA::inla, "scale": "Fixed (optional) scale parameters of the
#     precision ... default rep(1, n.data)").
#   - brms: NOT brms's own `weights()` addition term. brms's own
#     documentation (brmsformula.Rd, "Additional response information")
#     states that mechanism "multiplying the log-posterior values of each
#     observation by their corresponding weights" -- a LIKELIHOOD-POWER
#     weighting, algebraically a different quantity from precision
#     weighting for a Gaussian, confirmed via brms::make_stancode()'s
#     generated Stan (`target += weights[n] * normal_lpdf(...)`) and
#     empirically: fit against lme4::lmer(weights =) on identical
#     simulated data, that route's sigma diverged from lme4's by about
#     25%, stable across repeat fits -- not Monte Carlo noise. Instead,
#     flexyBayes lowers weights on brms's sigma DISTRIBUTIONAL parameter
#     (log link) as a known offset, log(sigma_i) = log(sigma_base) -
#     0.5 * log(w_i) -- see R/emit_brms.R's .FB_BRMS_WEIGHTS_OFFSET_COL
#     for the construction and R/priors_to_brms.R's
#     .brms_retarget_sigma_for_heterogeneous_residual() for why the
#     default sigma prior has to move with it. That route reproduces
#     lme4's sigma to within Monte Carlo error (grounded below).
#
# The metamorphic direction stated in the originating work package spec
# ("multiplying every weight by 2 ... halves the residual variance") does
# not hold under this -- or any -- precision-weighting implementation.
# Grounded independently three ways before this file was written:
#
#   1. lm(y ~ x, weights = w) vs lm(y ~ x, weights = 2 * w): fixed
#      effects (coef) identical; sigma^2 ratio (2w model / w model) = 2,
#      not 0.5. Analytically, sigma_hat^2 = sum(w_i * r_i^2) / (n - p);
#      residuals are unchanged by a common weight rescaling (weighted
#      least squares' normal equations are homogeneous of degree 0 in a
#      common weight scale), so the numerator -- and hence sigma_hat^2 --
#      scales LINEARLY with the rescaling constant.
#   2. lme4::lmer(y ~ x + (1|g), weights = w) vs weights = 2 * w) on the
#      same simulated random-intercept data: identical result, sigma^2
#      doubles.
#   3. Full flexyBayes fits on both engines (this file): sigma scales by
#      sqrt(2) between weights = w and weights = 2 * w, i.e. sigma^2
#      doubles, matching (1) and (2) on both INLA (scale =) and brms
#      (sigma offset).
#
# So this file tests the GROUNDED direction -- doubling weights doubles
# the residual variance -- documented here rather than the spec's stated
# one, per the charter's Independent Oracle Principle (an external claim
# about a modelling convention is verified by running it, not assumed)
# and its "never silently ignore an anomaly" rule. See the WP-C report,
# item C6, for the full trail (Stan code excerpts, brms doc quote, the
# lm()/lme4 spikes, and the full-fit numbers this file's tolerances are
# read off).
# =============================================================================

.wl_dat <- function(seed = 42L, n = 240L, n_group = 12L) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_group), each = n %/% n_group))
  x <- stats::rnorm(n)
  u <- stats::rnorm(n_group, 0, 0.6)
  w <- stats::runif(n, 0.5, 3)
  y <- 2 + 0.5 * x + u[as.integer(g)] + stats::rnorm(n, 0, 1) / sqrt(w)
  data.frame(y = y, x = x, g = g, w = w)
}

.wl_sigma <- function(fit) {
  vc <- summary(fit)$varcomp
  vc$estimate[vc$component == "sigma"]
}

test_that("INLA weighted fixed effects and sigma match lme4::lmer(weights=)", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wl_dat()
  lm4 <- lme4::lmer(y ~ x + (1 | g), data = d, weights = w)
  ref_fixef <- lme4::fixef(lm4)
  ref_sigma <- stats::sigma(lm4)

  fit <- flexybayes(
    y ~ x, random = ~g, data = d, weights = d$w, backend = "inla",
    verbose = FALSE
  )
  # INLA's Laplace approximation on a well-identified random-intercept
  # Gaussian model is close to REML on data this size -- 2% covers both
  # the approximation gap and simulation noise; the live figures this
  # file was grounded from ran at 0.2-0.4%.
  expect_equal(unname(coef(fit)["(Intercept)"]), unname(ref_fixef[["(Intercept)"]]),
    tolerance = 0.02)
  expect_equal(unname(coef(fit)["x"]), unname(ref_fixef[["x"]]), tolerance = 0.02)
  expect_equal(.wl_sigma(fit), ref_sigma, tolerance = 0.02)
})

test_that("brms weighted fixed effects and sigma match lme4::lmer(weights=)", {
  skip_if_not_installed("brms")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wl_dat()
  lm4 <- lme4::lmer(y ~ x + (1 | g), data = d, weights = w)
  ref_fixef <- lme4::fixef(lm4)
  ref_sigma <- stats::sigma(lm4)

  fit <- flexybayes(
    y ~ x, random = ~g, data = d, weights = d$w, backend = "brms",
    n_samples = 700L, warmup = 700L, chains = 2L, seed = 11L,
    verbose = FALSE
  )
  # A posterior mean against a REML point estimate, from an HMC fit at a
  # few hundred effective samples -- 3% covers the sampler's own Monte
  # Carlo error on top of the (small, here) posterior-mean/REML gap; the
  # live figures this file was grounded from ran at 0.1-0.3%.
  expect_equal(unname(coef(fit)["(Intercept)"]), unname(ref_fixef[["(Intercept)"]]),
    tolerance = 0.03)
  expect_equal(unname(coef(fit)["x"]), unname(ref_fixef[["x"]]), tolerance = 0.03)
  expect_equal(.wl_sigma(fit), ref_sigma, tolerance = 0.03)
})

test_that("doubling weights on INLA leaves fixed effects unchanged and doubles sigma^2", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wl_dat()
  fit_w <- flexybayes(
    y ~ x, random = ~g, data = d, weights = d$w, backend = "inla",
    verbose = FALSE
  )
  fit_2w <- flexybayes(
    y ~ x, random = ~g, data = d, weights = 2 * d$w, backend = "inla",
    verbose = FALSE
  )
  expect_equal(coef(fit_w), coef(fit_2w), tolerance = 1e-3)
  ratio <- (.wl_sigma(fit_2w))^2 / (.wl_sigma(fit_w))^2
  expect_equal(ratio, 2, tolerance = 0.02)
})

test_that("doubling weights on brms leaves fixed effects unchanged and doubles sigma^2", {
  skip_if_not_installed("brms")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wl_dat()
  # This is a quantitative metamorphic ratio (sigma^2 doubling), so a
  # Bulk-ESS-low warning here is a real small-sample signal about the
  # estimate itself, not something to muffle -- the fix is a bigger
  # sampler budget, verified by re-running: n_samples = 700 / warmup =
  # 700 / chains = 2 raised a Bulk-ESS warning (Min ESS bulk under the
  # 100 * n_chains = 200 threshold); n_samples = 1500 / warmup = 1000 /
  # chains = 2 does not (Min ESS bulk/tail both comfortably above 200
  # across two independent runs at this fixed seed), and the ratio and
  # coefficient-agreement assertions below hold at essentially the same
  # values either way (ratio 1.998, coef diff 0.005).
  fit_w <- flexybayes(
    y ~ x, random = ~g, data = d, weights = d$w, backend = "brms",
    n_samples = 1500L, warmup = 1000L, chains = 2L, seed = 21L,
    verbose = FALSE
  )
  fit_2w <- flexybayes(
    y ~ x, random = ~g, data = d, weights = 2 * d$w, backend = "brms",
    n_samples = 1500L, warmup = 1000L, chains = 2L, seed = 21L,
    verbose = FALSE
  )
  # Fixed effects: two independent HMC runs on the same data, so a few
  # percent of posterior spread is expected sampler noise, not drift --
  # 5% is loose on purpose (see the two engines' agreement instead for a
  # tight check: the INLA test above holds fixed effects to 0.1%).
  expect_equal(coef(fit_w), coef(fit_2w), tolerance = 0.05)
  ratio <- (.wl_sigma(fit_2w))^2 / (.wl_sigma(fit_w))^2
  expect_equal(ratio, 2, tolerance = 0.08)
})

test_that("backend = \"auto\" fits the Gaussian weighted model", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wl_dat(seed = 7L)
  fit <- flexybayes(
    y ~ x, random = ~g, data = d, weights = d$w, backend = "auto",
    verbose = FALSE
  )
  expect_true(inherits(fit, "flexybayes"))
  expect_true(is.finite(.wl_sigma(fit)))
})

test_that("brms's own weights() addition term is NOT what gets emitted", {
  # Regression guard for the mechanism itself: a naive `y | weights(w)`
  # addition term would put "weights" in the generated Stan code as an
  # addition-term data block. flexyBayes lowers weights as a sigma
  # offset instead (see this file's banner) -- the generated code must
  # show the offset, not brms's own weighted-likelihood construct.
  skip_if_not_installed("brms")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wl_dat(seed = 3L, n = 60L, n_group = 6L)
  code <- flexybayes(
    y ~ x, random = ~g, data = d, weights = d$w, backend = "brms",
    return_code = TRUE
  )
  stan_txt <- paste(code, collapse = "\n")
  expect_match(stan_txt, "offset")
  expect_false(grepl("weights[[]n[]] \\* normal_lpdf", stan_txt))
})
