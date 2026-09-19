# =============================================================================
# INLA containment (2026-09-13).
#
# Estimated correlation matrices -- the correlated random slope (x | g) and
# us(f):g -- are fitted on brms only (correlated_covariance_requires_brms),
# a recorded decision.
#
# The same change made a Gaussian fit on INLA refuse when the response
# standard deviation exceeded 150 (inla_response_scale_unvalidated). That
# refusal is removed: every Gaussian INLA fit now runs on a standardised
# response under scale-aware default priors. The second half of this file
# checks that a large-scale response fits, and fits the same model as the
# response on its original scale, on every route the refusal was written on.
#
# A refusal protects the route it is written on and no other, so every route
# that can put a model on INLA is exercised here: the gate, an explicit
# backend = "inla" fit, the fb_inla() pin, fb_plan() under "auto", the
# aggregated route and the streamed route, the last two of which do not pass
# through the gate.
# =============================================================================

.ct_rules <- function(res) {
  if (!inherits(res, "lgm_refusal")) {
    return(character(0))
  }
  vapply(res$failures, function(f) f$rule_id, character(1L))
}

# A response with an exact standard deviation, around a level far from zero.
.ct_scaled <- function(sd_target, n = 120L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(
    g = factor(rep(sprintf("g%02d", seq_len(12L)), length.out = n)),
    x = stats::rnorm(n)
  )
  z <- stats::rnorm(n)
  d$y <- 1000 + sd_target * (z - mean(z)) / stats::sd(z)
  d
}

.ct_sleepstudy <- function(multiplier = 1) {
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  d <- e$sleepstudy
  d$Reaction <- multiplier * d$Reaction
  d
}


# --- estimated correlation matrices ------------------------------------ #

test_that("the gate refuses both correlated structures and nothing near them", {
  d <- .ct_scaled(1)
  d$env <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  refused <- list(
    fb_from_brms(y ~ x + (x | g), data = d),
    fb_from_asreml(fixed = y ~ x, random = ~ us(env):g, data = d)
  )
  for (fb in refused) {
    res <- lgm_gate(fb)
    expect_identical(res$primary_rule, "correlated_covariance_requires_brms")
    expect_identical(
      res$failures[[1L]]$reason_code,
      "correlated_covariance_requires_brms"
    )
    expect_match(res$failures[[1L]]$reason, "backend = \"brms\"", fixed = TRUE)
  }

  # The neighbours that carry no estimated correlation stay on INLA: the
  # uncorrelated slope and one variance per level.
  admitted <- list(
    fb_from_brms(y ~ x + (x || g), data = d),
    fb_from_asreml(fixed = y ~ x, random = ~ diag(env):g, data = d)
  )
  for (fb in admitted) {
    expect_false(
      "correlated_covariance_requires_brms" %in% .ct_rules(lgm_gate(fb))
    )
  }
})

test_that("an explicit INLA request for either structure raises the typed code", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ct_scaled(1)
  d$env <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  calls <- list(
    function() flexybayes(y ~ x + (x | g), data = d, backend = "inla",
                          return_code = TRUE, verbose = FALSE),
    function() flexybayes(y ~ x, random = ~ us(env):g, data = d,
                          backend = "inla", return_code = TRUE,
                          verbose = FALSE),
    function() suppressWarnings(
      fb_inla(y ~ x + (x | g), data = d, return_code = TRUE, verbose = FALSE),
      classes = "flexybayes_deprecated"
    )
  )
  for (call in calls) {
    err <- tryCatch(suppressMessages(call()), error = function(e) e)
    expect_s3_class(
      err, "flexybayes_refusal_correlated_covariance_requires_brms"
    )
    expect_identical(err$reason_code, "correlated_covariance_requires_brms")
    expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
  }
})

test_that("auto routes both structures to brms before anything is fitted", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ct_scaled(1)
  d$env <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  p_slope <- fb_plan(y ~ x + (x | g), data = d, backend = "auto")
  expect_identical(p_slope$backend_chosen, "brms")
  expect_identical(p_slope$gate_primary_rule,
                   "correlated_covariance_requires_brms")

  p_us <- suppressMessages(flexybayes(
    y ~ x, random = ~ us(env):g, data = d, backend = "auto", plan = TRUE
  ))
  expect_identical(p_us$backend_chosen, "brms")
  expect_identical(p_us$gate_primary_rule,
                   "correlated_covariance_requires_brms")

  # An explicit INLA plan reports the same rule and does not claim a fit.
  p_inla <- fb_plan(y ~ x + (x | g), data = d, backend = "inla")
  expect_false(p_inla$will_fit)
  expect_identical(p_inla$gate_primary_rule,
                   "correlated_covariance_requires_brms")
})


# --- the response scale ------------------------------------------------ #

# Each route fits a response one thousand times its original scale and is
# compared with the same model on the original response, which is that
# response divided by the multiplier. Multiplying the response by a constant
# changes no model, so every fixed effect and every standard deviation,
# divided by the multiplier, is the same number. The expected values are
# therefore a second fit, never a number read off the code.
#
# Tolerances, from ten identical fits at each scale on each route (INLA
# 26.8.7), as the largest relative difference of any estimate, within one
# scale or between a large-scale fit and a reference fit:
#
#   fixed effects        at most 2.9e-15 on the per-row routes and 9.8e-12
#                        on the aggregated and streamed routes. Tolerance
#                        1e-9, the allowance test-inla-scale-aware-defaults.R
#                        makes for engine output that differs between
#                        platforms in its last digits.
#   standard deviations  at most 4.9e-6 on the per-row routes (sd_Subject)
#                        and 7.3e-7 on the aggregated and streamed routes.
#                        Tolerance 1.5e-5, about three times the larger.
#
# The block design has twenty blocks because a block standard deviation
# over four ran 2 per cent apart between identical aggregated fits at either
# scale; over twenty it holds to 7.3e-7. Before the response was
# standardised a thousandfold fit on sleepstudy returned ratios of the
# subject variance to the residual between 0.03 and 6.99 against 1.198,
# orders of magnitude outside either tolerance.
.ct_multiplier <- 1000
.ct_tol <- list(fixed = 1e-9, sd = 1.5e-5)

.ct_blocks <- function(multiplier = 1) {
  set.seed(20260913L)
  d <- expand.grid(
    rep = seq_len(4L),
    trt = factor(seq_len(5L)),
    blk = factor(seq_len(20L))
  )
  d$y <- multiplier * (10 + as.numeric(d$trt) * 0.5 +
                         stats::rnorm(nrow(d), 0, 0.3) +
                         stats::rnorm(20L, 0, 0.4)[as.integer(d$blk)])
  d
}

# A large-scale fit is a converged INLA fit of the same model as the
# reference: the same estimates, each divided by the multiplier.
.ct_expect_same_model <- function(big, reference) {
  expect_s3_class(big, "flexybayes_inla")
  expect_true(isTRUE(summary(big)$converge$ok))
  fixed_big <- coef(big) / .ct_multiplier
  fixed_ref <- coef(reference)
  expect_identical(names(fixed_big), names(fixed_ref))
  expect_lt(max(abs(fixed_big - fixed_ref) / abs(fixed_ref)), .ct_tol$fixed)

  vc_big <- summary(big)$varcomp
  vc_ref <- summary(reference)$varcomp
  expect_identical(vc_big$component, vc_ref$component)
  sd_big <- vc_big$estimate / .ct_multiplier
  expect_lt(max(abs(sd_big - vc_ref$estimate) / vc_ref$estimate), .ct_tol$sd)
}

test_that("the gate admits a Gaussian model at any response scale", {
  skip_if_not_installed("lme4")
  d <- .ct_sleepstudy(.ct_multiplier)
  expect_gt(stats::sd(d$Reaction), 150)

  expect_s3_class(
    lgm_gate(fb_from_brms(Reaction ~ Days + (1 | Subject), data = d)),
    "fb_terms"
  )
  expect_s3_class(
    lgm_gate(fb_from_asreml(fixed = Reaction ~ Days, random = ~ Subject,
                            data = d)),
    "fb_terms"
  )
  # A transformed response is gated as evaluated, and ten times a large one
  # is admitted too.
  expect_s3_class(
    lgm_gate(fb_from_asreml(fixed = I(10 * Reaction) ~ Days,
                            random = ~ Subject, data = d)),
    "fb_terms"
  )
})

test_that("a large-scale response fits on INLA through backend, pin and auto", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d_big <- .ct_sleepstudy(.ct_multiplier)
  d_ref <- d_big
  d_ref$Reaction <- d_ref$Reaction / .ct_multiplier

  reference <- suppressWarnings(flexybayes(
    Reaction ~ Days + (1 | Subject), data = d_ref, backend = "inla",
    verbose = FALSE
  ))

  explicit <- suppressWarnings(flexybayes(
    Reaction ~ Days + (1 | Subject), data = d_big, backend = "inla",
    verbose = FALSE
  ))
  .ct_expect_same_model(explicit, reference)

  pinned <- suppressWarnings(fb_inla(
    Reaction ~ Days + (1 | Subject), data = d_big, verbose = FALSE
  ))
  .ct_expect_same_model(pinned, reference)

  # Under "auto" the plan names INLA and the gate accepts, where the refusal
  # used to route this model to brms before anything was fitted.
  p <- fb_plan(Reaction ~ Days + (1 | Subject), data = d_big,
               backend = "auto")
  expect_identical(p$backend_chosen, "inla")
  expect_identical(p$gate_outcome, "accept")
  expect_true(p$will_fit)
  auto <- suppressWarnings(flexybayes(
    Reaction ~ Days + (1 | Subject), data = d_big, backend = "auto",
    verbose = FALSE
  ))
  .ct_expect_same_model(auto, reference)
})

test_that("a large-scale Gaussian response refuses by name on the aggregated and streamed routes", {
  # The aggregated and streamed Gaussian routes refuse whatever the response
  # scale, because they cannot carry the residual prior. The refusal is that
  # one, not a response-scale refusal; the per-row large-scale fits above
  # carry the scale check.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d_big <- .ct_blocks(.ct_multiplier)
  expect_gt(stats::sd(d_big$y), 150)

  expect_error(
    suppressWarnings(flexybayes(
      fixed = y ~ trt, random = ~ blk, data = d_big,
      backend = "inla", aggregate = TRUE, verbose = FALSE
    )),
    class = "flexybayes_refusal_aggregation_residual_prior_not_composed"
  )
  expect_error(
    flexybayes_stream(
      y ~ trt, random = ~ blk, source = d_big, backend = "inla",
      chunk_rows = 150L, verbose = FALSE
    ),
    class = "flexybayes_refusal_aggregation_residual_prior_not_composed"
  )
})
