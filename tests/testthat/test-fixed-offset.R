# =============================================================================
# offset() in the fixed formula.
#
# terms() files an offset apart from the term labels, and the fixed-formula
# parser read only the labels, so `Reaction ~ Days + offset(Days)` fitted
# `Reaction ~ Days` on INLA and on brms with no message: a Days coefficient
# of 10.47 where lme4 with the offset returns 9.47. The parser now records
# every offset on the model description and evaluates it against the data;
# a route that cannot carry one refuses by name.
# =============================================================================

.off_sleepstudy <- function() {
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  e$sleepstudy
}

.off_expect_refusal <- function(expr, code, pattern = NULL) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, paste0("flexybayes_refusal_", code))
  if (inherits(err, "condition") && !is.null(pattern)) {
    for (p in pattern) {
      expect_match(conditionMessage(err), p, fixed = TRUE)
    }
  }
  invisible(err)
}

# --- ingest ---------------------------------------------------------------- #

test_that("both grammars record the offsets on the model description", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  ir_a <- fb_from_asreml(
    fixed = Reaction ~ Days + offset(Days), random = ~ Subject, data = d
  )
  ir_b <- fb_from_brms(Reaction ~ Days + offset(Days) + (1 | Subject), data = d)
  expect_identical(ir_a$offsets, "offset(Days)")
  expect_identical(ir_b$offsets, "offset(Days)")
  # The offset is not a fixed term: it carries no coefficient.
  expect_identical(
    vapply(ir_a$fixed_terms, function(t) t$label, character(1L)),
    "Days"
  )
  expect_identical(ir_a$fixed_terms, ir_b$fixed_terms)
})

test_that("several offsets are recorded in formula order and add", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  ir <- fb_from_asreml(
    fixed = Reaction ~ offset(Days) + Days + offset(log(Days + 1)),
    random = ~ Subject,
    data = d
  )
  expect_identical(ir$offsets, c("offset(Days)", "offset(log(Days + 1))"))
  # Hand computation, not the helper's own arithmetic.
  expect_equal(
    .fb_offset_values(ir, d),
    d$Days + log(d$Days + 1),
    tolerance = 0
  )
  expect_null(.fb_offset_values(fb_from_asreml(Reaction ~ Days, data = d), d))
})

test_that("a column named offset is a covariate, not an offset", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  d$offset <- d$Days
  fb <- fb_from_asreml(fixed = Reaction ~ offset, random = ~ Subject, data = d)
  expect_identical(fb$fixed_terms[[1L]]$var, "offset")
  expect_identical(fb$fixed_terms[[1L]]$type, "continuous")
  expect_length(fb$offsets, 0L)
})

test_that("an offset that is not a finite number on every row refuses", {
  set.seed(3L)
  d <- data.frame(
    y = stats::rnorm(20),
    x = stats::rnorm(20),
    f = factor(rep(c("a", "b"), 10)),
    z = stats::rnorm(20)
  )
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(f), data = d),
    "fixed_offset_invalid",
    c("`offset(f)`", "factor")
  )
  d_na <- d
  d_na$z[c(4, 9)] <- NA
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(z), data = d_na),
    "fixed_offset_invalid",
    c("`offset(z)`", "row(s) 4, 9")
  )
  d_inf <- d
  d_inf$z[7] <- Inf
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(z), data = d_inf),
    "fixed_offset_invalid",
    c("not finite", "row(s) 7")
  )
  # A workspace object of the right length is not a column, and the
  # baseenv() enclosure keeps it from being read as one.
  zz <- stats::rnorm(20)
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(zz), data = d),
    "fixed_offset_invalid",
    "`offset(zz)`"
  )
  .off_expect_refusal(
    fb_from_brms(y ~ x + offset(zz), data = d),
    "fixed_offset_invalid"
  )
})

test_that("offsets in the random or residual formula keep their refusal", {
  set.seed(4L)
  d <- data.frame(
    y = stats::rnorm(20),
    g = factor(rep(1:4, 5)),
    z = stats::rnorm(20)
  )
  .off_expect_refusal(
    fb_from_asreml(y ~ 1, random = ~ g + offset(z), data = d),
    "asreml_function_not_recognised"
  )
  .off_expect_refusal(
    fb_from_asreml(y ~ 1, residual = ~ offset(z), data = d),
    "asreml_function_not_recognised"
  )
})

test_that("a design cell completed under augment needs its offset", {
  # The completed row would otherwise carry row 1's offset. `z` varies, so
  # its value at an absent cell cannot be determined and the completion
  # refuses, as it does for any other model variable.
  d <- expand.grid(row = factor(1:3), col = factor(1:3))
  d$y <- stats::rnorm(9)
  d$z <- seq_len(9) / 10
  d <- d[-5, ]
  fb <- fb_from_asreml(
    y ~ offset(z), random = ~ ar1(row):ar1(col), data = d
  )
  expect_true("z" %in% .fb_model_vars(fb))
  err <- tryCatch(
    .fb_apply_na_action(fb, d, "augment"),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "z", fixed = TRUE)
})

# --- brms ------------------------------------------------------------------ #

test_that("the brms formula carries each offset verbatim", {
  set.seed(5L)
  d <- data.frame(
    y = stats::rnorm(24),
    x = stats::rnorm(24),
    g = factor(rep(1:6, 4)),
    E = stats::runif(24, 1, 5)
  )
  ir <- fb_from_asreml(
    y ~ x + offset(log(E)) + offset(x), random = ~ g, data = d
  )
  expect_identical(
    deparse(.fb_to_brms_formula(ir)),
    "y ~ x + offset(log(E)) + offset(x) + (1 | g)"
  )
  ir0 <- fb_from_asreml(y ~ 0 + offset(x), random = ~ g, data = d)
  expect_identical(
    deparse(.fb_to_brms_formula(ir0)),
    "y ~ offset(x) - 1 + (1 | g)"
  )
})

test_that("a brms fit carries the offset into Stan and matches lme4", {
  # Four seeds of this call (2 chains, 1000 draws each) put the Days
  # posterior mean 0.034 or less from lme4's 9.467286 with the offset, a
  # spread of 0.067 across seeds (2026-09-13); without the offset the
  # coefficient is 10.47, a full unit away. The tolerance is 0.2.
  skip_if_not_installed("brms")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_sleepstudy()
  fit <- flexybayes(
    Reaction ~ Days + offset(Days),
    random = ~ Subject,
    data = d,
    backend = "brms",
    chains = 2L,
    warmup = 500L,
    n_samples = 1000L,
    seed = 3L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  expect_match(
    paste(deparse(fit$brms$formula$formula), collapse = ""),
    "offset(Days)",
    fixed = TRUE
  )
  expect_identical(
    as.numeric(brms::standata(fit$brms)$offsets),
    as.numeric(d$Days)
  )
  ref <- lme4::fixef(
    lme4::lmer(Reaction ~ Days + offset(Days) + (1 | Subject), data = d)
  )
  expect_lt(abs(coef(fit)[["Days"]] - ref[["Days"]]), 0.2)
  # The fitted values include the offset: posterior means of the intercept,
  # the slope and each subject's deviation, plus Days, by hand.
  re <- brms::ranef(fit$brms)$Subject[, "Estimate", "Intercept"]
  hand <- coef(fit)[["(Intercept)"]] + coef(fit)[["Days"]] * d$Days +
    re[as.character(d$Subject)] + d$Days
  expect_equal(unname(fitted(fit)), unname(hand), tolerance = 1e-8)

  # A population-level prediction on new data is linear in the draws, so its
  # posterior mean is the posterior-mean coefficients applied to newdata plus
  # the offset evaluated there.
  nd <- data.frame(Days = c(0, 3.5, 9), Subject = d$Subject[1:3])
  hand_nd <- coef(fit)[["(Intercept)"]] + coef(fit)[["Days"]] * nd$Days +
    nd$Days
  expect_equal(
    unname(predict(fit, newdata = nd, re_formula = NA)),
    hand_nd,
    tolerance = 1e-8
  )
  .off_expect_refusal(
    predict(fit, newdata = nd[, "Subject", drop = FALSE]),
    "newdata_offset_missing",
    "`Days`"
  )
})

# --- INLA, per row --------------------------------------------------------- #

.off_inla <- function(...) {
  suppressMessages(flexybayes(
    ...,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  ))
}

test_that("a per-row Gaussian INLA fit with an offset matches lme4", {
  # Five identical calls returned the same coefficients to 1.1e-13, and
  # they differed from lme4 with the offset by 8.4e-12 (intercept) and
  # 1.6e-12 (Days) (2026-09-13). The tolerance is 1e-6; dropping the offset
  # moves Days by 1.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_sleepstudy()
  fit <- .off_inla(Reaction ~ Days + offset(Days), random = ~ Subject, data = d)
  ref <- lme4::fixef(
    lme4::lmer(Reaction ~ Days + offset(Days) + (1 | Subject), data = d)
  )
  expect_lt(abs(coef(fit)[["(Intercept)"]] - ref[["(Intercept)"]]), 1e-6)
  expect_lt(abs(coef(fit)[["Days"]] - ref[["Days"]]), 1e-6)

  # offset(Days) is the model without it with the Days coefficient raised by
  # exactly one, so the two fits share every fitted value. The default
  # priors of the two differ slightly, since each reads the spread of its
  # own y - o. Measured over five pairs of these calls (2026-09-13): the
  # coefficients agree to 9.1e-14, while the fitted values, read from INLA's
  # per-row marginals, differ by 4.3e-3 to 4.5e-3 at INLA's default thread
  # setting and by 5.6e-5 on one thread. The fitted-value tolerance is 0.05;
  # the offset it has to carry reaches 9.
  plain <- .off_inla(Reaction ~ Days, random = ~ Subject, data = d)
  expect_lt(abs(coef(plain)[["Days"]] - coef(fit)[["Days"]] - 1), 1e-6)
  expect_lt(max(abs(fitted(fit) - fitted(plain))), 0.05)
  expect_lt(
    max(abs(residuals(fit) - (d$Reaction - fitted(plain)))),
    0.05
  )

  # Each posterior draw's Predictor row is the fixed part, the subject's
  # deviation and the offset, in that draw (measured to 5.7e-14).
  draw <- .fb_inla_posterior_sample(fit, 2L, num.threads = "1:1")[[1L]]$latent
  rn <- rownames(draw)
  eta <- draw["(Intercept):1", 1L] + draw["Days:1", 1L] * d$Days +
    draw[grep("^Subject:", rn), 1L][as.integer(d$Subject)]
  expect_lt(
    max(abs(draw[grep("^Predictor:", rn), 1L] - eta - d$Days)),
    1e-8
  )
})

test_that("a Poisson INLA fit with offset(log(E)) matches glmer", {
  # INLA against glmer's Laplace fit on this design: 3.6e-4 on the intercept
  # and 2.7e-7 on the slope, with five identical calls within 1.7e-9
  # (2026-09-13). Without the offset the intercept is 2.57 against 0.26.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(20260913)
  n <- 400
  g <- factor(sprintf("g%02d", rep(1:20, each = 20)))
  x <- stats::rnorm(n)
  E <- stats::runif(n, 0.5, 20)
  u <- stats::rnorm(20, 0, 0.4)
  d <- data.frame(
    y = stats::rpois(n, E * exp(0.3 + 0.5 * x + u[as.integer(g)])),
    x = x,
    g = g,
    E = E
  )
  ref <- lme4::fixef(lme4::glmer(
    y ~ x + offset(log(E)) + (1 | g),
    family = stats::poisson,
    data = d
  ))
  fit <- .off_inla(
    y ~ x + offset(log(E)), random = ~ g, family = "poisson", data = d
  )
  expect_lt(abs(coef(fit)[["(Intercept)"]] - ref[["(Intercept)"]]), 2e-3)
  expect_lt(abs(coef(fit)[["x"]] - ref[["x"]]), 1e-5)
  # The stored linear predictor carries log(E) (measured to 1.4e-4).
  eta <- coef(fit)[["(Intercept)"]] + coef(fit)[["x"]] * x +
    fit$inla$summary.random$g$mean[as.integer(g)] + log(E)
  expect_lt(max(abs(fit$inla$summary.linear.predictor$mean - eta)), 1e-3)
})

test_that("a Bernoulli INLA fit with an offset matches glmer", {
  # INLA against glmer: 3.6e-3 on the intercept, 4.9e-3 on the slope, two
  # identical calls identical (2026-09-13). Without the offset the slope is
  # 0.536 against 0.610.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(7)
  n <- 600
  g <- factor(rep(1:30, each = 20))
  x <- stats::rnorm(n)
  o <- stats::rnorm(n, 0, 0.8)
  u <- stats::rnorm(30, 0, 0.5)
  d <- data.frame(
    y = stats::rbinom(n, 1, stats::plogis(-0.3 + 0.7 * x + o + u[as.integer(g)])),
    x = x,
    g = g,
    o = o
  )
  ref <- lme4::fixef(lme4::glmer(
    y ~ x + offset(o) + (1 | g),
    family = stats::binomial,
    data = d
  ))
  fit <- .off_inla(y ~ x + offset(o), random = ~ g, family = "binomial", data = d)
  expect_lt(abs(coef(fit)[["(Intercept)"]] - ref[["(Intercept)"]]), 0.02)
  expect_lt(abs(coef(fit)[["x"]] - ref[["x"]]), 0.02)
})

test_that("an offset on a sectioned residual shifts the fit exactly", {
  # Two Gaussian likelihoods, one per section. offset(Days) against the fit
  # without it, five pairs at INLA's default thread setting (2026-09-13):
  # Days differs by 1 within 3.6e-5 and the intercept by up to 6.0e-5; under
  # one stated prior on one thread the gaps reached 8.0e-5 and 1.4e-4. The
  # tolerance is 1e-3; without the offset Days moves by 1.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_sleepstudy()
  d$half <- factor(ifelse(d$Days >= 5, "late", "early"))
  res <- ~ dsum(~ units | half)
  with_off <- .off_inla(
    Reaction ~ Days + offset(Days), random = ~ Subject, residual = res, data = d
  )
  plain <- .off_inla(Reaction ~ Days, random = ~ Subject, residual = res, data = d)
  expect_lt(abs(coef(plain)[["Days"]] - coef(with_off)[["Days"]] - 1), 1e-3)
  expect_lt(
    abs(coef(plain)[["(Intercept)"]] - coef(with_off)[["(Intercept)"]]),
    1e-3
  )
})

# --- standardisation and default priors ----------------------------------- #

.off_large <- function(sd_o, seed = 20260913L) {
  set.seed(seed)
  n_g <- 30L
  per <- 20L
  n <- n_g * per
  g <- factor(sprintf("g%02d", rep(seq_len(n_g), each = per)))
  x <- stats::rnorm(n)
  o <- stats::rnorm(n, 0, sd_o)
  u <- stats::rnorm(n_g, 0, 2)
  data.frame(
    y = o + 5 + 1.5 * x + u[as.integer(g)] + stats::rnorm(n, 0, 3),
    x = x,
    g = g,
    o = o
  )
}

test_that("the scale an offset model is fitted on is read from y - o", {
  # A 600-row design whose offset has standard deviation 1000: sd(y) = 999.3,
  # sd(y - o) = 4.11. Read from y, the response-scale bound of 150 refuses
  # the fit; read from y - o it fits, and three identical calls were 2.0e-6
  # (intercept), 1.0e-4 (slope), 0.103 (random-effect SD, a posterior mean
  # against REML, 4.3 per cent) and 5.7e-3 (residual SD) from lme4, the same
  # distances as at sd(o) = 0, 30 and 100 under either reading (2026-09-13).
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_large(1000)
  w <- d$y - d$o

  ir <- fb_from_asreml(y ~ x + offset(o), random = ~ g, data = d)
  expect_equal(ir$data_summary$response_sd, stats::sd(w), tolerance = 1e-12)

  fit <- .off_inla(y ~ x + offset(o), random = ~ g, data = d)
  spec <- fit$extras$response_standardisation
  expect_equal(spec$centre, mean(w), tolerance = 1e-12)
  expect_equal(spec$scale, stats::sd(w), tolerance = 1e-12)
  # The shared default bound, U = 5 sd, reads the same quantity; read from
  # y it would be 5 x 999.3.
  expect_equal(
    attr(fit$extras$fb_terms$priors, "fb_prior_default_scale"),
    5 * stats::sd(w),
    tolerance = 1e-12
  )

  ref_fit <- lme4::lmer(y ~ x + offset(o) + (1 | g), data = d)
  ref <- lme4::fixef(ref_fit)
  vc <- summary(fit)$varcomp
  expect_lt(abs(coef(fit)[["(Intercept)"]] - ref[["(Intercept)"]]), 1e-4)
  expect_lt(abs(coef(fit)[["x"]] - ref[["x"]]), 1e-3)
  expect_lt(
    abs(vc$estimate[vc$component == "sd_g"] -
      attr(lme4::VarCorr(ref_fit)$g, "stddev")[[1L]]),
    0.3
  )
  expect_lt(
    abs(vc$estimate[vc$component == "sigma"] - stats::sigma(ref_fit)),
    0.05
  )
})

test_that("an offset model is invariant to rescaling the response and offset", {
  # Response and offset multiplied by k in {0.01, 1, 1000}, the fit divided
  # back by k. Over three identical calls at each k (2026-09-13) the means
  # of the rescaled intercept, slope, random-effect SD, residual SD and first
  # fitted value ranged across k by 7.0e-10, 2.8e-8, 9.9e-7, 1.1e-7 and
  # 2.1e-7, inside the ranges of identical calls at one k (up to 1.3e-9,
  # 5.5e-8, 5.1e-6, 4.0e-7 and 4.0e-7). Relative tolerances: 1e-5 on the
  # coefficients and fitted values, 1e-3 on the standard deviations, whose
  # identical-call range is 1.9e-4 of their value.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(20260914L)
  n_g <- 30L
  per <- 20L
  n <- n_g * per
  g <- factor(sprintf("g%02d", rep(seq_len(n_g), each = per)))
  x <- stats::rnorm(n)
  o <- stats::rnorm(n)
  u <- stats::rnorm(n_g, 0, 0.02)
  base <- data.frame(
    y = o + 0.05 + 0.015 * x + u[as.integer(g)] + stats::rnorm(n, 0, 0.03),
    x = x,
    g = g,
    o = o
  )
  read <- function(k) {
    d <- base
    d$yk <- d$y * k
    d$ok <- d$o * k
    fit <- .off_inla(yk ~ x + offset(ok), random = ~ g, data = d)
    vc <- summary(fit)$varcomp
    list(
      coef = unname(coef(fit)) / k,
      sd = c(
        vc$estimate[vc$component == "sd_g"],
        vc$estimate[vc$component == "sigma"]
      ) / k,
      fitted = fitted(fit) / k
    )
  }
  unit <- read(1)
  for (k in c(0.01, 1000)) {
    scaled <- read(k)
    expect_equal(scaled$coef, unit$coef, tolerance = 1e-5)
    expect_equal(scaled$sd, unit$sd, tolerance = 1e-3)
    expect_equal(scaled$fitted, unit$fitted, tolerance = 1e-5)
  }
})

# --- aggregated and streamed routes ---------------------------------------- #

# A replicated factorial the aggregated route compresses: 480 rows in 48
# cells, an offset that varies inside every cell, and one INLA thread, as the
# aggregation equivalence tests are run.
.off_cells <- function() {
  set.seed(20260915L)
  n <- 480L
  d <- data.frame(
    f = factor(rep(sprintf("f%d", 1:4), length.out = n)),
    g = factor(rep(sprintf("g%02d", 1:12), each = n / 12))
  )
  fe <- c(f1 = 0, f2 = 0.5, f3 = -0.3, f4 = 0.8)
  u <- stats::rnorm(12, 0, 0.4)
  d$z <- stats::rnorm(n, 0, 2)
  d$y <- 3 + fe[as.character(d$f)] + u[as.integer(d$g)] + d$z +
    stats::rnorm(n, 0, 0.7)
  d$E <- stats::runif(n, 1, 10)
  d$count <- stats::rpois(
    n,
    d$E * exp(0.2 + fe[as.character(d$f)] + u[as.integer(d$g)])
  )
  d$bin <- stats::rbinom(n, 1, stats::plogis(-0.2 + fe[as.character(d$f)] + d$z / 2))
  d
}

.off_route <- function(route, fixed, family, d, ...) {
  prior <- fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 10))
  suppressMessages(switch(
    route,
    per_row = flexybayes(
      fixed, random = ~ g, data = d, family = family, prior = prior,
      backend = "inla", aggregate = FALSE, verbose = FALSE, ...
    ),
    aggregated = flexybayes(
      fixed, random = ~ g, data = d, family = family, prior = prior,
      backend = "inla", aggregate = TRUE, verbose = FALSE, ...
    ),
    stream = flexybayes_stream(
      fixed, random = ~ g, source = d, family = family, prior = prior,
      backend = "inla", chunk_rows = 100, verbose = FALSE, ...
    )
  ))
}

.off_fixed <- function(fit) {
  stats::setNames(fit$inla$summary.fixed$mean, rownames(fit$inla$summary.fixed))
}

.off_hyper <- function(fit) {
  stats::setNames(
    fit$inla$summary.hyperpar$mean,
    rownames(fit$inla$summary.hyperpar)
  )
}

test_that("aggregated and streamed Gaussian fits with an offset refuse by name", {
  # The cell statistics still subtract the offset exactly, which the plan
  # records. The aggregated and streamed Gaussian routes refuse, because they
  # cannot carry the residual prior; their agreement with the per-row fit
  # (measured 2026-09-13 at 6.6e-12 on the fixed effects) returns with that
  # composition.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_cells()

  plan <- .fb_aggregation_plan(
    fb_from_asreml(y ~ f + offset(z), random = ~ g, data = d),
    .fb_dataset(d)
  )
  expect_true(plan$eligible)
  expect_identical(plan$offset_route, "offset_subtracted_in_cell_statistics")

  for (route in c("aggregated", "stream")) {
    expect_error(
      .off_route(route, y ~ f + offset(z), "gaussian", d),
      class = "flexybayes_refusal_aggregation_residual_prior_not_composed"
    )
  }
})

test_that("a Poisson offset reaches the aggregated route as exposure", {
  # Against the per-row fit of count ~ f + offset(log(E)), three calls each
  # (2026-09-13): the aggregated and the streamed fit within 4.5e-10 (fixed)
  # and 5.7e-9 relative (hyperparameters), a second per-row call within
  # 1.5e-10 and 1.8e-9; the streamed fit written with exposure = "E" and no
  # offset within 9.6e-10 and 9.3e-9. Without the offset the aggregated
  # intercept is 1.95 against 0.26.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  withr::defer(INLA::inla.setOption(num.threads = threads))
  d <- .off_cells()

  plan <- .fb_aggregation_plan(
    fb_from_asreml(count ~ f + offset(log(E)), random = ~ g, data = d,
                   family = "poisson"),
    .fb_dataset(d)
  )
  expect_true(plan$eligible)
  expect_identical(plan$offset_route, "offset_folded_into_exposure")

  per_row <- .off_route("per_row", count ~ f + offset(log(E)), "poisson", d)
  check <- function(fit) {
    b <- .off_fixed(fit)
    expect_lt(max(abs(b[names(.off_fixed(per_row))] - .off_fixed(per_row))), 1e-7)
    h <- .off_hyper(fit)
    expect_equal(
      unname(h[names(.off_hyper(per_row))]),
      unname(.off_hyper(per_row)),
      tolerance = 1e-6
    )
  }
  agg <- .off_route("aggregated", count ~ f + offset(log(E)), "poisson", d)
  check(agg)
  check(.off_route("stream", count ~ f + offset(log(E)), "poisson", d))
  check(.off_route("stream", count ~ f, "poisson", d, exposure = "E"))

  X <- stats::model.matrix(~ f, d)
  expect_lt(
    max(abs(agg$glm$linear.predictors -
      (as.numeric(X %*% coef(agg)[colnames(X)]) + log(d$E)))),
    1e-10
  )
})

test_that("a binomial offset refuses the aggregated routes by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_cells()
  plan <- .fb_aggregation_plan(
    fb_from_asreml(bin ~ f + offset(z / 2), random = ~ g, data = d,
                   family = "binomial"),
    .fb_dataset(d)
  )
  expect_false(plan$eligible)
  expect_true("fixed_offset_unsupported" %in% plan$reason_codes)
  expect_identical(plan$offset_route, "fixed_offset_unsupported")

  err <- tryCatch(
    .off_route("aggregated", bin ~ f + offset(z / 2), "binomial", d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_true("fixed_offset_unsupported" %in% err$reason_codes)
  expect_match(
    conditionMessage(err),
    "binomial only when every row in the cell shares one success probability",
    fixed = TRUE
  )

  .off_expect_refusal(
    .off_route("stream", bin ~ f + offset(z / 2), "binomial", d),
    "fixed_offset_unsupported",
    "aggregated binomial route"
  )

  # aggregate = "auto" takes the per-row route, which carries the offset.
  auto <- suppressMessages(flexybayes(
    bin ~ f + offset(z / 2), random = ~ g, data = d, family = "binomial",
    backend = "inla", verbose = FALSE
  ))
  expect_identical(auto$exactness, "exact")
})

# --- prediction and the recorded model ------------------------------------ #

.off_pred_data <- function() {
  set.seed(20260915L)
  n <- 480L
  d <- data.frame(
    f = factor(rep(sprintf("f%d", 1:4), length.out = n)),
    g = factor(rep(sprintf("g%02d", 1:12), each = n / 12))
  )
  fe <- c(f1 = 0, f2 = 0.5, f3 = -0.3, f4 = 0.8)
  d$x <- stats::rnorm(n)
  d$z <- stats::rnorm(n, 0, 2)
  d$y <- 3 + fe[as.character(d$f)] + 0.4 * d$x +
    stats::rnorm(12, 0, 0.4)[as.integer(d$g)] + d$z + stats::rnorm(n, 0, 0.7)
  d
}

test_that("predict() on INLA carries the offset where it predicts", {
  # Against X b + z by hand the in-sample and new-data predictions and the
  # classify table matched to 0 (2026-09-13), so the tolerances are 1e-10.
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_pred_data()
  fit <- .off_inla(y ~ f + x + offset(z), random = ~ g, data = d)
  b <- coef(fit)

  X <- stats::model.matrix(~ f + x, d)
  expect_lt(
    max(abs(predict(fit) - (as.numeric(X %*% b[colnames(X)]) + d$z))),
    1e-10
  )

  nd <- data.frame(
    f = factor(c("f1", "f3", "f4"), levels = levels(d$f)),
    x = c(-1, 0, 2),
    z = c(10, -5, 0.25)
  )
  Xn <- stats::model.matrix(~ f + x, nd)
  expect_lt(
    max(abs(predict(fit, newdata = nd) -
      (as.numeric(Xn %*% b[colnames(Xn)]) + nd$z))),
    1e-10
  )
  err <- .off_expect_refusal(
    predict(fit, newdata = nd[, c("f", "x")]),
    "newdata_offset_missing",
    c("`offset(z)`", "`z`")
  )
  expect_identical(err$missing, "z")

  # emmeans evaluates the offset on its reference grid, where each covariate
  # sits at its mean: b0 + b_f + b_x mean(x) + mean(z) by hand. Without the
  # offset every mean would be 0.064 lower.
  cl <- predict(fit, classify = "f")
  hand <- b[["(Intercept)"]] + c(0, b[["ff2"]], b[["ff3"]], b[["ff4"]]) +
    b[["x"]] * mean(d$x) + mean(d$z)
  expect_lt(max(abs(cl$estimate - hand)), 1e-10)
})

test_that("predict() on brms refuses new data that lack an offset column", {
  # The guard runs before brms is asked for anything, so a fit object that
  # carries the model description is enough to exercise it.
  skip_if_not_installed("brms")
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  ir <- fb_from_asreml(
    Reaction ~ Days + offset(log(Days + 1)), random = ~ Subject, data = d
  )
  stub <- structure(
    list(brms = list(), extras = list(fb_terms = ir)),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  .off_expect_refusal(
    predict(stub, newdata = data.frame(Subject = d$Subject[1:2])),
    "newdata_offset_missing",
    c("`offset(log(Days + 1))`", "`Days`")
  )
})

test_that("the fingerprint and the printed surfaces record it", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_pred_data()
  fit <- .off_inla(y ~ f + x + offset(z), random = ~ g, data = d)
  plain <- .off_inla(y ~ f + x, random = ~ g, data = d)

  expect_identical(fit$extras$fingerprint$offset, "offset(z)")
  expect_identical(plain$extras$fingerprint$offset, "<none>")

  ps <- prior_summary(fit)
  expect_identical(ps$offsets, "offset(z)")
  expect_true(any(grepl(
    "Offset: offset(z)", utils::capture.output(print(ps)), fixed = TRUE
  )))

  expect_identical(summary(fit)$offset, "offset(z)")
  expect_false("offset" %in% names(summary(plain)))

  plan <- fb_plan(y ~ f + offset(z), random = ~ g, data = d)
  expect_identical(plan$offsets, "offset(z)")
  expect_identical(
    plan$aggregation$offset_route,
    "offset_subtracted_in_cell_statistics"
  )
  expect_true(any(grepl(
    "Offset:                  offset(z)",
    utils::capture.output(print(plan)),
    fixed = TRUE
  )))
})

# --- the two engines on one offset model ----------------------------------- #

test_that("brms and INLA agree on an offset model under one stated prior", {
  # lme4::sleepstudy, Reaction ~ Days + offset(Days) + (1 | Subject), uniform
  # priors on both standard deviations. brms at four seeds (4 chains, 2000
  # draws each) against INLA (2026-09-13): at most 0.77 apart on the
  # intercept, 0.031 on Days, 0.50 on the subject SD and 0.044 on the
  # residual SD, inside brms's own spread across the seeds (0.75, 0.038,
  # 0.39, 0.025). The same model without the offset gave 0.40, 0.007, 0.38
  # and 0.049. Tolerances 2, 0.1, 1.5 and 0.15; a dropped offset on either
  # engine moves Days by 1.
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_sleepstudy()
  prior <- fb_prior(
    sd(group = "Subject") ~ uniform(lower = 0, upper = 250),
    sigma ~ uniform(lower = 0, upper = 250)
  )
  read <- function(fit) {
    vc <- summary(fit)$varcomp
    c(
      coef(fit)[["(Intercept)"]],
      coef(fit)[["Days"]],
      vc$estimate[vc$component == "sd_Subject"],
      vc$estimate[vc$component == "sigma"]
    )
  }
  inla_fit <- .off_inla(
    Reaction ~ Days + offset(Days), random = ~ Subject, data = d, prior = prior
  )
  brms_fit <- suppressMessages(flexybayes(
    Reaction ~ Days + offset(Days),
    random = ~ Subject,
    data = d,
    prior = prior,
    backend = "brms",
    chains = 4L,
    warmup = 1000L,
    n_samples = 2000L,
    seed = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  gap <- abs(read(brms_fit) - read(inla_fit))
  expect_lt(gap[[1L]], 2)
  expect_lt(gap[[2L]], 0.1)
  expect_lt(gap[[3L]], 1.5)
  expect_lt(gap[[4L]], 0.15)
})
