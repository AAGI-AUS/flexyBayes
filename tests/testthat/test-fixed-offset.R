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
  expect_identical(fb$fixed[[1L]]$var, "offset")
  expect_identical(fb$fixed[[1L]]$type, "continuous")
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

# --- routes ---------------------------------------------------------------- #

test_that("a fitting route that does not carry an offset refuses by name", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  cells <- data.frame(
    f = factor(rep(sprintf("f%d", 1:4), length.out = 120)),
    g = factor(rep(sprintf("g%02d", 1:12), each = 10)),
    z = stats::rnorm(120)
  )
  cells$y <- stats::rnorm(120, 10) + cells$z
  .off_expect_refusal(
    flexybayes(
      y ~ f + offset(z),
      random = ~ g,
      data = cells,
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    ),
    "fixed_offset_unsupported",
    "`offset(z)`"
  )
  .off_expect_refusal(
    suppressMessages(flexybayes_stream(
      y ~ f + offset(z),
      random = ~ g,
      source = cells,
      backend = "inla",
      fit = FALSE,
      verbose = FALSE
    )),
    "fixed_offset_unsupported",
    "`offset(z)`"
  )
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
  # exactly one, so the two fits share every fitted value. Measured over six
  # pairs of identical calls (2026-09-13): the coefficients agree to 1.8e-13,
  # while the fitted values, read from INLA's per-row marginals, differ by
  # 2.3e-5 to 4.8e-5 between the two spellings and by up to 1.3e-5 between
  # identical calls of one spelling. The fitted-value tolerance is 1e-3;
  # the offset it has to carry reaches 9.
  plain <- .off_inla(Reaction ~ Days, random = ~ Subject, data = d)
  expect_lt(abs(coef(plain)[["Days"]] - coef(fit)[["Days"]] - 1), 1e-6)
  expect_lt(max(abs(fitted(fit) - fitted(plain))), 1e-3)
  expect_lt(
    max(abs(residuals(fit) - (d$Reaction - fitted(plain)))),
    1e-3
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
  # without it: Days differs by 1 - 3.9e-7 and the intercept by 6.7e-7
  # (2026-09-13). The tolerance is 1e-5.
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
  expect_lt(abs(coef(plain)[["Days"]] - coef(with_off)[["Days"]] - 1), 1e-5)
  expect_lt(
    abs(coef(plain)[["(Intercept)"]] - coef(with_off)[["(Intercept)"]]),
    1e-5
  )
})
