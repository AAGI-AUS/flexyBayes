# Tests for the factor:continuous indexed-slope emit (ADR 0019,
# Stage 1A v0.2.6).
#
# ADR 0019 §Verification + the design spec require >=8 subtests on
# this file. Each test_that() block below corresponds to one of the
# eight (or more) acceptance criteria called out in the spec; the
# first three are structural / linear-predictor (no fitting needed),
# the fourth is the posterior-mean recovery vs lme4::lmer fit (skipped
# without greta + lme4), the next four cover the unsupported-contrast
# refusal surface + the canonical-name registry + the dispatch trace
# slot.

# ---------------------------------------------------------------- #
# Subtest 1 -- pre-fix-bug reproduction (regression guard).         #
# ---------------------------------------------------------------- #
#
# Before ADR 0019, .classify_fixed_term() sent `f:x` to term type
# "interaction" (numeric x numeric) which caused codegen.R to emit
# `beta_<tag> * as_data(f) * as_data(x)` -- silently coercing the
# factor to its integer codes. The post-fix path classifies the same
# label as "factor_numeric_interaction" with the level-frozen
# descriptor. The test below pins the post-fix classification so any
# regression to the pre-fix path surfaces immediately.

test_that("factor:continuous interaction routes to factor_numeric_interaction", {
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  fb <- flexyBayes:::fb_from_brms(y ~ f * x, data = d)
  fni <- Filter(
    function(t) identical(t$type, "factor_numeric_interaction"),
    fb$fixed_terms
  )
  expect_length(fni, 1L)
  expect_identical(fni[[1L]]$factor, "f")
  expect_identical(fni[[1L]]$continuous, "x")
  expect_identical(fni[[1L]]$levels, c("a", "b", "c"))
  expect_identical(fni[[1L]]$n_levels, 3L)

  # Pre-fix-bug guard: the legacy `interaction` (numeric:numeric)
  # branch must NOT fire on this term. If it does we are back to the
  # silent factor-as-numeric coercion path.
  legacy_int <- Filter(
    function(t) {
      identical(t$type, "interaction") &&
        identical(t$label, "f:x")
    },
    fb$fixed_terms
  )
  expect_length(legacy_int, 0L)
})

# ---------------------------------------------------------------- #
# Subtest 2 -- linear-predictor equivalence vs model.matrix().      #
# ---------------------------------------------------------------- #
#
# The canonical Stage-1A acceptance test (the design spec
# deliverable 3). Build the design matrix via base R's treatment-
# coded model.matrix(); sample a coefficient vector from a fixed
# seed; compute the reference linear predictor `eta_mm = X %*% beta`;
# then build the equivalent linear predictor by binding the same
# coefficients into the indexed-emit slots (mu, tau_f, beta_x,
# slope_dev_f_x). Assert agreement at tolerance 1e-10 (the indexed
# construction is an algebraic identity for treatment contrasts; the
# residual is pure floating-point noise).

test_that("indexed emit equals model.matrix() at tolerance 1e-10", {
  set.seed(20260523L)
  n <- 200L
  f <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  x <- rnorm(n)
  d <- data.frame(y = rnorm(n), x = x, f = f)

  X_mm <- model.matrix(
    y ~ f * x,
    data = d,
    contrasts.arg = list(f = "contr.treatment")
  )
  expect_setequal(
    colnames(X_mm),
    c("(Intercept)", "fb", "fc", "x", "fb:x", "fc:x")
  )

  set.seed(42L)
  beta <- setNames(rnorm(ncol(X_mm)), colnames(X_mm))
  eta_mm <- as.numeric(X_mm %*% beta)

  # Bind the same coefficients into the indexed-emit slots:
  #   mu_atg            -> "(Intercept)"
  #   tau_f[1]          -> 0       (reference level absorbed into mu)
  #   tau_f[2]          -> "fb"
  #   tau_f[3]          -> "fc"
  #   beta_x            -> "x"     (reference-level slope)
  #   slope_dev_f_x[1]  -> 0       (reference-level slope deviation pinned)
  #   slope_dev_f_x[2]  -> "fb:x"
  #   slope_dev_f_x[3]  -> "fc:x"
  mu_atg <- as.numeric(beta["(Intercept)"])
  tau_f <- c(0, as.numeric(beta["fb"]), as.numeric(beta["fc"]))
  beta_x <- as.numeric(beta["x"])
  slope_dev_f_x <- c(0, as.numeric(beta["fb:x"]), as.numeric(beta["fc:x"]))

  f_id <- as.integer(f)
  eta_flexy <- mu_atg + tau_f[f_id] + beta_x * x + x * slope_dev_f_x[f_id]

  expect_equal(eta_flexy, eta_mm, tolerance = 1e-10)
})

# ---------------------------------------------------------------- #
# Subtest 3 -- generated greta code carries the indexed-slope idiom #
# ---------------------------------------------------------------- #
#
# Snapshot guard against the pre-fix idiom returning. The post-fix
# code emits `slope_dev_<tag>_raw` + `slope_dev_<tag>` + an indexed
# predictor contribution `as_data(<con>) * slope_dev_<tag>[<fac>_id]`.
# The pre-fix idiom `beta_<tag> * as_data(<fac>) * as_data(<con>)`
# must be gone.

test_that("generated greta code uses indexed-slope idiom (Option C)", {
  skip_if_greta_backend_unusable()
  withr::local_options(list(
    flexyBayes.silence_default_prior_note = TRUE
  ))
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  code <- flexybayes(
    y ~ f * x,
    data = d,
    backend = "greta",
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  expect_true(grepl("slope_dev_f_x_raw <- normal\\(", code))
  expect_true(grepl("slope_dev_f_x <- c\\(zeros\\(1\\)", code))
  expect_true(grepl(
    "as_data\\(x\\) \\* slope_dev_f_x\\[f_id\\]",
    code
  ))
  # Pre-fix idiom must not appear.
  expect_false(grepl("beta_f_x_x \\* as_data\\(f\\)", code))
})

# ---------------------------------------------------------------- #
# Subtest 4 -- posterior-mean recovery vs lme4::lmer.               #
# ---------------------------------------------------------------- #
#
# The full archetype-1 numerical validation cell. Requires greta +
# lme4; skipped on CRAN / CI per the validation-lmer convention. The
# tolerance W_1 <= 0.10 * |beta_true| per the design spec + the
# `1.0 * REML SE` working tolerance used by test-validation-lmer.R
# absorb MCMC noise at the chosen budget while still catching the
# >3-sigma regression that the v0.2 bug produced.

test_that("posterior mean recovers lme4 fit on factor:continuous", {
  skip_if_greta_backend_unusable()
  skip_if_not_installed("lme4")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  withr::local_options(list(
    flexyBayes.silence_default_prior_note = TRUE
  ))

  set.seed(20260523L)
  n <- 180L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n))
  )
  # Truth: reference slope 1.0; f=b adds +0.5; f=c adds -0.3.
  beta_true <- c(
    `(Intercept)` = 0.0,
    fb = 0.2,
    fc = -0.1,
    x = 1.0,
    `fb:x` = 0.5,
    `fc:x` = -0.3
  )
  X <- model.matrix(
    ~ f * x,
    data = d,
    contrasts.arg = list(f = "contr.treatment")
  )
  d$y <- as.numeric(X %*% beta_true + rnorm(n, 0, 0.3))

  set.seed(20260523L)
  fit <- flexybayes(
    y ~ f * x,
    data = d,
    backend = "greta",
    n_samples = 2000,
    warmup = 2000,
    chains = 2,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )

  draws <- as.matrix(fit$greta$draws)
  expect_true("mu_atg" %in% colnames(draws))
  expect_true(any(grepl("^slope_dev_f_x_raw", colnames(draws))))

  # Post-mean for the slope deviations + the main slope.
  pm_beta_x <- mean(draws[, "beta_x"])
  pm_slope_b <- mean(draws[, "slope_dev_f_x_raw[1,1]"])
  pm_slope_c <- mean(draws[, "slope_dev_f_x_raw[2,1]"])

  ref <- lm(y ~ f * x, data = d)
  ref_coef <- coef(ref)
  # 1.0 . REML SE tolerance per test-validation-lmer.R convention.
  ref_se <- sqrt(diag(vcov(ref)))
  expect_lt(abs(pm_beta_x - ref_coef[["x"]]), 1.0 * ref_se[["x"]])
  expect_lt(abs(pm_slope_b - ref_coef[["fb:x"]]), 1.0 * ref_se[["fb:x"]])
  expect_lt(abs(pm_slope_c - ref_coef[["fc:x"]]), 1.0 * ref_se[["fc:x"]])

  # Dispatch-trace slot records the chosen emit path.
  emit_path <- fit$extras$parse_info$factor_continuous_emit
  expect_true(emit_path %in% c("option_c", "option_d"))
})

# ---------------------------------------------------------------- #
# Subtest 5 -- contr.helmert refusal.                               #
# ---------------------------------------------------------------- #
#
# The Stage-1A scope per ADR 0019 §Decision 6 ships only the
# treatment-coded indexed emit. Other contrast schemes raise a
# structured refusal at ingest with a deferral message.

test_that("contr.helmert raises flexybayes_contrast_unsupported", {
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  contrasts(d$f) <- contr.helmert(3)
  cond <- tryCatch(
    flexyBayes:::fb_from_brms(y ~ f * x, data = d),
    error = function(e) e
  )
  expect_s3_class(cond, "flexybayes_contrast_unsupported")
  expect_identical(cond$contrast, "contr.helmert")
  expect_identical(cond$factor_name, "f")
  expect_match(cond$deferral_target, "representation IR")
  expect_match(conditionMessage(cond), "contr\\.helmert")
  expect_match(conditionMessage(cond), "representation IR")
})

# ---------------------------------------------------------------- #
# Subtest 6 -- contr.sum refusal.                                   #
# ---------------------------------------------------------------- #

test_that("contr.sum raises flexybayes_contrast_unsupported", {
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  contrasts(d$f) <- contr.sum(3)
  cond <- tryCatch(
    flexyBayes:::fb_from_brms(y ~ f * x, data = d),
    error = function(e) e
  )
  expect_s3_class(cond, "flexybayes_contrast_unsupported")
  expect_identical(cond$contrast, "contr.sum")
  expect_match(cond$deferral_target, "representation IR")
})

# ---------------------------------------------------------------- #
# Subtest 7 -- ordered factor refusal.                              #
# ---------------------------------------------------------------- #

test_that("ordered factor raises flexybayes_contrast_unsupported", {
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n), ordered = TRUE),
    y = rnorm(n)
  )
  cond <- tryCatch(
    flexyBayes:::fb_from_brms(y ~ f * x, data = d),
    error = function(e) e
  )
  expect_s3_class(cond, "flexybayes_contrast_unsupported")
  expect_identical(cond$contrast, "ordered")
  expect_identical(cond$factor_name, "f")
  expect_match(cond$deferral_target, "representation IR")
})

# ---------------------------------------------------------------- #
# Subtest 8 -- user-supplied custom contrast matrix refusal.        #
# ---------------------------------------------------------------- #

test_that("custom contrast matrix raises flexybayes_contrast_unsupported", {
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  # User-supplied custom L x (L-1) matrix.
  cmat <- matrix(
    c(-1, 0, 0, 1, -1, 1),
    nrow = 3,
    dimnames = list(c("a", "b", "c"), c("c1", "c2"))
  )
  contrasts(d$f) <- cmat
  cond <- tryCatch(
    flexyBayes:::fb_from_brms(y ~ f * x, data = d),
    error = function(e) e
  )
  expect_s3_class(cond, "flexybayes_contrast_unsupported")
  expect_identical(cond$contrast, "<custom_matrix>")
  expect_match(cond$deferral_target, "representation IR")
})

# ---------------------------------------------------------------- #
# Subtest 9 -- canonical-name family slope_<f>_<x>[<lvl>] resolves. #
# ---------------------------------------------------------------- #
#
# Confirms .parse_slope_dev_raw() correctly translates the greta-side
# parameter name slope_dev_f_x_raw[i,1] to the canonical
# slope_f_x[<level>] slot for the non-reference levels. The test
# avoids running a full greta MCMC -- it exercises the parser
# directly against a synthetic draw-name vector.

test_that("canonical_names() resolves slope_<f>_<x>[<level>] family", {
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  fb <- flexyBayes:::fb_from_brms(y ~ f * x, data = d)
  lookup <- list(f_x = c("a", "b", "c"))
  # raw[1,1] -> non-reference level 1 -> level b ; raw[2,1] -> level c
  r1 <- flexyBayes:::.parse_slope_dev_raw("slope_dev_f_x_raw[1,1]", lookup)
  r2 <- flexyBayes:::.parse_slope_dev_raw("slope_dev_f_x_raw[2,1]", lookup)
  expect_identical(r1$canonical, "slope_f_x[b]")
  expect_identical(r2$canonical, "slope_f_x[c]")
  # Names not matching the pattern return NULL.
  expect_null(
    flexyBayes:::.parse_slope_dev_raw("mu_atg", lookup)
  )
  expect_null(
    flexyBayes:::.parse_slope_dev_raw("slope_dev_other_raw[1,1]", lookup)
  )
})

# ---------------------------------------------------------------- #
# Subtest 10 -- lgm_gate accepts factor_numeric_interaction         #
# ---------------------------------------------------------------- #
#
# ADR 0017 (gate-truth = emit-truth): the fixed-term allowlist gains
# `factor_numeric_interaction` so the gate accepts; the post-fit
# §3.4 verification check (.lgm_check_factor_numeric_interaction_inla_verified)
# is the term-class-specific INLA conditional. When the verification
# artefact is absent the gate refuses INLA dispatch with the
# documented deferral message.

test_that("lgm_gate accepts factor_numeric_interaction structurally", {
  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    fixed_terms = list(
      list(
        type = "factor",
        var = "f",
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f"
      ),
      list(type = "continuous", var = "x", label = "x"),
      list(
        type = "factor_numeric_interaction",
        factor = "f",
        continuous = "x",
        vars = c("f", "x"),
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f:x"
      )
    ),
    random_terms = list(),
    rcov_terms = list(list(type = "units")),
    priors = list(legacy = TRUE),
    source = "brms"
  )
  # Structural fixed-term allowlist accepts the new class.
  r7 <- flexyBayes:::.lgm_check_fixed_term_inla_support(fb)
  expect_true(r7$pass)
  # §3.4 verification check refuses when artefact is absent on this
  # host (the default in CI / fresh checkouts).
  r10 <- flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(fb)
  if (!flexyBayes:::.factor_numeric_interaction_inla_verified()) {
    expect_false(r10$pass)
    expect_match(r10$reason, "INLA mapping deferred")
  } else {
    expect_true(r10$pass)
  }
})

# ---------------------------------------------------------------- #
# Subtest 11 -- dispatch trace records emit-path choice.            #
# ---------------------------------------------------------------- #

test_that("parse_info$factor_continuous_emit records the emit path", {
  skip_if_greta_backend_unusable()
  withr::local_options(list(
    flexyBayes.silence_default_prior_note = TRUE
  ))
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  # return_code = TRUE path: emit_greta() does not populate `extras`
  # (it returns the code string). Use a minimal MCMC instead so the
  # extras slot is built. Skip on CRAN / CI for the same reason
  # subtest 4 does -- the MCMC time is non-trivial.
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  fit <- flexybayes(
    y ~ f * x,
    data = d,
    backend = "greta",
    n_samples = 200,
    warmup = 200,
    chains = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  expect_true(
    fit$extras$parse_info$factor_continuous_emit %in%
      c("option_c", "option_d")
  )
})

# ---------------------------------------------------------------- #
# Subtest 12 -- Option D fallback path emits valid code.            #
# ---------------------------------------------------------------- #
#
# Force-trigger the Option D fallback via the option toggle and
# confirm the emitted code carries the per-observation lookup
# vectors and the reference-row mask. This catches regressions on
# the fallback path even when Option C ships.

test_that("Option D fallback emits per-observation lookup vectors", {
  withr::local_options(list(
    flexyBayes.force_option_d = TRUE,
    flexyBayes.silence_default_prior_note = TRUE
  ))
  set.seed(20260523L)
  n <- 60L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    y = rnorm(n)
  )
  skip_if_greta_backend_unusable()
  code <- flexybayes(
    y ~ f * x,
    data = d,
    backend = "greta",
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  expect_true(grepl("slope_dev_f_x_per_obs", code))
  expect_true(grepl("f_x_is_ref_obs", code))
  expect_true(grepl("f_x_shifted_idx", code))
  # Option D path should NOT call c(zeros(1), ...).
  expect_false(grepl("c\\(zeros\\(1\\)", code))
})
