# Tests for emit_inla() -- the INLA backend (deliverable 5 partial).
#
# All tests skip when INLA is not installed; INLA is in
# DESCRIPTION:Suggests with Additional_repositories pointing at
# inla.r-inla-download.org/R/stable. Tests that exercise only the
# pure helpers (.build_inla_formula, .resolve_inla_family) run
# without INLA.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

mk_inla_data <- function() {
  set.seed(2026L)
  n <- 30L
  data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = factor(rep(1:5, length.out = n))
  )
}

mk_fb_for_inla <- function(
  data,
  random_terms = list(list(type = "simple", var = "g"))
) {
  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    link = "identity",
    fixed_terms = list(list(type = "continuous", var = "x")),
    random_terms = random_terms,
    rcov_terms = list(list(type = "units")),
    source = "asreml"
  )
  flexyBayes:::lgm_gate(fb)
}

# ---------------------------------------------------------------- #
# Pure helpers (run without INLA)                                  #
# ---------------------------------------------------------------- #

test_that(".build_inla_formula() builds intercept + fixed + simple RE", {
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)
  form <- flexyBayes:::.build_inla_formula(fb)
  expect_s3_class(form, "formula")
  fstr <- deparse(form)
  expect_match(fstr, "y ~ ")
  expect_match(fstr, "1")
  expect_match(fstr, "x")
  expect_match(fstr, 'f\\(g, model = "iid"\\)')
})

# ADR 0017: .build_inla_formula()'s switch() defaults are now
# internal contract-violation assertions, not user-facing
# refusals. The gate (.lgm_check_random_term_inla_support()) owns
# the user-facing refusal. The assertion fires only when the gate
# is bypassed via the two-key override (force = "inla" +
# acknowledge_silent_bias_risk = TRUE) on a structurally-
# incompatible IR -- i.e., the path that a forced-INLA fit on a
# truly incompatible model takes. The assertion message names
# the gate / emit drift explicitly so a triggered assertion is
# immediately recognisable as a flexyBayes-side bug.

test_that(".build_inla_formula() fires the gate-broken-contract assertion on an unsupported random term type when the gate is bypassed", {
  fb_raw <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    link = "identity",
    fixed_terms = list(list(type = "continuous", var = "x")),
    random_terms = list(list(
      type = "fa_gxe",
      outer = "env",
      inner = "geno",
      k = 2L
    )),
    rcov_terms = list(list(type = "units")),
    source = "asreml"
  )
  forced <- flexyBayes:::lgm_gate(
    fb_raw,
    force = "inla",
    acknowledge_silent_bias_risk = TRUE,
    reason = "ADR 0017 internal-assertion test"
  )
  expect_s3_class(forced, "fb_terms")
  expect_true("lgm_force_overridden" %in% forced$capabilities)
  expect_error(
    flexyBayes:::.build_inla_formula(forced),
    "lgm_gate broken contract"
  )
})

test_that(".build_inla_formula() fires the gate-broken-contract assertion on an unsupported fixed term type when the gate is bypassed", {
  fb_raw <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    link = "identity",
    fixed_terms = list(list(type = "non_linear", expr = "a*exp(-b*x)")),
    random_terms = list(),
    rcov_terms = list(list(type = "units")),
    source = "asreml"
  )
  forced <- flexyBayes:::lgm_gate(
    fb_raw,
    force = "inla",
    acknowledge_silent_bias_risk = TRUE,
    reason = "ADR 0017 internal-assertion test"
  )
  expect_s3_class(forced, "fb_terms")
  expect_error(
    flexyBayes:::.build_inla_formula(forced),
    "lgm_gate broken contract.*fixed term type"
  )
})

test_that(".build_inla_formula() fires the gate-broken-contract assertion on an unsupported rcov term type when the gate is bypassed", {
  fb_raw <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    link = "identity",
    fixed_terms = list(list(type = "continuous", var = "x")),
    random_terms = list(),
    rcov_terms = list(list(type = "at_units", var = "env")),
    source = "asreml"
  )
  forced <- flexyBayes:::lgm_gate(
    fb_raw,
    force = "inla",
    acknowledge_silent_bias_risk = TRUE,
    reason = "ADR 0017 internal-assertion test"
  )
  expect_s3_class(forced, "fb_terms")
  expect_error(
    flexyBayes:::.build_inla_formula(forced),
    "lgm_gate broken contract.*rcov term type"
  )
})

test_that(".resolve_inla_family() maps common families", {
  mk <- function(fam) {
    flexyBayes:::new_fb_terms(response = "y", family = fam, source = "asreml")
  }
  expect_identical(
    flexyBayes:::.resolve_inla_family(mk("gaussian")),
    "gaussian"
  )
  expect_identical(
    flexyBayes:::.resolve_inla_family(mk("binomial")),
    "binomial"
  )
  expect_identical(flexyBayes:::.resolve_inla_family(mk("poisson")), "poisson")
  expect_identical(
    flexyBayes:::.resolve_inla_family(mk("negative_binomial")),
    "nbinomial"
  )
  expect_identical(flexyBayes:::.resolve_inla_family(mk("gamma")), "gamma")
})

test_that(".lgm_check_numerical() flags non-zero mode.status", {
  fake_fit <- list(mode = list(mode.status = 1L), mlik = matrix(-100, 1, 1))
  res <- flexyBayes:::.lgm_check_numerical(fake_fit)
  expect_false(res$pass)
  expect_match(paste(res$reasons, collapse = " "), "mode\\.status = 1")
})

test_that(".lgm_check_numerical() flags non-finite mlik", {
  fake_fit <- list(mode = list(mode.status = 0L), mlik = matrix(NaN, 1, 1))
  res <- flexyBayes:::.lgm_check_numerical(fake_fit)
  expect_false(res$pass)
  expect_match(paste(res$reasons, collapse = " "), "mlik|finite")
})

test_that(".lgm_check_numerical() passes a clean fit", {
  fake_fit <- list(mode = list(mode.status = 0L), mlik = matrix(-50, 1, 1))
  res <- flexyBayes:::.lgm_check_numerical(fake_fit)
  expect_true(res$pass)
  expect_length(res$reasons, 0L)
})

# ---------------------------------------------------------------- #
# emit_inla() requires fb_terms                                    #
# ---------------------------------------------------------------- #

test_that("emit_inla() rejects non-fb_terms input", {
  expect_error(
    flexyBayes:::emit_inla(fb = list(response = "y"), data = mk_inla_data()),
    "must be an fb_terms object"
  )
})

# ---------------------------------------------------------------- #
# Live fits -- skip when INLA is not installed                     #
# ---------------------------------------------------------------- #

test_that("emit_inla() fits a simple gaussian + RE model", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)

  fit <- flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE)

  expect_s3_class(fit, "flexybayes_inla")
  expect_true(!is.null(fit$inla))
  expect_true(!is.null(fit$extras$summary$fixed))
  expect_true(nrow(fit$extras$summary$fixed) >= 1L)
  expect_true(
    isTRUE(fit$num_check$pass) ||
      isTRUE(fit$num_check$pass == FALSE)
  )
})

test_that("emit_inla() return_code returns formula and family without fitting", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)

  res <- flexyBayes:::emit_inla(
    fb = fb,
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_named(res, c("formula", "family", "hyper", "control_family"))
  expect_s3_class(res$formula, "formula")
  expect_identical(res$family, "gaussian")
  # control_family is empty when the user supplies no sigma prior;
  # populated when fb_prior(sigma ~ pc(...)) is passed.
  expect_identical(res$control_family, list())
})

test_that("fb(... backend = 'inla') dispatches to emit_inla", {
  # ADR 0004 D1: v0.1's flexybayes() (and fb alias) is greta-only.
  # User-facing INLA dispatch is deferred to v0.2 with the brms-
  # format ingest path (fb_brms()). The internal emit_inla() pathway
  # is exercised by the test above ("emit_inla() returns formula +
  # family + hyper + control_family lists").
  skip("INLA dispatch from user API deferred to v0.2 (ADR 0004 D1).")
  skip_if_no_inla()
  d <- mk_inla_data()
  fit <- flexybayes(fixed = y ~ x, random = ~g, data = d, verbose = FALSE)
  expect_s3_class(fit, "flexybayes_inla")
})

test_that("print.flexybayes_inla emits a brief summary", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)
  fit <- flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE)
  out <- capture.output(print(fit))
  expect_true(any(grepl("flexybayes_inla", out)))
  expect_true(any(grepl("formula", out)))
  expect_true(any(grepl("numerical confirm", out)))
})

test_that("summary.flexybayes_inla prints fixed effects and hyperparameters", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)
  fit <- flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE)
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Fixed effects", out)))
  expect_true(any(grepl("Hyperparameters", out)))
})
