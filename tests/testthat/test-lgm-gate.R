# Tests for lgm_gate() — deliverable 2.
#
# Covers the structural-first 6 checks, the structured refusal
# class, and the two-key override path. Numerical-confirm gate
# (check 7) is post-fit and stubbed for v0.1; tests for it land
# alongside emit_inla() (deliverable 5+).

# ---------------------------------------------------------------- #
# Test fixtures                                                    #
# ---------------------------------------------------------------- #

mk_fb <- function(
  family = "gaussian",
  fixed_terms = list(),
  random_terms = list(),
  rcov_terms = list(list(type = "units")),
  addition_terms = list(),
  priors = list(legacy = TRUE),
  source = "asreml"
) {
  flexyBayes:::new_fb_terms(
    response = "y",
    family = family,
    fixed_terms = fixed_terms,
    random_terms = random_terms,
    rcov_terms = rcov_terms,
    addition_terms = addition_terms,
    priors = priors,
    source = source
  )
}

# ---------------------------------------------------------------- #
# Check 1 — family allowlist                                       #
# ---------------------------------------------------------------- #

test_that(".lgm_check_family() passes a gaussian model", {
  fb <- mk_fb(family = "gaussian")
  r <- flexyBayes:::.lgm_check_family(fb)
  expect_true(r$pass)
  expect_identical(r$rule_id, "family_allowlist")
})

test_that(".lgm_check_family() passes binomial / poisson / gamma", {
  for (fam in c("binomial", "poisson", "gamma", "beta", "lognormal")) {
    fb <- mk_fb(family = fam)
    expect_true(
      flexyBayes:::.lgm_check_family(fb)$pass,
      info = paste("family =", fam)
    )
  }
})

test_that(".lgm_check_family() refuses an unsupported family", {
  fb <- mk_fb(family = "weibull5p_hypothetical")
  r <- flexyBayes:::.lgm_check_family(fb)
  expect_false(r$pass)
  expect_match(r$reason, "INLA likelihood allowlist")
  expect_match(r$diagnostic, "weibull5p_hypothetical")
})

# ---------------------------------------------------------------- #
# Check 2 — predictor linearity                                    #
# ---------------------------------------------------------------- #

test_that(".lgm_check_predictor() passes a linear asreml model", {
  fb <- mk_fb(fixed_terms = list(list(type = "factor", var = "env")))
  expect_true(flexyBayes:::.lgm_check_predictor(fb)$pass)
})

test_that(".lgm_check_predictor() refuses an explicit non_linear term", {
  fb <- mk_fb(
    fixed_terms = list(list(type = "non_linear", expr = "a*exp(-b*x)"))
  )
  r <- flexyBayes:::.lgm_check_predictor(fb)
  expect_false(r$pass)
  expect_match(r$reason, "non-linear predictor")
})

test_that(".lgm_check_predictor() refuses term tagged nl = TRUE", {
  fb <- mk_fb(
    fixed_terms = list(list(type = "expression", nl = TRUE, expr = "x^a"))
  )
  expect_false(flexyBayes:::.lgm_check_predictor(fb)$pass)
})

# ---------------------------------------------------------------- #
# Check 3 — distributional regression                              #
# ---------------------------------------------------------------- #

test_that(".lgm_check_distributional() passes when no dpar terms", {
  fb <- mk_fb(addition_terms = list(list(type = "weights", values = 1:3)))
  expect_true(flexyBayes:::.lgm_check_distributional(fb)$pass)
})

test_that(".lgm_check_distributional() refuses non-trivial sigma RHS", {
  fb <- mk_fb(
    addition_terms = list(
      list(type = "dpar_sigma", is_intercept_only = FALSE, rhs = "~ x")
    )
  )
  r <- flexyBayes:::.lgm_check_distributional(fb)
  expect_false(r$pass)
  expect_match(r$reason, "auxiliary parameter")
})

test_that(".lgm_check_distributional() passes intercept-only sigma RHS", {
  fb <- mk_fb(
    addition_terms = list(
      list(type = "dpar_sigma", is_intercept_only = TRUE, rhs = "~ 1")
    )
  )
  expect_true(flexyBayes:::.lgm_check_distributional(fb)$pass)
})

# ---------------------------------------------------------------- #
# Check 4 — RE Gaussian prior                                      #
# ---------------------------------------------------------------- #

test_that(".lgm_check_re_prior() passes legacy asreml prior", {
  fb <- mk_fb(priors = list(legacy = TRUE, fixed_sd = 10, vc_sd = 1))
  expect_true(flexyBayes:::.lgm_check_re_prior(fb)$pass)
})

test_that(".lgm_check_re_prior() passes when priors is NULL", {
  fb <- mk_fb(priors = NULL)
  expect_true(flexyBayes:::.lgm_check_re_prior(fb)$pass)
})

test_that(".lgm_check_re_prior() refuses non-Gaussian RE prior in fb_prior", {
  fb_prior <- structure(
    list(re = list(list(family = "horseshoe", target = "geno"))),
    class = "fb_prior"
  )
  fb <- mk_fb(priors = fb_prior)
  r <- flexyBayes:::.lgm_check_re_prior(fb)
  expect_false(r$pass)
  expect_match(r$reason, "non-Gaussian random-effect prior")
})

# ---------------------------------------------------------------- #
# Check 5 — latent-class detection                                 #
# ---------------------------------------------------------------- #

test_that(".lgm_check_latent_class() passes a clean LMM", {
  fb <- mk_fb(
    fixed_terms = list(list(type = "factor", var = "env")),
    random_terms = list(list(type = "simple", var = "geno"))
  )
  expect_true(flexyBayes:::.lgm_check_latent_class(fb)$pass)
})

test_that(".lgm_check_latent_class() refuses a mixture term", {
  fb <- mk_fb(random_terms = list(list(type = "mixture", k = 2)))
  r <- flexyBayes:::.lgm_check_latent_class(fb)
  expect_false(r$pass)
  expect_match(r$reason, "latent-class structure")
})

test_that(".lgm_check_latent_class() refuses a mixture family", {
  fb <- mk_fb(family = "mixture_gaussian_gaussian")
  r <- flexyBayes:::.lgm_check_latent_class(fb)
  expect_false(r$pass)
  expect_match(r$reason, "latent-class")
})

# ---------------------------------------------------------------- #
# Check 6 — hyperparameter budget                                  #
# ---------------------------------------------------------------- #

test_that(".lgm_check_hyperparam_budget() passes a small model", {
  fb <- mk_fb(
    fixed_terms = list(list(type = "factor", var = "env")),
    random_terms = list(
      list(type = "simple", var = "geno"),
      list(type = "simple", var = "site")
    )
  )
  r <- flexyBayes:::.lgm_check_hyperparam_budget(fb)
  expect_true(r$pass)
  expect_null(r$warning)
})

test_that(".lgm_check_hyperparam_budget() warns above the soft limit", {
  # gaussian (1) + 10 simple REs (10) = 11 hypers => warn
  re <- replicate(10, list(type = "simple", var = "g"), simplify = FALSE)
  fb <- mk_fb(random_terms = re)
  r <- flexyBayes:::.lgm_check_hyperparam_budget(fb)
  expect_true(r$pass)
  expect_match(r$warning, "exceeds soft limit")
})

test_that(".lgm_check_hyperparam_budget() refuses above the hard limit", {
  # gaussian (1) + 16 simple REs (16) = 17 hypers => fail
  re <- replicate(16, list(type = "simple", var = "g"), simplify = FALSE)
  fb <- mk_fb(random_terms = re)
  r <- flexyBayes:::.lgm_check_hyperparam_budget(fb)
  expect_false(r$pass)
  expect_match(r$reason, "exceeds hard limit")
})

test_that(".lgm_count_hyperparams() handles fa_gxe via 2k contribution", {
  fb <- mk_fb(random_terms = list(list(type = "fa_gxe", k = 3L)))
  # gaussian (1) + 2*3 = 7
  expect_identical(flexyBayes:::.lgm_count_hyperparams(fb), 7L)
})

test_that(".lgm_count_hyperparams() handles us_gxe via n*(n+1)/2", {
  fb <- mk_fb(random_terms = list(list(type = "us_gxe", n_outer = 4L)))
  # gaussian (1) + 4*5/2 = 11
  expect_identical(flexyBayes:::.lgm_count_hyperparams(fb), 11L)
})

# ---------------------------------------------------------------- #
# Checks 7-9 — INLA emit-support allowlists (ADR 0017)             #
# ---------------------------------------------------------------- #

# Positive cases: each new check passes on its allowlist.

test_that(".lgm_check_fixed_term_inla_support() passes the allowlist", {
  fb <- mk_fb(
    fixed_terms = list(
      list(type = "factor", var = "env"),
      list(type = "continuous", var = "x"),
      list(type = "interaction", vars = c("a", "b")),
      list(type = "factor_interaction", vars = c("g", "e")),
      list(type = "expression", label = "I(x^2)")
    )
  )
  r <- flexyBayes:::.lgm_check_fixed_term_inla_support(fb)
  expect_true(r$pass)
  expect_identical(r$rule_id, "fixed_term_type_inla")
  expect_null(r$reason)
})

test_that(".lgm_check_random_term_inla_support() passes the allowlist", {
  fb <- mk_fb(
    random_terms = list(
      list(type = "simple", var = "geno"),
      list(type = "ide", var = "site"),
      list(type = "id", var = "block"),
      list(type = "spline", var = "x")
    )
  )
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_true(r$pass)
  expect_identical(r$rule_id, "random_term_type_inla")
  expect_null(r$reason)
})

test_that(".lgm_check_rcov_term_inla_support() passes the allowlist", {
  fb <- mk_fb(rcov_terms = list(list(type = "units")))
  r <- flexyBayes:::.lgm_check_rcov_term_inla_support(fb)
  expect_true(r$pass)
  expect_identical(r$rule_id, "rcov_term_type_inla")
  expect_null(r$reason)
})

# Negative cases: each structured-covariance random-term class
# triggers the random_term_type_inla rule with the class label in
# the diagnostic. Six classes covered: vm / ped / at / us / fa /
# ar1 — the asreml-style structured-covariance vocabulary.

test_that(".lgm_check_random_term_inla_support() refuses vm (GBLUP / kinship)", {
  fb <- mk_fb(random_terms = list(list(type = "vm", var = "geno")))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "random_term_type_inla")
  expect_match(r$reason, "\"vm\"")
  expect_match(r$reason, "variance-matrix")
  expect_match(r$reason, "backend = \"greta\"")
})

test_that(".lgm_check_random_term_inla_support() refuses ped (pedigree)", {
  fb <- mk_fb(random_terms = list(list(type = "ped", var = "anim")))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "random_term_type_inla")
  expect_match(r$reason, "\"ped\"")
  expect_match(r$reason, "pedigree")
})

test_that(".lgm_check_random_term_inla_support() refuses at (heterogeneous)", {
  fb <- mk_fb(random_terms = list(list(type = "at", var = "env")))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "random_term_type_inla")
  expect_match(r$reason, "\"at\"")
  expect_match(r$reason, "heterogeneous")
})

test_that(".lgm_check_random_term_inla_support() refuses us (unstructured)", {
  fb <- mk_fb(random_terms = list(list(type = "us", var = "env")))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "random_term_type_inla")
  expect_match(r$reason, "\"us\"")
  expect_match(r$reason, "unstructured")
})

test_that(".lgm_check_random_term_inla_support() refuses fa (factor-analytic)", {
  fb <- mk_fb(random_terms = list(list(type = "fa", var = "geno", k = 2L)))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "random_term_type_inla")
  expect_match(r$reason, "\"fa\"")
  expect_match(r$reason, "factor-analytic")
})

test_that(".lgm_check_random_term_inla_support() refuses ar1 (autoregressive)", {
  fb <- mk_fb(random_terms = list(list(type = "ar1", var = "time")))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "random_term_type_inla")
  expect_match(r$reason, "\"ar1\"")
  expect_match(r$reason, "autoregressive lag-1")
})

# Rcov negative case: at_units triggers rcov_term_type_inla.

test_that(".lgm_check_rcov_term_inla_support() refuses at_units (heterogeneous residual)", {
  fb <- mk_fb(rcov_terms = list(list(type = "at_units", var = "env")))
  r <- flexyBayes:::.lgm_check_rcov_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "rcov_term_type_inla")
  expect_match(r$reason, "\"at_units\"")
  expect_match(r$reason, "heterogeneous residual")
})

# Integration: lgm_gate() routes a structured-cov refusal through
# the new rule and the print method surfaces the new rule_id.

test_that("lgm_gate() refuses a vm random term via random_term_type_inla", {
  fb <- mk_fb(random_terms = list(list(type = "vm", var = "geno")))
  out <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(out, "lgm_refusal")
  expect_identical(out$primary_rule, "random_term_type_inla")
  text <- capture.output(print(out))
  expect_true(any(grepl("\\[random_term_type_inla\\]", text)))
  expect_true(any(grepl("variance-matrix", text)))
})

test_that("lgm_gate() refuses an at_units rcov term via rcov_term_type_inla", {
  fb <- mk_fb(rcov_terms = list(list(type = "at_units", var = "env")))
  out <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(out, "lgm_refusal")
  expect_identical(out$primary_rule, "rcov_term_type_inla")
})

# ---------------------------------------------------------------- #
# Integration — lgm_gate                                           #
# ---------------------------------------------------------------- #

test_that("lgm_gate() passes a typical asreml mixed model", {
  fb <- mk_fb(
    fixed_terms = list(list(type = "factor", var = "env")),
    random_terms = list(list(type = "simple", var = "geno", var_n = 5L))
  )
  out <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(out, "fb_terms")
  expect_true("lgm_compatible" %in% out$capabilities)
  expect_false(flexyBayes:::is_lgm_refusal(out))
})

test_that("lgm_gate() returns lgm_refusal on family failure", {
  fb <- mk_fb(family = "weibull5p_hypothetical")
  out <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(out, "lgm_refusal")
  expect_true(flexyBayes:::is_lgm_refusal(out))
  expect_identical(out$primary_rule, "family_allowlist")
  expect_identical(out$n_failures, 1L)
})

test_that("lgm_gate() returns lgm_refusal on multiple failures", {
  fb <- mk_fb(
    family = "weibull5p_hypothetical",
    random_terms = list(list(type = "mixture", k = 2))
  )
  out <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(out, "lgm_refusal")
  rules <- vapply(out$failures, function(f) f$rule_id, character(1))
  expect_true("family_allowlist" %in% rules)
  expect_true("latent_class" %in% rules)
})

test_that("lgm_gate() passes warnings into capabilities (soft-limit budget)", {
  re <- replicate(10, list(type = "simple", var = "g"), simplify = FALSE)
  fb <- mk_fb(random_terms = re)
  out <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(out, "fb_terms")
  expect_true("lgm_compatible" %in% out$capabilities)
  expect_true(any(grepl("^lgm_warning:hyperparam_budget", out$capabilities)))
})

# ---------------------------------------------------------------- #
# Override path — two-key armed                                    #
# ---------------------------------------------------------------- #

test_that("lgm_gate() overrides structural refusal when fully armed", {
  fb <- mk_fb(family = "weibull5p_hypothetical")
  out <- flexyBayes:::lgm_gate(
    fb,
    force = "inla",
    acknowledge_silent_bias_risk = TRUE,
    reason = "research validation against parametric AFT"
  )
  expect_s3_class(out, "fb_terms")
  expect_true("lgm_force_overridden" %in% out$capabilities)
  expect_true(any(grepl("^lgm_force_reason:", out$capabilities)))
  expect_true(any(grepl(
    "^lgm_force_bypassed:family_allowlist",
    out$capabilities
  )))
})

test_that("lgm_gate() rejects override without reason", {
  fb <- mk_fb(family = "weibull5p_hypothetical")
  expect_error(
    flexyBayes:::lgm_gate(
      fb,
      force = "inla",
      acknowledge_silent_bias_risk = TRUE
    ),
    "`reason`.*is required"
  )
})

test_that("lgm_gate() ignores override flag without acknowledgement", {
  fb <- mk_fb(family = "weibull5p_hypothetical")
  out <- flexyBayes:::lgm_gate(fb, force = "inla", reason = "x")
  expect_s3_class(out, "lgm_refusal")
})

test_that("lgm_gate() ignores acknowledgement without force", {
  fb <- mk_fb(family = "weibull5p_hypothetical")
  out <- flexyBayes:::lgm_gate(
    fb,
    acknowledge_silent_bias_risk = TRUE,
    reason = "x"
  )
  expect_s3_class(out, "lgm_refusal")
})

test_that("lgm_gate() rejects non-fb_terms input", {
  expect_error(
    flexyBayes:::lgm_gate(list(response = "y")),
    "must be an fb_terms object"
  )
})

# ---------------------------------------------------------------- #
# Refusal print                                                    #
# ---------------------------------------------------------------- #

test_that("print.lgm_refusal() emits the structured refusal template", {
  fb <- mk_fb(
    family = "weibull5p_hypothetical",
    random_terms = list(list(type = "mixture", k = 2))
  )
  out <- flexyBayes:::lgm_gate(fb)
  text <- capture.output(print(out))
  expect_true(any(grepl("flexyBayes: INLA backend refused", text)))
  expect_true(any(grepl("\\[family_allowlist\\]", text)))
  expect_true(any(grepl("\\[latent_class\\]", text)))
  expect_true(any(grepl("Re-route", text)))
  expect_true(any(grepl("Override", text)))
})
