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
  residual_terms = list(list(type = "units")),
  addition_terms = list(),
  priors = list(legacy = TRUE),
  source = "asreml"
) {
  flexyBayes:::new_fb_terms(
    response = "y",
    family = family,
    fixed_terms = fixed_terms,
    random_terms = random_terms,
    residual_terms = residual_terms,
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

test_that(".lgm_count_hyperparams() charges fa_gxe 2k hyperparameters", {
  fb <- mk_fb(random_terms = list(list(type = "fa_gxe", k = 3L)))
  # gaussian (1) + 2*3 = 7
  expect_identical(flexyBayes:::.lgm_count_hyperparams(fb), 7L)
})

test_that(".lgm_count_hyperparams() charges us_gxe n(n+1)/2 hyperparameters", {
  fb <- mk_fb(random_terms = list(list(type = "us_gxe", n_outer = 4L)))
  # gaussian (1) + 4*5/2 = 11
  expect_identical(flexyBayes:::.lgm_count_hyperparams(fb), 11L)
})

test_that(".inla_random_term_class_label() GxE labels use plain ASCII, not U+00D7", {
  # FB-16 (D3): the multiplication sign renders as U+FFFD under xelatex +
  # lmodern when these labels reach vignette/console output, and is
  # forbidden house style. Assert the literal ASCII spelling directly, not
  # just its absence, so a re-introduced Unicode glyph is caught by value.
  expect_identical(
    flexyBayes:::.inla_random_term_class_label("fa_gxe"), "factor-analytic GxE"
  )
  expect_identical(
    flexyBayes:::.inla_random_term_class_label("us_gxe"), "unstructured GxE"
  )
  expect_identical(
    flexyBayes:::.inla_random_term_class_label("corh_gxe"), "equicorrelated GxE"
  )
  for (lbl in c(
    flexyBayes:::.inla_random_term_class_label("fa_gxe"),
    flexyBayes:::.inla_random_term_class_label("us_gxe"),
    flexyBayes:::.inla_random_term_class_label("corh_gxe")
  )) {
    expect_false(grepl("×", lbl, fixed = TRUE))
  }
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

test_that(".lgm_check_residual_term_inla_support() passes the allowlist", {
  fb <- mk_fb(residual_terms = list(list(type = "units")))
  r <- flexyBayes:::.lgm_check_residual_term_inla_support(fb)
  expect_true(r$pass)
  expect_identical(r$rule_id, "residual_term_type_inla")
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
  # What the OTHER engine does is asked, not assumed. brms lowers a
  # dense-carrier vm() to (1 | gr(g, cov = K)), so the refusal has to say
  # so: it read "no other active backend represents it either" until the
  # 0.9.0 fidelity pass, which was false for every term brms had learned
  # to emit -- vm, us(f):g, and every nested random effect.
  expect_match(r$reason, "brms does represent this model", fixed = TRUE)
})

test_that(".lgm_check_random_term_inla_support() says so when brms cannot", {
  # 0.9.3: the message no longer names any specific engine as the
  # reason brms cannot represent this term either (a since-withdrawn
  # native engine used to be the exception the older wording carved
  # out) -- it states the fact plainly and engine-free: no active
  # backend can fit it.
  fb <- mk_fb(random_terms = list(list(type = "fa", var = "env", k = 2L)))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_false(r$pass)
  expect_match(r$reason, "brms cannot represent it either", fixed = TRUE)
  expect_match(r$reason, "no active backend can fit it faithfully", fixed = TRUE)
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

test_that(".lgm_check_random_term_inla_support() accepts ar1 (WP16 spatial field)", {
  # WP16 (2026-07-25): AR1 (1D) and separable AR1xAR1 spatial latent fields
  # are emitted on INLA (grouped-AR1 idiom); emit_inla() validates the
  # one-observation-per-node grid condition, so the gate accepts the class.
  fb <- mk_fb(random_terms = list(list(type = "ar1", var = "time")))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_true(r$pass)
})

test_that(".lgm_check_random_term_inla_support() accepts ar1_spatial (WP16)", {
  fb <- mk_fb(random_terms = list(list(
    type = "ar1_spatial", row_var = "row", col_var = "col", col_ar1 = TRUE
  )))
  r <- flexyBayes:::.lgm_check_random_term_inla_support(fb)
  expect_true(r$pass)
})

# Residual cases. at_units -- the sectioned residual, dsum(~units | f) --
# reaches INLA through the multiple-likelihood interface from 0.10.0.9000,
# so it is now an ACCEPTED type and the negative case moved to a residual
# shape INLA genuinely does not carry.

test_that(".lgm_check_residual_term_inla_support() accepts at_units", {
  fb <- mk_fb(residual_terms = list(list(type = "at_units", var = "env")))
  r <- flexyBayes:::.lgm_check_residual_term_inla_support(fb)
  expect_true(r$pass)
})

test_that(".lgm_check_residual_term_inla_support() refuses a structured residual", {
  fb <- mk_fb(residual_terms = list(list(
    type = "ar1_spatial", row_var = "row", col_var = "col"
  )))
  r <- flexyBayes:::.lgm_check_residual_term_inla_support(fb)
  expect_false(r$pass)
  expect_identical(r$rule_id, "residual_term_type_inla")
  expect_match(r$reason, "\"ar1_spatial\"")
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

test_that("lgm_gate() refuses a structured residual via residual_term_type_inla", {
  fb <- mk_fb(residual_terms = list(list(
    type = "ar1_spatial", row_var = "row", col_var = "col"
  )))
  out <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(out, "lgm_refusal")
  expect_identical(out$primary_rule, "residual_term_type_inla")
})

test_that("lgm_gate() accepts the sectioned residual it used to refuse", {
  fb <- mk_fb(residual_terms = list(list(type = "at_units", var = "env")))
  expect_false(inherits(flexyBayes:::lgm_gate(fb), "lgm_refusal"))
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

# ---------------------------------------------------------------- #
# Check 12 — joint shrinkage prior families                        #
# ---------------------------------------------------------------- #
#
# SPEC_bayesian_selection_v5_2026-09-12.md section 4.5. The check reads
# `fb$penalty_blocks`, which nothing in R/ writes: fb_select() fits its
# reference model on brms, so no public call reaches it. These tests set the
# field to exercise the check as an internal guard.

mk_penalty_block <- function(family, id = 1L) {
  list(
    id = id,
    family = family,
    target = "term(everything())",
    columns = 1:2,
    coef_names = c("x1", "x2"),
    protected = FALSE,
    hyperparameters = list(),
    hyperparameter_mode = "fixed",
    scaling = "none",
    route = if (family %in% c("regularized_horseshoe", "r2d2")) {
      "class_wide"
    } else {
      "coef_rows"
    }
  )
}

test_that(".lgm_check_prior_family_inla_support() passes with no blocks", {
  # Every model outside the selection module, and a selection whose penalty
  # is NULL. The default of NULL is what makes the check free for them.
  fb <- mk_fb()
  r <- flexyBayes:::.lgm_check_prior_family_inla_support(fb)
  expect_true(r$pass)
  expect_identical(r$rule_id, "prior_family_unsupported_inla")
  expect_null(r$reason)
})

test_that("a none() block is an opt-out of shrinkage and does not refuse", {
  fb <- mk_fb()
  fb$penalty_blocks <- list(mk_penalty_block("none"))
  expect_true(flexyBayes:::.lgm_check_prior_family_inla_support(fb)$pass)
})

test_that("every shrinkage family refuses at INLA, and the reason names it", {
  for (family in c("ridge", "laplace", "elastic_net",
                   "regularized_horseshoe", "r2d2")) {
    fb <- mk_fb()
    fb$penalty_blocks <- list(mk_penalty_block(family))
    r <- flexyBayes:::.lgm_check_prior_family_inla_support(fb)
    expect_false(r$pass, info = family)
    expect_match(r$reason, family, fixed = TRUE)
    expect_match(r$reason, "not fitted on INLA", fixed = TRUE)
    expect_match(r$diagnostic, family, fixed = TRUE)
  }
})

test_that("the check is wired into lgm_gate(), red on the block and green without", {
  # The injected defect is the block itself: the same model passes the gate
  # with the block removed, so the refusal is the prior's and not the
  # model's. Without this half, a check that refused every model would pass.
  fb <- mk_fb(random_terms = list(list(type = "simple", var = "g")))
  expect_false(flexyBayes:::is_lgm_refusal(flexyBayes:::lgm_gate(fb)))

  fb$penalty_blocks <- list(mk_penalty_block("regularized_horseshoe"))
  refused <- flexyBayes:::lgm_gate(fb)
  expect_true(flexyBayes:::is_lgm_refusal(refused))
  expect_true("prior_family_unsupported_inla" %in%
    vapply(refused$failures, function(f) f$rule_id, character(1L)))
})

test_that("the twelfth check is the twelfth: the gate runs it alongside the others", {
  # A check that is written but never called is a comment. The mixed case
  # asserts both: a model that fails on structure AND carries a block
  # reports both rule ids.
  fb <- mk_fb(family = "gaussian")
  fb$penalty_blocks <- list(mk_penalty_block("laplace"))
  refused <- flexyBayes:::lgm_gate(fb)
  expect_true(flexyBayes:::is_lgm_refusal(refused))
  ids <- vapply(refused$failures, function(f) f$rule_id, character(1L))
  expect_true("prior_family_unsupported_inla" %in% ids)
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
  expect_true(any(grepl("INLA cannot fit this model", text, fixed = TRUE)))
  expect_true(any(grepl("\\[family_allowlist\\]", text)))
  expect_true(any(grepl("\\[latent_class\\]", text)))
  # The refusal names one route onward and one place to read the rest.
  # It used to print four lines on what backend = "auto" already does, an
  # unconditional roadmap sentence about mixtures, and a two-argument
  # override incantation -- sixteen lines for a two-rule refusal, which is
  # what the R Journal review read.
  expect_true(any(grepl("backend = \"brms\"", text, fixed = TRUE)))
  expect_true(any(grepl("?lgm_gate", text, fixed = TRUE)))
  # The citation still belongs to the latent-class rule's own reason, which
  # this fixture triggers. What is gone is the unconditional advert that
  # printed it on every refusal, mixture or not.
  expect_false(any(grepl("on the roadmap", text, fixed = TRUE)))
  expect_lte(length(text), 8L)
})

# ---------------------------------------------------------------- #
# The gate's refusals carry a class (A5.5)                         #
# ---------------------------------------------------------------- #
#
# lgm_gate() returns an <lgm_refusal> value, and under an explicit
# backend = "inla" request dispatch turned it into a bare stop(). That
# left one refusal in the package a caller could not catch by name while
# every other refusal carried flexybayes_refusal_<code>. The value is
# unchanged -- rule ids, check-list and print method all stand -- and the
# throw site now wraps it in the shared condition hierarchy.

test_that("an explicit inla request refused by the gate raises a typed condition", {
  # The vehicle is a model the gate refuses, and it has to keep being one.
  # us(env):gen is refused on INLA at any order: estimated correlation
  # matrices are fitted on brms. Its rule names a refusal code of its own,
  # which the condition carries above the generic one.
  d <- data.frame(
    y = stats::rnorm(60L),
    x = stats::rnorm(60L),
    env = factor(rep(letters[1:6], each = 10L)),
    gen = factor(rep(seq_len(10L), times = 6L))
  )

  err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = y ~ x,
      random = ~ us(env):gen,
      data = d,
      backend = "inla",
      verbose = FALSE
    )),
    condition = function(e) e
  )

  expect_s3_class(err, "flexybayes_refusal_correlated_covariance_requires_brms")
  expect_identical(err$reason_code, "correlated_covariance_requires_brms")
  # The generic code stays beneath it, so a handler written against it
  # still catches this refusal.
  expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
  expect_s3_class(err, "flexybayes_refusal")
  expect_s3_class(err, "error")

  # Every existing gate message survives verbatim inside the condition.
  expect_match(conditionMessage(err), "INLA cannot fit this model", fixed = TRUE)
  expect_match(conditionMessage(err), "backend = \"brms\"", fixed = TRUE)
  expect_match(conditionMessage(err), "?lgm_gate", fixed = TRUE)

  # The check-list object itself rides on the condition, so a caller can
  # read the rules rather than parse the printed text.
  expect_s3_class(err$lgm_refusal, "lgm_refusal")
  expect_true(err$lgm_refusal$n_failures >= 1L)
  expect_identical(err$rule_id, err$lgm_refusal$primary_rule)
})

test_that("a gate refusal through the fb_inla() pin is typed and names its rule", {
  # Same throw site, reached through the engine pin rather than through
  # backend = "inla", and carrying a DIFFERENT rule from the test above so
  # the per-rule class is shown to track the rule rather than the site.
  #
  # The original vehicle here was the sectioned residual dsum(~units | f),
  # refused by residual_term_type_inla. That residual now fits. So does
  # every other residual spelling that reaches the gate: the structured
  # ones are refused earlier, at dispatch, and `units` and `at_units` both
  # pass. residual_term_type_inla is therefore no longer reachable from the
  # formula surface -- it remains as a defence, and the unit test above
  # exercises it directly on a hand-built IR.
  #
  # The second vehicle was a six-level us(), refused by the budget; us() is
  # now refused first by the correlated-covariance rule. One variance per
  # level over twenty levels is 21 hyperparameters, refused by the budget
  # alone.
  d <- data.frame(
    y = stats::rnorm(200L),
    x = stats::rnorm(200L),
    env = factor(rep(sprintf("e%02d", seq_len(20L)), each = 10L)),
    gen = factor(rep(seq_len(10L), times = 20L))
  )

  err <- tryCatch(
    suppressMessages(fb_inla(
      fixed = y ~ x,
      random = ~ diag(env):gen,
      data = d,
      verbose = FALSE
    )),
    condition = function(e) e
  )

  expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
  expect_identical(err$rule_id, "hyperparam_budget")
  # The per-rule class lets one rule be caught without catching the rest.
  expect_s3_class(err, "flexybayes_lgm_hyperparam_budget")
})

# ---------------------------------------------------------------- #
# D11 / WP-D: a vignette("...") reference in R/ must exist          #
# ---------------------------------------------------------------- #
# lgm_gate.R's own refusal-help banner used to name
# "flexyBayes-11-dispatch-and-refusals" after the 0.9.3 renumber moved
# that page to -10-; the string was stale and nothing caught it. This
# test scans every R/ source for a flexyBayes-<NN>-<slug> filename
# token and asserts each one names a real vignettes/*.Rmd.orig file, so
# a future renumber that misses one message fails CI instead of
# shipping a broken vignette("...") pointer.

.vignette_source_dir_d11 <- function() {
  cand <- testthat::test_path("..", "..", "vignettes")
  if (dir.exists(cand)) {
    return(normalizePath(cand))
  }
  NULL
}

test_that("every flexyBayes-<NN>-<slug> filename named in R/ exists under vignettes/", {
  vig_dir <- .vignette_source_dir_d11()
  skip_if(is.null(vig_dir), "vignette sources not present in this layout")

  r_dir <- testthat::test_path("..", "..", "R")
  skip_if(!dir.exists(r_dir), "R/ sources not present in this layout")
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  expect_gt(length(r_files), 0L)

  real_pages <- sub(
    "\\.Rmd\\.orig$", "",
    list.files(vig_dir, pattern = "^flexyBayes-[0-9]{2}-.*\\.Rmd\\.orig$")
  )
  expect_gt(length(real_pages), 0L)

  referenced <- character(0)
  for (f in r_files) {
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    hits <- regmatches(
      txt,
      gregexpr("flexyBayes-[0-9]{2}-[a-z-]+[a-z]", txt)
    )[[1L]]
    referenced <- c(referenced, hits)
  }
  referenced <- unique(referenced)
  expect_gt(length(referenced), 0L) # the gate banner itself, at minimum

  missing <- setdiff(referenced, real_pages)
  expect_identical(missing, character(0))
})
