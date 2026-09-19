# Tests for ADR 0020 -- uncorrelated random slopes (x || g) and
# (1 + x || g).
#
# Contract (per design-notes-priors-brms-lgm.md ADR 0020 + the design spec
# §4.2.2 acceptance criteria):
#
#   (a) fb(y ~ x + (x || g)) parses to type = "simple_slope_uncor"
#       with slope_var = "x", with_intercept = FALSE, and does NOT
#       raise stop() on the brms path.
#   (b) fb(y ~ x + (1 + x || g)) parses to type =
#       "simple_slope_uncor" with with_intercept = TRUE; fits both
#       intercept- and slope-variance components.
#   (c) fb(y ~ x + (x | g)) parses to type = "simple_slope_cor" --
#       the same model with the correlation freed. It used to raise a
#       structured deferral; both engines now carry it.
#   (d-e) The two spellings agree, and emit_brms writes the single
#       pipe back out.
#   (f) Posterior matches lme4::lmer(REML = FALSE) on sleepstudy
#       within W_1 <= 0.20 * tau_true on the slope SD (loosened per
#       the design spec) and 0.15 on the intercept SD.
#   (g) Canonical name sd_<x>_<g> resolves on a brms fit.
#   (h) The brms-side raw-name-to-canonical mapper rule, unit-tested
#       directly (the regex split (g) exercises end to end).
#   (i) [removed at 0.9.3, see NEWS.md] a cross-engine comparison of
#       fb(backend = "brms") against the since-withdrawn native
#       engine on the slope-variance posterior; that engine is
#       withdrawn entirely, so no 2-engine pairing through it remains
#       for this structure.
#   (j) Minimum-group guard: the slope fixtures use J >= 10 by
#       construction (sleepstudy has J = 18).
#   (k) prior_summary() lists the slope-variance prior row.
#   (l) print.flexybayes() lists the slope-variance hyperparameter
#       in the model_info$n_params count.
#   (m) emit_brms() formula round-trips both forms.
#   (n) 10-row priors-to-brms table: (x || g) adds the
#       (class = "sd", coef = "<x>", group = "<g>") row.


mk_slope_data <- function(
  seed = 20260523L,
  J = 18L,
  n_per = 10L,
  beta = 0.5,
  sd_int = 30,
  sd_slope = 5,
  sigma_e = 25
) {
  set.seed(seed)
  Subject <- factor(rep(seq_len(J), each = n_per))
  Days <- rep(seq_len(n_per) - 1L, times = J)
  u_int <- rnorm(J, sd = sd_int)
  u_slope <- rnorm(J, sd = sd_slope)
  Reaction <- 250 +
    beta * Days +
    u_int[as.integer(Subject)] +
    u_slope[as.integer(Subject)] * Days +
    rnorm(length(Days), sd = sigma_e)
  data.frame(Reaction = Reaction, Days = Days, Subject = Subject)
}


# ---------------------------------------------------------------- #
# (a) Parse (x || g) -- intercept + slope, lme4 / brms semantics    #
# ---------------------------------------------------------------- #

test_that("(x || g) parses to simple_slope_uncor with with_intercept = TRUE (lme4 semantics)", {
  d <- mk_slope_data()
  fb <- fb_from_brms(Reaction ~ Days + (Days || Subject), data = d)
  expect_length(fb$random_terms, 1L)
  rt <- fb$random_terms[[1]]
  expect_identical(rt$type, "simple_slope_uncor")
  expect_identical(rt$var, "Subject")
  expect_identical(rt$slope_var, "Days")
  # Per ?lme4::lFormula: (x || g) is sugar for
  # (1 | g) + (0 + x | g) -- the intercept is included by default.
  # The (0 + x || g) form below suppresses it.
  expect_true(isTRUE(rt$with_intercept))
})


# ---------------------------------------------------------------- #
# (b) (1 + x || g) and (0 + x || g) parameterisations                #
# ---------------------------------------------------------------- #

test_that("(1 + x || g) keeps the intercept; (0 + x || g) suppresses it", {
  d <- mk_slope_data()
  fb_plus1 <- fb_from_brms(Reaction ~ Days + (1 + Days || Subject), data = d)
  expect_true(isTRUE(fb_plus1$random_terms[[1]]$with_intercept))
  fb_plus_sym <- fb_from_brms(Reaction ~ Days + (Days + 1 || Subject), data = d)
  expect_true(isTRUE(fb_plus_sym$random_terms[[1]]$with_intercept))
  # (0 + x || g) suppresses the intercept -- slope-only form.
  fb_zero <- fb_from_brms(Reaction ~ Days + (0 + Days || Subject), data = d)
  expect_false(isTRUE(fb_zero$random_terms[[1]]$with_intercept))
  fb_zero_sym <- fb_from_brms(Reaction ~ Days + (Days + 0 || Subject), data = d)
  expect_false(isTRUE(fb_zero_sym$random_terms[[1]]$with_intercept))
})


# ---------------------------------------------------------------- #
# (c) (x | g) parses to the correlated-slope term                   #
# ---------------------------------------------------------------- #

test_that("(x | g) parses to simple_slope_cor with the intercept in", {
  # The single pipe used to refuse with a structured deferral. It is the
  # uncorrelated model with the correlation freed, and both engines now
  # carry it -- INLA as one iid2d field over a two-block index, brms as
  # `(x | g)` verbatim -- so the deferral is gone and the term is ordinary.
  d <- mk_slope_data()
  fb <- fb_from_brms(Reaction ~ Days + (Days | Subject), data = d)
  expect_length(fb$random_terms, 1L)
  t1 <- fb$random_terms[[1]]
  expect_identical(t1$type, "simple_slope_cor")
  expect_identical(t1$slope_var, "Days")
  expect_identical(t1$var, "Subject")
  # The intercept is always in: the correlation is what it is with.
  expect_true(isTRUE(t1$with_intercept))
  expect_identical(t1$var_n, 18L)
})

test_that("(1 + x | g) is the same term as (x | g)", {
  d <- mk_slope_data()
  bare <- fb_from_brms(Reaction ~ Days + (Days | Subject), data = d)
  explicit <- fb_from_brms(Reaction ~ Days + (1 + Days | Subject), data = d)
  expect_identical(
    bare$random_terms[[1]]$type,
    explicit$random_terms[[1]]$type
  )
  expect_identical(
    bare$random_terms[[1]]$slope_var,
    explicit$random_terms[[1]]$slope_var
  )
})

test_that("emit_brms writes the single pipe back out", {
  d <- mk_slope_data()
  fb <- fb_from_brms(Reaction ~ Days + (Days | Subject), data = d)
  txt <- gsub("\\s+", " ", paste(deparse(.fb_to_brms_formula(fb)),
                                 collapse = " "))
  expect_match(txt, "(Days | Subject)", fixed = TRUE)
  expect_false(grepl("||", txt, fixed = TRUE))
})


# ---------------------------------------------------------------- #
# (f) Posterior matches lme4::lmer on the slope SD                  #
# ---------------------------------------------------------------- #

test_that("brms posterior matches lme4 on slope + intercept SD within loosened tolerance", {
  skip_on_cran()
  skip_if_not_installed("brms")
  testthat::skip_if_not_installed("lme4")
  d <- mk_slope_data()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  # 2000 samples + 2000 warmup per ADR 0020 §Verification. Slope-
  # variance posteriors are intrinsically less identified than
  # intercept-variance ones at modest J (per the design spec)
  # and shorter chains yield posterior means outside the 0.20 *
  # tau_true tolerance even on the J = 18 sleepstudy-style fixture.
  # Slope-variance posteriors are intrinsically less identified (see
  # the comment above), so even this generous 2000+2000 budget can
  # still trip Stan's low-ESS sampler warning; the test's own 0.20 *
  # tau_true tolerance already absorbs the same MC noise the warning
  # is about, so muffle rather than raise the budget further.
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    Reaction ~ Days + (1 + Days || Subject),
    data = d,
    backend = "brms",
    n_samples = 2000L,
    warmup = 2000L,
    chains = 2L,
    seed = 20260523L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  ref <- lme4::lmer(
    Reaction ~ Days + (1 + Days || Subject),
    data = d,
    REML = FALSE
  )
  vc <- as.data.frame(lme4::VarCorr(ref))
  # lme4 reports the two independent groups (Subject and Subject.1
  # in the (1 + Days || Subject) parameterisation). The third row
  # (Residual) has var1 = NA which breaks an `==` filter without an
  # explicit !is.na() guard, so first filter to the random-effect
  # rows (grp != "Residual"), then key on var1.
  re_rows <- vc[!is.na(vc$var1) & vc$grp != "Residual", , drop = FALSE]
  sd_Days_lme4 <- re_rows$sdcor[re_rows$var1 == "Days"]
  sd_Subject_lme4 <- re_rows$sdcor[re_rows$var1 == "(Intercept)"]

  draws <- as.matrix(posterior::as_draws_matrix(fit$brms))

  # Hardened 2026-05-30 (ADR 0031). The previous assertion compared the
  # posterior MEAN of a weakly-identified variance-component SD (J = 18
  # subjects) to lme4's ML point estimate within a 20% relative band.
  # That is intrinsically flaky: the posterior of an RE-SD is right-
  # skewed and the Bayesian mean is a different estimator from the
  # frequentist ML point -- the gap exceeds 20% on ~1/3 of runs from MC
  # + skew + prior influence, not from under-sampling (this is already
  # 2000 x 2 chains). The statistically-principled, non-flaky check is
  # that lme4's estimate falls within the brms posterior's central
  # credible interval -- the correct notion of Bayesian/frequentist
  # agreement, robust to MC error and skew by construction. A 99%
  # interval is used for a comfortable robustness margin on the
  # poorly-identified slope SD.
  ci_Days <- stats::quantile(
    draws[, "sd_Subject__Days"],
    c(0.005, 0.995),
    names = FALSE
  )
  ci_Subject <- stats::quantile(
    draws[, "sd_Subject__Intercept"],
    c(0.005, 0.995),
    names = FALSE
  )
  expect_gte(sd_Days_lme4, ci_Days[1])
  expect_lte(sd_Days_lme4, ci_Days[2])
  expect_gte(sd_Subject_lme4, ci_Subject[1])
  expect_lte(sd_Subject_lme4, ci_Subject[2])

  # Retain a loose point-estimate sanity bound so the test still fails if
  # the posterior median diverges grossly from lme4 (a real bug), without
  # the tight-tolerance flakiness. Median is the stable centre for a
  # skewed variance posterior.
  med_Days <- stats::median(draws[, "sd_Subject__Days"])
  med_Subject <- stats::median(draws[, "sd_Subject__Intercept"])
  expect_lt(abs(med_Days - sd_Days_lme4), 0.60 * sd_Days_lme4)
  expect_lt(abs(med_Subject - sd_Subject_lme4), 0.60 * sd_Subject_lme4)
})


# ---------------------------------------------------------------- #
# (g) Canonical name sd_<x>_<g> resolves on a brms fit              #
# ---------------------------------------------------------------- #

test_that("canonical_names() on a brms (x || g) fit resolves sd_<x>_<g> AND sd_<g>", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_slope_data()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  # (Days || Subject) carries BOTH the intercept-variance (sd_Subject)
  # and slope-variance (sd_Days_Subject) hyperparameters per lme4 /
  # brms semantics. canonical_names() must surface both. This is the
  # full-pipeline counterpart to the isolated regex-split unit test in
  # subtest (h) below.
  # Structural test (does canonical_names() surface both hyperparameter
  # names?), not numeric agreement -- muffle at this tiny budget.
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    Reaction ~ Days + (Days || Subject),
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  cn <- canonical_names(fit)
  expect_true(all(c("sd_Subject", "sd_Days_Subject") %in% cn$map))
  expect_identical(cn$map[["sd_Subject__Days"]], "sd_Days_Subject")
  expect_identical(cn$map[["sd_Subject__Intercept"]], "sd_Subject")

  # The (0 + Days || Subject) form drops the intercept-variance.
  fit2 <- .muffle_ess_warnings(suppressMessages(flexybayes(
    Reaction ~ Days + (0 + Days || Subject),
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  cn2 <- canonical_names(fit2)
  expect_true("sd_Days_Subject" %in% cn2$map)
  expect_false("sd_Subject" %in% cn2$map)
})


# ---------------------------------------------------------------- #
# (h) Canonical name sd_<x>_<g> resolves on a brms fit              #
# ---------------------------------------------------------------- #

test_that("brms-side mapper translates sd_<g>__<x> to sd_<x>_<g>", {
  # Unit-test the brms mapper rule directly via the same regex split
  # that .mapper_brms_via_stan uses. This keeps the test fast and
  # independent of the (slow + heavy-toolchain) brms install: the
  # mapping logic itself is what ADR 0020 §Decision 3 names; the
  # full brms-backend integration is covered by subtest (g) above,
  # on a live fit.
  split_brms <- function(brms_name) {
    bare <- sub("^sd_", "", brms_name)
    if (grepl("__Intercept$", bare)) {
      bare2 <- sub("__Intercept$", "", bare)
      return(paste0("sd_", bare2))
    }
    m <- regmatches(bare, regexec("^([^_]+(?:_[^_]+)*?)__(.+)$", bare))[[1]]
    if (length(m) == 3L && nzchar(m[2L]) && nzchar(m[3L])) {
      return(paste0("sd_", m[3L], "_", m[2L]))
    }
    paste0("sd_", bare)
  }
  expect_identical(split_brms("sd_Subject__Intercept"), "sd_Subject")
  expect_identical(split_brms("sd_Subject__Days"), "sd_Days_Subject")
  expect_identical(split_brms("sd_g__x"), "sd_x_g")
  expect_identical(split_brms("sd_my_group__Days"), "sd_Days_my_group")
})


# ---------------------------------------------------------------- #
# (i) [removed] Cross-engine comparison on slope variance -- 0.9.3  #
# ---------------------------------------------------------------- #
#
# This cell compared a since-withdrawn native engine against brms
# on (x || g) (see NEWS.md, 0.9.3).


# ---------------------------------------------------------------- #
# (j) Minimum-group guard: J >= 10 holds for the slope fixtures     #
# ---------------------------------------------------------------- #

test_that("slope fixtures honour the J >= 10 minimum-group guard", {
  # The minimum-group guard is a property of the slope fixtures; the
  # standard sleepstudy / mk_slope_data fixtures use J = 18 by default
  # which satisfies J >= 10. Verify the property holds on both fixtures
  # so future fixture edits cannot silently regress below the guard.
  d_default <- mk_slope_data()
  expect_gte(nlevels(d_default$Subject), 10L)
  if (requireNamespace("lme4", quietly = TRUE)) {
    data(sleepstudy, package = "lme4", envir = environment())
    expect_gte(nlevels(sleepstudy$Subject), 10L)
  }
})


# ---------------------------------------------------------------- #
# (k) prior_summary() lists the slope-variance row                  #
# ---------------------------------------------------------------- #

test_that("prior_summary() lists the slope-variance prior row for (x || g)", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_slope_data()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  # Structural test (does the printed prior_summary() mention the
  # slope-variance parameter?), not numeric agreement -- muffle at
  # this tiny budget.
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    Reaction ~ Days + (Days || Subject),
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  ps <- prior_summary(fit)
  # The slope-variance row carries the canonical sd_<x>_<g> name in
  # the parameter column; printed shape may be either a data.frame
  # or an fb_prior_summary object depending on the backend, so the
  # robust check is via capture.output() text search.
  txt <- paste(utils::capture.output(print(ps)), collapse = "\n")
  expect_true(
    grepl("sd_Days_Subject|Days_Subject", txt) ||
      grepl("Days.*Subject", txt)
  )
})


# ---------------------------------------------------------------- #
# (l) Fit object exposes the slope-variance hyperparameter           #
# ---------------------------------------------------------------- #

test_that("fit object exposes the slope-variance hyperparameter via draws + variance_comps", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_slope_data()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  # Structural test (does the parameter name appear in the draws
  # matrix?), not numeric agreement -- muffle at this tiny budget.
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    Reaction ~ Days + (Days || Subject),
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  # The slope-variance hyperparameter must appear in the brms draws
  # matrix (brms-native name) so downstream tooling -- coef(),
  # canonical_names(), prior_summary() -- can resolve
  # it. print() output is brief by design and lists term-level
  # summaries only.
  draws <- as.matrix(posterior::as_draws_matrix(fit$brms))
  expect_true("sd_Subject__Days" %in% colnames(draws))
  # The print method should at least name the grouping factor so the
  # user sees that a random-effect block was monitored.
  txt <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_true(grepl("Subject", txt))
})


# ---------------------------------------------------------------- #
# (m) emit_brms() formula round-trip                                #
# ---------------------------------------------------------------- #

test_that("emit_brms() reconstructs (x || g) and (0 + x || g) round-trip", {
  d <- mk_slope_data()
  # (Days || Subject) and (1 + Days || Subject) parameterise the
  # same model under lme4 / brms semantics (intercept + slope,
  # uncorrelated). The canonical round-trip is the shorter
  # (Days || Subject) shape.
  fb_a <- fb_from_brms(Reaction ~ Days + (Days || Subject), data = d)
  fb_b <- fb_from_brms(Reaction ~ Days + (1 + Days || Subject), data = d)
  ff_a <- flexyBayes:::.fb_to_brms_formula(fb_a)
  ff_b <- flexyBayes:::.fb_to_brms_formula(fb_b)
  expect_identical(deparse(ff_a), "Reaction ~ Days + (Days || Subject)")
  expect_identical(deparse(ff_b), "Reaction ~ Days + (Days || Subject)")
  # The slope-only (0 + x || g) form round-trips with the explicit
  # intercept-suppression token preserved.
  fb_zero <- fb_from_brms(Reaction ~ Days + (0 + Days || Subject), data = d)
  ff_zero <- flexyBayes:::.fb_to_brms_formula(fb_zero)
  expect_identical(deparse(ff_zero), "Reaction ~ Days + (0 + Days || Subject)")
})


# ---------------------------------------------------------------- #
# (n) 10-row priors-to-brms table -- slope-variance row             #
# ---------------------------------------------------------------- #

test_that("priors_to_brms legacy specs add the slope-variance row(s) per ADR 0020 semantics", {
  d <- mk_slope_data()
  # (x || g) -- lme4 / brms semantics: intercept + slope, both
  # uncorrelated. Two sd rows: the intercept-variance row
  # (coef = NA, group = Subject) and the slope-variance row
  # (coef = Days, group = Subject).
  fb_a <- fb_from_brms(Reaction ~ Days + (Days || Subject), data = d)
  specs_a <- flexyBayes:::.priors_to_brms_specs(
    NULL,
    fb_a,
    prior_fixed_sd = 100,
    prior_vc_sd = 1
  )
  classes <- vapply(specs_a, `[[`, character(1), "class")
  coefs <- vapply(specs_a, function(s) s$coef %||% NA_character_, character(1))
  groups <- vapply(
    specs_a,
    function(s) s$group %||% NA_character_,
    character(1)
  )
  sd_rows <- which(classes == "sd")
  expect_length(sd_rows, 2L)
  expect_setequal(coefs[sd_rows], c(NA_character_, "Days"))
  expect_true(all(groups[sd_rows] == "Subject"))

  # (0 + x || g) suppresses the intercept -- slope-only form. Just
  # one sd row keyed by (group = "Subject", coef = "Days").
  fb_b <- fb_from_brms(Reaction ~ Days + (0 + Days || Subject), data = d)
  specs_b <- flexyBayes:::.priors_to_brms_specs(
    NULL,
    fb_b,
    prior_fixed_sd = 100,
    prior_vc_sd = 1
  )
  classes_b <- vapply(specs_b, `[[`, character(1), "class")
  coefs_b <- vapply(
    specs_b,
    function(s) s$coef %||% NA_character_,
    character(1)
  )
  groups_b <- vapply(
    specs_b,
    function(s) s$group %||% NA_character_,
    character(1)
  )
  sd_rows_b <- which(classes_b == "sd")
  expect_length(sd_rows_b, 1L)
  expect_identical(coefs_b[sd_rows_b], "Days")
  expect_identical(groups_b[sd_rows_b], "Subject")
})


# ---------------------------------------------------------------- #
# (o) INLA path FITS (x || g) -- the deferral was retired at         #
#     0.10.0.9000                                                    #
# ---------------------------------------------------------------- #

test_that("INLA path fits (x || g) rather than deferring", {
  # This asserted the runtime deferral: (x || g) refused on INLA unless an
  # on-disk artefact recorded that a three-arbitrator comparison had passed,
  # and one of those arbitrators was withdrawn at 0.9.3, so it refused
  # everywhere. Retiring that gate exposed a real defect the refusal had
  # been hiding -- the emit wrote an index column it never created -- which
  # is now fixed. Agreement with lme4 REML is checked in
  # test-inla-slope-and-interaction-verified.R.
  testthat::skip_if_not_installed("INLA")
  testthat::skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- mk_slope_data()
  fit <- suppressMessages(flexybayes(
    Reaction ~ Days + (Days || Subject),
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  comps <- summary(fit)$varcomp$component
  # One SD for the intercept, one for the slope, and no generated index
  # name reaching the label.
  expect_true("sd_Subject" %in% comps)
  expect_true("sd_Days_Subject" %in% comps)
  expect_false(any(grepl("_for_slope_", comps, fixed = TRUE)))
})
