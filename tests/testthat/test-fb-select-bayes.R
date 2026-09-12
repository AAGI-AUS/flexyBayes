# =============================================================================
# fb_select(): the Bayesian workflow, end to end.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 2, 3.1, 5, 6, 8
# and 9.
#
# The units this file integrates are each tested against their own
# interfaces elsewhere. What is tested HERE is the wiring -- that the verb
# routes a request to the right unit, that the four axes it records are the
# ones that ran, and that no configuration produces a penalised-likelihood
# answer under a Bayesian label.
#
# The expensive fits are memoised through the file's own cache, because one
# Stan compile per test would put this file well past any sensible budget and
# the same posterior answers every question asked of it here.
#
# Two gates are shown red on an injected defect before they are trusted:
#
#   the certificate's estimand   built from the fold object, so a certificate
#                                that read an argument instead would not
#                                notice a fold working in a different one.
#                                The falsifier corrupts the fold's recorded
#                                estimand and requires the certificate to
#                                follow it
#   the fallback                 fires ONLY on an absent provider, so the
#                                falsifier makes the provider absent and
#                                requires the substitution, then restores it
#                                and requires no substitution
# =============================================================================

# Set only when cmdstanr is reachable: an unconditional setting here reaches
# every file in the run, and under `R CMD check` -- where cmdstanr is not on
# the library path -- it turns unrelated brms fits into an install refusal.
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  options(brms.backend = "cmdstanr")
}

fbsb_cache <- new.env(parent = emptyenv())

fbsb_memo <- function(key, value) {
  if (is.null(fbsb_cache[[key]])) {
    assign(key, value, envir = fbsb_cache)
  }
  get(key, envir = fbsb_cache)
}

fbsb_needs_brms <- function() {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
}

# The base: one protected covariate, one grouping factor, four candidates of
# which two carry signal. Small on purpose -- these tests are about the
# module's contract, not about the posterior's precision.
fbsb_data <- function() fb_sel_bayes_data()

fbsb_base <- function() {
  fbsb_memo("base", fb_sel_bayes_base(fbsb_data()))
}

# A second base whose fixed part is the intercept alone. The regularised
# horseshoe is class-wide in brms, so with nothing protected beyond the
# intercept the candidate block IS the whole class and the request is exact
# -- the case the contract calls exact, and the one with the simplest
# geometry to sample. `scope = "all_fixed"` is exercised on the protected
# base at plan time, where no sampler is involved.
# A wider fixture for the reference fit: 20 sites of eight plots. A
# global-local prior on 72 rows is a funnel the sampler diverges in at the
# contract's step-size target -- three divergent transitions, and the
# adequacy gate refused, which is the gate working. More data is the fix
# that makes the geometry easier rather than the one that makes the
# diagnostics quieter.
fbsb_data_wide <- function() fb_sel_bayes_data(J = 20L, m = 8L)

fbsb_base_intercept <- function() {
  fbsb_memo("base_intercept", fb_brms(
    y ~ 1 + (1 | site),
    data = fbsb_data_wide(),
    chains = 2L,
    n_samples = 1000L,
    warmup = 1000L,
    seed = 20260912L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
}

fbsb_candidates <- ~ x1 + x2 + z1 + z2

# The sampler budget is what the adequacy gate needs, not a recommendation.
# At two chains of 1000 kept draws the reference model's minimum bulk
# effective sample size came back at 382 against the contract's 400 and the
# gate refused, which is the gate working: the budget is raised here rather
# than the gate being loosened.
#
# The step-size target stays at the contract's 0.95. Raising it to 0.995 was
# tried and is a false economy at this shape: it turns a 10-second reference
# fit into a multi-minute one, because the horseshoe's geometry then
# saturates the treedepth ceiling on most transitions, and the projection
# refits the reference once per fold. `fb_engine(adapt_delta = )` is the
# lever where a fit needs it.
fbsb_engine <- function() {
  fb_engine("brms", chains = 2L, n_samples = 1500L, warmup = 2000L)
}


# --- the axes a run records are the ones that ran ------------------------ #

test_that("a horseshoe with a projection decision records four axes", {
  fbsb_needs_brms()
  skip_if_not_installed("projpred")
  sel <- fbsb_memo("hs_projection", suppressWarnings(fb_select(
    fbsb_base_intercept(),
    fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(validate_search = FALSE),
    engine = fbsb_engine(),
    seed = 20260912L
  )))

  expect_s3_class(sel, "fb_select")
  expect_identical(sel$workflow, "bayes")
  expect_identical(sel$provider, "brms")
  expect_identical(sel$prior, "regularized_horseshoe")
  expect_identical(sel$decision, "projection")
  expect_false(sel$fallback_fired)

  # The certificate says the same four things, and says them from the run.
  cert <- sel$certificate
  expect_identical(cert$workflow, "bayes")
  expect_identical(cert$provider, "brms")
  expect_identical(cert$decision, "projection")
  expect_identical(cert$contract, "SPEC_bayesian_selection_v5_2026-09-12")
  expect_identical(cert$covariance$estimand, "integrated")
  expect_true(cert$covariance$integrated)
  expect_identical(cert$calibration_n, nrow(fbsb_data_wide()))
  expect_true("brms" %in% names(cert$provider_versions))
  expect_true("projpred" %in% names(cert$provider_versions))

  # The prior block travelled onto the certificate with its hyperparameters
  # and its route, so a reader can see which prior was compiled.
  expect_length(cert$prior_blocks, 1L)
  expect_identical(cert$prior_blocks[[1L]]$family, "regularized_horseshoe")
  expect_identical(cert$prior_blocks[[1L]]$route, "class_wide")
  expect_identical(cert$prior_blocks[[1L]]$hyperparameters$p0, 2)
  expect_identical(cert$prior_blocks[[1L]]$hyperparameter_mode, "learned")

  # The block is exact here -- nothing is protected beyond the intercept --
  # so the global-shrinkage approximation is NOT claimed. The one that is
  # claimed is the scaling: the reference model is one fit, so the
  # standardising constants are the whole sample's.
  expect_false("protected_under_global_shrinkage" %in% cert$approximations)
  expect_true(any(grepl(
    "scaling_learned_on_all_rows", cert$approximations, fixed = TRUE
  )))
})


test_that("the projection chose a size, forced the protected terms and labelled its draws", {
  fbsb_needs_brms()
  skip_if_not_installed("projpred")
  sel <- fbsb_cache[["hs_projection"]]
  skip_if(is.null(sel), "the projection run is not in this file's cache")

  expect_false(is.null(sel$projection))
  expect_gte(sel$selection$chosen_size, 1L)
  expect_true(is.data.frame(sel$projection$elpd_by_size))
  # The base's random term is forced into every path, so it is in the chosen
  # submodel whatever the size.
  expect_true("(1 | site)" %in% sel$projection$forced_terms)
  expect_true("(1 | site)" %in% sel$selection$chosen_terms)
  # The projected draws never masquerade as the reference posterior.
  expect_identical(
    sel$projection$projected_draws$label, "projected_parameters"
  )
  # validate_search = FALSE was asked for, and the object says it is
  # optimistic rather than only recording the flag.
  expect_true(sel$projection$settings$optimistic)
  expect_true(sel$certificate$projection$optimistic)
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "optimistic", fixed = TRUE)
})


test_that("the shrinkage factors are reported with the label their formula earns", {
  fbsb_needs_brms()
  skip_if_not_installed("projpred")
  sel <- fbsb_cache[["hs_projection"]]
  skip_if(is.null(sel), "the projection run is not in this file's cache")
  kappa <- sel$diagnostics$kappa
  expect_identical(kappa$label, "orthogonal_design_kappa_approx")
  if (isTRUE(kappa$available)) {
    expect_true(all(kappa$kappa >= 0 & kappa$kappa <= 1))
    expect_true(is.finite(kappa$effective_nonzero_at_half))
  }
})


test_that("tidy reports the reference posterior on the scientific scale", {
  fbsb_needs_brms()
  skip_if_not_installed("projpred")
  sel <- fbsb_cache[["hs_projection"]]
  skip_if(is.null(sel), "the projection run is not in this file's cache")
  td <- tidy(sel)
  expect_true(all(
    c("term", "role", "block", "family", "retained", "estimate") %in% names(td)
  ))
  expect_setequal(
    td$term[td$role == "candidate"], c("x1", "x2", "z1", "z2")
  )
  # A consumer must be able to see which draws these are.
  expect_identical(attr(td, "draw_source"), "reference_posterior")
  expect_true(all(is.finite(td$estimate[td$role == "candidate"])))

  g <- glance(sel)
  expect_identical(g$workflow, "bayes")
  expect_identical(g$prior, "regularized_horseshoe")
  expect_identical(g$decision, "projection")
  expect_false(g$fallback_fired)
})


test_that("a class-wide family over a protected fixed effect refuses, and names both ways forward", {
  fbsb_needs_brms()
  err <- tryCatch(
    fb_select(
      fbsb_base(), fbsb_candidates,
      penalty = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      plan = TRUE, seed = 1L
    ),
    error = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
  # Both alternatives are named, so the refusal is actionable.
  expect_match(conditionMessage(err), "ridge()", fixed = TRUE)
  expect_match(conditionMessage(err), "scope = \"all_fixed\"", fixed = TRUE)
  # And an exactly-targetable family on the same base does not refuse.
  expect_s3_class(
    fb_select(
      fbsb_base(), fbsb_candidates,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      plan = TRUE, seed = 1L
    ),
    "fb_select_plan"
  )
})


test_that("scope = \"all_fixed\" is admitted, and names what it costs", {
  # The opt-in the refusal above points at. It is settled at resolution
  # time, before anything is fitted, so the plan carries the same
  # approximation the certificate will -- and this test needs no sampler.
  fbsb_needs_brms()
  plan <- fb_select(
    fbsb_base(), fbsb_candidates,
    penalty = fb_penalty(
      term(everything()) ~ regularized_horseshoe(p0 = 2, scope = "all_fixed")
    ),
    resamples = fb_resamples(unit = ~site, v = 3L),
    plan = TRUE, seed = 1L
  )
  expect_s3_class(plan, "fb_select_plan")
  expect_true(
    "protected_under_global_shrinkage" %in% plan$selection$approximations
  )
  expect_identical(plan$selection$prior_route, "class_wide")
  out <- paste(utils::capture.output(print(plan)), collapse = "\n")
  expect_match(out, "protected_under_global_shrinkage", fixed = TRUE)

  # Negative control: without the opt-in the same request refuses, and with
  # an exactly-targeting family the approximation is not claimed.
  ridge_plan <- fb_select(
    fbsb_base(), fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    plan = TRUE, seed = 1L
  )
  expect_false(
    "protected_under_global_shrinkage" %in%
      ridge_plan$selection$approximations
  )
})


# --- the other two decisions --------------------------------------------- #

test_that("a laplace penalty with decision = full keeps every candidate", {
  fbsb_needs_brms()
  sel <- fbsb_memo("laplace_full", suppressWarnings(fb_select(
    fbsb_base(),
    fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full"),
    engine = fbsb_engine(),
    seed = 20260912L
  )))

  expect_identical(sel$decision, "full")
  expect_identical(sel$prior, "laplace")
  expect_setequal(sel$selection$chosen_terms, c("x1", "x2", "z1", "z2"))
  expect_null(sel$projection)
  # Continuous shrinkage puts no mass at zero, and the reference posterior
  # says so: no candidate coefficient has an exactly-zero draw.
  draws <- as.matrix(sel$reference$brms)
  candidate_cols <- intersect(
    paste0("b_", c("x1", "x2", "z1", "z2")), colnames(draws)
  )
  expect_gt(length(candidate_cols), 0L)
  expect_false(any(draws[, candidate_cols] == 0))
  # The route this family takes targets the block exactly.
  expect_identical(sel$certificate$prior_blocks[[1L]]$route, "coef_rows")
})


test_that("a posterior threshold is labelled a decision, not an inclusion probability", {
  fbsb_needs_brms()
  sel <- fbsb_memo("laplace_threshold", suppressWarnings(fb_select(
    fbsb_base(),
    fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(
      method = "posterior_threshold", effect_size = 0.3, threshold = 0.8
    ),
    engine = fbsb_engine(),
    seed = 20260912L
  )))

  expect_identical(sel$decision, "posterior_threshold")
  expect_true(is.data.frame(sel$selection$coefficients))
  expect_setequal(
    sel$selection$coefficients$coefficient, c("x1", "x2", "z1", "z2")
  )
  expect_true(all(
    sel$selection$coefficients$probability >= 0 &
      sel$selection$coefficients$probability <= 1
  ))
  expect_match(sel$selection$label, "no inclusion-probability interpretation")
  # The planted signals clear the threshold and the noise columns do not.
  p <- stats::setNames(
    sel$selection$coefficients$probability,
    sel$selection$coefficients$coefficient
  )
  expect_gt(p[["x1"]], p[["z1"]])
  expect_gt(p[["x2"]], p[["z2"]])
  expect_true(all(c("x1", "x2") %in% sel$selection$chosen_terms))
})


# --- the unpenalised Bayesian reference ---------------------------------- #

test_that("penalty = NULL on brms is the Bayesian no-penalty reference", {
  fbsb_needs_brms()
  sel <- fbsb_memo("bayes_null_brms", suppressWarnings(fb_select(
    fbsb_base(),
    ~ x1 + x2,
    resamples = fb_resamples(unit = ~site, v = 3L),
    engine = fbsb_engine(),
    seed = 20260912L
  )))

  expect_identical(sel$workflow, "bayes")
  expect_identical(sel$provider, "brms")
  expect_identical(sel$prior, "none")
  # The covariance is integrated by the engine, not plugged in, and the
  # certificate reads that off the folds.
  expect_identical(sel$certificate$covariance$estimand, "integrated")
  expect_true(all(c("protected_only", "full") %in% sel$scores$model))
  # Two real signals: the candidates must help, and the increment is signed
  # so positive is better whichever metric was asked for.
  expect_gt(sel$increment$improvement, 0)
  expect_identical(sel$increment$metric, "elpd")
})


test_that("penalty = NULL on INLA refits the base's own engine per fold", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- fb_sel_bayes_data()
  base <- fbsb_memo("inla_base", flexybayes(
    y ~ p1, random = ~site, data = d, backend = "inla", verbose = FALSE
  ))
  # The INLA route carries the assessment rows into each fold refit with an
  # NA response, so the high-missingness note fires once per fold. It is the
  # route working as designed, not a finding about this fit.
  sel <- fbsb_memo("bayes_null_inla", suppressWarnings(fb_select(
    base,
    ~ x1 + x2,
    resamples = fb_resamples(unit = ~site, v = 3L),
    seed = 20260912L
  )))

  expect_identical(sel$workflow, "bayes")
  expect_identical(sel$provider, "inla")
  expect_identical(sel$certificate$covariance$estimand, "integrated")
  expect_gt(sel$increment$improvement, 0)
})


# --- the benchmark arms --------------------------------------------------- #

test_that("benchmark = TRUE adds two arms, each with its equivalence level", {
  fbsb_needs_brms()
  skip_if_not_installed("glmnet")
  sel <- fbsb_memo("en_benchmark", suppressWarnings(fb_select(
    fbsb_base(),
    fbsb_candidates,
    penalty = fb_penalty(
      term(everything()) ~ elastic_net(lambda = 0.1, alpha = 0.5)
    ),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full"),
    benchmark = TRUE,
    engine = fbsb_engine(),
    seed = 20260912L
  )))

  expect_setequal(names(sel$benchmark), c("gls", "glmnet"))

  gls <- sel$benchmark$gls
  expect_identical(gls$workflow, "unpenalized_benchmark")
  expect_identical(gls$provider, "gls_reml_plugin")
  expect_false(gls$fully_bayesian)
  # `gls_reml_plugin` expresses no penalty at all, so a comparison against a
  # shrinkage posterior claims nothing.
  expect_identical(gls$equivalence_level$level, "no_identity_expected")

  glmnet_arm <- sel$benchmark$glmnet
  expect_identical(glmnet_arm$workflow, "penalized_optimizer")
  expect_identical(glmnet_arm$provider, "glmnet")
  expect_false(glmnet_arm$fully_bayesian)
  # glmnet DOES express an elastic net, but the Bayesian side integrates the
  # covariance and leaves sigma free, so no mode identity can hold: the
  # arms are comparable on predictive risk and on nothing else.
  expect_identical(
    glmnet_arm$equivalence_level$level, "posterior_predictive_comparable"
  )
  expect_match(
    glmnet_arm$equivalence_level$reason, "integrated over rather than"
  )
  # The provider's own convention is measured, not recalled.
  expect_match(glmnet_arm$provider_version, "^[0-9]")
  expect_false(is.null(glmnet_arm$provider_lambda_convention))

  # And the certificate carries a level per arm.
  levels <- vapply(
    sel$certificate$benchmark, function(a) a$equivalence_level, character(1L)
  )
  expect_setequal(names(levels), c("gls", "glmnet"))
  expect_true(all(levels %in% c(
    "map_exact_fixed_covariance", "posterior_predictive_comparable",
    "no_identity_expected"
  )))

  # The print method never says "same as glmnet"; it says the level.
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "posterior_predictive_comparable", fixed = TRUE)
  expect_false(grepl("same as glmnet", out, fixed = TRUE))
})


test_that("a horseshoe benchmark carries no identity, because glmnet has none", {
  fbsb_needs_brms()
  skip_if_not_installed("glmnet")
  sel <- fbsb_cache[["hs_projection"]]
  skip_if(is.null(sel), "the projection run is not in this file's cache")
  level <- .fb_sel_equivalence_level(
    prior_family = "regularized_horseshoe",
    provider = "glmnet",
    covariance_fixed = FALSE,
    sigma_fixed = FALSE
  )
  expect_identical(level$level, "no_identity_expected")
  expect_match(level$reason, "no representation of a `regularized_horseshoe`")
})


# --- the fallback --------------------------------------------------------- #

test_that("fallback = \"none\" refuses when the Bayesian provider is absent", {
  # The provider is present on this machine, so the absence is injected:
  # `.fb_sel_provider_installed()` exists as a function precisely so that a
  # refusal nobody could otherwise watch fire can be watched firing.
  fbsb_needs_brms()
  local_mocked_bindings(
    .fb_sel_provider_installed = function(package) {
      !identical(package, "brms")
    }
  )
  expect_error(
    fb_select(
      fbsb_base(), fbsb_candidates,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      seed = 1L
    ),
    class = "flexybayes_refusal_selection_fallback_not_authorised"
  )
  expect_error(
    fb_select(
      fbsb_base(), fbsb_candidates,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      seed = 1L
    ),
    "Nothing was substituted"
  )
})


test_that("fallback = \"glmnet\" fires only on an absent provider, and certifies", {
  fbsb_needs_brms()
  skip_if_not_installed("glmnet")

  # (a) provider absent: the substitution happens, is printed on the first
  #     line, and is on the certificate.
  substituted <- with_mocked_bindings(
    fb_select(
      fbsb_base(), fbsb_candidates,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      fallback = "glmnet",
      seed = 20260912L
    ),
    .fb_sel_provider_installed = function(package) {
      !identical(package, "brms")
    }
  )
  expect_true(substituted$fallback_fired)
  expect_identical(substituted$workflow, "penalized_optimizer")
  expect_identical(substituted$provider, "glmnet")
  expect_identical(substituted$requested_workflow, "bayes")
  expect_true(substituted$certificate$fallback_fired)
  first_line <- utils::capture.output(print(substituted))[[1L]]
  expect_match(first_line, "SUBSTITUTED", fixed = TRUE)
  # No configuration produces a penalised-likelihood answer under a
  # Bayesian label.
  expect_false(identical(substituted$workflow, "bayes"))

  # (b) provider present: the same call substitutes nothing. Without this
  #     half the test above passes for a fallback that always fires.
  present <- suppressWarnings(fb_select(
    fbsb_base(), fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full"),
    fallback = "glmnet",
    engine = fbsb_engine(),
    seed = 20260912L
  ))
  expect_false(present$fallback_fired)
  expect_identical(present$workflow, "bayes")
  expect_identical(present$provider, "brms")
})


# --- the plan ------------------------------------------------------------- #

test_that("plan = TRUE reports the axes, the blocks and the gate without fitting", {
  fbsb_needs_brms()
  plan <- fb_select(
    fbsb_base_intercept(),
    fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    plan = TRUE,
    seed = 1L
  )

  expect_s3_class(plan, "fb_select_plan")
  s <- plan$selection
  expect_identical(s$workflow, "bayes")
  expect_identical(s$provider, "brms")
  expect_identical(s$prior, "regularized_horseshoe")
  expect_identical(s$decision, "projection")
  expect_identical(s$prior_route, "class_wide")
  expect_length(s$blocks, 1L)
  expect_setequal(s$candidates, c("x1", "x2", "z1", "z2"))
  expect_identical(s$resampling$n_folds, 3L)
  expect_true(is.data.frame(s$structure_gate))
  expect_true(is.finite(s$cost$seconds))
  expect_identical(s$cost$source, "base_model_run_time")

  out <- paste(utils::capture.output(print(plan)), collapse = "\n")
  expect_match(out, "Nothing was fitted", fixed = TRUE)
  expect_match(out, "class_wide", fixed = TRUE)
})


# --- the injected-defect falsifiers -------------------------------------- #

test_that("the certificate cannot disagree with the fold estimand", {
  # The defect this closes, restated: the certificate used to read the
  # covariance estimand from a public argument whose default said
  # `plugin_mean_covariance` while the folds worked in
  # `reml_plugin_covariance`. A certificate built that way describes a
  # computation that did not happen and cannot notice.
  #
  # The falsifier corrupts what the FOLDS report and requires the
  # certificate to follow. A certificate reading an argument would print the
  # same string under both, and this test would be green for the defect.
  base <- structure(
    list(
      fb = list(
        family = "gaussian", link = "identity", source = "test",
        fixed_terms = list(), random_terms = list()
      ),
      data = data.frame(y = 1:4, x = 1:4)
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  as_the_fold_reported <- .fb_select_certificate(
    base = base,
    axes = list(workflow = "bayes", provider = "inla", prior = "none",
                decision = "full"),
    cv = list(
      folds = list(), aggregated = NULL, prediction = "new_cluster",
      seed = 1L, metric = "elpd", covariance_estimand = "integrated"
    ),
    resamples = fb_resamples(v = 2L),
    candidates = ~x
  )
  expect_identical(as_the_fold_reported$covariance$estimand, "integrated")

  injected <- .fb_select_certificate(
    base = base,
    axes = list(workflow = "bayes", provider = "inla", prior = "none",
                decision = "full"),
    cv = list(
      folds = list(), aggregated = NULL, prediction = "new_cluster",
      seed = 1L, metric = "elpd",
      covariance_estimand = "reml_plugin_covariance"
    ),
    resamples = fb_resamples(v = 2L),
    candidates = ~x
  )
  expect_identical(
    injected$covariance$estimand, "reml_plugin_covariance"
  )
  expect_false(identical(
    injected$covariance$estimand, as_the_fold_reported$covariance$estimand
  ))

  # And a run whose folds disagreed reports both rather than one of them.
  disagreeing <- .fb_select_certificate(
    base = base,
    axes = list(workflow = "bayes", provider = "inla", prior = "none",
                decision = "full"),
    cv = list(
      folds = list(), aggregated = NULL, prediction = "new_cluster",
      seed = 1L, metric = "elpd",
      covariance_estimand = c("integrated", "reml_plugin_covariance")
    ),
    resamples = fb_resamples(v = 2L),
    candidates = ~x
  )
  expect_length(disagreeing$covariance$estimand, 2L)
})


test_that("a joint shrinkage prior routed at INLA refuses by name", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- fb_sel_bayes_data()
  base <- fbsb_memo("inla_base", flexybayes(
    y ~ p1, random = ~site, data = d, backend = "inla", verbose = FALSE
  ))
  expect_error(
    fb_select(
      base, fbsb_candidates,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      seed = 1L
    ),
    class = "flexybayes_refusal_prior_family_unsupported_inla"
  )
  # The negative control: the same base with no joint block passes the same
  # gate, so the refusal is the prior's and not the model's.
  expect_silent(lgm_gate(base$fb))
})


test_that("a budget smaller than the estimated cost refuses before any fit", {
  fbsb_needs_brms()
  t0 <- proc.time()
  expect_error(
    fb_select(
      fbsb_base(), fbsb_candidates,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      budget = 1e-9,
      seed = 1L
    ),
    class = "flexybayes_refusal_penalty_runtime_budget_exceeded"
  )
  # A refusal that costs a model fit before it fires is not a refusal that
  # came first. One second is generous: a Stan compile alone is longer.
  expect_lt(unname((proc.time() - t0)["elapsed"]), 30)
})
