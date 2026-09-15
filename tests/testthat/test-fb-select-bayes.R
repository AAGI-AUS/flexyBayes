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
  # 20 sites of eight plots: the prior scales are calibrated against the 20
  # sites, not against the 160 rows, and the certificate says which.
  expect_identical(cert$calibration_n, 20L)
  expect_identical(cert$n_rows, nrow(fbsb_data_wide()))
  expect_identical(cert$calibration_basis, "resampling_units")
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
  # claimed is the projection's own search, whose k-fold refits reuse the
  # reference fit's whole-sample preprocessing. Whole-sample scaling is no
  # longer claimed for the out-of-fold scores, which learn it per fold.
  expect_false("protected_under_global_shrinkage" %in% cert$approximations)
  expect_true(any(grepl(
    "projection_search_folds_reuse_whole_sample_preprocessing",
    cert$approximations, fixed = TRUE
  )))
  expect_false(any(grepl(
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


test_that("the projection's own site refuses a non-elpd metric by class", {
  # The late plain stop() inside `.fb_sel_project()` is retired: a caller that
  # reaches the projection directly, past the verb's early check, meets the
  # same typed refusal rather than a bare error.
  fbsb_needs_brms()
  skip_if_not_installed("projpred")
  sel <- fbsb_cache[["hs_projection"]]
  skip_if(is.null(sel), "the projection run is not in this file's cache")
  expect_error(
    .fb_sel_project(
      sel$reference,
      folds = sel$folds,
      decision = fb_decision(method = "projection", metric = "rmse"),
      seed = 1L
    ),
    class = "flexybayes_refusal_projection_metric_unsupported"
  )
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


test_that("tidy retains a factor and an interaction by term, not by contrast column", {
  # The defect this closes. Retention was `map$column %in% chosen_terms`.
  # `chosen_terms` are formula term labels -- `f`, `x1:f` -- and `map$column`
  # is model-matrix columns -- `fB`, `fC`, `x1:fB`, `x1:fC`. A projection
  # that kept `x1` and `f` therefore reported every factor and interaction
  # column as dropped, on exactly the formulas this package is written for.
  #
  # No sampler is involved: the column map and the decision are the whole
  # of what `tidy()` reads for the retained flag.
  set.seed(4021L)
  d <- data.frame(
    y = stats::rnorm(60),
    p1 = stats::rnorm(60),
    x1 = stats::rnorm(60),
    f = factor(rep(c("A", "B", "C"), length.out = 60))
  )
  map <- .fb_sel_term_map(~ p1 + x1 + f + x1:f, d)
  part <- .fb_sel_partition_design(map, protected = "p1")
  column_map <- .fb_select_column_map(part, list(), "p1", map)

  # The design really does carry contrast columns whose names are not the
  # term labels, which is the precondition the defect needed.
  expect_true(all(c("fB", "fC", "x1:fB", "x1:fC") %in% column_map$column))
  expect_false(any(c("f", "x1:f") %in% column_map$column))

  chosen <- c("x1", "f", "x1:f")
  sel <- structure(
    list(
      workflow = "bayes", provider = "brms", prior = "ridge",
      decision = "projection",
      column_map = column_map,
      blocks = NULL,
      reference = NULL,
      selection = list(
        chosen_terms = chosen, chosen_vocabulary = "term", chosen_size = 3L
      )
    ),
    class = "fb_select"
  )
  td <- tidy(sel)
  retained <- stats::setNames(td$retained, td$term)

  # The injected defect, shown red on this same table: the column-string
  # rule marks every contrast column of a chosen factor as dropped.
  by_column <- column_map$column %in% chosen | column_map$role == "protected"
  expect_false(by_column[column_map$column == "fB"])
  expect_false(by_column[column_map$column == "x1:fB"])

  # And the term rule follows the term for every contrast column of it.
  expect_true(retained[["x1"]])
  expect_true(all(retained[c("fB", "fC")]))
  expect_true(all(retained[c("x1:fB", "x1:fC")]))
  expect_true(retained[["p1"]])
  expect_identical(td$source_term[td$term == "x1:fB"], "x1:f")

  # The negative control, and the other half of "follows the term": a term
  # the decision did not keep is dropped at every one of its columns, so
  # the rule is neither retaining everything nor splitting a term in two.
  dropped <- structure(
    list(
      workflow = "bayes", provider = "brms", prior = "ridge",
      decision = "projection",
      column_map = column_map, blocks = NULL, reference = NULL,
      selection = list(
        chosen_terms = c("x1", "f"), chosen_vocabulary = "term",
        chosen_size = 2L
      )
    ),
    class = "fb_select"
  )
  td2 <- tidy(dropped)
  expect_true(all(td2$retained[td2$source_term == "f"]))
  expect_false(any(td2$retained[td2$source_term == "x1:f"]))
  expect_true(td2$retained[td2$term == "x1"])

  # A posterior threshold is the one decision written in coefficients, and
  # it is still matched on the column.
  coefwise <- structure(
    list(
      workflow = "bayes", provider = "brms", prior = "laplace",
      decision = "posterior_threshold",
      column_map = column_map, blocks = NULL, reference = NULL,
      selection = list(
        chosen_terms = "fB", chosen_vocabulary = "coefficient",
        chosen_size = 1L
      )
    ),
    class = "fb_select"
  )
  td3 <- tidy(coefwise)
  expect_true(td3$retained[td3$term == "fB"])
  expect_false(td3$retained[td3$term == "fC"])
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


# --- the primary arm is scored, or the comparison is not claimed ---------- #

test_that("the stage tables say why they are absent, and substitute nothing", {
  # The defect this closes, at the unit that decides it. The Bayesian
  # penalised route set `aggregated = NULL`, so there was no out-of-fold
  # score on the primary arm at all -- and the glmnet arm beside it was
  # still labelled `posterior_predictive_comparable`, whose own definition
  # is "scored on the same folds, target, loss, unit and candidate design".
  # The contract itself is tested in test-fb-select-score-contract.R; this
  # is the wiring that decides which arms reach each stage.
  base <- structure(
    list(
      fb = list(response = "y", family = "gaussian", link = "identity"),
      data = data.frame(y = stats::rnorm(4L))
    ),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  arms <- list(gls = list())

  # No primary table: both stages are absent and say what would produce one.
  none <- .fb_select_stage_tables(
    base,
    axes = list(workflow = "bayes", prior = "ridge", decision = "full"),
    cv = list(clusters = NULL, metric = "elpd"),
    benchmark = list(arms = arms),
    metric = "nlpd"
  )
  expect_false(none$reference$available)
  expect_match(none$reference$reason, "benchmark = TRUE", fixed = TRUE)
  expect_false(none$decision$available)

  # A posterior threshold defines no predictive distribution: the decision
  # stage names that, and the reference is not put in its place.
  threshold <- .fb_select_stage_tables(
    base,
    axes = list(
      workflow = "bayes", prior = "laplace", decision = "posterior_threshold"
    ),
    cv = list(clusters = NULL, metric = "elpd"),
    benchmark = NULL,
    metric = "nlpd"
  )
  expect_false(threshold$decision$available)
  expect_match(threshold$decision$reason, "posterior threshold", fixed = TRUE)
  expect_null(threshold$decision$arms)

  # With every candidate protected there is no decision at all.
  protected <- .fb_select_stage_tables(
    base,
    axes = list(workflow = "bayes", prior = "none", decision = "full"),
    cv = list(clusters = NULL, metric = "elpd"),
    benchmark = NULL,
    metric = "nlpd"
  )
  expect_match(
    protected$decision$reason, "every candidate is protected", fixed = TRUE
  )
})


test_that("a ridge benchmark scores both arms on one loss and leaves the mode oracle empty", {
  fbsb_needs_brms()
  skip_if_not_installed("glmnet")
  sel <- fbsb_memo("ridge_benchmark", suppressWarnings(fb_select(
    fbsb_base(),
    fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 2L),
    decision = fb_decision(method = "full"),
    benchmark = TRUE,
    engine = fbsb_engine(),
    seed = 20260912L
  )))

  # The primary arm now writes the same unit-level out-of-fold table the
  # `penalty = NULL` route writes: protected baseline and reference, scored
  # from posterior predictive draws at the resampling unit.
  expect_true(is.data.frame(sel$scores))
  expect_true(all(c("protected_only", "full") %in% sel$scores$model))
  expect_false(is.null(sel$increment))
  expect_true(all(
    c("elpd", "rmse", "nlpd") %in% names(as.data.frame(sel$clusters))
  ))

  ps <- sel$certificate$reference_predictive
  expect_true(ps$available)
  expect_true(all(c("bayes", "gls", "glmnet") %in% ps$arms$arm))
  expect_true(all(ps$arms$stage == "reference"))
  # One contract, one loss, one unit, and the same number of them on every
  # arm.
  expect_length(unique(ps$arms$contract_hash), 1L)
  expect_identical(length(unique(ps$arms$n_units)), 1L)
  expect_identical(length(unique(ps$arms$unit)), 1L)

  # The decision stage places the models each procedure chose: decision =
  # "full" keeps every candidate, so the Bayesian choice is the reference,
  # beside the optimiser's own chosen lambda -- under a contract that
  # differs from the reference stage's in the stage and nothing else.
  dp <- sel$certificate$decision_predictive
  expect_true(dp$available)
  expect_setequal(dp$arms$arm, c("bayes", "glmnet"))
  expect_identical(dp$arms$model[dp$arms$arm == "bayes"], "full")
  expect_match(dp$arms$model[dp$arms$arm == "glmnet"], "^benchmark_")
  expect_true(all(dp$arms$stage == "decision"))
  expect_identical(
    .fb_select_contract_differences(ps$contract, dp$contract), "stage"
  )

  # The fold-local constants are on the certificate, one row per fold.
  fp <- sel$certificate$fold_preprocessing
  expect_identical(fp$learned_on, "analysis_rows")
  expect_identical(nrow(fp$folds), length(sel$folds))

  # The comparison is now one that happened, so the middle level is earned.
  expect_identical(
    sel$benchmark$glmnet$equivalence_level$level,
    "posterior_predictive_comparable"
  )
  expect_true(
    sel$benchmark$glmnet$equivalence_level$shared_predictive_table
  )

  # And the Tier-1 stanza is empty, because a sampled route integrates the
  # covariance and samples sigma. It is a stanza rather than an absent
  # field: "no mode oracle ran" is what a reader has to be told.
  mo <- sel$certificate$mode_oracle
  expect_false(mo$available)
  expect_length(mo$arms, 0L)
  expect_match(mo$reason, "held fixed", fixed = TRUE)

  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "reference stage", fixed = TRUE)
  expect_match(out, "decision stage", fixed = TRUE)
  expect_match(out, "mode oracle: none", fixed = TRUE)
  cert_out <- paste(
    utils::capture.output(print(sel$certificate)), collapse = "\n"
  )
  expect_match(cert_out, "mode oracle: none", fixed = TRUE)
  expect_match(cert_out, "reference stage", fixed = TRUE)
  expect_match(cert_out, ps$contract$hash, fixed = TRUE)
  expect_match(cert_out, "fold prep", fixed = TRUE)
})


test_that("a horseshoe benchmark claims no identity against glmnet", {
  fbsb_needs_brms()
  skip_if_not_installed("glmnet")
  sel <- fbsb_memo("hs_benchmark", suppressWarnings(fb_select(
    fbsb_base_intercept(),
    fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
    resamples = fb_resamples(unit = ~site, v = 2L),
    decision = fb_decision(method = "full"),
    benchmark = TRUE,
    engine = fbsb_engine(),
    seed = 20260912L
  )))

  # glmnet has no representation of a scale mixture, so the optimiser arm
  # does not run at all and the refusal is carried by name rather than the
  # arm being reported with a level it did not earn.
  expect_false("glmnet" %in% names(sel$benchmark))
  expect_true(any(grepl(
    "regularized_horseshoe", sel$certificate$refusals, fixed = TRUE
  )))
  # The level the family would carry, at every setting of the two
  # conditions and with a shared table available.
  for (shared in c(TRUE, FALSE)) {
    expect_identical(
      .fb_sel_equivalence_level(
        "regularized_horseshoe", "glmnet", FALSE, FALSE,
        shared_predictive_table = shared
      )$level,
      "no_identity_expected"
    )
  }

  # The unpenalised arm and the Bayesian arm DO share a contract at the
  # reference stage, and are reported side by side on it.
  ps <- sel$certificate$reference_predictive
  expect_true(ps$available)
  expect_setequal(ps$arms$arm, c("bayes", "gls"))
  expect_identical(length(unique(ps$arms$n_units)), 1L)
  expect_false(sel$certificate$mode_oracle$available)
  # With no optimiser arm, the decision stage has one arm and no comparison.
  expect_false(sel$certificate$decision_predictive$available)
  expect_null(sel$certificate$decision_predictive$arms)
})


test_that("a projection benchmark scores the chosen submodel at the decision stage", {
  # The defect this closes. The one side-by-side table placed the Bayesian
  # reference beside the optimiser's SELECTED model, and the projection's
  # chosen submodel was not in it at all. Both stages are now scored under
  # one contract each, on a live run.
  fbsb_needs_brms()
  skip_if_not_installed("projpred")
  skip_if_not_installed("glmnet")
  sel <- fbsb_memo("laplace_projection_benchmark", suppressWarnings(fb_select(
    fbsb_base(),
    fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "projection", validate_search = FALSE),
    benchmark = TRUE,
    engine = fbsb_engine(),
    seed = 20260912L
  )))
  cert <- sel$certificate
  ref <- cert$reference_predictive
  dec <- cert$decision_predictive

  expect_true(ref$available)
  expect_true(dec$available)
  expect_true(all(ref$arms$stage == "reference"))
  expect_true(all(dec$arms$stage == "decision"))
  expect_identical(ref$arms$model[ref$arms$arm == "bayes"], "full")
  expect_identical(
    dec$arms$model[dec$arms$arm == "bayes"], "projected_submodel"
  )
  expect_match(dec$arms$model[dec$arms$arm == "glmnet"], "^benchmark_")
  expect_length(unique(ref$arms$contract_hash), 1L)
  expect_length(unique(dec$arms$contract_hash), 1L)
  expect_identical(
    .fb_select_contract_differences(ref$contract, dec$contract), "stage"
  )
  # The projected submodel is scored on its own draws: its score is not the
  # reference's under another name.
  expect_false(isTRUE(all.equal(
    dec$arms$score[dec$arms$arm == "bayes"],
    ref$arms$score[ref$arms$arm == "bayes"]
  )))
  expect_true("projected_submodel" %in% sel$clusters$model)
  # The chosen submodel's terms were searched again on each fold's training
  # rows, and the certificate says where.
  expect_identical(
    cert$projection$decision_stage_search, "each_fold_training_rows"
  )
  expect_length(sel$projection$decision_stage_terms, 3L)
  expect_match(
    dec$arms$label[dec$arms$arm == "bayes"], "each fold's training rows",
    fixed = TRUE
  )

  # The out-of-fold score learned its constants per fold, and the search
  # inside the provider is named as the place it did not.
  expect_identical(cert$fold_preprocessing$learned_on, "analysis_rows")
  expect_identical(nrow(cert$fold_preprocessing$folds), 3L)
  expect_true(any(grepl(
    "projection_search_folds_reuse_whole_sample_preprocessing",
    cert$approximations, fixed = TRUE
  )))

  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "reference stage", fixed = TRUE)
  expect_match(out, "projected_submodel", fixed = TRUE)
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


# --- the structure gate leads, and it resolves the decision --------------- #

fbsb_genomic <- function() {
  set.seed(7L)
  d <- data.frame(
    site = factor(rep(sprintf("S%d", seq_len(8L)), each = 6L)),
    gen = factor(rep(sprintf("G%d", seq_len(6L)), times = 8L))
  )
  d$p1 <- stats::rnorm(nrow(d))
  d$x1 <- stats::rnorm(nrow(d))
  d$x2 <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$p1 + stats::rnorm(nrow(d))
  K <- diag(6L)
  dimnames(K) <- list(levels(d$gen), levels(d$gen))
  list(data = d, K = K)
}


test_that("the plan and the printed result lead with the structure-gate table", {
  fbsb_needs_brms()
  plan <- fb_select(
    fbsb_base(), fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    plan = TRUE, seed = 1L
  )
  out <- utils::capture.output(print(plan))
  # The table, not a one-word verdict, and before anything else.
  expect_match(out[[2L]], "Projection structure gate", fixed = TRUE)
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "response family", fixed = TRUE)
  expect_match(joined, "(1 | site)", fixed = TRUE)
  expect_match(joined, "probes C1-C3", fixed = TRUE)
  # One row per structural feature the model carries, with a verdict.
  gate <- plan$selection$structure_gate
  expect_true(is.data.frame(gate))
  expect_true(all(c("feature", "verdict", "reason") %in% names(gate)))
  expect_true(all(gate$verdict == "supported"))

  # The clustered-draw counts are on the plan too, not only on the object.
  expect_match(joined, "20 clustered draws", fixed = TRUE)
  expect_match(joined, "evaluation on 100", fixed = TRUE)
})


test_that("decision = \"auto\" resolves against the gate and records the resolution", {
  fbsb_needs_brms()
  skip_if_not_installed("projpred")

  # Green gate: the projection runs, and the object says the method was
  # resolved rather than asked for.
  green <- fb_select(
    fbsb_base(), fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    plan = TRUE, seed = 1L
  )
  expect_identical(green$selection$decision, "projection")
  expect_identical(green$selection$decision_requested, "auto")
  expect_identical(green$selection$decision_resolved_from, "auto")

  # Red gate: a known-covariance genomic term lowers to
  # `(1 | gr(gen, cov = K))`, which projpred refuses (probe C5). The
  # decision falls to `full` BEFORE a reference model is sampled -- which
  # is the whole point: it used to error after the fit.
  bits <- fbsb_genomic()
  base <- fbsb_memo("genomic_base", suppressWarnings(fb_brms(
    y ~ p1, random = ~ vm(gen, K), data = bits$data,
    known_matrices = list(K = bits$K), chains = 2L, n_samples = 600L,
    warmup = 600L, seed = 11L, verbose = FALSE, mcmc_verbose = FALSE
  )))
  red <- fb_select(
    base, ~ x1 + x2,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~gen, v = 3L),
    plan = TRUE, seed = 1L
  )
  expect_identical(red$selection$decision, "full")
  expect_identical(red$selection$decision_resolved_from, "auto")
  expect_match(
    red$selection$decision_resolution_reason, "gate is red", fixed = TRUE
  )
  red_out <- paste(utils::capture.output(print(red)), collapse = "\n")
  expect_match(red_out, "gr(gen, cov = K)", fixed = TRUE)
  expect_match(red_out, "probe C5", fixed = TRUE)

  # And naming the projection explicitly on the same base still refuses --
  # before the reference model is sampled, not several minutes into it.
  t0 <- proc.time()
  err <- tryCatch(
    fb_select(
      base, ~ x1 + x2,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
      resamples = fb_resamples(unit = ~gen, v = 3L),
      decision = fb_decision(method = "projection"),
      engine = fbsb_engine(), seed = 1L
    ),
    error = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_projection_structure_unsupported"
  )
  expect_match(conditionMessage(err), "gr(gen, cov = K)", fixed = TRUE)
  expect_match(conditionMessage(err), "decision = \"full\"", fixed = TRUE)
  # A refusal that costs a reference fit before it fires is not a refusal
  # that came first. A Stan compile alone is longer than this.
  expect_lt(unname((proc.time() - t0)["elapsed"]), 30)
})


test_that("a projection result prints its clustered-draw counts on the first line", {
  fbsb_needs_brms()
  skip_if_not_installed("projpred")
  sel <- fbsb_cache[["hs_projection"]]
  skip_if(is.null(sel), "the projection run is not in this file's cache")
  first_line <- utils::capture.output(print(sel))[[1L]]
  expect_match(first_line, "decision projection", fixed = TRUE)
  expect_match(first_line, "20 clustered draws", fixed = TRUE)
  expect_match(first_line, "evaluated on 100", fixed = TRUE)
  # The counts are the ones the object actually ran at, not a recital.
  expect_identical(
    sel$projection$settings$nclusters, .FB_SEL_PROJECTION_NCLUSTERS
  )
  expect_identical(
    sel$projection$settings$nclusters_pred,
    .FB_SEL_PROJECTION_NCLUSTERS_PRED
  )
  # And the gate leads the printed result.
  out <- utils::capture.output(print(sel))
  expect_match(out[[2L]], "Projection structure gate", fixed = TRUE)
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


# --- what the prior scales are calibrated against ------------------------ #

test_that("calibration_n counts resampling units, and the certificate says so", {
  # The defect this closes. `calibration_n` was `nrow(data)`, and the row
  # count went into brms's horseshoe calibration and into the elastic net's
  # canonical `n`. In a clustered design that overstates the information
  # about a coefficient by roughly the cluster size, so the global scale is
  # too small and the prior over-shrinks -- and the certificate reported the
  # row count, so nothing said by how much.
  set.seed(2026L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(20L)), each = 10L)))
  d$x1 <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d))
  expect_identical(nrow(d), 200L)

  at_unit <- .fb_select_calibration_n(fb_resamples(unit = ~site, v = 3L), d)
  expect_identical(at_unit$n, 20L)
  expect_identical(at_unit$basis, "resampling_units")
  expect_identical(at_unit$n_rows, 200L)

  # The injected defect, shown red on the same data: the row count is ten
  # times the unit count, and it is what the earlier code passed on.
  expect_identical(nrow(d), 10L * at_unit$n)

  # And where no unit was declared a row IS a unit, which is correct rather
  # than a fallback -- the basis says which.
  at_rows <- .fb_select_calibration_n(fb_resamples(v = 3L), d)
  expect_identical(at_rows$n, 200L)
  expect_identical(at_rows$basis, "rows")

  # The certificate prints the number, its basis and the row count beside
  # it, so a reader can see the ratio without going back to the data.
  base <- structure(
    list(
      fb = list(
        family = "gaussian", link = "identity", source = "test",
        fixed_terms = list(), random_terms = list()
      ),
      data = d
    ),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  cert <- .fb_select_certificate(
    base = base,
    axes = list(workflow = "bayes", provider = "brms",
                prior = "regularized_horseshoe", decision = "full"),
    cv = list(folds = list(), aggregated = NULL, prediction = "new_cluster",
              seed = 1L, metric = "elpd", covariance_estimand = "integrated"),
    resamples = fb_resamples(unit = ~site, v = 3L),
    candidates = ~x1,
    calibration = at_unit
  )
  expect_identical(cert$calibration_n, 20L)
  expect_identical(cert$n_rows, 200L)
  out <- paste(utils::capture.output(print(cert)), collapse = "\n")
  expect_match(out, "calibration_n: 20", fixed = TRUE)
  expect_match(out, "against 200 rows", fixed = TRUE)
})


test_that("the horseshoe global scale is moved to the unit, not left at the row", {
  # brms derives its global scale from `par_ratio` and the number of FITTED
  # ROWS, and it offers no way to say otherwise, so the calibration is moved
  # by moving `par_ratio`. The relationship was measured against
  # brms::make_standata() rather than read off the documentation:
  # hs_scale_global = par_ratio / sqrt(N).
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(2026L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(20L)), each = 10L)))
  for (v in paste0("x", 1:4)) d[[v]] <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d))

  map <- .fb_sel_term_map(~ x1 + x2 + x3 + x4, d)
  blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
    map,
    protected = character(0)
  )
  canonical <- blocks[[1L]]$hyperparameters$par_ratio
  expect_equal(canonical, 2 / (4 - 2))

  moved <- .fb_pen_calibrate_n(blocks, calibration_n = 20L, n_rows = 200L)
  h <- moved[[1L]]$hyperparameters
  expect_equal(h$par_ratio_canonical, canonical)
  expect_identical(h$calibration_n, 20)
  # The global scale brms will derive is the canonical ratio over the square
  # root of the UNIT count, which is the whole point of moving the ratio.
  expect_equal(h$par_ratio / sqrt(200), canonical / sqrt(20))

  # Measured against brms itself rather than against the algebra above.
  compiled <- .fb_pen_compile_brms(
    moved, .fb_pen_brms_coef_names(colnames(map$X)[map$assign != 0L]),
    n_train = 20L
  )
  sdata <- brms::make_standata(
    y ~ x1 + x2 + x3 + x4, d, prior = compiled$prior
  )
  expect_equal(sdata$hs_scale_global, canonical / sqrt(20), tolerance = 1e-10)

  # The injected defect, shown red: left at the row count the same request
  # gives a global scale smaller by the square root of the cluster size.
  at_rows <- .fb_pen_calibrate_n(blocks, calibration_n = 200L, n_rows = 200L)
  compiled_rows <- .fb_pen_compile_brms(
    at_rows, .fb_pen_brms_coef_names(colnames(map$X)[map$assign != 0L]),
    n_train = 200L
  )
  sdata_rows <- brms::make_standata(
    y ~ x1 + x2 + x3 + x4, d, prior = compiled_rows$prior
  )
  expect_lt(sdata_rows$hs_scale_global, sdata$hs_scale_global)
  expect_equal(
    sdata$hs_scale_global / sdata_rows$hs_scale_global, sqrt(10),
    tolerance = 1e-8
  )
})


test_that("a class-wide family refuses an unstandardised candidate column", {
  # Finding 6.3. The reference fit is a re-emission of the formula, so a
  # factor contrast has no data column to centre and divide and is left at
  # centre 0, scale 1. brms puts ONE global-local prior on the whole
  # population-level vector, so that 0/1 dummy then sits under the same
  # local scale as a numeric already at unit standard deviation. The
  # coefficient-row families state a scale per column and keep the
  # approximation they already carried.
  fbsb_needs_brms()
  set.seed(4021L)
  d <- fbsb_data_wide()
  d$f <- factor(rep(c("A", "B", "C"), length.out = nrow(d)))

  base <- fbsb_memo("base_intercept_factor", suppressWarnings(fb_brms(
    y ~ 1 + (1 | site), data = d, chains = 2L, n_samples = 500L,
    warmup = 500L, seed = 20260912L, verbose = FALSE, mcmc_verbose = FALSE
  )))

  err <- tryCatch(
    fb_select(
      base, ~ x1 + x2 + f,
      penalty = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
      resamples = fb_resamples(unit = ~site, v = 3L),
      plan = TRUE, seed = 1L
    ),
    error = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_penalty_candidate_columns_unscaled"
  )
  # The columns are named, and so are the two ways out.
  expect_match(conditionMessage(err), "fB", fixed = TRUE)
  expect_match(conditionMessage(err), "fC", fixed = TRUE)
  expect_match(conditionMessage(err), "ridge()", fixed = TRUE)
  expect_match(conditionMessage(err), "numeric columns", fixed = TRUE)
  expect_setequal(err$columns, c("fB", "fC"))

  # The negative control, twice over. A coefficient-row family on the same
  # block plans, and keeps the approximation rather than refusing.
  ridge_plan <- fb_select(
    base, ~ x1 + x2 + f,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full"),
    plan = TRUE, seed = 1L
  )
  expect_s3_class(ridge_plan, "fb_select_plan")
  expect_true(any(grepl(
    "candidate_columns_not_standardised",
    ridge_plan$selection$approximations, fixed = TRUE
  )))

  # And a numeric-only candidate block under the same class-wide family
  # plans, so the refusal is the unstandardised column's and not the
  # family's.
  hs_plan <- fb_select(
    base, ~ x1 + x2 + z1,
    penalty = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 1)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    plan = TRUE, seed = 1L
  )
  expect_s3_class(hs_plan, "fb_select_plan")
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


test_that("an INLA-fitted base has its reference re-emitted on brms, and says so", {
  # The defect this closes. `fb_select()` always fits the reference model
  # through brms, but it ran lgm_gate()'s twelfth check on the BASE's IR
  # whenever the base was INLA -- so an INLA user asking for a Laplace prior
  # was told the prior was impossible on INLA, about a fit that was never
  # going to INLA. The check belongs to the engine that will fit the
  # reference.
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- fb_sel_bayes_data()
  base <- fbsb_memo("inla_base", flexybayes(
    y ~ p1, random = ~site, data = d, backend = "inla", verbose = FALSE
  ))

  plan <- fb_select(
    base, fbsb_candidates,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full"),
    plan = TRUE, seed = 1L
  )
  expect_s3_class(plan, "fb_select_plan")
  expect_identical(plan$selection$provider, "brms")
  expect_identical(plan$selection$reference_engine, "brms")
  expect_identical(plan$selection$base_engine, "inla")
  expect_true(any(grepl(
    "reference_fitted_on_brms", plan$selection$approximations, fixed = TRUE
  )))
  out <- paste(utils::capture.output(print(plan)), collapse = "\n")
  expect_match(out, "the base was fitted by inla", fixed = TRUE)

  # The check still fires where it means something: a genuine INLA fit
  # request carrying a joint block on the IR is refused by the gate itself,
  # which is where the refusal belongs.
  fb_gate <- base$fb
  fb_gate$penalty_blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    .fb_sel_term_map(~ p1 + x1 + x2 + z1 + z2, d),
    protected = "p1"
  )
  gated <- lgm_gate(fb_gate)
  expect_true(is_lgm_refusal(gated))
  expect_true("prior_family_unsupported_inla" %in%
    vapply(gated$failures, function(f) f$rule_id, character(1L)))
  # The negative control: the same base with no joint block passes the same
  # gate, so the refusal is the prior's and not the model's.
  expect_silent(lgm_gate(base$fb))
})


test_that("an INLA base with a joint prior samples on brms", {
  # The other half of the same finding, executed rather than planned. A
  # small design on purpose: what is under test is the route, not the
  # posterior.
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  d <- fb_sel_bayes_data(J = 15L, m = 6L)
  base <- fbsb_memo("inla_base_wide", flexybayes(
    y ~ p1, random = ~site, data = d, backend = "inla", verbose = FALSE
  ))
  sel <- fbsb_memo("inla_base_laplace", suppressWarnings(fb_select(
    base, ~ x1 + x2,
    penalty = fb_penalty(term(everything()) ~ laplace(0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full"),
    engine = fb_engine("brms", chains = 4L, n_samples = 1500L,
                       warmup = 1500L),
    seed = 20260912L
  )))

  expect_identical(sel$workflow, "bayes")
  expect_identical(sel$provider, "brms")
  expect_identical(sel$prior, "laplace")
  expect_identical(sel$certificate$reference_engine, "brms")
  # It really did sample: the reference carries a brms posterior with the
  # candidate coefficients in it.
  expect_s3_class(sel$reference$brms, "brmsfit")
  expect_true(all(
    c("b_x1", "b_x2") %in% colnames(as.matrix(sel$reference$brms))
  ))
  # And the object says the reference was fitted by a different engine from
  # the base, on its own first block rather than only on the certificate.
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "fitted by brms", fixed = TRUE)
  expect_match(out, "the base was fitted by inla", fixed = TRUE)
  expect_true(any(grepl(
    "reference_fitted_on_brms", sel$certificate$approximations, fixed = TRUE
  )))
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
