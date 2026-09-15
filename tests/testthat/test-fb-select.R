# =============================================================================
# fb_select(penalty = NULL): the B0 deliverable.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 2, 3.1, 11.
#
# The partition principle is the thing under test. The base model, the
# extended base model and a penalised selection are one computation under
# three partitions of the design, so `penalty = NULL` is not a degenerate
# case kept for tidiness -- it is the lambda = 0 corner, it exercises every
# part of the substrate except the provider adapter, and it is what lets the
# eventual penalty be a small diff against machinery already proved.
#
# What these tests hold it to: the gate refuses what it should before
# anything runs, the increment is computed against a protected-only baseline
# from the SAME folds, and the certificate records what was done including
# where it was approximate.
#
# The route these tests exercise now has a name. SPEC v5 section 2 makes the
# closed-form GLS-on-a-plugged-in-fold-covariance computation
# `workflow = "unpenalized_benchmark"` and reserves the default,
# `workflow = "bayes"`, for the per-fold engine refit -- a different
# estimator that integrates the covariance rather than plugging it in. Every
# call below therefore names the workflow it was written for. The Bayesian
# route has its own file, test-fb-select-bayes.R.
#
# One thing these tests used to assert is now a defect they would have
# hidden. The certificate's covariance estimand was read from the public
# `covariance` argument, which defaulted to `plugin_mean_covariance` while
# the folds used `reml_plugin_covariance`. SPEC v5 section 9.2 closes that by
# reading the field from the fold object; the assertions below read it the
# same way.
# =============================================================================


fb_sel_verb_data <- function(seed = 20260912L, J = 24L, m = 6L) {
  set.seed(seed)
  d <- data.frame(
    site = factor(rep(sprintf("S%02d", seq_len(J)), each = m))
  )
  per_site <- function() {
    stats::setNames(stats::rnorm(J), levels(d$site))
  }
  clay <- per_site()
  rain <- per_site()
  noise <- per_site()
  d$clay <- as.numeric(clay[as.character(d$site)])
  d$rain <- as.numeric(rain[as.character(d$site)])
  d$noise <- as.numeric(noise[as.character(d$site)])
  # A within-site numeric fixed effect. With only factor predictors the fit
  # takes the aggregated route, which carries no parsed model and no data.
  d$plot_z <- stats::rnorm(nrow(d))
  u <- stats::setNames(stats::rnorm(J, sd = 0.6), levels(d$site))
  d$yield <- 3 + 0.8 * d$clay + 0.5 * d$rain + 0.2 * d$plot_z +
    as.numeric(u[as.character(d$site)]) + stats::rnorm(nrow(d), sd = 0.7)
  d
}

fb_sel_verb_base <- function(d) {
  skip_on_cran()
  skip_if_not_installed("INLA")
  flexybayes(yield ~ plot_z, random = ~ site, data = d,
             backend = "inla", verbose = FALSE)
}


# --- the happy path ---------------------------------------------------- #

test_that("fb_select returns an increment over a protected-only baseline", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d),
    candidates = ~ clay + rain + noise,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    workflow = "unpenalized_benchmark",
    seed = 20260912L
  )

  expect_s3_class(sel, "fb_select")
  expect_true(all(c("protected_only", "full") %in% sel$scores$model))
  # The baseline is scored on the SAME folds, not on a separate run, which is
  # what makes the difference an increment rather than two numbers.
  expect_identical(nrow(sel$scores), 2L)
  expect_identical(
    sel$increment$improvement,
    sel$increment$protected_only - sel$increment$full
  )
  # Three real signals against one noise column: the candidates must help.
  expect_gt(sel$increment$improvement, 0)
})


test_that("a candidate block of pure noise does not earn a large increment", {
  # The negative control. Without it the test above passes for a function
  # that always reports an improvement.
  d <- fb_sel_verb_data()
  real <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    workflow = "unpenalized_benchmark", seed = 20260912L
  )
  noise_only <- fb_select(
    fb_sel_verb_base(d), candidates = ~ noise,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    workflow = "unpenalized_benchmark", seed = 20260912L
  )
  expect_gt(real$increment$improvement, noise_only$increment$improvement)
})


test_that("the same seed gives the same folds and the same answer", {
  d <- fb_sel_verb_data()
  base <- fb_sel_verb_base(d)
  args <- list(
    base = base, candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    workflow = "unpenalized_benchmark", seed = 4242L
  )
  a <- do.call(fb_select, args)
  b <- do.call(fb_select, args)
  expect_equal(a$increment$improvement, b$increment$improvement,
               tolerance = 1e-12)
})


# --- the gate ----------------------------------------------------------- #

test_that("a plain list of tuning arguments is not a penalty", {
  # The earlier release's shape. It now refuses because a penalty is an
  # <fb_penalty> that validated at construction, not because penalties are
  # unimplemented -- and the message says which.
  d <- fb_sel_verb_data()
  expect_error(
    fb_select(fb_sel_verb_base(d), ~ clay, penalty = list(alpha = 1),
              resamples = fb_resamples(unit = ~ site, v = 6L), seed = 1L),
    class = "flexybayes_refusal_selection_penalty_spec_malformed"
  )
  expect_error(
    fb_select(fb_sel_verb_base(d), ~ clay, penalty = list(alpha = 1),
              resamples = fb_resamples(unit = ~ site, v = 6L), seed = 1L),
    "must be NULL or an <fb_penalty>",
    fixed = TRUE
  )
})


test_that("a penalty on the unpenalised benchmark arm refuses by name", {
  d <- fb_sel_verb_data()
  expect_error(
    fb_select(
      fb_sel_verb_base(d), ~ clay,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 1)),
      workflow = "unpenalized_benchmark",
      resamples = fb_resamples(unit = ~ site, v = 6L), seed = 1L
    ),
    "has no penalty term"
  )
})


test_that("a non-flexyBayes object refuses before anything runs", {
  expect_error(fb_select(structure(list(), class = "lm"), ~ x),
               "must be a fitted flexyBayes model")
})


test_that("an aggregated-route fit refuses by name, not by falling over", {
  # An intercept-only fixed model takes the aggregated route and carries no
  # parsed model and no data. It must say so rather than fail somewhere
  # downstream with a message about row counts.
  skip_on_cran()
  skip_if_not_installed("INLA")
  # A Gaussian model is fitted per row, so the aggregated route is taken by
  # a count response.
  d <- fb_sel_verb_data()
  d$count <- stats::rpois(nrow(d), 3)
  agg <- flexybayes(count ~ 1, random = ~ site, data = d, family = "poisson",
                    backend = "inla", aggregate = TRUE, verbose = FALSE)
  expect_s3_class(agg, "flexybayes_aggregated")
  expect_error(
    fb_select(agg, ~ clay, resamples = fb_resamples(unit = ~ site, v = 6L),
              workflow = "unpenalized_benchmark", seed = 1L),
    "does not carry the parsed model"
  )
})


# --- the correlated random slope, and why the engine decides it --------- #
#
# Finding S01's defect D1. `(x | g)` on INLA lowers to an iid2d field whose
# slope variance comes back 32 per cent above REML's at the posterior mode,
# and whose intercept-slope correlation moves from 0.02 to 0.969 when the
# response is multiplied by ten -- a transformation that changes no model.
# The selection substrate whitens by the base's own covariance, so a base
# carrying that block sets the rotation every score is read through.
#
# The refusal is engine-specific, and both halves are asserted. A gate that
# refused the structure everywhere would be simpler and wrong: brms writes
# `(x | g)` natively and the measurement says nothing against it.

fb_sel_slope_data <- function(seed = 20260913L, J = 12L, m = 6L) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(sprintf("G%02d", seq_len(J)), each = m)))
  d$x <- stats::rnorm(nrow(d))
  u <- stats::rnorm(J, sd = 0.8)
  v <- stats::rnorm(J, sd = 0.4)
  gi <- as.integer(d$g)
  d$clay <- stats::rnorm(nrow(d))
  d$y <- 2 + 0.6 * d$x + 0.4 * d$clay + u[gi] + v[gi] * d$x +
    stats::rnorm(nrow(d), sd = 0.5)
  d
}


test_that("a correlated random slope cannot reach selection through an INLA base", {
  # INLA refuses the structure outright since 2026-09-13, so no new INLA
  # base carries it.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- fb_sel_slope_data()
  expect_error(
    suppressWarnings(suppressMessages(flexybayes(
      y ~ x + (x | g), data = d, backend = "inla", verbose = FALSE
    ))),
    class = "flexybayes_refusal_correlated_covariance_requires_brms"
  )

  # The selection gate keeps its own refusal for a base fitted on INLA
  # before that, which a saved object can still be. The stub carries what
  # the gate reads: the class, the parsed model and the data.
  fb <- fb_from_brms(y ~ x + (x | g), data = d, family = "gaussian")
  fb$link <- fb$link %||% "identity"
  base <- structure(
    list(fb = fb, data = d),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  expect_error(
    flexyBayes:::.fb_select_gate(base),
    class = "flexybayes_refusal_selection_base_not_allowlisted"
  )
  # The refusal carries the measurement that put the structure out, not a
  # bare verdict, and it names the route that is admitted.
  expect_error(
    flexyBayes:::.fb_select_gate(base),
    "46.4 at the posterior mode against 35.08 in ASReml",
    fixed = TRUE
  )
  expect_error(flexyBayes:::.fb_select_gate(base), "0.02 to 0.969",
               fixed = TRUE)
  expect_error(flexyBayes:::.fb_select_gate(base), "backend = \"brms\"",
               fixed = TRUE)
})


test_that("a brms base carrying the same term stays admitted", {
  # The negative control for the refusal above. Without it the gate passes
  # for an implementation that refuses `(x | g)` on every engine, which
  # would withdraw a capability the measurement does not touch.
  skip_on_cran()
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- fb_sel_slope_data()
  base <- suppressWarnings(suppressMessages(fb_brms(
    y ~ x + (x | g), data = d, chains = 1L, n_samples = 200L,
    warmup = 200L, seed = 20260913L, verbose = FALSE, mcmc_verbose = FALSE
  )))
  expect_true("simple_slope_cor" %in%
                vapply(base$fb$random_terms, function(t) t$type, character(1L)))
  expect_true(flexyBayes:::.fb_select_gate(base))
})


# --- the certificate ----------------------------------------------------- #

test_that("the certificate records the covariance, the estimand and the unit", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    workflow = "unpenalized_benchmark", seed = 20260912L
  )
  cert <- sel$certificate
  expect_s3_class(cert, "fb_select_certificate")
  # Read off the fold object, not off an argument. The public argument used
  # to default to `plugin_mean_covariance` while the folds worked in
  # `reml_plugin_covariance`, so the certificate described a computation
  # that had not happened. It now cannot: the field is whatever the
  # substrates carry, and a run whose folds disagreed reports both.
  expect_identical(cert$covariance$estimand, "reml_plugin_covariance")
  expect_identical(cert$covariance$builder, "block_cholesky")
  expect_identical(cert$workflow, "unpenalized_benchmark")
  expect_identical(cert$provider, "gls_reml_plugin")
  expect_identical(cert$resampling$prediction, "new_cluster")
  expect_identical(cert$resampling$seed, 20260912L)
  # Convergence is the tri-state verdict, not a recomputed local rule.
  expect_identical(cert$base$convergence, .fb_convergence_verdict(sel$base)$ok)
  expect_true(is.character(cert$approximations))
})


test_that("glance reports the increment and how it was scored", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    workflow = "unpenalized_benchmark", seed = 20260912L
  )
  g <- glance(sel)
  expect_identical(nrow(g), 1L)
  expect_true(all(c(
    "metric", "protected_only", "full", "improvement", "n_units",
    "se_unit", "prediction", "covariance", "builder"
  ) %in% names(g)))
  expect_identical(g$n_units, 24L)
})


test_that("print states the estimand rather than only the number", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    workflow = "unpenalized_benchmark", seed = 20260912L
  )
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "reml_plugin_covariance", fixed = TRUE)
  expect_match(out, "block_cholesky", fixed = TRUE)
  expect_match(out, "new_cluster", fixed = TRUE)
  expect_match(out, "independent", fixed = TRUE)
  # And it never calls this route Bayesian.
  expect_match(out, "unpenalized_benchmark", fixed = TRUE)
  expect_match(out, "gls_reml_plugin", fixed = TRUE)
})
