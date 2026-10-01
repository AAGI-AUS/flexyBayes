# =============================================================================
# The unpenalised Bayesian reference: a per-fold engine refit scored on
# posterior predictive draws.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 2, 3.4, 7, 9.2,
# 10 and 12 (the "unpenalised Bayesian reference" bullet).
#
# What is under test, and why each gate is here rather than somewhere cheaper.
#
#   The route is genuinely Bayesian. The covariance is integrated over by the
#   engine, not estimated once per fold and plugged in, and the estimand that
#   says so is read off the fold object rather than off an argument.
#   Section 9.2 asks for that because the argument and the fold disagreed once
#   already: `fb_select()`'s public `covariance = "plugin_mean_covariance"`
#   was copied into the certificate while the folds ran on
#   `reml_plugin_covariance`. The last block of this file holds a checker
#   green on the fold-derived value and red on the argument-derived one, so
#   the old failure mode is a test rather than a memory.
#
#   The fold boundary holds. The B0 leakage falsifier is run against THIS
#   route's builder, and an injected leak -- the assessment responses left in
#   the vector handed to the engine -- must turn it red.
#
#   The prediction target is real. Under `new_cluster` the held-out unit's
#   random effect must not be conditioned on, which on a design where the
#   folds SPLIT a unit is a property of the data handed over rather than an
#   accident of the fold. The `known_cluster` arm must therefore score
#   strictly better on a dataset with a large cluster variance, and the test
#   below is the only place that difference is visible.
#
# Every engine test is `skip_on_cran()` and gated on its provider. The INLA
# arm fits in seconds; the brms arm compiles once through CmdStan.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

# A clustered Gaussian mixed model: one protected covariate, two candidates of
# which one carries a real effect and one is null.
#
# The default cluster and residual scales are 1.0 and 0.6. That ratio matters
# to the arithmetic and not only to the story: under `new_cluster` the
# held-out unit's random effect is integrated by Monte Carlo, so only the
# draws landing within a residual standard deviation of the observation
# contribute to a row's predictive density. Measured on this fixture at the
# larger 2.5 / 0.4 ratio, four repeated 2,000-draw scorings of one cluster
# spread over 3.6 nats; at 8,000 draws, 1.4. The tests that assert on a score
# difference therefore run at a ratio where the estimator is well conditioned,
# and the one test that needs the extreme ratio states its measured tolerance.
fb_bcv_data <- function(
  seed = 20260912L,
  n_cluster = 12L,
  per = 10L,
  sd_u = 1.0,
  sd_e = 0.6,
  beta_z1 = 1.5
) {
  .fb_sel_with_seed(seed, {
    n <- n_cluster * per
    site <- factor(rep(sprintf("S%02d", seq_len(n_cluster)), each = per))
    x1 <- stats::rnorm(n)
    z1 <- stats::rnorm(n)
    z2 <- stats::rnorm(n)
    u <- stats::rnorm(n_cluster, sd = sd_u)
    y <- 2 +
      0.8 * x1 +
      beta_z1 * z1 +
      u[as.integer(site)] +
      stats::rnorm(n, sd = sd_e)
    data.frame(y = y, x1 = x1, z1 = z1, z2 = z2, site = site)
  })
}


# Fits are the expensive part of this file, so each one is built once and
# reused. The cache is lazy: a fixture whose tests all skip is never built.
fb_bcv_cache <- new.env(parent = emptyenv())

fb_bcv_cached <- function(key, expr) {
  if (!exists(key, envir = fb_bcv_cache, inherits = FALSE)) {
    assign(key, expr, envir = fb_bcv_cache)
  }
  get(key, envir = fb_bcv_cache, inherits = FALSE)
}

fb_bcv_inla_base <- function() {
  fb_bcv_cached(
    "inla_base",
    flexybayes(
      y ~ x1,
      random = ~site,
      data = fb_bcv_data(),
      backend = "inla",
      verbose = FALSE
    )
  )
}

# `adapt_delta` is raised at the base and travels into every fold refit
# through `call_info$control`: at brms's 0.8 default this model reported a
# divergent transition, and a divergence in a fold refit is a fold whose
# posterior is not the one being scored.
fb_bcv_brms_base <- function() {
  fb_bcv_cached(
    "brms_base",
    flexybayes(
      y ~ x1,
      random = ~site,
      data = fb_bcv_data(),
      backend = "brms",
      chains = 2L,
      n_samples = 500L,
      warmup = 500L,
      control = list(adapt_delta = 0.95),
      seed = 4041L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  )
}

# The INLA arm carries the assessment rows into the fit with an NA response,
# so a v-fold run hands the engine 1/v of its rows unobserved and the
# package's high-missingness warning is a property of the ROUTE rather than of
# any model under test. It is asserted on its own below, once, and silenced in
# the engine tests so those report only what is new to them.
fb_bcv_quiet_missingness <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    .local_envir = envir
  )
}


# The builder the leakage gate calls, wrapped to the (data, split) signature.
fb_bcv_builder <- function(data, split, prediction = "new_cluster") {
  .fb_sel_bayes_fold_inputs(
    data,
    split,
    response = "y",
    protected = ~x1,
    candidates = ~ z1 + z2,
    random = "site",
    prediction = prediction
  )
}


# --- (i) the brms arm ----------------------------------------------------- #

test_that("the Bayesian route runs on brms over three folds", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_long_tier()
  withr::local_options(brms.backend = "cmdstanr")

  # A brms fit records its dataset as metadata only, so the data the folds are
  # cut from is passed in. It is the same frame the base was fitted to.
  cv <- fb_bcv_cached(
    "brms_cv",
    .fb_sel_bayes_cv(
      fb_bcv_brms_base(),
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      prediction = "new_cluster",
      metric = "elpd",
      engine = fb_engine("brms", chains = 2L, n_samples = 500L, warmup = 500L),
      seed = 20260912L,
      data = fb_bcv_data(),
      n_draws = 4000L,
      probe_seconds = 5
    )
  )

  # The four axes of section 2, resolved from what ran.
  expect_identical(cv$workflow, "bayes")
  expect_identical(cv$provider, "brms")
  expect_identical(cv$covariance_estimand, "integrated")
  expect_identical(cv$fits_per_fold, 2L)
  expect_identical(cv$timing$n_fits, 6L)
  expect_length(cv$folds, 3L)
  expect_setequal(cv$aggregated$model, c("protected_only", "full"))
  expect_true(all(is.finite(cv$aggregated$score)))
  expect_true(all(cv$clusters$n > 0L))
  # elpd is a log density and rmse a distance: neither may come back as the
  # other's default fill value.
  expect_true(all(is.finite(cv$clusters$elpd)))
  expect_true(all(cv$clusters$rmse > 0))
})


test_that("the brms predictive draws agree with brms's own response sampler", {
  # An independent oracle on the construction this module uses. The module
  # scores on the per-draw predictive MEAN and the per-draw residual scale,
  # because a mean and a density taken from two separate calls to
  # posterior_predict() would describe two different posterior draws. That
  # construction is only legitimate if it reproduces what posterior_predict()
  # itself returns, so it is measured against it rather than assumed.
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_long_tier()
  withr::local_options(brms.backend = "cmdstanr")

  d <- fb_bcv_data()
  base <- fb_bcv_brms_base()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), d)[[1L]]
  inputs <- fb_bcv_builder(d, split)

  fit <- fb_bcv_cached(
    "brms_fold_fit",
    .fb_sel_bayes_refit(
      .fb_sel_bayes_normalise_base(base, d),
      inputs$train_data,
      .fb_sel_bayes_fixed("y", "x1", c("z1", "z2")),
      "brms",
      engine = fb_engine("brms", chains = 2L, n_samples = 500L, warmup = 500L),
      seed = 909L
    )
  )

  predictive <- .fb_sel_posterior_predict_fold(
    fit,
    seq_len(nrow(inputs$assess_data)),
    inputs$assess_data,
    prediction = "new_cluster",
    n_draws = 4000L,
    seed = 909L
  )
  set.seed(909L)
  own <- brms::posterior_predict(
    fit$brms,
    newdata = inputs$assess_data,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )

  # Column means agree because the response sampler adds mean-zero noise; the
  # tolerance is Monte-Carlo, on a residual scale of about 0.6 against a
  # cluster scale of 1.5.
  expect_lt(max(abs(colMeans(own) - colMeans(predictive$mu))), 0.5)
  own_sd <- apply(own, 2L, stats::sd)
  ours_sd <- apply(
    predictive$mu +
      matrix(
        stats::rnorm(
          length(predictive$mu),
          sd = rep(
            predictive$sigma,
            ncol(predictive$mu)
          )
        ),
        nrow = nrow(predictive$mu)
      ),
    2L,
    stats::sd
  )
  expect_lt(max(abs(own_sd - ours_sd)), 0.5)
})


# --- (ii) the INLA arm ---------------------------------------------------- #

test_that("the Bayesian route runs on INLA over three folds", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  fb_bcv_quiet_missingness()

  cv <- fb_bcv_cached(
    "inla_cv",
    .fb_sel_bayes_cv(
      fb_bcv_inla_base(),
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      prediction = "new_cluster",
      metric = "elpd",
      seed = 20260912L,
      n_draws = 4000L,
      probe_seconds = 5
    )
  )

  expect_identical(cv$workflow, "bayes")
  expect_identical(cv$provider, "inla")
  expect_identical(cv$covariance_estimand, "integrated")
  expect_identical(cv$timing$n_fits, 6L)
  expect_length(cv$folds, 3L)
  expect_true(all(is.finite(cv$aggregated$score)))
  # Every fold's substrate declares the same estimand, and it is the fit's
  # rather than a default carried down from an argument.
  expect_true(all(vapply(
    cv$substrates,
    function(z) identical(z$covariance_estimand, "integrated"),
    logical(1L)
  )))
})


test_that("the INLA arm's held-out sites draw the marginal random effect", {
  # The mechanism the new_cluster estimand rests on, measured rather than
  # assumed: a site whose rows all carry an NA response has nothing in the
  # likelihood to inform its effect, so its predictive spread must be the
  # random-effect scale and not the much smaller conditional one.
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()

  d <- fb_bcv_data(n_cluster = 8L, per = 8L, sd_u = 2.5, sd_e = 0.4)
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  inputs <- fb_bcv_builder(d, split)
  fit <- .fb_sel_bayes_refit(
    fb_bcv_cached(
      "inla_small_base",
      flexybayes(
        y ~ x1,
        random = ~site,
        data = d,
        backend = "inla",
        verbose = FALSE
      )
    ),
    inputs$engine_data,
    .fb_sel_bayes_fixed("y", "x1", c("z1", "z2")),
    "inla",
    na_action = "augment"
  )
  predictive <- .fb_sel_posterior_predict_fold(
    fit,
    split$assess_rows,
    inputs$engine_data,
    prediction = "new_cluster",
    n_draws = 4000L,
    seed = 5L
  )

  expect_identical(dim(predictive$mu), c(4000L, length(split$assess_rows)))
  expect_length(predictive$sigma, 4000L)
  expect_true(all(predictive$sigma > 0))
  # The held-out rows' predictive means move with the random effect, which on
  # this design has a standard deviation of 2.5 against a residual 0.4.
  expect_gt(mean(apply(predictive$mu, 2L, stats::sd)), 1)
})


test_that("the INLA fold refit reports its own missing fraction", {
  # A property of the route, asserted once so the engine tests can silence it.
  # The assessment rows enter the training refit with an NA response, so a
  # v-fold run hands the engine 1/v of its rows unobserved; at v = 3 that is
  # 33%, above the threshold at which the package warns that a variance
  # component can sit near a boundary. The warning is correct and the route
  # does not suppress it -- a reader choosing v is choosing that fraction.
  skip_on_cran()
  skip_if_not_installed("INLA")
  # The note is emit-once for the session, so a sibling file that exercises
  # this route first would consume it and leave this assertion reading an
  # empty warning list. The latch is reset for the scope of this test and
  # restored on exit.
  local_clean_emit_state()

  d <- fb_bcv_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), d)[[1L]]
  inputs <- fb_bcv_builder(d, split)
  expect_warning(
    .fb_sel_bayes_refit(
      fb_bcv_inla_base(),
      inputs$engine_data,
      .fb_sel_bayes_fixed("y", "x1", c("z1", "z2")),
      "inla",
      na_action = "augment"
    ),
    "carry no observed response",
    fixed = TRUE
  )
})


# --- (iii) the increment -------------------------------------------------- #

test_that("the increment is the protected-only score less the full score", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  fb_bcv_quiet_missingness()

  cv <- fb_bcv_cached(
    "inla_cv",
    .fb_sel_bayes_cv(
      fb_bcv_inla_base(),
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      prediction = "new_cluster",
      metric = "elpd",
      seed = 20260912L,
      n_draws = 4000L,
      probe_seconds = 5
    )
  )

  agg <- cv$aggregated
  protected <- agg$score[agg$model == "protected_only"]
  full <- agg$score[agg$model == "full"]
  increment <- protected - full

  # Both arms were scored on the SAME folds -- one fold list, one loop -- so
  # the difference is an increment rather than two numbers from two runs.
  expect_identical(
    sort(unique(cv$clusters$fold[cv$clusters$model == "protected_only"])),
    sort(unique(cv$clusters$fold[cv$clusters$model == "full"]))
  )
  expect_identical(increment, protected - full)
  # `score` is a loss on every metric -- elpd is carried negated -- so a
  # positive increment means the candidates helped. z1 carries an effect of
  # 1.5 against a residual scale of 0.6, which on ten-row clusters is worth
  # several nats per unit: measured at 3.2 here and 3.6 on a second fold
  # assignment, against a null-candidate control of -0.1. The threshold is set
  # well inside that so it is not a restatement of one run's number.
  expect_identical(cv$score_orientation, "loss")
  expect_gt(increment, 1)
  # And the raw per-unit elpd is kept unnegated, so a reader can see the
  # quantity itself rather than only its loss form.
  expect_lt(
    sum(cv$clusters$elpd[cv$clusters$model == "protected_only"]),
    sum(cv$clusters$elpd[cv$clusters$model == "full"])
  )
})


test_that("a null candidate block does not earn the increment a real one does", {
  # The negative control. Without it the test above passes for a route that
  # always reports an improvement.
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  fb_bcv_quiet_missingness()

  real <- fb_bcv_cached(
    "inla_cv",
    .fb_sel_bayes_cv(
      fb_bcv_inla_base(),
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      prediction = "new_cluster",
      metric = "elpd",
      seed = 20260912L,
      n_draws = 4000L,
      probe_seconds = 5
    )
  )
  null_only <- fb_bcv_cached(
    "inla_cv_null",
    .fb_sel_bayes_cv(
      fb_bcv_inla_base(),
      candidates = ~z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      prediction = "new_cluster",
      metric = "elpd",
      seed = 20260912L,
      n_draws = 4000L,
      probe_seconds = 5
    )
  )

  gain <- function(cv) {
    agg <- cv$aggregated
    agg$score[agg$model == "protected_only"] - agg$score[agg$model == "full"]
  }
  expect_gt(gain(real), gain(null_only))
})


# --- (iv) the leakage falsifier ------------------------------------------- #

test_that("the Bayesian fold builder passes the B0 leakage gate", {
  d <- fb_bcv_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  report <- .fb_sel_leakage_gate(fb_bcv_builder, d, split)

  expect_true(report$clean)
  expect_length(report$leaked, 0L)
  # Nothing inert: every declared quantity moved when the TRAINING rows moved,
  # so the gate is looking at something rather than at a set of constants.
  expect_length(report$inert, 0L)
  expect_setequal(
    report$checked,
    c("engine_response", "train_response", "train_xtx", "train_xty")
  )
  expect_silent(.fb_sel_certify_no_leakage(fb_bcv_builder, d, split))
})


test_that("FALSIFIER: assessment rows left in training turn the gate red", {
  # The injected defect is the one this route is most exposed to. INLA has no
  # newdata surface, so the assessment rows are carried INTO the fit with the
  # response set to NA; a builder that forgot the NA would hand the engine the
  # answers and every score in the run would be optimistic with nothing else
  # to show for it. The gate must see that at the exact object the engine
  # reads.
  d <- fb_bcv_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]

  leaky <- function(data, split) {
    inputs <- fb_bcv_builder(data, split)
    all_rows <- seq_len(nrow(data))
    y <- as.numeric(data$y)
    x <- cbind(
      .fb_sel_design(~x1, data)$X,
      .fb_sel_design(~ z1 + z2, data)$X[, c("z1", "z2"), drop = FALSE]
    )
    inputs$learned$engine_response <- y
    inputs$learned$train_response <- y[all_rows]
    inputs$learned$train_xtx <- crossprod(x)
    inputs$learned$train_xty <- crossprod(x, y)
    inputs
  }

  report <- .fb_sel_leakage_gate(leaky, d, split)
  expect_false(report$clean)
  expect_setequal(
    report$leaked,
    c("engine_response", "train_response", "train_xtx", "train_xty")
  )
  expect_error(
    .fb_sel_certify_no_leakage(leaky, d, split),
    "were not learned on its training rows",
    fixed = TRUE
  )
})


test_that("FALSIFIER: the candidate design read off all rows turns it red", {
  # A narrower leak than the one above and a more plausible one: the response
  # partition is right, but the design the fit is built on was read from every
  # row rather than from the training rows.
  d <- fb_bcv_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]

  leaky <- function(data, split) {
    inputs <- fb_bcv_builder(data, split)
    x <- cbind(
      .fb_sel_design(~x1, data)$X,
      .fb_sel_design(~ z1 + z2, data)$X[, c("z1", "z2"), drop = FALSE]
    )
    inputs$learned$train_xtx <- crossprod(x)
    inputs
  }

  report <- .fb_sel_leakage_gate(leaky, d, split)
  expect_false(report$clean)
  expect_setequal(report$leaked, "train_xtx")
})


test_that("FALSIFIER: an inert declared quantity is reported, not passed", {
  # A gate that only asks "did it move with the assessment rows" stays green
  # on a hard-coded constant, which is how a gate comes to pass vacuously.
  d <- fb_bcv_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]

  inert <- function(data, split) {
    inputs <- fb_bcv_builder(data, split)
    inputs$learned$train_response <- rep(0, length(split$train_rows))
    inputs
  }

  report <- .fb_sel_leakage_gate(inert, d, split)
  expect_false(report$clean)
  expect_length(report$leaked, 0L)
  expect_identical(report$inert, "train_response")
  expect_error(
    .fb_sel_certify_no_leakage(inert, d, split),
    "the gate cannot see them",
    fixed = TRUE
  )
})


test_that("the engine data carries NA at every assessment row and only there", {
  # The claim the gate rests on, asserted directly as well: the response
  # handed to the engine is missing exactly where the fold is scored.
  d <- fb_bcv_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  inputs <- fb_bcv_builder(d, split)

  expect_identical(which(is.na(inputs$engine_data$y)), split$assess_rows)
  expect_identical(
    as.numeric(inputs$engine_data$y[split$train_rows]),
    as.numeric(d$y[split$train_rows])
  )
  # And under new_cluster the assessment rows' grouping level is one no
  # training row carries, so the engine has nothing to condition on.
  expect_length(
    intersect(
      as.character(inputs$engine_data$site[split$assess_rows]),
      as.character(inputs$engine_data$site[split$train_rows])
    ),
    0L
  )
  # Under known_cluster the labels are left exactly as the data had them.
  known <- fb_bcv_builder(d, split, prediction = "known_cluster")
  expect_identical(
    as.character(known$engine_data$site),
    as.character(d$site)
  )
})


# --- (v) new_cluster against known_cluster -------------------------------- #

test_that("known_cluster scores no worse than new_cluster when clusters split", {
  # The estimand's whole content. The folds here split each site rather than
  # holding it out, so `known_cluster` has training rows from the very cluster
  # it is scoring and `new_cluster` must not be allowed to use them. On a
  # design whose cluster standard deviation is 2.5 against a residual 0.4,
  # conditioning is worth most of the predictive density, and the two arms
  # would be indistinguishable if the relabelling were a no-op.
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  fb_bcv_quiet_missingness()

  d <- fb_bcv_data(n_cluster = 8L, per = 8L, sd_u = 2.5, sd_e = 0.4)
  base <- fb_bcv_cached(
    "inla_small_base",
    flexybayes(
      y ~ x1,
      random = ~site,
      data = d,
      backend = "inla",
      verbose = FALSE
    )
  )

  run <- function(prediction) {
    .fb_sel_bayes_cv(
      base,
      candidates = ~ z1 + z2,
      # No `unit`: the folds fall on rows, so each site has rows on both
      # sides of every split.
      resamples = fb_resamples(v = 2L),
      prediction = prediction,
      metric = "nlpd",
      seed = 31L,
      n_draws = 4000L,
      candidate_sets = list(full = 1:2),
      cluster = "site",
      probe_seconds = 5
    )
  }
  new_cluster <- run("new_cluster")
  known_cluster <- run("known_cluster")

  expect_lte(known_cluster$aggregated$score, new_cluster$aggregated$score)
  # Not merely no worse -- markedly better, which is what rules out the two
  # arms having computed the same thing.
  expect_lt(known_cluster$aggregated$score, new_cluster$aggregated$score - 0.2)
})


# --- (vi) determinism under a seed ---------------------------------------- #

test_that("the brms arm reproduces its draws under the same seed and chains", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_long_tier()
  withr::local_options(brms.backend = "cmdstanr")

  d <- fb_bcv_data()
  base <- fb_bcv_brms_base()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), d)[[1L]]
  inputs <- fb_bcv_builder(d, split)
  fit <- fb_bcv_cached(
    "brms_fold_fit",
    .fb_sel_bayes_refit(
      .fb_sel_bayes_normalise_base(base, d),
      inputs$train_data,
      .fb_sel_bayes_fixed("y", "x1", c("z1", "z2")),
      "brms",
      engine = fb_engine("brms", chains = 2L, n_samples = 500L, warmup = 500L),
      seed = 909L
    )
  )

  draw <- function() {
    .fb_sel_posterior_predict_fold(
      fit,
      seq_len(nrow(inputs$assess_data)),
      inputs$assess_data,
      prediction = "new_cluster",
      n_draws = 4000L,
      seed = 909L
    )
  }
  a <- draw()
  b <- draw()
  expect_identical(a$mu, b$mu)
  expect_identical(a$sigma, b$sigma)

  # A refit under the same seed and chains returns the same posterior, so the
  # whole arm is reproducible and not only its sampling step.
  again <- .fb_sel_bayes_refit(
    .fb_sel_bayes_normalise_base(base, d),
    inputs$train_data,
    .fb_sel_bayes_fixed("y", "x1", c("z1", "z2")),
    "brms",
    engine = fb_engine("brms", chains = 2L, n_samples = 500L, warmup = 500L),
    seed = 909L
  )
  expect_identical(as.matrix(fit$brms), as.matrix(again$brms))

  # A different seed moves the new-level draws, so the equality above is a
  # property of the seed rather than of a degenerate predictive.
  different <- .fb_sel_posterior_predict_fold(
    fit,
    seq_len(nrow(inputs$assess_data)),
    inputs$assess_data,
    prediction = "new_cluster",
    n_draws = 4000L,
    seed = 910L
  )
  expect_false(identical(a$mu, different$mu))
})


test_that("the INLA arm reproduces its score, and its draws do not", {
  # Measured, not assumed, and the negative half is the point.
  # `INLA::inla.posterior.sample()` on INLA 26.8.7 does not return the same
  # draws twice for the same `seed`, whether `num.threads = "1:1"`,
  # `parallel.configs = FALSE`, `skew.corr = FALSE`, `inla.mode = "classic"`
  # or an R-side `set.seed()` with `seed = 0L`; the underlying
  # `INLA::inla.qsample(seed = )` IS bit-reproducible, so the loss is above
  # it. The arm's reproducibility is therefore a Monte-Carlo statement about
  # the score, and this test states it as one rather than claiming an equality
  # the instrument cannot deliver.
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()

  d <- fb_bcv_data(n_cluster = 8L, per = 8L, sd_u = 2.5, sd_e = 0.4)
  base <- fb_bcv_cached(
    "inla_small_base",
    flexybayes(
      y ~ x1,
      random = ~site,
      data = d,
      backend = "inla",
      verbose = FALSE
    )
  )
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  inputs <- fb_bcv_builder(d, split)
  fit <- .fb_sel_bayes_refit(
    base,
    inputs$engine_data,
    .fb_sel_bayes_fixed("y", "x1", c("z1", "z2")),
    "inla",
    na_action = "augment"
  )

  draw <- function(seed) {
    .fb_sel_posterior_predict_fold(
      fit,
      split$assess_rows,
      inputs$engine_data,
      prediction = "new_cluster",
      n_draws = 8000L,
      seed = seed
    )
  }
  a <- draw(17L)
  b <- draw(17L)
  # The negative finding, recorded as an assertion so a future INLA that fixes
  # it turns this red and the comment above gets revisited.
  expect_false(identical(a$mu, b$mu))

  # What does hold: the scored quantity agrees to Monte-Carlo error. The
  # tolerance is measured rather than chosen -- four repeated scorings of this
  # fold spread over 3.6 nats at 2,000 draws and 1.4 at 8,000, on a per-cluster
  # elpd of about -43 -- and it is set at roughly twice the measured spread.
  # This design is the worst case on purpose: a cluster scale of 2.5 against a
  # residual 0.4 is where a marginal predictive integrated by Monte Carlo is
  # least well conditioned.
  score <- function(p) {
    .fb_sel_bayes_score_fold(
      inputs$assessment$y,
      p$mu,
      p$sigma,
      .fb_sel_key(d, "site")[split$assess_rows]
    )$clusters$elpd
  }
  expect_lt(max(abs(score(a) - score(b))), 3)
})


# --- (vii) the runtime budget --------------------------------------------- #

test_that("the budget refusal fires before a single fold is fitted", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  base <- fb_bcv_inla_base()
  fitted_folds <- 0L
  counting_refit <- function(fold_data, fixed, fit_seed) {
    fitted_folds <<- fitted_folds + 1L
    stop("the budget gate let a fold fit start.", call. = FALSE)
  }

  expect_error(
    .fb_sel_bayes_cv(
      base,
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      seed = 1L,
      budget = 0.5,
      # The per-fit cost is supplied rather than probed, so the refusal is
      # provably ahead of EVERY fit and not merely ahead of the loop.
      probe_seconds = 8,
      refit_fn = counting_refit
    ),
    class = "flexybayes_refusal_penalty_runtime_budget_exceeded"
  )
  expect_identical(fitted_folds, 0L)
})


test_that("the budget estimate is one timed fit times the fits the folds imply", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  cost <- .fb_sel_cost_estimate(
    fb_bcv_inla_base(),
    fb_resamples(unit = ~site, v = 5L),
    fits_per_fold = 2L,
    probe_seconds = 3
  )
  expect_identical(cost$n_folds, 5L)
  expect_identical(cost$n_fits, 10L)
  expect_identical(cost$seconds, 30)
  expect_identical(cost$source, "supplied")

  # Green: a budget the estimate fits inside returns the cost invisibly.
  expect_silent(.fb_sel_bayes_check_budget(cost, 60))
  expect_identical(.fb_sel_bayes_check_budget(cost, 60)$seconds, 30)
  # Red: the same estimate against a smaller budget refuses by name.
  expect_error(
    .fb_sel_bayes_check_budget(cost, 10),
    class = "flexybayes_refusal_penalty_runtime_budget_exceeded"
  )
  expect_error(
    .fb_sel_bayes_check_budget(cost, 10),
    "Nothing was fitted",
    fixed = TRUE
  )
  # And no budget is no gate.
  expect_silent(.fb_sel_bayes_check_budget(cost, NULL))
})


test_that("a probed estimate times a real fit rather than guessing", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  probes <- 0L
  cost <- .fb_sel_cost_estimate(
    fb_bcv_inla_base(),
    fb_resamples(unit = ~site, v = 4L),
    fits_per_fold = 2L,
    probe_fn = function() {
      probes <<- probes + 1L
      Sys.sleep(0.05)
      NULL
    }
  )
  expect_identical(probes, 1L)
  expect_identical(cost$source, "timed_probe")
  expect_gt(cost$per_fit_seconds, 0.04)
  expect_identical(cost$n_fits, 8L)
})


# --- the absent provider --------------------------------------------------- #

test_that("an absent provider refuses by name and substitutes nothing", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  base <- fb_bcv_inla_base()
  # The provider IS installed here, so the only way to see the refusal fire is
  # to substitute the one answer it turns on.
  testthat::local_mocked_bindings(
    .fb_sel_provider_installed = function(package) FALSE
  )
  expect_error(
    .fb_sel_bayes_provider(base),
    class = "flexybayes_refusal_selection_provider_absent"
  )
  expect_error(
    .fb_sel_bayes_provider(base),
    "No substitute provider is chosen here",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_bayes_cv(
      base,
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      seed = 1L
    ),
    class = "flexybayes_refusal_selection_provider_absent"
  )
})


test_that("the provider is read from the fitted class, not from a record", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  base <- fb_bcv_inla_base()
  # `call_info$backend` holds the routing POLICY the call asked for. Rewriting
  # it must not change which engine the route says fitted the base.
  base$extras$call_info$backend <- "brms"
  expect_identical(.fb_sel_bayes_provider(base), "inla")
})


# --- (viii) the certificate falsifier -------------------------------------- #

# The check the certificate has to pass: what a run DECLARES about its
# covariance estimand must be what the folds actually worked in. It is written
# as a function so it can be shown green on the fold-derived value and red on
# the argument-derived one.
fb_bcv_estimand_agrees <- function(cv, declared) {
  fold_estimand <- cv$covariance_estimand
  if (!identical(as.character(declared), as.character(fold_estimand))) {
    stop(
      "the certificate declares covariance estimand \"",
      declared,
      "\" and the folds ran on \"",
      paste(fold_estimand, collapse = ", "),
      "\".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


test_that("the benchmark route's estimand comes off the fold, not an argument", {
  d <- fb_bcv_data()
  cv <- .fb_sel_cv(
    data = d,
    response = "y",
    protected = ~x1,
    candidates = ~ z1 + z2,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 3L),
    metric = "rmse"
  )

  # The fold object says what it worked in, and `.fb_sel_cv()` has no
  # `covariance` argument at all, so there is nothing for it to have copied.
  expect_identical(cv$covariance_estimand, "reml_plugin_covariance")
  expect_identical(cv$workflow, "unpenalized_benchmark")
  expect_identical(cv$provider, "gls_reml_plugin")
  expect_false("covariance" %in% names(formals(.fb_sel_cv)))
  expect_true(all(vapply(
    cv$substrates,
    function(z) identical(z$covariance_estimand, "reml_plugin_covariance"),
    logical(1L)
  )))

  # GREEN: a certificate built from the fold object agrees with it.
  expect_silent(fb_bcv_estimand_agrees(cv, cv$covariance_estimand))
})


test_that("FALSIFIER: the old argument-fed certificate turns the check red", {
  # The defect the review found: `fb_select()`'s public argument defaulted to
  # "plugin_mean_covariance" and was copied straight into the certificate,
  # while every fold ran the independent REML plugin. Nothing in the object
  # disagreed out loud, because nothing compared the two. This is that
  # comparison, and it is red on the old value.
  d <- fb_bcv_data()
  cv <- .fb_sel_cv(
    data = d,
    response = "y",
    protected = ~x1,
    candidates = ~ z1 + z2,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 3L),
    metric = "rmse"
  )

  expect_error(
    fb_bcv_estimand_agrees(cv, "plugin_mean_covariance"),
    "the folds ran on \"reml_plugin_covariance\"",
    fixed = TRUE
  )
  # Any other argument value is caught the same way: the check reads the fold,
  # so no argument can talk it round.
  expect_error(
    fb_bcv_estimand_agrees(cv, "integrated"),
    "the certificate declares covariance estimand \"integrated\"",
    fixed = TRUE
  )

  # And an injected disagreement BETWEEN folds is reported as two estimands
  # rather than silently reduced to the first one.
  mixed <- cv
  mixed$substrates[[1L]]$covariance_estimand <- "plugin_mean_covariance"
  mixed$covariance_estimand <- unique(vapply(
    mixed$substrates,
    function(z) z$covariance_estimand,
    character(1L)
  ))
  expect_length(mixed$covariance_estimand, 2L)
  expect_error(
    fb_bcv_estimand_agrees(mixed, "reml_plugin_covariance"),
    "the folds ran on \"plugin_mean_covariance, reml_plugin_covariance\"",
    fixed = TRUE
  )
})


test_that("the Bayesian route's estimand is integrated, and says so", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fb_bcv_quiet_missingness()

  cv <- fb_bcv_cached(
    "inla_cv",
    .fb_sel_bayes_cv(
      fb_bcv_inla_base(),
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      prediction = "new_cluster",
      metric = "elpd",
      seed = 20260912L,
      n_draws = 4000L,
      probe_seconds = 5
    )
  )
  expect_silent(fb_bcv_estimand_agrees(cv, "integrated"))
  expect_error(
    fb_bcv_estimand_agrees(cv, "reml_plugin_covariance"),
    "the folds ran on \"integrated\"",
    fixed = TRUE
  )
})


# --- scoring and aggregation, without an engine ---------------------------- #

test_that("the log predictive density is a log-mean-exp over draws", {
  # An arithmetic oracle: with a fixed mean and scale on every draw, the
  # log-mean-exp collapses to the Gaussian log density itself, which is
  # computable in closed form and independent of this module.
  y <- c(1.2, -0.4, 2.0)
  mu <- matrix(rep(c(1.0, 0.0, 1.5), each = 50L), nrow = 50L)
  sigma <- rep(0.8, 50L)
  scored <- .fb_sel_bayes_score_fold(y, mu, sigma, rep("A", 3L))

  expect_equal(
    scored$rows$lpd_row,
    stats::dnorm(y, mean = c(1.0, 0.0, 1.5), sd = 0.8, log = TRUE),
    tolerance = 1e-12
  )
  # With no draw-to-draw variation the rows are conditionally independent, so
  # the joint density is the product and `nlpd` is the row mean of -lpd.
  expect_equal(
    scored$clusters$nlpd,
    -mean(scored$rows$lpd_row),
    tolerance = 1e-12
  )
  expect_equal(
    scored$clusters$elpd,
    sum(scored$rows$lpd_row),
    tolerance = 1e-12
  )
  expect_equal(
    scored$clusters$rmse,
    sqrt(mean((y - c(1.0, 0.0, 1.5))^2)),
    tolerance = 1e-12
  )
})


test_that("the log-mean-exp helper matches its own definition", {
  x <- c(-800, -801, -799.5)
  expect_equal(
    .fb_sel_logmeanexp(x),
    log(mean(exp(x - max(x)))) + max(x),
    tolerance = 1e-12
  )
  # And it survives a scale that overflows the naive expression.
  expect_true(is.finite(.fb_sel_logmeanexp(c(-1e5, -1e5 + 1))))
  expect_false(is.finite(log(mean(exp(c(-1e5, -1e5 + 1))))))
})


test_that("elpd aggregates at the unit and is reported as a loss", {
  clusters <- data.table::data.table(
    model = "full",
    fold = rep(c("1", "2"), each = 3L),
    repeat_id = 1L,
    cluster = c("A", "B", "C", "D", "E", "F"),
    n = 4L,
    mse = 1,
    rmse = 1,
    elpd = c(-10, -12, -8, -11, -9, -13),
    nlpd = 2
  )
  summary <- .fb_sel_aggregate_bayes_scores(clusters, metric = "elpd")

  expect_identical(summary$n_units, 6L)
  expect_equal(
    summary$score,
    -mean(c(-10, -12, -8, -11, -9, -13)),
    tolerance = 1e-12
  )
  expect_equal(
    summary$se,
    stats::sd(clusters$elpd) / sqrt(6),
    tolerance = 1e-12
  )
  expect_identical(summary$orientation, "loss")
  expect_equal(summary$elpd_per_unit[["A"]], -10, tolerance = 1e-12)

  # The three loss-oriented metrics go through the benchmark route's own
  # aggregator, so the two workflows cannot drift apart on how a unit is
  # averaged.
  rmse <- .fb_sel_aggregate_bayes_scores(clusters, metric = "rmse")
  expect_identical(
    rmse[c("metric", "unit", "score", "se", "n_units")],
    .fb_sel_aggregate_scores(clusters, metric = "rmse")[
      c("metric", "unit", "score", "se", "n_units")
    ]
  )
})


test_that("a single independent unit refuses rather than reporting an SE", {
  clusters <- data.table::data.table(
    model = "full",
    fold = "1",
    repeat_id = 1L,
    cluster = "A",
    n = 4L,
    mse = 1,
    rmse = 1,
    elpd = -10,
    nlpd = 2
  )
  expect_error(
    .fb_sel_aggregate_bayes_scores(clusters, metric = "elpd"),
    class = "flexybayes_refusal_selection_scoring_unit_undeclared"
  )
})


test_that("mismatched draw dimensions are an error, not a recycled score", {
  expect_error(
    .fb_sel_bayes_score_fold(
      c(1, 2, 3),
      matrix(0, nrow = 10L, ncol = 2L),
      rep(1, 10L),
      rep("A", 3L)
    ),
    "10 x 2",
    fixed = TRUE
  )
})


# --- shapes the route refuses ---------------------------------------------- #

test_that("a grouping column that also carries a fixed effect refuses", {
  d <- fb_bcv_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  expect_error(
    .fb_sel_bayes_fold_inputs(
      d,
      split,
      response = "y",
      protected = ~ x1 + site,
      candidates = ~ z1 + z2,
      random = "site",
      prediction = "new_cluster"
    ),
    "would change the fixed design",
    fixed = TRUE
  )
  # It is only new_cluster that relabels, so known_cluster is unaffected.
  expect_s3_class(
    .fb_sel_bayes_fold_inputs(
      d,
      split,
      response = "y",
      protected = ~x1,
      candidates = ~ z1 + z2,
      random = "site",
      prediction = "known_cluster"
    ),
    "fb_sel_bayes_inputs"
  )
})


test_that("a fit class with no predictive route says so rather than guessing", {
  # Dispatch is a plain inherits() branch over the two engines rather than an
  # S3 generic, so this is the branch's own refusal and not a default method.
  expect_error(
    .fb_sel_posterior_predict_fold(
      structure(list(), class = "not_a_flexybayes_fit"),
      1L,
      data.frame(y = 1)
    ),
    "no posterior predictive route",
    fixed = TRUE
  )
})


test_that("the INLA arm samples the scored rows' predictors only", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fb_bcv_quiet_missingness()

  d <- fb_bcv_data(n_cluster = 8L, per = 8L)
  base <- fb_bcv_cached(
    "inla_default_small_base",
    flexybayes(y ~ x1, random = ~site, data = d, backend = "inla",
               verbose = FALSE)
  )
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  inputs <- fb_bcv_builder(d, split)
  fit <- .fb_sel_bayes_refit(
    base,
    inputs$engine_data,
    .fb_sel_bayes_fixed("y", "x1", c("z1", "z2")),
    "inla",
    na_action = "augment"
  )
  seen <- new.env()
  real <- .fb_inla_posterior_sample
  testthat::local_mocked_bindings(
    .fb_inla_posterior_sample = function(object, n, ...) {
      seen$args <- list(...)
      draws <- real(object, n, ...)
      seen$latent_rows <- rownames(draws[[1L]]$latent)
      draws
    }
  )
  predictive <- .fb_sel_posterior_predict_fold(
    fit,
    split$assess_rows,
    inputs$engine_data,
    prediction = "new_cluster",
    n_draws = 300L,
    seed = 3L
  )
  expect_identical(
    seen$args$selection,
    list(Predictor = as.integer(split$assess_rows))
  )
  expect_identical(
    seen$latent_rows,
    paste0("Predictor:", split$assess_rows)
  )
  expect_identical(dim(predictive$mu), c(300L, length(split$assess_rows)))
  # Each column is its own row's predictor: its mean over the draws is the
  # fit's posterior mean for that row, to Monte-Carlo error.
  expected <- fit$inla$summary.linear.predictor$mean[split$assess_rows]
  mc_error <- max(apply(predictive$mu, 2L, stats::sd)) / sqrt(300)
  expect_lt(max(abs(colMeans(predictive$mu) - expected)), 4 * mc_error)
})


test_that("fold scoring defaults to the draw count set from Monte-Carlo error", {
  expect_identical(eval(formals(.fb_sel_bayes_cv)$n_draws), 8000L)
  expect_identical(
    eval(formals(.fb_sel_posterior_predict_fold)$n_draws),
    8000L
  )
})


test_that("the default draw count moves the score by less than its Monte-Carlo error at 4,000", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  fb_bcv_quiet_missingness()
  # 800 rows over 20 sites with a weak candidate effect. Twelve replicate
  # scorings at 4,000 draws gave a standard deviation of 0.164 nats on the
  # aggregated score, and six at 8,000 gave 0.070.
  d <- .fb_sel_with_seed(20261001L, {
    n <- 800L
    out <- data.frame(site = factor(rep(sprintf("S%02d", 1:20), each = 40L)))
    for (v in c("x", "z1", "z2", "z3")) out[[v]] <- stats::rnorm(n)
    u <- stats::rnorm(20L, sd = 0.8)
    out$y <- 2 + 0.5 * out$x + 0.15 * out$z1 + u[as.integer(out$site)] +
      stats::rnorm(n, sd = 0.5)
    out
  })
  base <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                     verbose = FALSE)
  spec <- .fb_sel_bayes_base_spec(base)
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), spec$data,
                         seed = 1L)
  key <- .fb_sel_key(spec$data, "site")
  fold_fits <- lapply(folds, function(split) {
    inputs <- .fb_sel_bayes_fold_inputs(
      spec$data, split, response = "y", protected = ~x,
      candidates = ~ z1 + z2 + z3, random = "site", prediction = "new_cluster"
    )
    fit <- .fb_sel_bayes_refit(
      base, inputs$engine_data, y ~ x + z1 + z2 + z3, "inla",
      na_action = "augment"
    )
    list(fit = fit, inputs = inputs, split = split)
  })
  score <- function(n_draws) {
    clusters <- lapply(fold_fits, function(z) {
      p <- .fb_sel_posterior_predict_fold(
        z$fit, z$split$assess_rows, z$inputs$engine_data,
        prediction = "new_cluster", n_draws = n_draws, seed = 11L
      )
      .fb_sel_bayes_score_fold(
        z$inputs$assessment$y, p$mu, p$sigma, key[z$split$assess_rows],
        label = z$split$label
      )$clusters
    })
    .fb_sel_aggregate_bayes_scores(
      data.table::rbindlist(clusters), metric = "elpd"
    )$score
  }
  # Three standard deviations of the difference of two independent scorings.
  expect_lt(abs(score(8000L) - score(4000L)), 3 * sqrt(0.164^2 + 0.070^2))
})
