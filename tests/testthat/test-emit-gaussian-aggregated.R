# test-emit-gaussian-aggregated.R -- Stage 3A backend wiring tests.
# Covers the ten ADR 0022 §Consequences gates (a) through (j) for the
# aggregated-gaussian emit path. (a)/(b)/(c) -- bit-exact log-likelihood
# equivalence -- live in test-aggregate-gaussian.R (foundation). The
# gates exercised here:
#   (e) plan eligibility gates dispatch
#   (f) out-of-scope IRs fall through to per-row under aggregate = "auto"
#   (g) heterogeneous residual at_units cell-key validity
#   (h) exactness field on both per-row and aggregated fits
#   (i) compression-ratio line gating (ratio >= 2)
#   (j) refusal taxonomy (greta + aggregate=TRUE; brms + aggregate=TRUE;
#       smooth + aggregate=TRUE)
# Plus posterior-equivalence smoke tests for the INLA aggregated path
# (greta posterior-equivalence is deferred to v0.3.3 with the greta
# emit; see emit_gaussian_aggregated.R header).
#
# Reference for the algebraic identity gate: day-1 spikes at
# `spikes_stage3a/spike_inla.R` (matched-prior bit-exact agreement) and
# `spikes_stage3a/spike_greta.R` (algebraic viability of the cell-mean
# + gamma-on-WSS pattern; viable for v0.3.3).

suppressPackageStartupMessages({
  library(testthat)
})

# Common silenced options across the file.
old_opts <- options(
  flexyBayes.silence_default_prior_note = TRUE,
  flexyBayes.silence_uniform_inla_approx = TRUE,
  flexyBayes.silence_auto_fallback_note = TRUE,
  flexyBayes.silence_auto_inla_missing_note = TRUE
)
on.exit(options(old_opts), add = TRUE)


# ---------------------------------------------------------------- #
# Test fixture -- in-scope Stage 3A IR with productive compression  #
# ---------------------------------------------------------------- #
mk_agg_data <- function(seed = 1L, N = 500L, J = 15L) {
  set.seed(seed)
  dat <- data.frame(
    f1 = factor(sample(letters[1:3], N, replace = TRUE)),
    f2 = factor(sample(letters[1:3], N, replace = TRUE)),
    g = factor(sample.int(J, N, replace = TRUE))
  )
  dat$y <- 1.0 +
    0.5 * as.integer(dat$f1) -
    0.3 * as.integer(dat$f2) +
    rnorm(J, 0, 0.3)[as.integer(dat$g)] +
    rnorm(N, 0, 0.7)
  dat
}


# ---------------------------------------------------------------- #
# Gate (h): exactness field on every flexybayes fit                 #
# ---------------------------------------------------------------- #

test_that("(h) per-row fit carries exactness = 'exact'", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = FALSE,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_identical(fit$exactness, "exact")
})

test_that("(h) aggregated fit carries exactness = 'aggregated_exact'", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_identical(fit$exactness, "aggregated_exact")
})


# ---------------------------------------------------------------- #
# Gate (e): plan eligibility gates dispatch                         #
# ---------------------------------------------------------------- #

test_that("(e) eligible plan routes to aggregated_gaussian path", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "inla")
  expect_identical(bd$path, "aggregated_gaussian")
  expect_true(grepl("aggregation plan eligible", bd$reason))
})

test_that("(e) backend_decision reason carries N, K, ratio", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  reason <- backend_decision(fit)$reason
  expect_true(grepl("N = \\d+", reason))
  expect_true(grepl("K = \\d+", reason))
  expect_true(grepl("ratio = \\d+\\.\\d+:1", reason))
})


# ---------------------------------------------------------------- #
# Gate (f): out-of-scope IRs fall through cleanly                   #
# ---------------------------------------------------------------- #

test_that("(f) continuous fixed term falls through to per-row", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  dat$x_cont <- rnorm(nrow(dat))
  fit <- suppressMessages(fb(
    y ~ x_cont + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_identical(fit$exactness, "exact")
  bd <- backend_decision(fit)
  expect_true(bd$path %in% c("auto_accept", "explicit_inla_accept"))
})

test_that("(f) a Bernoulli binomial auto-aggregates exactly", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  dat$y_bin <- as.integer(dat$y > median(dat$y))
  fit <- suppressMessages(fb(
    y_bin ~ f1 + (1 | g),
    data = dat,
    family = "binomial",
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  # Binomial / poisson now share the gaussian exact-aggregation path
  # (count sufficient statistics). A 0/1 response with a factor cell key
  # routes to the aggregated count emit.
  expect_identical(fit$exactness, "aggregated_exact")
  expect_identical(backend_decision(fit)$path, "aggregated_count")
})


# ---------------------------------------------------------------- #
# Gate (j): refusal taxonomy                                        #
# ---------------------------------------------------------------- #

test_that("(j) aggregate = TRUE on out-of-scope IR raises refusal", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  dat$x_cont <- rnorm(nrow(dat))
  err <- tryCatch(
    suppressMessages(fb(
      y ~ x_cont + (1 | g),
      data = dat,
      backend = "inla",
      aggregate = TRUE,
      prior_vc_sd = 1,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_true(any(grepl(
    "continuous_cell_key|compression_unproductive",
    err$reason_codes
  )))
})

test_that("(j) aggregate = TRUE on backend = 'greta' now succeeds (v0.3.3)", {
  skip_if_greta_backend_unusable()
  # Fits via greta -> needs a usable TensorFlow backend, which CI does
  # not set up (absent on Windows even with the greta R package present).
  # Matches the skip_on_ci() guard on the suite's other greta-fit tests.
  testthat::skip_on_ci()
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + (1 | g),
    data = dat,
    backend = "greta",
    aggregate = TRUE,
    n_samples = 200,
    warmup = 100,
    chains = 1,
    prior_fixed_sd = 10,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$exactness, "aggregated_exact")
  expect_identical(backend_decision(fit)$backend, "greta")
  expect_identical(backend_decision(fit)$path, "aggregated_gaussian")
})

test_that("(j) aggregate = TRUE on backend = 'brms' refuses", {
  dat <- mk_agg_data()
  err <- tryCatch(
    suppressMessages(fb(
      y ~ f1 + (1 | g),
      data = dat,
      backend = "brms",
      aggregate = TRUE,
      n_samples = 100,
      warmup = 50,
      chains = 1,
      prior_vc_sd = 1,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    error = function(e) e
  )
  expect_true(grepl("aggregate = TRUE", conditionMessage(err)))
  expect_true(grepl("brms|greta and inla", conditionMessage(err)))
})

test_that("(j) aggregate = 'foo' (invalid) raises clear validation error", {
  dat <- mk_agg_data()
  err <- tryCatch(
    suppressMessages(fb(
      y ~ f1 + (1 | g),
      data = dat,
      backend = "inla",
      aggregate = "foo",
      prior_vc_sd = 1,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    error = function(e) e
  )
  expect_true(grepl("must be TRUE, FALSE, or \"auto\"", conditionMessage(err)))
})


# ---------------------------------------------------------------- #
# Gate (g): heterogeneous residual at_units cell-key validity       #
# ---------------------------------------------------------------- #
# Stage 3A allows at(f):units when f is in the cell key (so each cell
# has a single residual sigma). When f is NOT a cell key, the per-row
# / per-cell algebraic identity breaks and the emit refuses.

test_that("(g) at(f):units refuses when f is not a cell key", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  # f3 is NOT in the formula's cell key but referenced in rcov.
  dat$f3 <- factor(sample(letters[1:2], nrow(dat), replace = TRUE))
  # NOTE: this exercises the dispatch.R refusal for unsupported INLA
  # rcov shapes at v0.3.2 (heterogeneous at_units on INLA is deferred);
  # the at-not-in-cell-key refusal lives in the greta path (also
  # deferred). The test confirms the IR doesn't silently produce an
  # invalid aggregated fit.
  err <- tryCatch(
    suppressMessages(flexybayes(
      yield = y,
      fixed = ~ f1 + f2,
      random = ~g,
      rcov = ~ at(f3):units,
      data = dat,
      backend = "inla",
      aggregate = TRUE,
      prior_vc_sd = 1,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "error")
})


# ---------------------------------------------------------------- #
# Gate (i): compression-ratio line gating                           #
# ---------------------------------------------------------------- #

test_that("(i) summary() prints compression line on aggregated fit", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  out <- capture.output(summary(fit))
  expect_true(any(grepl("aggregation:.*->.*cells", out)))
  expect_true(any(grepl("ratio.*:1", out)))
})

test_that("(i) per-row summary does NOT print compression line", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = FALSE,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  out <- capture.output(summary(fit))
  expect_false(any(grepl("aggregation:.*->.*cells", out)))
})

test_that("(i) print() shows Exact. line on per-row and aggregated", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit_row <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = FALSE,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_agg <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  # Per-row INLA goes through print.flexybayes_inla which does NOT
  # currently surface exactness; aggregated goes through
  # print.flexybayes_aggregated which does.
  agg_out <- capture.output(print(fit_agg))
  expect_true(any(grepl("aggregated_exact", agg_out)))
  expect_true(any(grepl("aggregated-gaussian", agg_out)))
})


# ---------------------------------------------------------------- #
# Posterior-equivalence smoke (INLA path; matched explicit priors)  #
# ---------------------------------------------------------------- #
# Algebraic identity verified by the day-1 spike on the matched-prior
# reference. Here we confirm the per-row vs aggregated INLA fits agree
# on the fixed-effect posterior under the same effective prior path.

test_that("INLA per-row + aggregated agree on beta within 5e-3 under matched defaults", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data(seed = 7L)
  fit_row <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = FALSE,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_agg <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  beta_row <- fit_row$inla$summary.fixed$mean
  beta_agg <- fit_agg$extras$summary$beta_means
  # The fixed-effect posteriors are bit-exactly equivalent under the
  # algebraic identity (independent of the RI / residual prior tail).
  # v0.3.9 Phase B' amendment: the 1e-3 threshold was theoretical
  # (the day-1 spike at R/emit_gaussian_aggregated.R:30-41 reported
  # differences <= 1e-4 on a 1000-row example, and isolated runs of
  # this fixture deliver ~5e-4). In full-suite ordering, INLA's
  # internal numerical engine inherits BLAS reduction-order drift +
  # accumulated session state from prior fits, lifting the observed
  # delta to ~3.4e-3 on this fixture (documented in the v0.3.8 Phase
  # D triage that flagged this site, state_v038_released.md). The
  # 5e-3 band reflects the empirical worst-case at this fixture size
  # while keeping a 5x guard against gross algebraic regression. The
  # stress-gated counterpart below pins the tighter spike threshold
  # for the local-developer rehearsal probe. See ledger cairn
  # 2026-05-27-aggregated-inla-tolerance-amendment.cairn.md and the
  # v0.3.5 BLAS-non-associativity amendment on ADR 0023 for the
  # parallel precedent.
  expect_lt(max(abs(beta_row - beta_agg)), 5e-3)
})

test_that("intercept-only aggregated INLA fit summarises (1-coef beta_vcov guard)", {
  testthat::skip_if_not_installed("INLA")
  # Regression guard: diag(sd^2) on a length-one vector (intercept-only)
  # was read as diag(n) and returned a 0 x 0 matrix, crashing
  # .agg_inla_summarise() with a dimnames-extent mismatch.
  dat <- mk_agg_data(seed = 7L)
  fit <- suppressMessages(fb(
    y ~ 1 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = TRUE,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  bv <- fit$extras$summary$beta_vcov
  expect_identical(dim(bv), c(1L, 1L))
  expect_true(is.finite(bv[1, 1]))
})

test_that("INLA per-row + aggregated agree on beta within 1e-3 (stress, strict)", {
  testthat::skip_if_not_installed("INLA")
  if (!identical(Sys.getenv("FLEXYBAYES_RUN_STRESS"), "true")) {
    skip(
      "[stress] strict 1e-3 algebraic-identity band; opt-in via FLEXYBAYES_RUN_STRESS=true"
    )
  }
  dat <- mk_agg_data(seed = 7L)
  fit_row <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = FALSE,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_agg <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  beta_row <- fit_row$inla$summary.fixed$mean
  beta_agg <- fit_agg$extras$summary$beta_means
  # Spike-anchored algebraic-identity floor.
  expect_lt(max(abs(beta_row - beta_agg)), 1e-3)
})


# ---------------------------------------------------------------- #
# Aggregation_meta slot shape                                       #
# ---------------------------------------------------------------- #

test_that("aggregation_meta carries N, K, compression, residual kind", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  am <- fit$extras$aggregation_meta
  expect_true(is.list(am))
  expect_equal(am$N, nrow(dat))
  expect_true(is.integer(am$K) && am$K >= 1L)
  expect_equal(am$compression, am$K / am$N)
  expect_identical(am$residual, "homogeneous")
})


# ---------------------------------------------------------------- #
# Class + structure invariants                                      #
# ---------------------------------------------------------------- #

test_that("flexybayes_aggregated inherits from flexybayes", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
  expect_s3_class(fit, "flexybayes")
  # Backend identity in the class vector lets S3 generics that have
  # no .flexybayes_aggregated method (e.g., fb_as_draws_simple) fall
  # through to the correct backend-specific method. Without this,
  # triangulate(fit_g, fit_i) on aggregated INLA fits dispatched to
  # the greta extractor and failed.
  expect_s3_class(fit, "flexybayes_inla")
  expect_identical(
    class(fit),
    c("flexybayes_aggregated", "flexybayes_inla", "flexybayes")
  )
  expect_true(!is.null(fit$inla))
  expect_true(!is.null(fit$glm))
  expect_true(!is.null(fit$extras$summary))
})

test_that("flexybayes_aggregated on greta has no flexybayes_inla label", {
  # backend = "greta", aggregate = TRUE forces the aggregated path on
  # greta. (aggregate = "auto" defers to per-row on greta per ADR 0022,
  # which would give class "flexybayes" -- a different shape that
  # neither needs nor receives the inla label).
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "greta",
    aggregate = TRUE,
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
  expect_s3_class(fit, "flexybayes")
  expect_false(inherits(fit, "flexybayes_inla"))
  expect_identical(
    class(fit),
    c("flexybayes_aggregated", "flexybayes")
  )
})

test_that("aggregation_meta carries prior_parametrization, surfaced by canonical_names", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  # Default prior -> per-row-equivalent (the matched-prior guarantee).
  fit_def <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_identical(
    fit_def$extras$aggregation_meta$prior_parametrization,
    "per_row_equivalent"
  )
  expect_identical(
    canonical_names(fit_def)$prior_parametrization,
    "per_row_equivalent"
  )
  # Explicit prior -> custom (equivalence against a default-prior per-row
  # fit no longer holds; flagged so the user does not misread agreement).
  fit_cus <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior = fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 5)),
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_identical(
    fit_cus$extras$aggregation_meta$prior_parametrization,
    "custom"
  )
  expect_identical(
    canonical_names(fit_cus)$prior_parametrization,
    "custom"
  )
  # The label reaches the user-facing print surface.
  out <- utils::capture.output(print(fit_def))
  expect_true(any(grepl("per-row-equivalent", out)))
})

test_that("greta-aggregated fit exposes $greta$draws for triangulate / extractors", {
  # Regression guard for the v0.4.0 fix: the greta-aggregated path
  # stored the coda draws bare at fit$greta, leaving fit$greta$draws
  # NULL, so fb_as_draws_simple() / canonical_names() / confint() /
  # plot() / predict() -- and hence triangulate() -- all failed on
  # greta-aggregated fits with "fit$greta$draws is missing". Surfaced
  # by adding the live aggregated-path triangulate test; the fix wraps
  # the draws as list(draws = ...) to mirror the per-row emit_greta()
  # slot shape.
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  dat <- mk_agg_data(seed = 3L, N = 200L, J = 6L)
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "greta",
    aggregate = TRUE,
    n_samples = 60L,
    warmup = 60L,
    chains = 1L,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_false(is.null(fit$greta$draws))
  draws <- flexyBayes:::fb_as_draws_simple(fit)
  expect_true(is.list(draws) && length(draws) >= 1L)
  expect_true(all(vapply(draws, is.numeric, logical(1))))
  # Self-triangulation is exact: every parameter common, zero drift.
  tri <- triangulate(fit, fit)
  expect_identical(tri$n_common, length(draws))
  expect_equal(max(abs(tri$metrics$mean_diff)), 0)
})

test_that("$glm shim carries per-row reconstructed fitted values", {
  testthat::skip_if_not_installed("INLA")
  dat <- mk_agg_data()
  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_equal(length(fit$glm$fitted.values), nrow(dat))
  expect_equal(length(fit$glm$y), nrow(dat))
  expect_equal(length(fit$glm$residuals), nrow(dat))
})


# ---------------------------------------------------------------- #
# Greta posterior-equivalence (v0.3.3 activation)                    #
# ---------------------------------------------------------------- #
# Aggregated greta backend-fidelity check: on a fixed-seed dataset
# generated from a known (beta, sigma, tau) triple, the aggregated
# greta posterior means recover the truth within MC tolerance. The
# bit-exact algebraic identity between per-row and aggregated forms
# is already tested in `test-aggregate-gaussian.R` (pure-R, no
# backend). This test is the backend-fidelity check: did greta
# correctly translate the cell-mean weighted gaussian + gamma-on-WSS
# pattern into a TF graph + MCMC that recovers truth? Spike gates
# from `spikes_stage3a/spike_greta.R`: beta <= 0.10, sigma <= 0.10,
# tau <= 0.15 (tau gets a looser bound because the J=20 RI factor
# has only 20 levels, so its posterior is wider than beta/sigma's).
#
# Stress-gated by FLEXYBAYES_RUN_STRESS=true (one greta MCMC fit ~30-
# 60 sec on the dev laptop, TF graph compile included).
#
# A direct per-row vs aggregated cross-path comparison is impractical
# here because the per-row greta path uses an over-parameterised
# fixed-effect contrast scheme (4 columns for a 3-level factor: the
# intercept + all 3 levels, vs treatment contrasts in the aggregated
# path producing 3 columns). The two posteriors live on different
# parametrisations. The pure-R log-likelihood equivalence test in
# the foundation file covers the algebraic identity; this test
# covers the backend-side translation only.

test_that("(d) greta aggregated posterior recovers truth (MC tol)", {
  skip_if_greta_backend_unusable()
  if (!identical(Sys.getenv("FLEXYBAYES_RUN_STRESS"), "true")) {
    testthat::skip(
      "FLEXYBAYES_RUN_STRESS != \"true\"; skip greta backend-fidelity check."
    )
  }

  # Generation matches spike_greta.R: known beta on treatment-contrast
  # columns, J=20 random-intercept levels, residual sigma_true.
  set.seed(2026L)
  N <- 800L
  J <- 20L
  beta_true <- c(2.0, -1.0, 0.5, 0.8, -0.3)
  sigma_true <- 0.7
  tau_true <- 0.4
  dat <- data.frame(
    f1 = factor(sample(letters[1:3], N, replace = TRUE)),
    f2 = factor(sample(letters[1:3], N, replace = TRUE)),
    g = factor(sample.int(J, N, replace = TRUE))
  )
  u_true <- rnorm(J, 0, tau_true)
  Xmm <- model.matrix(~ f1 + f2, data = dat)
  mu <- as.numeric(Xmm %*% beta_true) + u_true[as.integer(dat$g)]
  dat$y <- rnorm(N, mu, sigma_true)

  fit <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = dat,
    backend = "greta",
    aggregate = TRUE,
    n_samples = 1000L,
    warmup = 500L,
    chains = 2L,
    prior_fixed_sd = 10,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  expect_identical(fit$exactness, "aggregated_exact")
  expect_identical(backend_decision(fit)$backend, "greta")
  expect_identical(backend_decision(fit)$path, "aggregated_gaussian")

  beta_post <- summary(fit)$beta_means
  sigma_post <- as.numeric(summary(fit)$sigma_means)
  tau_post <- as.numeric(summary(fit)$tau_means)

  expect_true(
    all(abs(beta_post - beta_true) <= 0.10),
    info = sprintf(
      "max |beta_post - beta_true| = %.4f (truth: %s)",
      max(abs(beta_post - beta_true)),
      paste(sprintf("%.3f", beta_true), collapse = ", ")
    )
  )
  expect_true(
    abs(sigma_post - sigma_true) <= 0.10,
    info = sprintf(
      "sigma_post = %.4f, sigma_true = %.4f",
      sigma_post,
      sigma_true
    )
  )
  expect_true(
    abs(tau_post - tau_true) <= 0.15,
    info = sprintf("tau_post = %.4f, tau_true = %.4f", tau_post, tau_true)
  )
})
