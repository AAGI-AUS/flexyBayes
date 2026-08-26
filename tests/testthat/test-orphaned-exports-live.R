# test-orphaned-exports-live.R -- the three exports with no live coverage.
#
# `prior_summary()`, `fb_met_summary()` and `fb_log_posterior()` were the
# package's only exports whose behavioural tests were all gated by a skip
# on the engine withdrawn entirely in 0.9.3 (see NEWS.md), which skipped
# unconditionally once that engine was quarantined -- or which pinned a
# stub's class rather than a fit's behaviour. Each is exercised here on a
# live fit of an active engine.
#
# prior_summary() also needed a correctness fix before a test was worth
# writing: on a fit with a distributional residual it reported
# `sigma ~ uniform(0, U)` for a model with no scalar sigma, and it omitted
# every term the default-prior walker does not reach. Both are asserted
# below.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

mk_orphan_data <- function(n_group = 12L, n_rep = 5L) {
  set.seed(20260815L)
  g <- factor(rep(seq_len(n_group), each = n_rep))
  x <- stats::rnorm(n_group * n_rep)
  u <- stats::rnorm(n_group, sd = 0.8)
  y <- 1.2 + 0.6 * x + u[as.integer(g)] +
    stats::rnorm(n_group * n_rep, sd = 0.5)
  data.frame(y = y, x = x, g = g)
}

mk_orphan_met_data <- function() {
  set.seed(20260816L)
  n_env <- 3L
  n_gen <- 8L
  n_rep <- 3L
  d <- expand.grid(
    rep = factor(seq_len(n_rep)),
    gen = factor(seq_len(n_gen)),
    env = factor(paste0("E", seq_len(n_env)))
  )
  env_mean <- c(E1 = 4, E2 = 6, E3 = 5)
  g_eff <- stats::rnorm(n_gen, sd = 1)
  ge_eff <- stats::rnorm(n_gen * n_env, sd = 0.6)
  sd_by_env <- c(E1 = 0.4, E2 = 0.9, E3 = 1.6)
  d$yield <- env_mean[as.character(d$env)] +
    g_eff[as.integer(d$gen)] +
    ge_eff[(as.integer(d$env) - 1L) * n_gen + as.integer(d$gen)] +
    stats::rnorm(nrow(d), sd = sd_by_env[as.character(d$env)])
  d
}

# ---------------------------------------------------------------- #
# prior_summary() on both active engines                            #
# ---------------------------------------------------------------- #

test_that("prior_summary() reports the auto-default on a live INLA fit", {
  skip_if_no_inla()
  d <- mk_orphan_data()

  fit <- suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
  ps <- prior_summary(fit)

  expect_s3_class(ps, "prior_summary_flexybayes")
  expect_identical(ps$kind, "fb_prior")
  expect_identical(ps$backend, "inla")
  expect_identical(ps$default_origin, "auto")
  expect_identical(ps$default_basis, "identity_sd_uniform")
  expect_true(is.numeric(ps$default_scale) && ps$default_scale > 0)

  txt <- utils::capture.output(print(ps))
  expect_true(any(grepl("backend = inla", txt, fixed = TRUE)))
  expect_true(any(grepl("auto-default bounded uniform", txt, fixed = TRUE)))
  # A simple random-intercept model has nothing outside the walker, so
  # there is no engine-default block to print.
  expect_length(ps$engine_default, 0L)
  expect_null(ps$residual_lowered_to)
})

test_that("prior_summary() reports what reached Stan on a live brms fit", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_orphan_data()

  fit <- suppressWarnings(suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "brms",
    n_samples = 300L,
    warmup = 300L,
    chains = 2L,
    seed = 20260815L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  ps <- prior_summary(fit)

  expect_s3_class(ps, "prior_summary_flexybayes")
  expect_identical(ps$backend, "brms")
  expect_identical(ps$default_origin, "auto")

  # brms is the authority on what reached Stan, so its own table is
  # carried rather than reconstructed from the declaration.
  expect_s3_class(ps$engine_prior_table, "brmsprior")
  txt <- utils::capture.output(print(ps))
  expect_true(any(grepl("As it reached Stan", txt, fixed = TRUE)))
})

test_that("prior_summary() names the retarget on a distributional fit", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_orphan_met_data()

  fit <- suppressWarnings(suppressMessages(flexybayes(
    fixed = yield ~ env,
    random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "brms",
    n_samples = 400L,
    warmup = 400L,
    chains = 2L,
    seed = 20260816L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  ps <- prior_summary(fit)
  txt <- utils::capture.output(print(ps))

  # The defect this test exists for: the model has no scalar sigma, and
  # the declared uniform on the SD scale was retargeted onto the
  # per-level log-sigma coefficients. Printing the declaration alone
  # named a parameter the model does not contain.
  expect_identical(ps$residual_lowered_to, "env")
  expect_true(any(grepl("no scalar `sigma`", txt, fixed = TRUE)))
  expect_true(any(grepl("retargeted", txt, fixed = TRUE)))

  # The Stan-side table is the check on the claim: there is a b_sigma
  # class and no scalar sigma row with a uniform on it.
  tab <- ps$engine_prior_table
  expect_true(any(tab$dpar == "sigma"))
  expect_false(any(tab$class == "sigma" & nzchar(tab$prior)))
})

test_that("prior_summary() lists parameters left to the engine's default", {
  # A term outside the default-prior walker must appear by name with its
  # reason, not vanish from the summary. The record is the same one
  # triangulate()'s matched-prior gate reads.
  d <- mk_orphan_met_data()
  d$t <- as.integer(d$rep)
  fb <- flexyBayes:::.build_ir_polymorphic(
    fixed = yield ~ env,
    random = ~ gen + ar1(t),
    residual = NULL,
    data = d,
    family = "gaussian",
    link = NULL,
    weights = NULL,
    known_matrices = list(),
    prior = NULL,
    prior_fixed_sd = 100,
    prior_vc_sd = 1,
    syntax = "asreml"
  )
  rec <- flexyBayes:::.fb_prior_record(fb)

  # The AR1 field sits outside the default-prior walker, so it must be
  # named as engine-default rather than dropped from the record.
  expect_true(length(rec$engine_default) > 0L)
  expect_true(all(nzchar(rec$engine_default)))
  expect_false(any(names(rec$engine_default) %in% names(rec$recorded)))
})

# ---------------------------------------------------------------- #
# fb_met_summary() -- typed abstention on both active engines       #
# ---------------------------------------------------------------- #

test_that("fb_met_summary() abstains with a typed condition on INLA", {
  skip_if_no_inla()
  d <- mk_orphan_data()

  fit <- suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))

  err <- tryCatch(fb_met_summary(fit), condition = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_met_summary_not_available")
  expect_s3_class(err, "flexybayes_refusal")
  # 0.9.3: the refusal is unconditional (no per-engine `$backend` field or
  # wording any more -- the reason is identical on every active engine, so
  # it is stated once rather than repeated per engine). Confirmed here on
  # an actual INLA fit, not just a hand-built stub.
  expect_null(err$backend)
  expect_match(conditionMessage(err), "INLA")
  expect_match(conditionMessage(err), "brms")
  # The pointer names something that works for an active engine. It used
  # to name fb_structured_cov(), which abstains for every structure an
  # active engine can fit.
  expect_match(conditionMessage(err), "VarCorr")
  expect_false(grepl("fb_structured_cov", conditionMessage(err), fixed = TRUE))
})

test_that("fb_met_summary() abstains with a typed condition on brms", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_orphan_data()

  fit <- suppressWarnings(suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "brms",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    seed = 20260817L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))

  err <- tryCatch(fb_met_summary(fit), condition = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_met_summary_not_available")
  # Same unconditional refusal as the INLA fit above -- confirms an
  # actual brms fit reaches the identical typed condition, not a
  # brms-specific variant.
  expect_null(err$backend)
  expect_match(conditionMessage(err), "brms")
})

# ---------------------------------------------------------------- #
# fb_log_posterior() -- typed abstention, on fits not on stubs      #
# ---------------------------------------------------------------- #

test_that("fb_log_posterior() abstains typed on a live INLA fit", {
  skip_if_no_inla()
  d <- mk_orphan_data()

  fit <- suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))

  err <- tryCatch(fb_log_posterior(fit), condition = function(e) e)
  expect_s3_class(err, "fb_c4_unavailable")
  expect_match(conditionMessage(err), "INLA backend")
})

test_that("fb_log_posterior() abstains typed on a live brms fit", {
  skip_on_cran()
  skip_if_not_installed("brms")
  d <- mk_orphan_data()

  fit <- suppressWarnings(suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "brms",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    seed = 20260818L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))

  err <- tryCatch(fb_log_posterior(fit), condition = function(e) e)
  expect_s3_class(err, "fb_c4_unavailable")
  expect_match(conditionMessage(err), "brms backend")
  expect_false(grepl("honest", conditionMessage(err), fixed = TRUE))
})
