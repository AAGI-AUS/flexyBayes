# Tests for prior_summary() -- user-facing accessor for the
# resolved prior on a flexyBayes fit.
#
# Contract:
#   - S3 generic + methods for flexybayes / flexybayes_inla, plus a
#     default that refuses with a structured message.
#   - Returns a prior_summary_flexybayes object carrying `kind`
#     (one of "fb_prior" / "legacy_scalar" / "no_prior_recorded" /
#     "unknown_shape"), `backend`, and (when applicable) the
#     fb_prior + auto-default origin metadata.
#
# 0.9.3: this file used to carry four cases gated behind a skip on the
# native engine withdrawn entirely in 0.9.3 (see NEWS.md) -- the
# auto-default / user-origin / legacy-scalar-bridge tracking, and a
# declaration-only flag specific to that engine's direct-model entry
# point. All four routes are gone along with the engine (fitting a
# `backend` naming it now refuses with `unknown_backend` before a fit
# object exists, so a fit to call prior_summary() on can no longer be
# built). The live contract on both active engines is covered without
# a skip in test-orphaned-exports-live.R. This file keeps only the
# one case that was never engine-specific.

test_that("prior_summary() default method refuses unknown classes", {
  expect_error(
    prior_summary(lm(mpg ~ wt, data = mtcars)),
    "prior_summary"
  )
})


# =============================================================================
# The prior column on a heterogeneous route (ASReml oracle suite, finding O01a).
#
# diag() / idh() / at() emit one iid field per level of the outer factor, and
# each carries the hyper block keyed by the inner group. The canonical names
# are `sd_<group>__<outer><level>`, which never equals the prior spec's key
# `sd_<group>`, so every per-level row fell through to the closing
# "engine default" branch while the emitted INLA formula carried the prior
# verbatim.
#
# The cost of that was not cosmetic. Reading the column was how the oracle
# suite's first pass concluded that the prior layer did not reach diag() at
# all -- a false finding, written up before the emitted formula was read. A
# column that misreports the fit is what the next reader reasons from.
# =============================================================================

.po_data <- function() {
  set.seed(20260911L)
  s1 <- 96L
  d <- data.frame(
    env = factor(rep(paste0("E", 1:4), each = 24L)),
    gen = factor(rep(paste0("g", 1:12), times = 8L)),
    rep = factor(rep(rep(c("r1", "r2"), each = 12L), times = 4L))
  )
  d$yield <- 10 + rnorm(4L, 0, 1)[as.integer(d$env)] +
    rnorm(12L, 0, 1.5)[as.integer(d$gen)] + rnorm(s1, 0, 1)
  d
}

test_that("a user prior on a diag() term is reported on every per-level row", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .po_data()
  fit <- suppressWarnings(flexybayes(
    fixed = yield ~ env, random = ~ diag(env):gen, data = d,
    backend = "inla", prior = fb_prior(sd(gen) ~ uniform(lower = 0, upper = 40)),
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  rows <- v[startsWith(v$component, "sd_gen__"), ]
  expect_gt(nrow(rows), 1L)
  expect_true(all(rows$prior == "uniform(lower=0, upper=40)"))
  expect_false(any(rows$prior == "engine default"))
})

test_that("the package default is reported on every per-level row too", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .po_data()
  fit <- suppressWarnings(flexybayes(
    fixed = yield ~ env, random = ~ diag(env):gen, data = d,
    backend = "inla", verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  rows <- v[startsWith(v$component, "sd_gen__"), ]
  expect_true(all(startsWith(rows$prior, "uniform(lower=0")))
  # The cell must agree with the residual row, which shares the default.
  expect_equal(unique(rows$prior), v$prior[v$component == "sigma"])
})

test_that("a component with no slot for the prior still reports the engine default", {
  # The falsifier. The fix widens a lookup, and a widened lookup that matched
  # everything would pass both tests above while destroying the column's
  # meaning. `not_applied` is tested ahead of the widened branch, so a
  # structure the prior cannot reach still says so.
  expect_true(is.function(flexyBayes:::.fb_varcomp_priors))
  body_txt <- paste(deparse(flexyBayes:::.fb_varcomp_priors), collapse = " ")
  expect_match(body_txt, "not_applied", fixed = TRUE)
  expect_lt(
    regexpr("not_applied", body_txt, fixed = TRUE),
    regexpr("base_key", body_txt, fixed = TRUE)
  )
})
