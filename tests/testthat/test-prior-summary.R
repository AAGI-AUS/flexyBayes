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
