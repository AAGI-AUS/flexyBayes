# test-emit-state-isolation.R -- Regression suite for the v0.3.9
# migration of the emit-once message latches out of the options()
# namespace into a package-internal env (R/emit_state.R).
#
# The pre-v0.3.9 latch was keyed on
# options(flexyBayes._default_prior_note_emitted) and
# options(flexyBayes._uniform_inla_approx_emitted); any caller that
# touched those names could consume the flag before the test
# intended to re-emit, producing suite-order flakes. These tests pin
# the new contract: the live store is .flexybayes_emit_state, the
# public-API reset is local_clean_emit_state(), and stray options()
# writes have no effect on the latch decision.

# Local fixtures (parallel testthat: each test file runs in its own
# R process, so cross-file helpers are not in scope).
mk_emit_data <- function() {
  set.seed(42)
  n <- 30L
  data.frame(
    yield = rnorm(n),
    env = factor(rep(1:3, length.out = n)),
    geno = factor(rep(1:5, length.out = n))
  )
}

test_that(".flexybayes_emit_state initialises to FALSE/FALSE on load", {
  # Reset then read; both keys must be FALSE after a fresh init.
  flexyBayes:::.reset_emit_state_for_test()
  expect_false(flexyBayes:::.emit_state_get("default_prior_note"))
  expect_false(flexyBayes:::.emit_state_get("uniform_inla_approx"))
})

test_that(".emit_state_get returns FALSE for unknown keys", {
  flexyBayes:::.reset_emit_state_for_test()
  expect_false(flexyBayes:::.emit_state_get("no_such_key"))
})

test_that(".emit_state_set toggles the once-flag", {
  flexyBayes:::.reset_emit_state_for_test()
  flexyBayes:::.emit_state_set("default_prior_note", TRUE)
  expect_true(flexyBayes:::.emit_state_get("default_prior_note"))
  flexyBayes:::.emit_state_set("default_prior_note", FALSE)
  expect_false(flexyBayes:::.emit_state_get("default_prior_note"))
})

test_that("options() writes do not influence the latch decision", {
  # Set the legacy option-keyed flag to TRUE; the new env-keyed latch
  # must still report FALSE since the option no longer drives the
  # decision in R/fb_prior.R.
  flexyBayes:::.reset_emit_state_for_test()
  withr::local_options(
    flexyBayes._default_prior_note_emitted = TRUE,
    flexyBayes._uniform_inla_approx_emitted = TRUE
  )
  expect_false(flexyBayes:::.emit_state_get("default_prior_note"))
  expect_false(flexyBayes:::.emit_state_get("uniform_inla_approx"))
})

test_that("local_clean_emit_state() resets on entry + restores on exit", {
  # Pre-seed a known state, then exercise the scope wrapper.
  flexyBayes:::.emit_state_set("default_prior_note", TRUE)
  flexyBayes:::.emit_state_set("uniform_inla_approx", TRUE)
  local({
    local_clean_emit_state()
    # Inside the scope: latch is fresh.
    expect_false(flexyBayes:::.emit_state_get("default_prior_note"))
    expect_false(flexyBayes:::.emit_state_get("uniform_inla_approx"))
  })
  # Outside the scope: prior state restored.
  expect_true(flexyBayes:::.emit_state_get("default_prior_note"))
  expect_true(flexyBayes:::.emit_state_get("uniform_inla_approx"))
  # Tidy up.
  flexyBayes:::.reset_emit_state_for_test()
})

test_that("default-prior-note re-emits in a back-to-back call pattern", {
  # The v0.3.8 flake reproduced when an earlier caller (e.g.
  # fb_plan(plan = TRUE)) flipped the option-keyed flag and the test's
  # reset was no-op'd by an intervening evaluation. Under the v0.3.9
  # env-keyed latch the reset is authoritative.
  skip_if_no_greta()
  d <- mk_emit_data()

  # First call: emit, set the latch.
  local({
    local_clean_emit_state()
    expect_message(
      flexybayes(
        yield ~ env,
        random = ~geno,
        data = d,
        verbose = FALSE,
        return_code = TRUE
      ),
      regexp = "uniform\\(0,"
    )
  })

  # Second call after a deliberate option-namespace poke: still emits
  # because the env-keyed latch reset is invisible to options().
  withr::local_options(
    flexyBayes._default_prior_note_emitted = TRUE
  )
  local({
    local_clean_emit_state()
    expect_message(
      flexybayes(
        yield ~ env,
        random = ~geno,
        data = d,
        verbose = FALSE,
        return_code = TRUE
      ),
      regexp = "uniform\\(0,"
    )
  })
})

test_that("silence option short-circuits the latch", {
  # The public-API silence toggle remains options()-keyed and must
  # still suppress the message even when the env-latch is FALSE.
  skip_if_no_greta()
  d <- mk_emit_data()
  local_clean_emit_state()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  expect_no_message(
    flexybayes(
      yield ~ env,
      random = ~geno,
      data = d,
      verbose = FALSE,
      return_code = TRUE
    )
  )
  # Latch was not flipped because the silence path returns early.
  expect_false(flexyBayes:::.emit_state_get("default_prior_note"))
})
