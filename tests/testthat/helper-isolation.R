# helper-isolation.R -- Umbrella test-isolation contract.
#
# The v0.3.9 test-isolation pass added two helpers covering two
# distinct leak classes that turned out to share a root cause
# (state held outside the test_that frame, surviving across files
# in tools/tally.R's single-process loop):
#
#   * helper-emit-state.R -- the once-per-session message latches
#     (default-prior note, uniform-on-INLA approx note) sit in a
#     package-private env (R/emit_state.R). local_clean_emit_state()
#     scopes a snapshot + reset + restore around a test_that body
#     so a downstream test's expected emission is not consumed by
#     an upstream caller. The package's public-API silence toggles
#     (options(flexyBayes.silence_*)) are unchanged; only the
#     internal once-flag store is encapsulated.
#
#   * helper-rng.R -- TensorFlow's session-internal RNG state is
#     not reachable from R's set.seed(). Any greta::mcmc() call
#     earlier in the session shifts TF state; a downstream
#     posterior depends on the cumulative call history.
#     local_tf_seed(seed) reseeds both R and TF at scope entry;
#     deterministic_ops = TRUE additionally enables
#     tf$config$experimental$enable_op_determinism() (opt-in;
#     not on by default because of process-global TF state +
#     throughput cost).
#
# local_flexybayes_clean_state() is the recommended one-call
# umbrella for tests that exercise both surfaces (a greta-MCMC
# fit that also triggers the default-prior emission). It is a
# thin composer; individual helpers remain the right entry
# point for tests that only need one surface.
#
# Cross-helper discipline rules (apply to any test in this
# directory):
#
#   * Never call bare `options(flexyBayes.* = ...)` -- always
#     wrap in withr::local_options(...) so the option is
#     restored on scope exit. tally.R sources files in a single
#     R process; un-restored options() leak downstream.
#   * Greta-MCMC tests that assert numerical agreement with a
#     deterministic reference must call local_tf_seed() at
#     entry. set.seed() alone is insufficient -- R and TF have
#     independent RNG streams.
#   * Tests that depend on a deterministic emit-state must call
#     local_clean_emit_state() at entry. The package's once-flag
#     latch persists across test_that blocks in the same session.
#
# See:
#   ledger/r-package/2026-05-27-emit-state-latch-env-migration.cairn.md
#   ledger/r-package/2026-05-27-tf-rng-containment-hybrid-sleepstudy.cairn.md
#   feedback_bare_options_in_r_tests.md (cross-project memory)
#
# The umbrella ties the two surfaces together; future helpers
# (e.g., INLA-side session-state isolation if a process-isolated
# harness lands) compose by adding another argument to this
# function.

local_flexybayes_clean_state <- function(
  seed = NULL,
  .local_envir = parent.frame(),
  deterministic_ops = FALSE
) {
  local_clean_emit_state(.local_envir = .local_envir)
  if (!is.null(seed)) {
    local_tf_seed(
      seed,
      .local_envir = .local_envir,
      deterministic_ops = deterministic_ops
    )
  }
  invisible(NULL)
}
