# helper-rng.R -- Per-test isolation of TensorFlow's internal RNG
# state for greta-backed assertions.
#
# greta builds models on top of TensorFlow, which maintains its own
# session-internal random state that R's set.seed() does not reach.
# Any greta::mcmc() call earlier in the test session shifts that
# state, so a downstream test's posterior depends on the cumulative
# call history even with a fixed R seed. local_tf_seed() narrows the
# leak surface by reseeding both R and TF at scope entry. Perfect
# determinism additionally requires
# tensorflow::tf$config$experimental$enable_op_determinism() which
# carries a ~5-15% throughput penalty and changes process-global
# state; the helper does not enable it by default (opt-in via
# deterministic_ops = TRUE for tests that genuinely need bit-exact
# repeatability).
#
# Discipline rule for new tests:
#   Any test that asserts numerical agreement with a deterministic
#   reference (lme4, INLA point estimate, BLUP, REML) on a
#   greta-MCMC posterior must call local_tf_seed() at test entry.
#   Bare set.seed() alone is insufficient -- R's RNG and TF's RNG
#   are independent streams. The discipline pairs with the
#   bare-options() rule documented in helper-emit-state.R: tests
#   leak across files in tally.R's single-process loop, so all
#   reset-style side effects must be scoped to the test.

local_tf_seed <- function(
  seed,
  .local_envir = parent.frame(),
  deterministic_ops = FALSE
) {
  set.seed(seed)
  if (requireNamespace("tensorflow", quietly = TRUE)) {
    tf <- tensorflow::tf
    tryCatch(tf$random$set_seed(as.integer(seed)), error = function(e) {
      invisible(NULL)
    })
    if (isTRUE(deterministic_ops)) {
      tryCatch(
        tf$config$experimental$enable_op_determinism(),
        error = function(e) invisible(NULL)
      )
    }
  }
  invisible(NULL)
}
