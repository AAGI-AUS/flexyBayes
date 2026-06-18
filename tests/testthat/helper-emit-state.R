# helper-emit-state.R -- Per-test isolation of the flexyBayes
# emit-once message latches.
#
# The package keeps "emit-once" announcements (default-prior note,
# uniform-on-INLA approx note) in a package-internal environment
# exposed only to internal accessors. Tests that depend on a
# deterministic emission state must reset that env before the call
# that should re-emit; the helpers here are the supported entry
# points. local_clean_emit_state() snapshots the env, resets it to
# the fresh shape, and restores the snapshot on scope exit so the
# isolation does not leak into sibling tests in the same file.
#
# Discipline rule for new tests:
#   Never call bare `options(flexyBayes.* = ...)` -- always wrap in
#   withr::local_options(...) so the option is restored on scope
#   exit. The tally.R full-suite runner sources tests in a single
#   R process, so an un-restored options() in one file leaks into
#   downstream files. Bare options(flexyBayes.silence_*) was the
#   upstream cause of the v0.3.8 suite-order flake at
#   test-smooth.R:88 (six leak sites repaired during v0.3.9).
#
# See R/emit_state.R for the live store and
# ledger/r-package/2026-05-27-emit-state-latch-env-migration.cairn.md
# for the rationale behind the v0.3.9 migration out of options().

reset_emit_state <- function() {
  flexyBayes:::.reset_emit_state_for_test()
}

local_clean_emit_state <- function(.local_envir = parent.frame()) {
  snapshot <- as.list(flexyBayes:::.flexybayes_emit_state)
  flexyBayes:::.reset_emit_state_for_test()
  withr::defer(
    {
      flexyBayes:::.reset_emit_state_for_test()
      for (k in names(snapshot)) {
        flexyBayes:::.emit_state_set(k, snapshot[[k]])
      }
    },
    envir = .local_envir
  )
  invisible(NULL)
}

set_emit_state <- function(key, value = TRUE) {
  flexyBayes:::.emit_state_set(key, value)
}
