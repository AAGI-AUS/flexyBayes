# emit_state.R -- Package-internal store for emit-once message latches.
#
# Holds the boolean flags that gate one-time announcements such as the
# default variance-component prior note and the uniform-on-INLA approx
# note. The store sits in a package-private environment so it cannot be
# read or mutated via the public options() namespace; test code reaches
# it through the .reset_emit_state_for_test() shim and the
# helper-emit-state.R local_clean_emit_state() wrapper. The public
# silence-mechanism remains options(flexyBayes.silence_*); only the
# internal once-flag store is encapsulated here.

.flexybayes_emit_state <- new.env(parent = emptyenv())

# Initialise the store on package load. Called from .onLoad() in zzz.R.
.init_emit_state <- function() {
  .flexybayes_emit_state$default_prior_note <- FALSE
  .flexybayes_emit_state$uniform_inla_approx <- FALSE
  # With backend = "auto" now the default, the
  # auto-fallback notes (lgm_gate refused -> greta; INLA not installed ->
  # greta) fire once per session rather than on every default call.
  .flexybayes_emit_state$auto_fallback_note <- FALSE
  .flexybayes_emit_state$auto_inla_missing_note <- FALSE
  # One-time note when the auto path's INLA fit
  # fails numerically and falls back to greta.
  .flexybayes_emit_state$auto_inla_numerical_fallback_note <- FALSE
  invisible(NULL)
}

# Read a once-flag. Returns FALSE for any key that has not been
# initialised, so callers can rely on a boolean without an existence
# check.
.emit_state_get <- function(key) {
  if (!exists(key, envir = .flexybayes_emit_state, inherits = FALSE)) {
    return(FALSE)
  }
  isTRUE(.flexybayes_emit_state[[key]])
}

# Set a once-flag.
.emit_state_set <- function(key, value = TRUE) {
  .flexybayes_emit_state[[key]] <- isTRUE(value)
  invisible(NULL)
}

# Reset the store to the freshly-loaded shape. Intended for test code
# (via flexyBayes:::.reset_emit_state_for_test() from
# helper-emit-state.R); not for user-facing call sites, which rely on
# the once-per-session semantics or the public
# options(flexyBayes.silence_*) toggles.
.reset_emit_state_for_test <- function() {
  .init_emit_state()
}
