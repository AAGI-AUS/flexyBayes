# Process-isolation refusals: registry population for
# R/process_isolation.R.
#
# Codes register as their throwing sites land, never in advance. Both codes
# below are thrown by R/process_isolation.R and triggered in
# tests/testthat/test-process-isolation.R: the first with callr's absence
# simulated, the second by a signal sent to the child's own process ID while
# it fits.
#
# Why a death needs a name of its own. A condition raised inside the child (a
# refusal, an engine error, an untyped failure) travels back and is
# re-raised as itself, so every existing code keeps its name under
# isolation. A child that the operating system kills, or that crashes in
# compiled code, sends nothing back: all that is left is its exit status or
# the signal that ended it. No existing code names a fit whose process ended
# without returning.

.populate_refusal_registry_process_isolation <- function() {
  .register_refusal(
    reason_code = "isolation_requires_callr",
    description = paste0(
      "options(flexyBayes.isolation = \"process\") runs each engine fit, ",
      "fb_isolated_lapply() each job, and option flexyBayes.inla_best_of ",
      "above 1 the further INLA fits, in a separate R process started ",
      "through callr, which is in Suggests and not installed. The request ",
      "is refused before anything is fitted, so no fit runs without the ",
      "protection that was asked for."
    ),
    message_template = paste0(
      "%s: options(flexyBayes.isolation = \"process\") needs the callr ",
      "package, which is not installed."
    ),
    registered_in_adr = "0.14.0 process isolation",
    plan_field = NA_character_,
    since_version = "0.14.0"
  )

  .register_refusal(
    reason_code = "isolated_process_died",
    description = paste0(
      "A separate R process running an engine fit or a job of ",
      "fb_isolated_lapply() ",
      "ended before it returned: a signal ended it (SIGKILL is sent from ",
      "outside the process, SIGSEGV or SIGABRT comes from a crash in ",
      "compiled code), or it exited with a non-zero status. The refusal ",
      "names the engine, the entry point, the fold where there is one, the ",
      "exit status and the signal. The session that asked for the fit is ",
      "unaffected, and nothing from the process is returned."
    ),
    message_template = paste0(
      "%s: the separate R process fitting with %s ended before it ",
      "returned: %s."
    ),
    registered_in_adr = "0.14.0 process isolation",
    plan_field = NA_character_,
    since_version = "0.14.0"
  )

  invisible(NULL)
}
