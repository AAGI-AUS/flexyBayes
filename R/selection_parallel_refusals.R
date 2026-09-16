# Fold-concurrency refusals --- registry population for the worker loop.
#
# Codes register as their throwing sites land, never in advance. A registered
# code with no triggering test is a defect: a refusal nobody has seen fire is
# not a refusal, it is a string. The one code below is thrown by
# R/selection_parallel.R and triggered in
# tests/testthat/test-selection-parallel-envelope.R, where a fold function
# kills its own worker.
#
# Why the failure needs a name of its own. A fold that refuses, or errors,
# inside a live worker sends its condition back and it is re-raised as
# itself, so every existing refusal keeps its own name on the concurrent
# path. A worker the operating system kills -- reclaiming memory, or losing
# the process with an engine crash -- sends nothing back at all: measured
# here, `parallel::mclapply()` returns a bare NULL in that element, with no
# condition attached. There is no existing code for a fold that neither
# completed nor refused, and aggregating the folds that did come back would
# report a cross-validation over a subset of its folds under the name of the
# whole.

.populate_refusal_registry_selection_parallel <- function() {
  .register_refusal(
    reason_code = "selection_parallel_worker_failed",
    description = paste0(
      "fb_select(): a worker running one fold did not return a result and ",
      "left no error to report, which is what a child killed from outside ",
      "looks like -- the operating system reclaiming memory under the ",
      "aggregate demand of concurrent fold refits, or an engine crash ",
      "taking the process with it. The scores are not aggregated from the ",
      "folds that survived: a cross-validation missing a fold is not the ",
      "cross-validation that was asked for. A fold that refuses inside a ",
      "live worker is unaffected, and is re-raised under its own reason ",
      "code."
    ),
    message_template = paste0(
      "fb_select(): the worker running fold %s did not return a result and ",
      "left no error to report."
    ),
    registered_in_adr = "PLAN_green_0.12.1_to_0.15.0_2026-09-16 section 6 (P4)",
    plan_field = "selection_plan",
    since_version = "0.13.0"
  )

  invisible(NULL)
}
