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
# path. A forked worker that ends without returning sends nothing back:
# `parallel::mclapply()` returns a bare NULL in that element, with no
# condition and no exit status attached. There is no existing code for a fold that neither
# completed nor refused, and aggregating the folds that did come back would
# report a cross-validation over a subset of its folds under the name of the
# whole.

.populate_refusal_registry_selection_parallel <- function() {
  .register_refusal(
    reason_code = "selection_parallel_worker_failed",
    description = paste0(
      "fb_select(): under options(flexyBayes.isolation = \"fork\"), a forked ",
      "worker running one fold ended without returning a result and left ",
      "no error; a forked worker's exit status is not reported, so memory ",
      "reclaimed by the operating system and a crash in compiled code are ",
      "both possible causes. The scores are not aggregated from the ",
      "folds that survived: a cross-validation missing a fold is not the ",
      "cross-validation that was asked for. A fold that refuses inside a ",
      "live worker is unaffected, and is re-raised under its own reason ",
      "code."
    ),
    message_template = paste0(
      "fb_select(): the forked worker running fold %s ended without ",
      "returning a result and left no error."
    ),
    registered_in_adr = "PLAN_green_0.12.1_to_0.15.0_2026-09-16 section 6 (P4)",
    plan_field = "selection_plan",
    since_version = "0.13.0"
  )

  invisible(NULL)
}
