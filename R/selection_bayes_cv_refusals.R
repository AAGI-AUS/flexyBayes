# Bayesian selection refusals --- registry population for the per-fold
# engine-refit route.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 3.3, 7 and 10.
#
# Two codes land with this tranche, and each is thrown by
# R/selection_bayes_cv.R and triggered by a test in
# tests/testthat/test-selection-bayes-cv.R. A registered code with no
# triggering test is a string, not a refusal.
#
# `selection_provider_absent` is registered here because this is the first
# module that needs it, and it is described for the whole provider set the
# contract names -- brms with its Stan toolchain, projpred, glmnet -- rather
# than for the one arm that happens to raise it first, so the projection and
# optimiser units can raise it without re-describing it.

.populate_refusal_registry_selection_bayes_cv <- function() {
  .register_refusal(
    reason_code = "selection_provider_absent",
    description = paste0(
      "fb_select(): the provider the requested workflow needs -- the ",
      "Bayesian engine that fitted the base and its toolchain, the ",
      "projection package, or the penalised optimiser -- is not installed. ",
      "No substitute provider is chosen: a different engine is a different ",
      "answer under the same name, so the refusal names the package and the ",
      "install line instead."
    ),
    message_template = "fb_select(): the %s workflow needs package `%s`, which is not installed.",
    registered_in_adr = "SPEC v5 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "penalty_runtime_budget_exceeded",
    description = paste0(
      "fb_select(): the estimated cost of the requested route, measured as ",
      "one timed engine fit times the number of fits the folds imply, ",
      "exceeds the declared `budget`. The estimate is formed and compared ",
      "before the fold loop starts, so nothing is fitted and no partial ",
      "result is returned."
    ),
    message_template = "fb_select(): the Bayesian route needs %d engine refits, an estimated %ss against a budget of %ss.",
    registered_in_adr = "SPEC v5 section 10",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  invisible(NULL)
}
