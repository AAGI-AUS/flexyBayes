# Optimiser-provider refusals --- registry population for the penalised
# optimiser workflow.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 10.
#
# Two codes registered here (the absent-provider and withdrawn-family codes
# are shared with the Bayesian route and the grammar),
# each with a triggering test in
# `tests/testthat/test-selection-glmnet.R` or `test-selection-bridges.R`. A
# registered code nobody has seen fire is a string, not a refusal, so none of
# these is registered ahead of the site that raises it.

.populate_refusal_registry_selection_glmnet <- function() {
  # `selection_provider_absent` is registered once, by the Bayesian route
  # (R/selection_bayes_cv_refusals.R), with a description covering every
  # provider; the optimiser raises it without re-registering it.

  # `selection_penalty_family_withdrawn` is registered once, by the grammar
  # (R/penalty_refusals.R); the optimiser raises it without re-registering it.

  .register_refusal(
    reason_code = "selection_prox_mode_not_converged",
    description = paste0(
      "fb_select(): the internal proximal-gradient mode finished above its ",
      "KKT tolerance, so no point is returned. Probe B4 found an optimiser ",
      "reporting convergence 0.16 away from a lasso optimum; the certificate ",
      "here is the residual measured at the returned point, and a point that ",
      "fails it is withheld rather than labelled a mode."
    ),
    message_template = "fb_select(): the proximal-gradient mode stopped with a KKT residual of %s, above the tolerance %s.",
    registered_in_adr = "SPEC v5 section 4.3",
    plan_field = NA_character_,
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "selection_penalty_blocks_not_expressible",
    description = paste0(
      "fb_select(): the resolved penalty blocks cannot be written as one ",
      "`glmnet` problem. Its `alpha` is scalar and its `penalty.factor` is ",
      "renormalised to sum to the column count, so two blocks with different ",
      "families -- or the same family at different hyperparameters -- are not ",
      "the problem the printed lambda names. The Bayesian workflow carries ",
      "block-specific priors exactly; the optimiser benchmark does not."
    ),
    message_template = "fb_select(): the penalty blocks cannot be expressed as one glmnet problem: %s.",
    registered_in_adr = "SPEC v3 sections 5.2 and 6",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  invisible(NULL)
}
