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
    message_template = paste0(
      "fb_select(): the proximal-gradient mode stopped with a KKT residual of ",
      "%s, above the tolerance %s."
    ),
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
    message_template = paste0(
      "fb_select(): the penalty blocks cannot be expressed as one glmnet ",
      "problem: %s."
    ),
    registered_in_adr = "SPEC v3 sections 5.2 and 6",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "selection_penalty_group_structure_absent",
    description = paste0(
      "fb_select(): a grouped penalty family was declared over a column set ",
      "in which the declared grouping makes every group a singleton, so the ",
      "grouping does nothing: a group lasso would be a lasso and an ",
      "exclusive lasso a ridge. The family's answer would then be a ",
      "different penalty than the one its name claims, which is the silent ",
      "substitution this module refuses everywhere else. Declare the ",
      "grouping the question needs -- `group = \"term\"` moves a factor's ",
      "contrast columns together, `group = \"block\"` makes the ",
      "specification's own targets one set of competitors."
    ),
    message_template = paste0(
      "fb_select(): the declared grouping names no group of more than one ",
      "column: %s."
    ),
    registered_in_adr = "final specification B9",
    plan_field = "selection_plan",
    since_version = "0.13.0"
  )

  .register_refusal(
    reason_code = "selection_penalty_family_optimizer_only",
    description = paste0(
      "fb_penalty(): a grouped family -- group lasso, sparse group lasso, ",
      "exclusive lasso or the group-exclusive composite -- was asked for on ",
      "the Bayesian route. These families have a penalised-optimiser ",
      "provider and no brms one: each is a penalty carrying a proximal ",
      "operator, and brms takes a density. Writing the nearest available ",
      "prior under the requested family's name would be the substitution ",
      "the equivalence levels exist to prevent."
    ),
    message_template =
      "fb_penalty(): %s has an optimiser provider and no Bayesian one.",
    registered_in_adr = "final specification B9",
    plan_field = "selection_plan",
    since_version = "0.13.0"
  )

  invisible(NULL)
}
