# The public verb's own refusals --- registry population for fb_select(),
# fb_decision() and the INLA prior gate.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 3.1, 3.3, 4.5
# and 10.
#
# Three codes land with the integration tranche, and each is thrown by a
# named site and fired by a named test:
#
#   selection_fallback_not_authorised   R/fb_select.R, the provider check
#                                       test-fb-select-bayes.R
#   prior_family_unsupported_inla       R/lgm_gate.R check 12, an internal
#                                       guard no public call reaches
#                                       test-lgm-gate.R, test-fb-select-bayes.R
#   selection_decision_invalid          R/fb_decision.R, at construction
#                                       test-fb-decision.R
#
# One code lands with the score contract that identifies values (0.13.0):
#
#   selection_score_contract_unverifiable
#                                       R/fb_select.R,
#                                       .fb_select_verify_contract()
#                                       test-fb-select-score-contract.R
#
# Two codes land with the Bayesian refit of each fold's frozen support
# (0.14.0):
#
#   selection_refit_support_undefined   R/selection_value.R, the route and
#                                       the support at the chosen point
#                                       test-selection-value.R
#   selection_value_contract_mismatch   R/selection_value.R,
#                                       .fb_sel_value_pair()
#                                       test-selection-value.R
#
# A registered code with no triggering test is a string, not a refusal.

.populate_refusal_registry_selection_verb <- function() {
  .register_refusal(
    reason_code = "selection_fallback_not_authorised",
    description = paste0(
      "fb_select(): the Bayesian workflow's provider is not installed and ",
      "`fallback = \"none\"`, which is the default. Nothing is substituted. ",
      "A penalised optimiser at a mode is a different estimator from a ",
      "posterior under a shrinkage prior, and there is no configuration ",
      "under which one is reported under the other's label; ",
      "`fallback = \"glmnet\"` runs the optimiser under its own name, ",
      "relabels the workflow, prints the substitution and certifies it."
    ),
    message_template = "fb_select(): the Bayesian workflow needs %s, which is not installed, and fallback = \"none\".",
    registered_in_adr = "SPEC v5 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "prior_family_unsupported_inla",
    description = paste0(
      "A joint shrinkage prior over the candidate fixed effects reached ",
      "lgm_gate() on a model routed at INLA. A Laplace, regularised ",
      "horseshoe or R2D2 prior is Gaussian given one local scale per ",
      "coefficient, so the model can be written as a latent Gaussian model, ",
      "but it then carries a hyperparameter for every coefficient in the ",
      "block, and INLA's numerical integration is built for models with a ",
      "few hyperparameters. Ridge and elastic-net blocks have no INLA ",
      "implementation. No public call reaches this refusal: fb_select() ",
      "fits its reference model on brms whatever engine fitted the base. ",
      "The check is kept as an internal guard."
    ),
    message_template = "the joint shrinkage prior(s) %s cannot be fitted by INLA.",
    registered_in_adr = "SPEC v5 section 4.5",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "selection_decision_invalid",
    description = paste0(
      "fb_decision(): the decision rule is missing the number it needs, or ",
      "carries one it does not read. A rule that has no tolerance is not a ",
      "tolerance rule and a posterior-threshold decision with no practical ",
      "effect size retains every candidate whatever the data say, so both ",
      "refuse at construction where the message can name the missing part ",
      "rather than several minutes into a sampling run."
    ),
    message_template = "fb_decision(): %s",
    registered_in_adr = "SPEC v5 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "selection_score_contract_unverifiable",
    description = paste0(
      "fb_select(): an arm's score contract cannot be verified, so its ",
      "out-of-fold score is not compared with another arm's. A contract is ",
      "verified when it carries this release's version field, an MD5 hash of ",
      "the canonical encoding of the values the arm was scored on, and a ",
      "hash equal to the one recomputed from its fields. A contract written ",
      "before 0.13.0 carries no version field and a polynomial checksum over ",
      "column names, which cannot tell a rescaled covariate from the original."
    ),
    message_template = "fb_select(): the score contract of the `%s` arm cannot be verified: %s.",
    registered_in_adr = "U-SCORE-1 (DEC-HASH-01)",
    plan_field = "selection_plan",
    since_version = "0.13.0"
  )

  .register_refusal(
    reason_code = "selection_refit_support_undefined",
    description = paste0(
      "fb_select(): refit = list(method = \"bayes\") refits the support the ",
      "penalised optimiser froze on each fold, and the request has none to ",
      "refit: the workflow is not the penalised optimiser, which alone ",
      "chooses a support on each fold's training rows, or the penalty ",
      "family, such as ridge(), shrinks every coefficient without setting ",
      "one to zero. The refusal comes before any refit."
    ),
    message_template = "fb_select(): refit = list(method = \"bayes\") has no frozen support to refit: %s.",
    registered_in_adr = "fb_select(refit = ), 0.14.0",
    plan_field = "selection_plan",
    since_version = "0.14.0"
  )

  .register_refusal(
    reason_code = "selection_value_contract_mismatch",
    description = paste0(
      "glance() or augment() on an fb_select result was asked to difference ",
      "its valuation against another, and the two value contracts differ: ",
      "the rows each fold holds out, the scoring unit, the prediction target ",
      "or the response values are not the same, or one result carries no ",
      "valuation. Units that were not held out together do not pair, so no ",
      "difference is formed. The fields that differ are named."
    ),
    message_template = "fb_select(): the two valuations' value contracts differ in: %s.",
    registered_in_adr = "fb_select(refit = ), 0.14.0",
    plan_field = "selection_plan",
    since_version = "0.14.0"
  )

  invisible(NULL)
}
