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
# One code lands with 0.14.0, for a missing response on the routes that
# estimate a fold's covariance by REML:
#
#   selection_missing_response_not_carried
#                                       R/fb_select.R,
#                                       .fb_select_require_observed_response()
#                                       test-fb-select-missing-response.R
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
      "lgm_gate() on a model routed at INLA. flexyBayes fits joint ",
      "shrinkage priors on brms only and has no INLA route for them. No ",
      "public call reaches this refusal: fb_select() fits its reference ",
      "model on brms whatever engine fitted the base. The check is kept as ",
      "an internal guard."
    ),
    message_template = paste0(
      "joint shrinkage priors are fitted on brms only (%s)."
    ),
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
    reason_code = "selection_missing_response_not_carried",
    description = paste0(
      "fb_select(): the base's data have a missing response, and the ",
      "unpenalised benchmark, the penalised optimiser or the comparison ",
      "`benchmark = TRUE` adds was asked for. Each reads the response of ",
      "every training row to estimate a fold's covariance and of every ",
      "held-out row to score it."
    ),
    message_template = paste0(
      "fb_select(): the base's data have a missing response in %s, and %s ",
      "does not carry missing responses."
    ),
    registered_in_adr = "0.14.0 release plan",
    plan_field = "selection_plan",
    since_version = "0.14.0"
  )

  invisible(NULL)
}
