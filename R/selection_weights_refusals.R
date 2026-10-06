# Selection weights refusals --- registry population for section 5.5
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 5.5 and 10.2.
#
# Section 10.2 lists the selection module's refusal vocabulary. The codes
# below are the ones the weight, offset and missing-response absorption
# layer can throw. They are absent from that table because it was written
# before the absorption point existed as a named object. Each is
# registered here and each has a test that fires it --- a registered code
# with no triggering test is a defect by the same standard the module
# applies to everything else.
#
# Nothing here invents a second vocabulary for a fact the package already
# names. Missing responses reuse `missing_response_refused` and the
# augment / omit / fail estimand names from R/na_action.R, and the weight
# semantics is the package's own Var(y_i) = sigma^2 / w_i (see
# .refuse_unsupported_weights() in R/dispatch.R). Only the absorption
# accounting is new, because only it is new.

.populate_refusal_registry_selection_weights <- function() {
  .register_selection_refusal(
    reason_code = "selection_weights_not_positive",
    description = paste0(
      "fb_select(): an observation weight is not a finite positive ",
      "number. A zero weight deletes a row and a missing weight is not ",
      "a weight, so neither is silently absorbed."
    ),
    message_template = paste0(
      "fb_select(): observation weight(s) at %s are not finite and ",
      "positive."
    ),
    registered_in_adr = "SPEC v3 section 5.5",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_weights_absorbed_twice",
    description = paste0(
      "fb_select(): the weighting the absorbed object actually carries ",
      "is not the declared weighting applied exactly once. Weights live ",
      "in R, and therefore in V, or in the provider loss -- never both."
    ),
    message_template = paste0(
      "fb_select(): observation weights are absorbed %s, not exactly ",
      "once."
    ),
    registered_in_adr = "SPEC v3 section 5.5",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_weights_loss_route_requires_residual_only",
    description = paste0(
      "fb_select(): the provider-loss route reproduces weighted GLS ",
      "only when V has no random part. With random effects present the ",
      "two routes are different estimators, so the loss route refuses ",
      "rather than quietly answering a different question."
    ),
    message_template = paste0(
      "fb_select(): the provider-loss weight route needs a ",
      "residual-only V; %s."
    ),
    registered_in_adr = "SPEC v3 section 5.5",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_weights_not_carried",
    description = paste0(
      "fb_select(): the base model was fitted with observation weights, and ",
      "the requested route does not carry them into its fold fits and ",
      "scores. Running it would select on an unweighted model under the ",
      "weighted base's name."
    ),
    message_template = paste0(
      "fb_select(): the base carries observation weights, and %s does not ",
      "carry them."
    ),
    registered_in_adr = "PLAN 0.14.2 unit W0",
    plan_field = "selection_plan",
    since_version = "0.14.2"
  )

  .register_selection_refusal(
    reason_code = "selection_weights_by_unweighted_base",
    description = paste0(
      "fb_select(): `fb_resamples(weights_by = )` re-estimates the base's ",
      "observation weights in each fold, and the base was fitted without ",
      "weights."
    ),
    message_template = paste0(
      "fb_select(): `weights_by` needs a base fitted with `weights`."
    ),
    registered_in_adr = "PLAN 0.14.2 unit W3",
    plan_field = "selection_plan",
    since_version = "0.14.2"
  )

  .register_selection_refusal(
    reason_code = "selection_offset_absorbed_twice",
    description = paste0(
      "fb_select(): the Gaussian offset does not round-trip. It is ",
      "subtracted once before whitening and added back once to ",
      "predictions, and the object does not show that."
    ),
    message_template = "fb_select(): the offset does not round-trip; %s.",
    registered_in_adr = "SPEC v3 section 5.5",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_missing_response_estimand_unnamed",
    description = paste0(
      "fb_select(): rows whose response is missing may enter prediction ",
      "or covariance blocks only under a named estimand, and the name ",
      "given is not one the package carries."
    ),
    message_template = paste0(
      "fb_select(): `%s` is not a named missing-response estimand; the ",
      "named ones are %s."
    ),
    registered_in_adr = "SPEC v3 section 5.5",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_missing_response_in_objective",
    description = paste0(
      "fb_select(): a row carrying no observed outcome reached the ",
      "selection objective. Only observed outcomes enter the objective, ",
      "under every estimand."
    ),
    message_template = paste0(
      "fb_select(): %d row(s) with no observed response reached the ",
      "objective."
    ),
    registered_in_adr = "SPEC v3 section 5.5",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_preprocessing_leak",
    description = paste0(
      "fb_select(): a transformation the object records as learned on ",
      "the training rows does not reproduce from those rows, so it saw ",
      "rows it was not entitled to."
    ),
    message_template = paste0(
      "fb_select(): recorded preprocessing does not reproduce from its ",
      "own rows: %s."
    ),
    registered_in_adr = "SPEC v3 sections 4.1 and 5.5",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  invisible(NULL)
}
