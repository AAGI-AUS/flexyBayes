# Penalty refusals --- registry population for the joint-prior layer
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 10, Unit A.
#
# The same posture as the rest of the selection module: a code is registered
# when its throwing site lands, and every code here is fired by a test in
# tests/testthat/test-penalty-grammar.R, test-penalty-compile.R or
# test-penalty-stancode.R. A registered code nobody has watched fire is a
# string, not a refusal.
#
# Four codes the penalty grammar needs already exist and are reused rather
# than respelled: grammar_argument_duplicated, grammar_named_argument_unknown,
# grammar_arity_too_many and grammar_argument_empty. SPEC v3 section 6 names
# the first of those as the constructor-level refusal, and the defect class is
# identical -- a call whose arguments do not bind.

.populate_refusal_registry_penalty <- function() {
  .register_selection_refusal(
    reason_code = "selection_penalty_spec_malformed",
    description = paste0(
      "fb_penalty(): an argument is not a two-sided `target ~ family(...)` ",
      "formula, or no specification was written at all."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 3.2",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "selection_penalty_target_unsupported",
    description = paste0(
      "fb_penalty(): the left-hand side is not one of the four supported ",
      "term selectors, or a selector was written with the wrong shape. ",
      "Call-shaped pseudo-selectors are refused rather than deparsed, ",
      "because deparsing a call head-first is what let `b(c(\"x\", \"z\"))` ",
      "become three coefficient names in the prior DSL."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 3.2",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "selection_penalty_target_not_in_model",
    description = paste0(
      "fb_penalty(): a named target is not a candidate term of the model ",
      "being selected over -- it is absent, or it is protected."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 3.2",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "selection_penalty_target_ambiguous",
    description = paste0(
      "fb_penalty(): two equally broad selectors name the same coordinate, ",
      "so the later-overrides-earlier rule cannot say which family applies."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 3.2",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "selection_penalty_family_unknown",
    description = paste0(
      "fb_penalty(): the right-hand side names no shrinkage family this ",
      "release compiles."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 3.2",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "selection_penalty_family_withdrawn",
    description = paste0(
      "fb_penalty(): a group, exclusive or adaptive family of the v3 ",
      "grammar was named. Their providers are not shipped in this release, ",
      "and a family without a provider is a promise rather than a method."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 8.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "penalty_hyperparameter_invalid",
    description = paste0(
      "fb_penalty(): a hyperparameter is outside the family table's domain, ",
      "or `p0` was given to a family whose grammar has no such quantity."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 sections 3.2 and 4.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "penalty_block_overlap",
    description = paste0(
      "fb_penalty(): one candidate design column ended up in two prior ",
      "blocks. Every resolved column carries exactly one role, and a column ",
      "under two densities is two priors on one coefficient."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 3.2",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "penalty_factor_contrast_ambiguous",
    description = paste0(
      "fb_penalty(): a coordinate-wise family reaches one contrast column ",
      "of a multi-column term at resolution = \"coefficient\". The support ",
      "would move under a relevel that changes no model."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 10",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "penalty_block_targeting_unsupported",
    description = paste0(
      "fb_penalty(): a class-wide family cannot be aimed at a subset of the ",
      "population-level coefficients. brms applies it to every one of them ",
      "except the intercept, so a protected or untargeted coefficient would ",
      "be shrunk without anyone asking for it."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 4.4",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "penalty_block_part_name_taken",
    description = paste0(
      "fb_penalty(): scope = \"block\" fits the model as the brms ",
      "non-linear parts `a` and `b`, and the model uses one of those names ",
      "as a variable, which brms would read as the part."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "Block-scoped class-wide prior (2026-10-01)",
    plan_field = "selection_plan",
    since_version = "0.15.0"
  )

  .register_selection_refusal(
    reason_code = "penalty_scaling_column_constant",
    description = paste0(
      "fb_penalty(): a candidate column has no spread inside the training ",
      "fold, so the train-fold standardisation has no scale to divide by."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_selection_refusal(
    reason_code = "penalty_candidate_columns_unscaled",
    description = paste0(
      "fb_penalty(): a class-wide global-local family -- the regularised ",
      "horseshoe or R2D2 -- was aimed at a candidate block containing a ",
      "design column the reference fit cannot standardise. The reference is ",
      "a re-emission of the formula rather than a design-matrix fit, so a ",
      "factor contrast or an interaction product has no data column to ",
      "centre and divide, and brms puts ONE global-local prior on the whole ",
      "population-level vector. A shared local scale then treats a 0/1 ",
      "dummy as comparable to a numeric already at unit standard deviation, ",
      "and the shrinkage applied is a property of the contrast coding. The ",
      "coefficient-row families state a scale per column and are not ",
      "affected; for them the departure stays an approximation on the ",
      "certificate."
    ),
    message_template = "fb_penalty(): %s",
    registered_in_adr = "SPEC v5 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  invisible(NULL)
}
