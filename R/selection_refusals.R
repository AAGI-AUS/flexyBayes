# Selection refusals --- registry population for the selection module
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 10.2.
#
# Codes are registered as their throwing sites land, never in advance. A
# registered code with no triggering test is a defect by the same standard
# the module applies to everything else: a refusal nobody has seen fire is
# not a refusal, it is a string.

.populate_refusal_registry_selection <- function() {
  .register_refusal(
    reason_code = "selection_covariance_not_buildable",
    description = paste0(
      "fb_select(): no admitted covariance builder applies to the base ",
      "fit, or the fit does not report what a builder needs."
    ),
    message_template = "fb_select(): %s",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_covariance_not_positive_definite",
    description = paste0(
      "fb_select(): a covariance block has no Cholesky factor. No jitter ",
      "is applied, because a silent repair here is misspecification."
    ),
    message_template = "fb_select(): block %d of the covariance is not positive definite.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_dense_reference_too_large",
    description = paste0(
      "fb_select(): the dense covariance reference was asked for above its ",
      "size cap. It is an oracle for small cases, not a production path."
    ),
    message_template = "fb_select(): the dense covariance reference is %d by %d, above the cap of %d.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_reduction_term_not_in_fit",
    description = paste0(
      "fb_select(): the random term the reduction was asked for is not in ",
      "the base fit, or its levels cannot be recovered from the data."
    ),
    message_template = "fb_select(): the base fit carries no random term `%s`.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_reduction_level_mismatch",
    description = paste0(
      "fb_select(): a random term reports a level that is absent from the ",
      "data, so effects and covariates cannot be matched by name."
    ),
    message_template = "fb_select(): random term `%s` reports level(s) absent from the data: %s.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_candidates_not_constant_within_term",
    description = paste0(
      "fb_select(): a candidate is not admissible for a term's reduction -- ",
      "it is not constant within the term's levels, or it does not vary ",
      "across them."
    ),
    message_template = "fb_select(): these candidates cannot enter the reduction for `%s`: %s.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_base_shape_unreadable",
    description = paste0(
      "fb_select(): the base fit does not carry the parsed model and data ",
      "the selection builders read, as on the aggregated route."
    ),
    message_template = "fb_select(): this base fit does not carry the parsed model and data the reduction reads.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  invisible(NULL)
}
