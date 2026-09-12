# Reference-fit and projection refusals --- registry population
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 10.
#
# Four codes, all of them gates that run before a decision is reported:
# two on the reference model (its sampler and its fit to the data) and two
# on the projection provider (the structure it was handed and the family).
# Registered here as their throwing sites land, never in advance -- a
# registered code with no triggering test is a string, not a refusal.
#
# The structure codes carry the section 6.3 matrix's own discipline: a row
# that has not been executed refuses exactly as a row that was executed and
# came back red. Unrun is not a permission.

.populate_refusal_registry_selection_projection <- function() {
  .register_refusal(
    reason_code = "penalty_sampler_diagnostics_failed",
    description = paste0(
      "The reference model's convergence contract returned `fail`. A ",
      "selection decision read off a posterior the sampler did not explore ",
      "is a decision about the sampler."
    ),
    message_template = "fb_select(): the reference model did not converge: %s.",
    registered_in_adr = "SPEC v5 section 5",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "penalty_reference_model_inadequate",
    description = paste0(
      "The reference model's posterior predictive check puts the observed ",
      "mean or standard deviation in the tail of the replicated data. A ",
      "projection inherits the reference model's misfit."
    ),
    message_template = paste0(
      "fb_select(): the reference model does not describe the data: %s."
    ),
    registered_in_adr = "SPEC v5 section 5",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "projection_structure_unsupported",
    description = paste0(
      "The reference model carries a structure the projection provider ",
      "refuses, or one whose gate has not been executed. Both refuse: an ",
      "unrun cell of the structure matrix is not a supported cell."
    ),
    message_template = "fb_select(): projection does not support %s.",
    registered_in_adr = "SPEC v5 section 6.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "projection_family_unsupported",
    description = paste0(
      "The reference model's response family is not Gaussian with an ",
      "identity link, which is the only family this release projects."
    ),
    message_template = paste0(
      "fb_select(): projection is Gaussian-only in this release, and this ",
      "reference is %s/%s."
    ),
    registered_in_adr = "SPEC v5 section 6.3",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  invisible(NULL)
}
