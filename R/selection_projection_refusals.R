# Reference-fit and projection refusals --- registry population
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 10.
#
# Six codes, all of them gates that run before a decision is reported: two
# on the reference model (its sampler and its fit to the data), two on the
# projection provider (the structure it was handed and the family), one on
# the scoring unit the elpd is aggregated at, and one on the search method.
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

  .register_refusal(
    reason_code = "selection_scoring_unit_mismatch",
    description = paste0(
      "The scoring unit key and the provider's pointwise log predictive ",
      "densities name different numbers of observations. The aggregation ",
      "used to fall back to row numbers and carry on, which reports a ",
      "standard error at the row rather than at the resampling unit -- too ",
      "small by roughly the square root of the cluster size, and silently, ",
      "on exactly the clustered models the unit was declared for. Both ",
      "lengths are named so the caller can see which side moved."
    ),
    message_template = paste0(
      "fb_select(): the scoring unit key names %s observation(s) and the ",
      "provider returned %s pointwise log predictive densities."
    ),
    registered_in_adr = "SPEC v5 section 6.1",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  .register_refusal(
    reason_code = "projection_search_unsupported",
    description = paste0(
      "A projection search method with no executed gate was asked for. ",
      "`forward` is the search this release runs and validates. `L1` is a ",
      "penalised path over the projection, and this module forces the ",
      "protected and random terms into every submodel through a ",
      "`search_terms` construction written for a forward search; no gate ",
      "has run L1 against that construction and read the ranking back, so ",
      "whether the forced block survives it is unmeasured. An unrun search ",
      "refuses exactly as an unrun structure does."
    ),
    message_template = "fb_decision(): search = \"%s\" has no executed gate.",
    registered_in_adr = "SPEC v5 section 6.1",
    plan_field = "selection_plan",
    since_version = "0.11.0"
  )

  invisible(NULL)
}
