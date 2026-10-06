# Selection design refusals --- registry population for oracles 5 and 6
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 10.2.
#
# Two codes, both from the contract's table, both with a triggering test in
# test-selection-design.R. A registered code with no triggering test is a
# defect by the module's own standard, so nothing is registered here in
# advance of the site that throws it.

.populate_refusal_registry_selection_design <- function() {
  .register_selection_refusal(
    reason_code = "selection_rank_deficient_protected",
    description = paste0(
      "fb_select(): the protected design is rank deficient after whitening. ",
      "The estimable space and the aliased columns are both named, because ",
      "a rank drop that is not reported is indistinguishable from a model ",
      "the user did not write."
    ),
    message_template = paste0(
      "fb_select(): the protected design is rank deficient after whitening ",
      "-- rank %d over %d columns."
    ),
    registered_in_adr = "SPEC v3 section 5.2",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_penalty_coordinate_basis_dependent",
    description = paste0(
      "fb_select(): a coordinate-wise L1 component was asked for over a ",
      "target whose coordinates are contrast columns. This release does ",
      "not promise contrast-invariant within-factor selection, so the ",
      "combination refuses rather than returning a support that a relevel ",
      "would change."
    ),
    message_template = paste0(
      "fb_select(): an L1 component (%s) acts one contrast column at a ",
      "time, and these targets are basis dependent: %s."
    ),
    registered_in_adr = "SPEC v3 section 5.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  invisible(NULL)
}
