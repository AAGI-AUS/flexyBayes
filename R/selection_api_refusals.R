# Selection fold-data refusals: registry population for R/selection_api.R.
#
# The code is thrown by fb_select(keep_fold_data = TRUE) on the Bayesian
# workflow and by fb_select_fold_data(), fb_select_path() and
# fb_select_refit() on a run that kept no fold data, and triggered in
# tests/testthat/test-selection-api.R.

.populate_refusal_registry_selection_api <- function() {
  .register_refusal(
    reason_code = "selection_fold_data_not_kept",
    description = paste0(
      "A selection run's fold data was asked for and the run kept none: it ",
      "ran without keep_fold_data = TRUE, or on the Bayesian workflow, ",
      "whose folds are model refits with no whitened data."
    ),
    message_template = "%s",
    registered_in_adr = "PLAN_v10_two_stage_study_01Oct2026 (U-SEL-API)",
    plan_field = NA_character_,
    since_version = "0.15.0"
  )
}
