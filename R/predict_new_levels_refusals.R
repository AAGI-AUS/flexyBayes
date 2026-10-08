# New-level prediction refusals --- registry population for
# R/predict_new_levels.R.
#
# Thrown by predict(allow_new_levels = TRUE) on either engine, and triggered
# in tests/testthat/test-predict-new-levels.R.

.populate_refusal_registry_predict_new_levels <- function() {
  .register_refusal(
    reason_code = "predict_new_levels_structure_not_supported",
    description = paste0(
      "predict(allow_new_levels = TRUE) on a model carrying a random or ",
      "residual structure, observation weights or a known variance for ",
      "which the draw at a new level is not built. The refusal names them."
    ),
    message_template = "%s",
    registered_in_adr = "PLAN_0.16.0 C1",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
  .register_refusal(
    reason_code = "predict_new_levels_no_variance",
    description = paste0(
      "predict(allow_new_levels = TRUE) at a new level the fit has no ",
      "variance for: a new level of f in diag(f):g, a new section of ",
      "dsum(~ units | f), or a row or column off the lattice of an AR1 x ",
      "AR1 field."
    ),
    message_template = "%s",
    registered_in_adr = "PLAN_0.16.0 C1",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
  .register_refusal(
    reason_code = "predict_new_levels_family_not_supported",
    description = paste0(
      "predict(allow_new_levels = TRUE) on a fit whose response is not ",
      "Gaussian on the identity link."
    ),
    message_template = "%s",
    registered_in_adr = "PLAN_0.16.0 C1",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
}
