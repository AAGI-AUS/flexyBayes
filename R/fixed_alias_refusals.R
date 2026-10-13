# Fixed-alias refusal: registry population for R/fixed_alias.R.
#
# Thrown by flexybayes() and triggered in tests/testthat/test-fixed-alias.R.

.populate_refusal_registry_fixed_alias <- function() {
  .register_refusal(
    reason_code = "fixed_effects_aliased",
    description = paste0(
      "flexybayes(): fixed-effect columns are linear combinations of the ",
      "columns before them on the rows with a response, as in a factorial ",
      "with empty cells, so the data give no estimate for them and an ",
      "engine would report them at their prior. The refusal names the ",
      "columns, fixed_alias = \"drop\" drops them as asreml does."
    ),
    message_template = paste0(
      "flexybayes(): %s fixed-effect columns have no estimate on the rows ",
      "with a response: %s."
    ),
    registered_in_adr = "0.14.1 aliased fixed effects",
    plan_field = NA_character_,
    since_version = "0.14.1"
  )
}
