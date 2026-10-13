# fb_refit() refusals: registry population for R/fb_refit.R.
#
# Thrown by fb_refit() and triggered in tests/testthat/test-fb-refit.R.

.populate_refusal_registry_refit <- function() {
  .register_refusal(
    reason_code = "refit_not_reconstructable",
    description = paste0(
      "fb_refit(): the refit cannot be built from what the fit records. ",
      "The fit keeps no fixed formula or no data, or it carries ",
      "observation weights and new data arrived without `weights`, so ",
      "the recorded weights cannot be matched to the new rows."
    ),
    message_template = "%s",
    registered_in_adr = "0.18.0 refit",
    plan_field = NA_character_,
    since_version = "0.15.0.9000"
  )
  .register_refusal(
    reason_code = "refit_hooks_require_brms",
    description = paste0(
      "fb_refit(): `extra_prior`, `stanvars` and `formula_hook` are brms ",
      "objects and reach only a brms fit. The fit was not fitted with ",
      "brms and no `backend = \"brms\"` was given."
    ),
    message_template = "%s",
    registered_in_adr = "0.18.0 refit",
    plan_field = NA_character_,
    since_version = "0.15.0.9000"
  )
}
