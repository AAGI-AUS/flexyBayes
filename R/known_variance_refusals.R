# Known-variance refusals --- registry population for R/known_variance.R.
#
# Both codes are thrown by R/known_variance.R and by the fb_select() routes
# that cannot carry a known variance, and triggered in
# tests/testthat/test-known-variance.R.

.populate_refusal_registry_known_variance <- function() {
  .register_refusal(
    reason_code = "known_variance_invalid",
    description = paste0(
      "`known_variance` is not one finite, non-negative variance for each ",
      "row of `data`: it is not numeric, its length differs from the row ",
      "count, or a value is missing, infinite or negative."
    ),
    message_template = paste0(
      "`known_variance` must be a numeric vector with one finite, ",
      "non-negative variance for each row of `data`, and %s."
    ),
    registered_in_adr = "PLAN_AND_UNIT_BRIEFS_0.14.0_2026-09-19 (U-STAGE2-VAR)",
    plan_field = NA_character_,
    since_version = "0.14.0"
  )

  .register_refusal(
    reason_code = "known_variance_not_supported",
    description = paste0(
      "A known variance per row was asked of a model or a route that cannot ",
      "carry it: a family other than the Gaussian on the identity link, ",
      "the aggregated route, or an fb_select() route whose fold refit, ",
      "covariance or predictive score would leave it out. The message ",
      "names the route that carries it."
    ),
    message_template = "%s",
    registered_in_adr = "PLAN_AND_UNIT_BRIEFS_0.14.0_2026-09-19 (U-STAGE2-VAR)",
    plan_field = NA_character_,
    since_version = "0.14.0"
  )
}
