# Marginal-covariance refusal: registry population for R/marginal_cov.R.
#
# Thrown by fb_marginal_cov() and its as.matrix() method, and triggered in
# tests/testthat/test-marginal-cov.R.

.populate_refusal_registry_marginal_cov <- function() {
  .register_refusal(
    reason_code = "marginal_cov_not_buildable",
    description = paste0(
      "fb_marginal_cov(): the model carries a random or residual structure ",
      "the marginal covariance is not built for, its variance-component ",
      "table lacks a variance the structure needs, or as.matrix() was asked ",
      "of a structure with no variances. The refusal names the terms."
    ),
    message_template = paste0(
      "fb_marginal_cov(): the marginal covariance cannot be built for %s."
    ),
    registered_in_adr = "0.18.0 marginal covariance",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
}
