# INLA configuration refusals --- registry population for the configuration
# Gaussian in R/methods_inla.R.
#
# Thrown by vcov(), predict(se.fit = TRUE), predict(classify = ) and the
# posterior draws of an INLA fit, and triggered in
# tests/testthat/test-inla-config-constrained.R.

.populate_refusal_registry_inla_config <- function() {
  .register_refusal(
    reason_code = "inla_config_precision_not_positive_definite",
    description = paste0(
      "The latent-field precision INLA stored for a hyperparameter ",
      "configuration is not positive definite, and the fit's linear ",
      "constraints do not make it so: the posterior of the latent field ",
      "is improper in some direction."
    ),
    message_template = "%s",
    registered_in_adr = "PLAN_0.16.0 DEF-205",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
}
