# Breeder-facing MET summary (G2). A factor-analytic G x E fit
# fa(env, k):gen answers the two questions a plant breeder asks of every
# multi-environment trial:
#
#   * overall performance (OP): which genotypes are best on average across
#     environments? -- the across-environment mean of the realised
#     genotype-by-environment effects.
#   * stability: which genotypes hold their performance across
#     environments, and which trade mean for sensitivity? -- the
#     across-environment spread of those effects.
#
# plus the environment genetic-correlation matrix (which environments rank
# genotypes alike, and which reverse them -- the crossover structure a
# reverse-U / G x E paradox investigation reads directly), and the
# genotype-by-environment BLUPs themselves.
#
# These would be computed from the *realised* effects g_mat = F Lambda' +
# delta (monitored by the fa codegen), which are identified -- rotation-
# and sign-invariant -- unlike the raw loadings. No active backend fits an
# fa() term (INLA and brms both refuse it as a structured-covariance term
# outside their representable set -- see .capability_inla() / lgm_gate()
# and .capability_brms()), so this function always abstains.

# --- the summary -------------------------------------------------- #

#' Breeder summary of a factor-analytic multi-environment-trial fit
#'
#' For a `fa(env, k):gen` factor-analytic G x E fit, summarise the
#' quantities a plant breeder acts on: each genotype's overall performance
#' (the across-environment mean of its realised effects) and stability (the
#' across-environment spread), the genotype-by-environment BLUPs, and the
#' environment genetic-correlation matrix (the crossover structure). The
#' realised effects are identified -- invariant to the rotation and sign
#' ambiguity of the raw loadings -- so their posterior summaries would be
#' interpretable.
#'
#' @param fit A flexybayes fit.
#' @param genotype_levels,environment_levels Optional character labels for
#'   the inner (genotype) and outer (environment) factors; unused (see
#'   Lifecycle).
#'
#' @return Does not return: raises the classed `met_summary_not_available`
#'   refusal.
#'
#' @seealso [fb_structured_cov()] for the identified environment covariance
#'   and its convergence diagnostic.
#' @section Lifecycle:
#' No active engine emits an `fa(env, k):gen` term -- both INLA and brms
#' refuse a factor-analytic structured-covariance term before a fit object
#' exists -- so this function always abstains (`met_summary_not_available`).
#' What an active engine does report for a multi-environment trial is the
#' variance components, through [summary()], and on brms the
#' genotype-by-environment covariance of a `diag()` or `us()` term through
#' `brms::VarCorr()`.
#'
#' @keywords internal
fb_met_summary <- function(
  fit,
  genotype_levels = NULL,
  environment_levels = NULL
) {
  .check_flexybayes_fit(fit, "`fit` must be a flexybayes object.")
  stop(.fb_refusal_condition(
    reason_code = "met_summary_not_available",
    message = paste0(
      "fb_met_summary() is not available: breeder summaries are computed ",
      "from realised factor-analytic effects, and no active backend emits ",
      "an fa(env, k):gen term (INLA and brms both refuse a factor-",
      "analytic structured-covariance term before a fit object exists). ",
      "What an active engine does report for a multi-environment trial is ",
      "the variance components, through summary(), and on brms the ",
      "genotype-by-environment covariance of a diag() or us() term ",
      "through brms::VarCorr()."
    )
  ))
}
