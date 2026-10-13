# Breeder-facing MET summary (G2). A factor-analytic G x E fit
# fa(env, k):gen answers the two questions a plant breeder asks of every
# multi-environment trial:
#
#   * overall performance (OP): which genotypes are best on average across
#     environments?: the across-environment mean of the realised
#     genotype-by-environment effects.
#   * stability: which genotypes hold their performance across
#     environments, and which trade mean for sensitivity? The
#     across-environment spread of those effects.
#
# plus the environment genetic-correlation matrix (which environments rank
# genotypes alike, and which reverse them, the crossover structure a
# reverse-U / G x E paradox investigation reads directly), and the
# genotype-by-environment BLUPs themselves.
#
# These would be computed from the *realised* effects g_mat = F Lambda' +
# delta, which are identified, rotation- and sign-invariant, unlike the
# raw loadings. No active backend fits an fa() term, so this function always
# abstains.
#
# Between 0.10.0.9000 and 0.12.1 that last sentence was false: INLA fitted
# fa(f, k):g and this function refused the breeder summary of a model the
# package had just fitted. It holds again from 0.13.0 because the structure
# is refused on both engines by factor_analytic_not_repeatable
# (R/dispatch.R), INLA fits of simulated trials missed the interval
# coverage required, and not because no engine can represent it.

# --- the summary -------------------------------------------------- #

#' Breeder summary of a factor-analytic multi-environment-trial fit
#'
#' For a `fa(env, k):gen` factor-analytic G x E fit, summarise the
#' quantities a plant breeder acts on: each genotype's overall performance
#' (the across-environment mean of its realised effects) and stability (the
#' across-environment spread), the genotype-by-environment BLUPs, and the
#' environment genetic-correlation matrix (the crossover structure). The
#' realised effects are identified, invariant to the rotation and sign
#' ambiguity of the raw loadings, so their posterior summaries would be
#' interpretable.
#'
#' @param fit A flexybayes fit.
#' @param genotype_levels,environment_levels Optional character labels for
#'   the inner (genotype) and outer (environment) factors. Unused (see
#'   Lifecycle).
#'
#' @return Does not return: raises the classed `met_summary_not_available`
#'   refusal.
#'
#' @seealso [vpredict()] for posterior functionals of the variance
#'   components, including genetic correlations.
#' @section Lifecycle:
#' No active engine fits an `fa(env, k):gen` term: both refuse it with
#' `factor_analytic_not_repeatable`. This function therefore always abstains
#' (`met_summary_not_available`).
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
      "from realised factor-analytic effects, and no active backend fits ",
      "an fa(env, k):gen term, both refuse it with ",
      "factor_analytic_not_repeatable. ",
      "What an active engine does report for a multi-environment trial is ",
      "the variance components, through summary(), and on brms the ",
      "genotype-by-environment covariance of a diag() or us() term ",
      "through brms::VarCorr()."
    )
  ))
}
