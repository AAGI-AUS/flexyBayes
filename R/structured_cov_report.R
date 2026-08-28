# Identified-quantity reporting for factor-analytic structured covariance.
#
# A factor-analytic term fa(outer, k) decomposes the outer-factor
# covariance as G = Lambda Lambda' + diag(psi), with Lambda an
# (n_outer x k) loadings matrix. The loadings themselves are only
# identified up to an orthogonal rotation and a per-column sign flip, so
# the marginal posterior of an individual Lambda entry is multimodal and
# its Rhat is meaningless -- a large Rhat on a raw loading is expected and
# is not evidence of non-convergence. The implied covariance G (and the
# correlation derived from it) is invariant to that rotation and sign, so
# G would be the identified quantity to report.
#
# No active backend fits an fa() term: INLA and brms both refuse a
# factor-analytic structured-covariance term before a fit object exists
# (.capability_inla() / lgm_gate() and .capability_brms()), so no surviving
# fit ever carries a fa_gxe entry in parse_info$random -- the reconstruction
# below is therefore never reached in practice, and this function reports
# that there is nothing to summarise.

#' Identified covariance for factor-analytic structured-covariance terms
#'
#' For each `fa(outer, k)` term in a fit, reconstruct the implied
#' outer-factor covariance \eqn{G = \Lambda\Lambda^\top + \mathrm{diag}(\psi)}
#' from the posterior draws and summarise it. Unlike the raw loadings
#' \eqn{\Lambda} -- which are identified only up to rotation and sign, so
#' their per-entry Rhat is meaningless -- the covariance \eqn{G} and the
#' correlation derived from it are rotation- and sign-invariant.
#'
#' @param fit A flexybayes fit.
#'
#' @return An empty list (with a message): no active backend fits an
#'   `fa()` term, so no fit ever carries one to reconstruct. Non-factor-
#'   analytic structured terms (`us`, `ar1`) are reported as
#'   not-yet-reconstructed.
#'
#' @export
fb_structured_cov <- function(fit) {
  .check_flexybayes_fit(fit, "`fit` must be a flexybayes object.")
  rt <- fit$extras$parse_info$random %||% list()
  other_struct <- Filter(
    function(t) (t$type %||% "") %in% c("us_gxe", "ar1_spatial"),
    rt
  )

  # No active backend fits an fa() term (see file header), so a
  # surviving fit carries none and this message-and-return path is the
  # only reachable outcome.
  if (length(other_struct)) {
    message(
      "flexyBayes: fb_structured_cov() reconstructs the identified ",
      "covariance for factor-analytic fa() terms; reconstruction for ",
      "us()/ar1() terms is not yet implemented."
    )
  } else {
    message(
      "flexyBayes: this fit carries no factor-analytic structured-",
      "covariance term; nothing to report."
    )
  }
  invisible(list())
}
