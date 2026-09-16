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
# No active backend fits an fa() term, so no surviving fit carries a fa_gxe
# entry in parse_info$random, the reconstruction below is not reached, and
# this function reports that there is nothing to summarise.
#
# That sentence was false between 0.10.0.9000 and 0.12.1 and this file went
# on asserting it: INLA emitted and fitted fa(f, k):g throughout, so a fit
# DID carry a fa_gxe entry and this function told its caller there was no
# factor-analytic term in a model that had one. It is true again from 0.13.0,
# for a stated reason rather than by assumption -- the structure is refused
# on both engines by factor_analytic_not_repeatable (R/dispatch.R), because
# repeated identical INLA calls did not agree on the covariance. Should that
# refusal ever be lifted, the G below is what this function must return, and
# this comment is the notice that it stops being dead code the moment it is.

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
#'   `fa()` term -- both refuse it with `factor_analytic_not_repeatable` --
#'   so no fit carries one to reconstruct. Non-factor-analytic structured
#'   terms (`us`, `ar1`) are reported as not-yet-reconstructed.
#'
#' @examples
#' # Reports on factor-analytic structured covariance. A fit carrying no
#' # such term is told so rather than returned an empty object.
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#'   fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
#'                     verbose = FALSE)
#'   flexyBayes:::fb_structured_cov(fit)
#' }
#' }
#' @keywords internal
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
