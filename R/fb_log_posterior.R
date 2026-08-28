# fb_log_posterior() -- exposes a fitted model's log-posterior as a callable.
#
# This is the one inference-result outflow from flexyBayes: a fitted
# flexyBayes object is turned into a vectorised, domain-safe, unnormalised
# log-posterior callable that a downstream tool (for example the proxymix
# package, via `proxymix::from_fb_posterior()`) can compress into a
# closed-form Gaussian-mixture proxy. flexyBayes owns the producer; the
# consumer never depends on flexyBayes at run time (an acyclic dependency
# invariant), so this file emits the *callable*, read by the consumer as a
# bare callable carrying its metadata as attributes.
#
# No active backend currently produces a real C4 log-density. brms's
# log-density (`rstan::log_prob`) lives on the Stan unconstrained scale, and
# the parameter-name to Stan-unconstrained-index mapping is version-fragile
# (transformed parameters, declaration order); a wrong mapping would yield a
# wrong log-density silently, so abstaining beats emitting a
# plausible-but-wrong producer (Independent Oracle Principle). INLA's
# posterior is a deterministic Laplace/grid approximation, not a sampling
# log-density at arbitrary natural-scale points, so it abstains too. Both
# raise an informative, classed condition rather than guessing.

# ---------------------------------------------------------------- #
# The generic                                                      #
# ---------------------------------------------------------------- #

#' Emit a flexyBayes posterior as a log-density producer
#'
#' Turns a fitted flexyBayes object into a log-posterior producer:
#' a vectorised, domain-safe, unnormalised log-posterior callable that
#' `proxymix::from_fb_posterior()` compresses into a closed-form
#' Gaussian-mixture proxy. It is the single inference-result outflow from
#' flexyBayes; the contract is the *log-density*, not the draws, so the
#' returned object is addressed purely through its callable.
#'
#' The returned value is a **bare callable**
#' `function(theta_matrix) -> numeric`. Its input is a numeric matrix whose
#' rows index independent parameter draws and whose columns index
#' parameters, in `attr(., "parameter_names")` order, on the natural
#' (constrained) scale. Its output is a length-`nrow(theta_matrix)` numeric
#' vector of `log p(theta | data) + const` (unnormalised). The callable is
#' vectorised, side-effect free, and domain-safe: a row outside the
#' parameters' support returns `-Inf` rather than raising an error (the
#' consumer probes it at construction). It carries, as attributes:
#'
#' \describe{
#'   \item{`parameter_names`}{Character vector naming the parameters; its
#'     length fixes the proxy's ambient dimension. Vector-valued targets
#'     are flattened in column-major order with index suffixes (e.g.
#'     `beta[1,1]`, `beta[2,1]`).}
#'   \item{`log_normalizer`}{The additive correction that would normalise
#'     the density, i.e. `-log Z`. For a posterior the marginal likelihood
#'     is generally unknown, so this is `NA_real_` and the consumer
#'     reports a shifted (not absolute) divergence.}
#'   \item{`support_lower`, `support_upper`}{Length-`n_dim` numeric support
#'     bounds taken from the model's parameter constraints (`NA` for an
#'     unbounded coordinate). A variance / scale parameter is bounded below
#'     by zero, for instance. Used only to centre and scale the consumer's
#'     default importance proposal.}
#'   \item{`draws`}{An `n` by `n_dim` numeric matrix of the fit's posterior
#'     draws on the natural scale, column-aligned to `parameter_names`.
#'     Used only to seed the consumer's default proposal; never required.}
#' }
#'
#' Backend support. No active backend produces a real C4 log-density today.
#' The **brms** and **INLA** backends both abstain with an informative
#' condition -- brms's log-density lives on the Stan unconstrained scale
#' with a version-fragile name mapping, and INLA's posterior is a
#' deterministic approximation, not a sampling log-density. Abstaining is
#' preferred to emitting a plausible-but-wrong log-density.
#'
#' Acyclic note. A consumer such as proxymix uses this callable without
#' depending on flexyBayes; flexyBayes does not list proxymix in `Imports`
#' or `Suggests`. The cross-package demonstration lives in a separate
#' integration harness, not in this package, preserving the acyclic
#' dependency graph.
#'
#' @param fit A fitted flexyBayes object. Every current backend class
#'   abstains; the generic and its methods are retained so a consumer can
#'   dispatch on the abstention rather than discovering the gap by a
#'   missing method.
#' @param ... Reserved for future producer options; currently unused.
#'
#' @return Does not return on any current backend: it raises a classed
#'   `fb_c4_unavailable` condition naming the backend and the reason.
#'
#' @family flexyBayes interop
#' @seealso `proxymix::from_fb_posterior()` for the consumer (compresses
#'   the returned callable into a Gaussian-mixture proxy).
#' @section Lifecycle:
#' No active backend produces a real C4 log-density; both `flexybayes_brms`
#' and `flexybayes_inla` fits abstain with the classed condition
#' `fb_c4_unavailable`, for the reasons given above. The generic and its
#' methods are retained so a consumer can dispatch on the abstention rather
#' than discovering the gap by a missing method.
#'
#' @keywords internal
fb_log_posterior <- function(fit, ...) {
  UseMethod("fb_log_posterior")
}

# ---------------------------------------------------------------- #
# Abstain: every current backend                                   #
# ---------------------------------------------------------------- #

#' @rdname fb_log_posterior
#' @export
fb_log_posterior.default <- function(fit, ...) {
  .fb_c4_abstain(
    backend = "this object class",
    detail = paste0(
      "fb_log_posterior() produces a C4 log-density only for a fitted ",
      "flexyBayes object. Got class: ",
      paste(class(fit), collapse = "/"),
      "."
    )
  )
}

#' @rdname fb_log_posterior
#' @export
fb_log_posterior.flexybayes_brms <- function(fit, ...) {
  .fb_c4_abstain(
    backend = "the brms backend",
    detail = paste0(
      "C4 log-density producer not available for the brms backend in ",
      "this version. brms's log-density (rstan::log_prob) is defined on ",
      "the Stan unconstrained scale, and the mapping from flexyBayes ",
      "parameter names to Stan unconstrained indices is version-fragile ",
      "(transformed parameters, declaration order). Abstaining is ",
      "preferred to emitting a plausible-but-wrong log-density. No ",
      "active backend produces a C4 log-density."
    )
  )
}

#' @rdname fb_log_posterior
#' @export
fb_log_posterior.flexybayes_inla <- function(fit, ...) {
  .fb_c4_abstain(
    backend = "the INLA backend",
    detail = paste0(
      "C4 log-density producer not available for the INLA backend. INLA's ",
      "posterior is a deterministic Laplace / grid approximation, not a ",
      "sampling log-density evaluable at arbitrary natural-scale points, ",
      "so there is no faithful unnormalised log-posterior to emit. No ",
      "active backend produces one."
    )
  )
}

# Raise the classed abstention condition. Never returns.
.fb_c4_abstain <- function(backend, detail) {
  cond <- structure(
    class = c("fb_c4_unavailable", "error", "condition"),
    list(
      message = paste0(
        "No C4 log-density producer is available for ", backend, ".\n",
        detail
      ),
      call = sys.call(-1L)
    )
  )
  stop(cond)
}
