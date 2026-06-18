# fb_brms() -- the brms (Stan) engine pin.
#
# v0.5.0 backend-axis recovery. `fb_brms()` is one
# of the three engine pins (alongside `fb_inla()` and `fb_greta()`): sugar
# over the universal entry with the engine fixed --- `fb_brms(...)` ==
# `fb(..., backend = "brms")`. It fits via Stan (through brms) and accepts
# every grammar the universal entry does (an ASReml `fixed` / `random` /
# `rcov` specification or a brms-style bar-grouped formula); the grammar
# is detected from the call shape, exactly as on `fb()`.
#
# Before v0.5.0 `fb_brms()` was the brms-GRAMMAR verb and carried a
# `backend` argument spanning greta / inla / brms / auto. That
# multi-backend role moved to the universal entry `fb()` / `flexybayes()`
# (which now detects brms grammar and reaches every backend); the name is
# reused here as the brms-ENGINE pin. The brms engine cannot represent an
# ASReml structured-covariance term (vm / ped / fa / us / ar1) or a
# low_rank approximation; those raise the registry's structured refusals
# (the brms capability_predicate in R/backend_registry.R).
#
# The pin forwards the user's call verbatim via .fb_engine_pin() (the
# match.call() rewrite seam in R/fb_inla.R), so missing()-based defaults
# and the data-name capture behave exactly as on a direct fb() call. A
# conflicting `backend` argument raises engine_pin_backend_conflict; the
# redundant self-pin fb_brms(backend = "brms") is accepted (the v0.4.1
# deprecation notice promised it would survive the rename).

#' Fit a flexyBayes model via the brms (Stan) engine
#'
#' Engine pin: fits the model with Stan through brms only. This is sugar
#' for [flexybayes()]`(..., backend = "brms")` and accepts the same
#' arguments and grammars --- an ASReml `fixed` / `random` / `rcov`
#' specification or a brms-style bar-grouped formula (see [flexybayes()]
#' for the full argument list). flexyBayes builds the intermediate
#' representation, translates the prior, calls `brms::brm()`, and wraps
#' the result; the live `brmsfit` is available on the `$brms` slot for
#' brms's own posterior tooling (`loo()`, `posterior_predict()`,
#' `bayes_factor()`).
#'
#' The brms / Stan engine cannot represent an ASReml structured-covariance
#' term (`vm`, `ped`, `fa`, `us`, `ar1`) or a `low_rank` smooth
#' approximation; such a model raises a structured refusal naming the
#' offending construct. Re-fit with [fb_greta()] (full MCMC) or, when the
#' model is latent-Gaussian feasible, [fb_inla()].
#'
#' @param ... Arguments passed to [flexybayes()] (e.g. `formula` / `fixed`,
#'   `random`, `rcov`, `data`, `family`, `prior`, `syntax`). The `backend`
#'   argument is pinned to `"brms"`; a conflicting `backend` value raises a
#'   structured refusal (the redundant `backend = "brms"` is accepted).
#'   The pre-v0.5.0 `formula = ` argument is remapped to the universal
#'   entry's model-spec slot for call-compatibility.
#'
#' @return An object of class `"flexybayes_brms"` (a subclass of
#'   `"flexybayes"`) carrying the live `brmsfit` on `$brms`; see
#'   [flexybayes()] for the shared structure.
#'
#' @family flexyBayes engine pins
#' @seealso [flexybayes()] and [fb()] for the universal entry that picks a
#'   backend; [fb_greta()] / [fb_inla()] for the other engine pins;
#'   [fb_from_brms()] for building a brms-grammar IR.
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE) &&
#'     requireNamespace("lme4", quietly = TRUE)) {
#'   data(sleepstudy, package = "lme4")
#'   fit <- fb_brms(Reaction ~ Days + (1 | Subject), data = sleepstudy,
#'                  chains = 1)
#'   coef(fit)
#' }
#' }
#' @export
fb_brms <- function(...) {
  cl <- match.call()
  # Call-compatibility: the pre-v0.5.0 `formula = ` argument names the
  # same model-spec slot as the universal entry's `fixed`.
  nm <- names(cl)
  if (!is.null(nm) && "formula" %in% nm) {
    if ("fixed" %in% nm) {
      stop(
        "fb_brms(): pass the model once -- `formula` and `fixed` name ",
        "the same model-spec slot.",
        call. = FALSE
      )
    }
    nm[nm == "formula"] <- "fixed"
    names(cl) <- nm
  }
  .fb_engine_pin("brms", cl, parent.frame())
}
