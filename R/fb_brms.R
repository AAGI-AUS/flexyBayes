# fb_brms(): the brms (Stan) engine pin.
#
# v0.5.0 backend-axis recovery. `fb_brms()` is one
# of the two active engine pins (alongside `fb_inla()`): sugar
# over the universal entry with the engine fixed, `fb_brms(...)` ==
# `fb(..., backend = "brms")`. It fits via Stan (through brms) and accepts
# every grammar the universal entry does (an ASReml `fixed` / `random` /
# `residual` specification or a brms-style bar-grouped formula). The grammar
# is detected from the call shape, exactly as on `fb()`.
#
# The multi-backend role belongs to the universal entry `fb()` /
# `flexybayes()`, which detects brms grammar and reaches every backend.
# `fb_brms()` is the brms-ENGINE pin. The brms engine cannot represent an
# ASReml structured-covariance term (vm / ped / fa / us / ar1) or a
# low_rank approximation. Those raise the registry's structured refusals
# (the brms capability_predicate in R/backend_registry.R).
#
# The pin forwards the user's call verbatim via .fb_engine_pin() (the
# match.call() rewrite seam in R/fb_inla.R), so missing()-based defaults
# and the data-name capture behave exactly as on a direct fb() call. A
# conflicting `backend` argument raises engine_pin_backend_conflict. The
# redundant self-pin fb_brms(backend = "brms") is accepted.

#' Fit a flexyBayes model via the brms (Stan) engine
#'
#' Engine pin: fits the model with Stan through brms only. This is sugar
#' for [flexybayes()]`(..., backend = "brms")` and accepts the same
#' arguments and grammars, an ASReml `fixed` / `random` / `residual`
#' specification or a brms-style bar-grouped formula (see [flexybayes()]
#' for the full argument list). flexyBayes builds the intermediate
#' representation, translates the prior, calls `brms::brm()`, and wraps
#' the result. The live `brmsfit` is available on the `$brms` slot for
#' brms's own posterior tooling (`loo()`, `posterior_predict()`,
#' `bayes_factor()`).
#'
#' flexyBayes has no brms lowering for some ASReml structures: `fa()`, the
#' autoregressive fields, a blocks or low-rank `vm()` / `ped()` carrier, and
#' a `low_rank` smooth approximation. Such a model raises a structured
#' refusal naming the term. Several of these fit on INLA through
#' [fb_inla()]. The capability table in the README lists what each engine
#' fits.
#'
#' @param ... Arguments passed to [flexybayes()] (e.g. `formula` / `fixed`,
#'   `random`, `residual`, `data`, `family`, `prior`, `syntax`). The `backend`
#'   argument is pinned to `"brms"`. A conflicting `backend` value raises a
#'   structured refusal (the redundant `backend = "brms"` is accepted).
#'   The pre-v0.5.0 `formula = ` argument is remapped to the universal
#'   entry's model-spec slot for call-compatibility.
#'
#' @return An object of class `"flexybayes_brms"` (a subclass of
#'   `"flexybayes"`) carrying the live `brmsfit` on `$brms`. See
#'   [flexybayes()] for the shared structure.
#'
#' @section Deprecated:
#' `fb_brms()` is deprecated as of 0.13.0. Use
#' `flexybayes(..., backend = "brms")`.
#'
#' @family flexyBayes engine pins
#' @seealso [flexybayes()] and [fb()] for the universal entry that picks a
#'   backend. [fb_inla()] for the other active engine pin and
#'   [fb_from_brms()] for building a brms-grammar IR.
#' @examples
#' # Held back from the default example run because the brms route compiles
#' # a Stan program before it samples, and the compile alone takes far
#' # longer than an example should. The sampler settings are the smallest
#' # that reach an effective sample size the fit does not warn about, and
#' # the seed makes that reproducible. On a single chain they cost about
#' # half a second next to the compile.
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE) &&
#'     requireNamespace("lme4", quietly = TRUE)) {
#'   data(sleepstudy, package = "lme4")
#'   fit <- flexybayes(
#'     Reaction ~ Days + (1 | Subject), data = sleepstudy,
#'     backend = "brms", chains = 1, n_samples = 4000, warmup = 1000,
#'     seed = 1
#'   )
#'   coef(fit)
#' }
#' }
#' @export
fb_brms <- function(...) {
  .fb_warn_deprecated_export(
    "fb_brms",
    "Use `flexybayes(..., backend = \"brms\")`."
  )
  cl <- match.call()
  # Call-compatibility: the pre-v0.5.0 `formula = ` argument names the
  # same model-spec slot as the universal entry's `fixed`.
  nm <- names(cl)
  if (!is.null(nm) && "formula" %in% nm) {
    if ("fixed" %in% nm) {
      stop(
        "fb_brms(): pass the model once: `formula` and `fixed` name ",
        "the same model-spec slot.",
        call. = FALSE
      )
    }
    nm[nm == "formula"] <- "fixed"
    names(cl) <- nm
  }
  .fb_engine_pin("brms", cl, parent.frame())
}
