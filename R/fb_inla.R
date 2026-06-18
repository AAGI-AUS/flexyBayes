# fb_inla -- INLA engine pin.
#
# An engine pin is sugar over the universal entry with `backend` fixed:
# fb_inla(...) == flexybayes(..., backend = "inla"). It fits via INLA
# only -- the model must be latent-Gaussian feasible, else the universal
# entry's lgm_gate() raises the usual structured refusal.
#
# The pin forwards the user's call VERBATIM (a match.call() rewrite)
# rather than re-listing flexybayes()'s ~20 arguments. This is what makes
# it true sugar: arguments the user did not supply stay unsupplied, so
# flexybayes()'s missing()-based defaults (the default-prior note) and its
# deparse(substitute(data)) data-name capture behave exactly as on a
# direct call. fb_greta() and fb_brms() are recast onto the same
# .fb_engine_pin() seam at v0.5.0.

# .fb_engine_pin() -- rewrite a pin verb's own call into a flexybayes()
# call with `backend` pinned, and evaluate it in the caller's frame so
# the data argument and the default-detection semantics resolve against
# the user's environment. Namespaced target (flexyBayes::flexybayes) so
# it resolves whether or not the package is attached.
#
# Backend-argument discipline. A pin fixes the
# engine, so an explicit `backend` is a contradiction --- EXCEPT the
# redundant self-pin (e.g. fb_brms(backend = "brms")), which the v0.4.1
# deprecation notice promised would remain forward-compatible. A
# conflicting `backend` raises the structured engine_pin_backend_conflict
# refusal; a matching one is accepted and overwritten with the same value.
.fb_engine_pin <- function(engine, call, env) {
  if ("backend" %in% names(call)) {
    supplied <- tryCatch(
      as.character(eval(call$backend, env)),
      error = function(e) as.character(call$backend)
    )
    if (!identical(supplied, engine)) {
      stop(.fb_refusal_condition(
        reason_code = "engine_pin_backend_conflict",
        message = paste0(
          "fb_",
          engine,
          "() is an engine pin: it fits via the ",
          engine,
          " engine only, so it cannot take backend = \"",
          paste(supplied, collapse = "\", \""),
          "\". Drop `backend`, or ",
          "use fb() / flexybayes() to choose a backend."
        ),
        backend = supplied
      ))
    }
  }
  call[[1L]] <- quote(flexyBayes::flexybayes)
  call$backend <- engine
  eval(call, env)
}

#' Fit a flexyBayes model via the INLA engine
#'
#' Engine pin: fits the model with INLA (integrated nested Laplace
#' approximation) only. This is sugar for
#' [flexybayes()]`(..., backend = "inla")` and accepts the same arguments
#' and grammars (an ASReml `fixed` / `random` / `rcov` specification or a
#' brms-style bar-grouped formula -- see [flexybayes()] for the full
#' argument list). The model must be latent-Gaussian feasible; if it is
#' not, the shared `lgm_gate()` raises a structured refusal naming the
#' offending term, exactly as `flexybayes(backend = "inla")` does.
#'
#' Sampling-control arguments (`n_samples`, `warmup`, `chains`,
#' `mcmc_verbose`) are accepted for call-compatibility with the other
#' engine pins but are inert under INLA's deterministic Laplace
#' approximation.
#'
#' @param ... Arguments passed to [flexybayes()] (e.g. `fixed`, `random`,
#'   `rcov`, `data`, `family`, `prior`, `syntax`). The `backend` argument
#'   is pinned to `"inla"` and must not be supplied.
#'
#' @return An object of class `"flexybayes"` (specifically a
#'   `flexybayes_inla` fit); see [flexybayes()] for the structure.
#'
#' @family flexyBayes engine pins
#' @seealso [flexybayes()] and [fb()] for the universal entry that picks a
#'   backend; [fb_from_asreml()] / [fb_from_brms()] for the ingest layer.
#' @examples
#' df <- data.frame(
#'   yield = rnorm(40),
#'   geno  = factor(rep(letters[1:8], 5)),
#'   env   = factor(rep(c("a", "b"), 20))
#' )
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   fit <- fb_inla(yield ~ env, random = ~ geno, data = df)
#' }
#' }
#' @export
fb_inla <- function(...) {
  .fb_engine_pin("inla", match.call(), parent.frame())
}
