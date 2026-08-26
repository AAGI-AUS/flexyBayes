# validate_approximation() --- v0.4.0 exported verb.
#
# The sixth and architecturally-final exported verb. Every approximate
# fit carries, on the registry side, a validation procedure that
# surfaces the realised approximation error against a pass threshold;
# validate_approximation() is the user-facing entry to that procedure.
# Dispatch is on the fit's registered approximation scheme: the fit is
# transiently tagged with its scheme class so a standard S3 UseMethod()
# routes to validate_approximation.<scheme> (e.g. .low_rank_smooth),
# whose body is the per-scheme validation_fn held in
# .approximation_registry. A fit with no recognised approximation falls
# to validate_approximation.default(), which refuses rather than
# returning a vacuous pass.

# --- scheme resolution -------------------------------------------- #

# .fit_approximation_scheme() --- the registered approximation
# scheme(s) carried by `fit`, read from the emit-time slot
# fit$extras$parse_info$approx (a named list keyed by smooth variable,
# each entry tagged with $scheme). Returns the single scheme string
# when the fit carries one scheme, NULL when it carries none, and ---
# defensively, since only one scheme ships at v0.4.0 --- the first
# scheme when a fit somehow mixes schemes.
.fit_approximation_scheme <- function(fit) {
  approx <- fit$extras$parse_info$approx
  if (is.null(approx) || length(approx) == 0L) {
    return(NULL)
  }
  schemes <- unique(vapply(
    approx,
    function(a) a$scheme %||% NA_character_,
    character(1)
  ))
  schemes <- schemes[!is.na(schemes)]
  if (length(schemes) == 0L) {
    return(NULL)
  }
  schemes[[1L]]
}


# --- validation-result object ------------------------------------- #

# .new_approximation_validation() --- constructor for the
# <fb_approximation_validation> result object returned by
# validate_approximation(). Carries the scheme, the per-smooth result
# rows (each with realised Frobenius capture, bias bound, threshold,
# and pass flag), the overall verdict, the threshold, and the
# registry's fallback hint.
.new_approximation_validation <- function(
  scheme,
  per_smooth,
  pass,
  threshold,
  fallback_hint
) {
  structure(
    list(
      scheme = scheme,
      pass = isTRUE(pass),
      threshold = threshold,
      per_smooth = per_smooth,
      fallback_hint = fallback_hint
    ),
    class = c("fb_approximation_validation", "list")
  )
}

#' @export
print.fb_approximation_validation <- function(x, ...) {
  verdict <- if (isTRUE(x$pass)) "PASS" else "FAIL"
  cat("<fb_approximation_validation>\n")
  cat("  scheme:    ", x$scheme, "\n", sep = "")
  cat("  threshold: Frobenius capture >= ", format(x$threshold), "\n", sep = "")
  cat("  verdict:   ", verdict, "\n", sep = "")
  for (r in x$per_smooth) {
    mark <- if (isTRUE(r$pass)) "ok " else "XX "
    cat(
      "    ",
      mark,
      "s(",
      r$smooth,
      "): rank ",
      r$rank,
      "/",
      r$k,
      "  capture ",
      format(round(r$frobenius_capture, 4L)),
      "  (bias bound ",
      format(signif(r$bias_bound, 3L)),
      ")\n",
      sep = ""
    )
  }
  if (!isTRUE(x$pass) && !is.null(x$fallback_hint)) {
    cat("  fallback:  ", x$fallback_hint, "\n", sep = "")
  }
  invisible(x)
}


# --- exported generic + dispatch ---------------------------------- #

#' Validate an approximate model fit against its bias bound
#'
#' `validate_approximation()` reports how much of a fitted model's
#' structure was lost to its approximation scheme, measured against the
#' scheme's declared pass threshold. It is the user-facing entry to the
#' per-scheme validation procedure registered for every approximate
#' route; the contract surfaces the realised error number while the
#' user keeps the accept / re-fit judgement.
#'
#' Dispatch is on the fit's registered approximation scheme. At present
#' the only registered scheme is `low_rank_smooth` (the rank-K
#' principal-component truncation of an `s()` smooth basis): for such a
#' fit, the procedure reports the realised Frobenius capture
#' \eqn{\sum_{i \le K} d_i^2 / \sum_i d_i^2} of each truncated smooth
#' against the default pass threshold of `0.99`, where \eqn{d_i} are
#' the singular values of the full smooth basis.
#'
#' A fit carrying no recognised approximation (an exact fit) is
#' refused rather than returned as a vacuous pass.
#'
#' @param fit A fitted `flexybayes` object.
#' @param ... Passed to the per-scheme validation procedure (e.g.
#'   `threshold` for `low_rank_smooth`).
#'
#' @return An `<fb_approximation_validation>` object: the scheme, the
#'   overall pass / fail verdict, the pass threshold, one result row
#'   per approximated smooth (realised capture, bias bound, per-smooth
#'   pass flag), and the registry's fallback hint.
#'
#' @seealso [flexybayes()] for fitting; the approximation registry
#'   records each scheme's bias bound and fallback.
#' @examples
#' # An exact fit (no registered approximation) refuses rather than
#' # returning a vacuous pass. \donttest{} + an INLA guard: this needs
#' # a live fit and does not fire on a machine without INLA installed.
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#'   fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
#'                      verbose = FALSE)
#'   tryCatch(validate_approximation(fit), error = function(e) {
#'     message(conditionMessage(e)) # "carries no recognised approximation"
#'   })
#' }
#' }
#'
#' # The low_rank_smooth verdict shape, on a plan-only object carrying
#' # the scheme (no active engine reaches low_rank_smooth on the
#' # current backends -- it was built for a since-withdrawn engine and
#' # now refuses at dispatch; see R/emit_smooth_low_rank.R). This part
#' # needs no engine and always runs.
#' demo_fit <- structure(
#'   list(
#'     exactness = "approximate_low_rank_smooth",
#'     extras = list(parse_info = list(approx = list(
#'       x = list(
#'         scheme = "low_rank_smooth", rank = 4L, k = 9L,
#'         frobenius_capture = 0.995,
#'         V_K = matrix(0, 9L, 4L), singular_values = rep(1, 9L)
#'       )
#'     )))
#'   ),
#'   class = c("flexybayes", "list")
#' )
#' v <- validate_approximation(demo_fit)
#' print(v)
#' v$pass
#' @export
validate_approximation <- function(fit, ...) {
  # Dispatch on the fit's *registered approximation scheme* (data on
  # the fit), not on the fit's own class. The scheme is a
  # closed-vocabulary key, so .approximation_registry is itself the
  # dispatch table: each entry carries the per-scheme validation_fn.
  # An absent or unregistered scheme falls to
  # validate_approximation.default(), which refuses. The
  # scheme-specific S3 methods (validate_approximation.low_rank_smooth)
  # wrap the same procedure bodies so direct method calls and the
  # exported-method inventory stay consistent with the registry.
  scheme <- .fit_approximation_scheme(fit)
  if (is.null(scheme)) {
    return(validate_approximation.default(fit, ...))
  }
  entry <- tryCatch(.lookup_approximation(scheme), error = function(e) NULL)
  if (is.null(entry)) {
    return(validate_approximation.default(fit, ...))
  }
  entry$validation_fn(fit, ...)
}

#' @exportS3Method validate_approximation low_rank_smooth
validate_approximation.low_rank_smooth <- function(fit, ...) {
  .validate_low_rank_smooth(fit, ...)
}

#' @exportS3Method validate_approximation default
validate_approximation.default <- function(fit, ...) {
  exactness <- fit$exactness %||% "exact"
  stop(.fb_refusal_condition(
    reason_code = "approximation_scheme_unknown",
    message = paste0(
      "validate_approximation(): this fit carries no recognised ",
      "approximation to validate (exactness = '",
      exactness,
      "'). ",
      "validate_approximation() applies only to approximate fits ",
      "(those routed through a registered approximation scheme such ",
      "as low_rank_smooth); an exact fit has no bias bound to ",
      "report."
    ),
    family_class = "flexybayes_approximation_scheme_unknown"
  ))
}
