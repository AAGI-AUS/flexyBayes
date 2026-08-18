# validators.R -- Shared input checks used across the flexyBayes surface.
#
# Each helper carries one guard that the entry points had repeated inline,
# extracted verbatim. The caller supplies the message, so a refusal reads
# exactly as it read before the extraction: these helpers move the test,
# never the wording. Structured refusals stay where they are raised --
# only the plain-guard duplication lives here.

#' Require an optional package before using it
#'
#' Stops with the caller's message when `pkg` is not installed. This is
#' the guard every backend, ecosystem shim, and file-format path
#' repeated inline, each with its own installation advice.
#'
#' @param pkg A length-one character string naming the package to test
#'   for, passed to [requireNamespace()].
#' @param ... The message pieces handed to [stop()] when the package is
#'   missing, written exactly as the call site would have written them.
#' @returns Invisibly, `TRUE` when the package is available. Stops
#'   otherwise, without a call context.
#'
#' @noRd
#' @keywords internal
.check_installed <- function(pkg, ...) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(..., call. = FALSE)
  }
  invisible(TRUE)
}


#' Require an `fb_terms` intermediate representation
#'
#' The type guard every gate and every emit repeated on its first
#' argument. Each caller names itself and its own re-route in the
#' message, so the wording travels with the call rather than with the
#' test.
#'
#' @param fb The object that should be the parsed IR.
#' @param ... The message pieces handed to [stop()] when it is not,
#'   written exactly as the call site would have written them.
#' @returns Invisibly, `TRUE` when `fb` is an `fb_terms` object. Stops
#'   otherwise.
#'
#' @noRd
#' @keywords internal
.check_fb_terms <- function(fb, ...) {
  if (!is_fb_terms(fb)) {
    stop(..., call. = FALSE)
  }
  invisible(TRUE)
}


#' Require an `fb_dataset` wrapper
#'
#' The companion guard to `.check_fb_terms()` on the preflight and
#' aggregation entry points, which take the IR and the dataset wrapper
#' as a pair.
#'
#' @param x The object that should be the dataset wrapper.
#' @param ... The message pieces handed to [stop()] when it is not.
#' @returns Invisibly, `TRUE` when `x` is an `fb_dataset` object. Stops
#'   otherwise.
#'
#' @noRd
#' @keywords internal
.check_fb_dataset <- function(x, ...) {
  if (!inherits(x, "fb_dataset")) {
    stop(..., call. = FALSE)
  }
  invisible(TRUE)
}


#' Require an `fb_aggregated` sufficient-statistics carrier
#'
#' The guard both aggregated emits repeat on the carrier built by the
#' aggregation layer.
#'
#' @param x The object that should be the aggregation carrier.
#' @param ... The message pieces handed to [stop()] when it is not.
#' @returns Invisibly, `TRUE` when `x` is an `fb_aggregated` object.
#'   Stops otherwise.
#'
#' @noRd
#' @keywords internal
.check_fb_aggregated <- function(x, ...) {
  if (!inherits(x, "fb_aggregated")) {
    stop(..., call. = FALSE)
  }
  invisible(TRUE)
}


#' Require a non-empty single string
#'
#' The scalar-character guard the registries and the specification
#' constructors repeat on a name or a scheme token.
#'
#' @param x The object that should be a length-one, non-empty character
#'   vector.
#' @param ... The message pieces handed to [stop()] when it is not.
#' @returns Invisibly, `TRUE` when `x` is a non-empty single string.
#'   Stops otherwise.
#'
#' @noRd
#' @keywords internal
.check_string_scalar <- function(x, ...) {
  if (!is.character(x) || length(x) != 1L || !nzchar(x)) {
    stop(..., call. = FALSE)
  }
  invisible(TRUE)
}


#' Require an `fb_prior` object
#'
#' The guard the prior compilers repeat before reading `$specs`.
#'
#' @param x The object that should be the prior specification.
#' @param ... The message pieces handed to [stop()] when it is not.
#' @returns Invisibly, `TRUE` when `x` is an `fb_prior` object. Stops
#'   otherwise.
#'
#' @noRd
#' @keywords internal
.check_fb_prior <- function(x, ...) {
  if (!inherits(x, "fb_prior")) {
    stop(..., call. = FALSE)
  }
  invisible(TRUE)
}


#' Require a fitted `flexybayes` object
#'
#' The guard the post-fit reporting functions repeat before reaching
#' into a fit's slots.
#'
#' @param fit The object that should be a fitted model.
#' @param ... The message pieces handed to [stop()] when it is not.
#' @returns Invisibly, `TRUE` when `fit` inherits from `flexybayes`.
#'   Stops otherwise.
#'
#' @noRd
#' @keywords internal
.check_flexybayes_fit <- function(fit, ...) {
  if (!inherits(fit, "flexybayes")) {
    stop(..., call. = FALSE)
  }
  invisible(TRUE)
}


#' Is the posterior package available?
#'
#' A named predicate rather than an inline [requireNamespace()] call, so
#' the guard that consults it can be reached from a test without
#' uninstalling the package. Same shape as `.fb_emmeans_available()` on
#' the classify path.
#'
#' @returns `TRUE` when \pkg{posterior} can be loaded, `FALSE` otherwise.
#'
#' @noRd
#' @keywords internal
.fb_posterior_available <- function() {
  requireNamespace("posterior", quietly = TRUE)
}


#' Is the brms package available?
#'
#' The companion predicate for the two door methods that delegate to
#' \pkg{brms} on a brms-engine fit: `loo()` and `pp_check()`. A fit
#' carrying a `brmsfit` normally travels with the package that built it,
#' but a saved fit read back on a machine without \pkg{brms} reaches the
#' delegation with nothing to delegate to.
#'
#' @returns `TRUE` when \pkg{brms} can be loaded, `FALSE` otherwise.
#'
#' @noRd
#' @keywords internal
.fb_brms_available <- function() {
  requireNamespace("brms", quietly = TRUE)
}


#' Is the bayesplot package available?
#'
#' The third of the door's package predicates. `pp_check()` builds its
#' display with \pkg{bayesplot}, and the guard is reached through
#' `plot(fit, type = "pp_check")`, which does not load the generic's home
#' package on the way in.
#'
#' @returns `TRUE` when \pkg{bayesplot} can be loaded, `FALSE` otherwise.
#'
#' @noRd
#' @keywords internal
.fb_bayesplot_available <- function() {
  requireNamespace("bayesplot", quietly = TRUE)
}
