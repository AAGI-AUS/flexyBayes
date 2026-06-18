# fb_approx() -- constructor noun for approximation schemes.
#
# One of the four flexyBayes constructor nouns. `fb_approx(scheme, ...)`
# names an approximation scheme from the locked approximation registry
# and carries its scheme-specific tuning (e.g. `rank` for
# `low_rank_smooth`). It is consumed wherever an approximate
# representation is accepted, most directly as
# `s(x, representation = fb_approx("low_rank_smooth", rank = 5L))`.
#
# The returned object is a classed list whose `scheme` and tuning
# kwargs are plain list elements, so the consuming parse path reads
# `spec$scheme` / `spec$rank` directly; the bias-bound description is an
# attribute surfaced by the print method.

#' Construct an approximation-scheme specification
#'
#' Names an approximation scheme and its tuning, validated against the
#' package's locked approximation registry. The result is consumed where
#' an approximate representation is accepted -- most directly in a
#' smooth term, `s(x, representation = fb_approx("low_rank_smooth",
#' rank = 5L))`.
#'
#' The only scheme with a smooth-basis fitting path in this release is
#' `"low_rank_smooth"`, a rank-`K` principal-component truncation of the
#' smooth basis. Its bias is the relative squared Frobenius error of the
#' truncation (Wood, 2017, chapter 5); [validate_approximation()]
#' reports the realised capture against the threshold for a fitted
#' model.
#'
#' @param scheme Character(1): the approximation-scheme name. Must be a
#'   scheme registered in the approximation registry; an unregistered
#'   name is refused with the supported vocabulary.
#' @param ... Scheme-specific tuning carried on the object, for example
#'   `rank` for `"low_rank_smooth"`.
#' @return An `fb_approx` object: a classed list with the `scheme` and
#'   the tuning kwargs as elements, and a `bias_bound_promise` attribute
#'   describing the scheme's bias bound.
#' @seealso [validate_approximation()], [fb_engine()]
#' @examples
#' a <- fb_approx("low_rank_smooth", rank = 5L)
#' a$scheme
#' a$rank
#' @export
fb_approx <- function(scheme, ...) {
  if (!is.character(scheme) || length(scheme) != 1L || !nzchar(scheme)) {
    stop(
      "fb_approx(): `scheme` must be a non-empty single string.",
      call. = FALSE
    )
  }

  # Validate against the locked registry; an unknown scheme is a
  # catchable structured refusal naming the supported vocabulary.
  entry <- tryCatch(.lookup_approximation(scheme), error = function(e) NULL)
  if (is.null(entry)) {
    known <- sort(ls(envir = .approximation_registry, all.names = FALSE))
    stop(.fb_refusal_condition(
      reason_code = "approximation_scheme_unknown",
      message = paste0(
        "fb_approx(): '",
        scheme,
        "' is not a registered approximation ",
        "scheme. Supported scheme",
        if (length(known) != 1L) "s" else "",
        ": ",
        if (length(known)) {
          paste(known, collapse = ", ")
        } else {
          "(none registered yet)"
        },
        "."
      ),
      family_class = "flexybayes_approximation_scheme_unknown"
    ))
  }

  kwargs <- list(...)
  out <- c(list(scheme = scheme), kwargs)
  structure(
    out,
    class = c("fb_approx", "list"),
    bias_bound_promise = .fb_bias_bound_promise(entry)
  )
}

#' Test whether an object is an `fb_approx` specification
#'
#' @param x An object.
#' @return `TRUE` if `x` is an `fb_approx` object.
#' @export
is_fb_approx <- function(x) inherits(x, "fb_approx")

# Human-readable, literature-referenced bias-bound description built
# from a registry entry. Provenance metadata (the internal design
# record a scheme traces to) is deliberately not surfaced here.
.fb_bias_bound_promise <- function(entry) {
  bb <- entry$bias_bound
  if (!is.list(bb) || is.null(bb$interpretation)) {
    return(NA_character_)
  }
  ref <- switch(
    bb$reference %||% "",
    wood_2017_chapter_5 = "Wood (2017), chapter 5",
    bb$reference %||% ""
  )
  if (nzchar(ref)) {
    paste0(bb$interpretation, " (", ref, ")")
  } else {
    bb$interpretation
  }
}

#' Print an `fb_approx` specification
#'
#' @param x An `fb_approx` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.fb_approx <- function(x, ...) {
  kw <- x[setdiff(names(x), "scheme")]
  cat("<fb_approx> scheme = \"", x$scheme, "\"", sep = "")
  if (length(kw)) {
    parts <- vapply(
      seq_along(kw),
      function(i) {
        paste0(names(kw)[i], " = ", format(kw[[i]]))
      },
      character(1L)
    )
    cat(" | ", paste(parts, collapse = ", "), sep = "")
  }
  cat("\n")
  bb <- attr(x, "bias_bound_promise")
  if (!is.null(bb) && !is.na(bb)) {
    cat("  bias bound: ", bb, "\n", sep = "")
  }
  invisible(x)
}
