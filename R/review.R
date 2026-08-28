# <flexybayes_review> -- deferred-execution token returned when a
# user passes `review_code = TRUE` to `flexybayes()` (and, when
# that entry lands, `fb_brms()`). The object
# carries the generated backend code together with the IR
# (intermediate representation), resolved prior, captured call,
# RNG snapshot, and the arguments needed to advance the deferred
# fit -- so the user can inspect the code first, then `proceed()`
# into a fit without re-issuing the original call.
#
# Storage is `environment` so that the proceed-cache contract
# (`proceed(rev)` a second time returns the cached fit) holds
# without forcing the caller to write `rev <- proceed(rev)`. A
# list-based class would reset to `proceeded = FALSE` on every
# call because R's copy-on-modify semantics would not propagate
# back to the caller's frame.
#
# Internal note: closest published precedent for the inspect-
# then-fit pattern is brms's `make_stancode()` / `make_standata()`
# pair plus the `chains = 0` idiom. The `<flexybayes_review>`
# surface subsumes both: `cat_code(rev)` reproduces
# `make_stancode()`-style display, and `proceed(rev)` runs the
# deferred fit with the captured RNG state.

# ---------------------------------------------------------------- #
# Constructor (not exported)                                       #
# ---------------------------------------------------------------- #

.new_flexybayes_review <- function(
  code,
  backend,
  ir,
  prior,
  data_name,
  call,
  seed,
  proceed_args,
  preflight = NULL
) {
  if (!is.character(code) || length(code) != 1L) {
    stop("`code` must be a length-1 character string.", call. = FALSE)
  }
  if (
    !is.character(backend) ||
      length(backend) != 1L ||
      !backend %in% c("stan_via_brms", "inla")
  ) {
    stop(
      "`backend` must be one of \"stan_via_brms\", ",
      "\"inla\".",
      call. = FALSE
    )
  }
  if (!inherits(ir, "fb_terms")) {
    stop("`ir` must be an `fb_terms` object.", call. = FALSE)
  }
  if (!is.list(proceed_args)) {
    stop("`proceed_args` must be a list of backend arguments.", call. = FALSE)
  }
  if (!is.null(preflight) && !inherits(preflight, "fb_preflight")) {
    stop(
      "`preflight` must be NULL or an `<fb_preflight>` object.",
      call. = FALSE
    )
  }

  env <- new.env(parent = emptyenv())
  env$code <- code
  env$backend <- backend
  env$ir <- ir
  env$prior <- prior
  env$data_name <- data_name
  env$call <- call
  env$seed <- seed
  env$proceed_args <- proceed_args
  env$proceeded <- FALSE
  env$fit <- NULL
  env$preflight <- preflight # NULL below 1e5 rows
  class(env) <- c("flexybayes_review", "environment")
  env
}


# ---------------------------------------------------------------- #
# print method                                                     #
# ---------------------------------------------------------------- #

#' Print method for a deferred-execution review object
#'
#' Two-line summary of a `<flexybayes_review>` object returned when
#' `review_code = TRUE` is passed to [flexybayes()]. The first line
#' reports the backend and the IR (intermediate representation)
#' dimensions; the second is a prompt directing the user to
#' [cat_code()] for inspection and [proceed()] to advance the fit.
#'
#' @param x   A `<flexybayes_review>` object as returned by
#'   `flexybayes(review_code = TRUE)`.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the two summary lines it
#'   prints.
#' @keywords internal
#' @export
print.flexybayes_review <- function(x, ...) {
  n_obs <- x$ir$data_summary$n
  n_fix <- length(x$ir$fixed_terms) + as.integer(isTRUE(x$ir$intercept))
  n_rand <- length(x$ir$random_terms)
  cat(sprintf(
    "<flexybayes_review> backend=%s  n=%d  fixed=%d  random=%d\n",
    x$backend,
    n_obs,
    n_fix,
    n_rand
  ))

  # Show the preflight summary above the cat_code()
  # prompt when the IR carries n_rows >= 1e5 (the threshold the
  # dispatch / review-entry preflight call gates on).
  if (inherits(x$preflight, "fb_preflight")) {
    cat("\nPreflight summary:\n")
    fmt_bytes <- function(b) {
      format(round(b), big.mark = " ", scientific = FALSE)
    }
    for (nm in names(x$preflight$per_term_estimate)) {
      e <- x$preflight$per_term_estimate[[nm]]
      cat(sprintf(
        "  %-28s  %14s B  %-26s\n",
        nm,
        fmt_bytes(e$design_memory_bytes),
        e$representation_class
      ))
    }
    cat(sprintf(
      "  total = %s B   ceiling = %s B   (accepted)\n",
      fmt_bytes(x$preflight$total_estimate_bytes),
      fmt_bytes(x$preflight$ceiling_bytes)
    ))
    cat("\n")
  }

  if (isTRUE(x$proceeded)) {
    cat("--- already proceeded; see x$fit\n")
  } else {
    cat(
      "Run cat_code(x) to view the generated code; ",
      "proceed(x) to fit.\n",
      sep = ""
    )
  }
  invisible(x)
}


# ---------------------------------------------------------------- #
# cat_code -- generic + method                                     #
# ---------------------------------------------------------------- #

#' Emit the generated backend code for a deferred review object
#'
#' Writes the backend code carried by a `<flexybayes_review>`
#' object (Stan code via [brms::stancode()] for
#' [flexybayes()] / `fb_brms()`) to a connection.
#'
#' @param x   A `<flexybayes_review>` object carrying the generated
#'   backend code.
#' @param ... Method-specific arguments. The `flexybayes_review` method
#'   accepts `file`, a connection defaulting to `stdout()`.
#' @returns Invisibly, the code as a character vector. Called for the
#'   code it writes to `file`.
#' @examples
#' # review_code = TRUE builds the engine code without fitting, so this
#' # runs without INLA or brms installed.
#' set.seed(1)
#' d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#' rv <- flexybayes(y ~ x + (1 | g), data = d, backend = "brms",
#'                  review_code = TRUE)
#' cat_code(rv)
#' @export
cat_code <- function(x, ...) {
  UseMethod("cat_code")
}

#' @rdname cat_code
#' @param file A connection to write the code to, defaulting to
#'   `stdout()`.
#' @export
cat_code.flexybayes_review <- function(x, file = stdout(), ...) {
  writeLines(x$code, con = file)
  invisible(x$code)
}


# ---------------------------------------------------------------- #
# proceed -- generic + method                                      #
# ---------------------------------------------------------------- #

#' Advance a deferred-execution object into its fit
#'
#' Generic for the inspect-then-fit pattern. The
#' `<flexybayes_review>` method restores the RNG snapshot
#' captured at review-object construction, runs the deferred fit
#' via the backend driver, caches the result in-place, and
#' returns the fit. A second call returns the cached fit.
#'
#' @param x   A `<flexybayes_review>` object holding the deferred fit.
#' @param ... Reserved for future deferred-execution classes, such as a
#'   deferred triangulation.
#' @returns The fit object the originating call would have returned, of
#'   class `flexybayes`. A second call returns the cached fit.
#' @examples
#' # The review object holds the model back until proceed() is called.
#' set.seed(1)
#' d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#' rv <- flexybayes(y ~ x + (1 | g), data = d, backend = "brms",
#'                  review_code = TRUE)
#' rv
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit <- proceed(rv)
#'   class(fit)
#' }
#' }
#' @export
proceed <- function(x, ...) {
  UseMethod("proceed")
}

#' @rdname proceed
#' @export
proceed.flexybayes_review <- function(x, ...) {
  if (isTRUE(x$proceeded)) {
    if (!isTRUE(getOption("flexyBayes.silence_review_cached_note", FALSE))) {
      message(
        "flexyBayes: returning cached fit from prior ",
        "proceed(). Silence this note via ",
        "options(flexyBayes.silence_review_cached_note ",
        "= TRUE)."
      )
    }
    return(x$fit)
  }

  # Restore the .Random.seed snapshot captured at review-object
  # construction so that proceed(rev) reproduces the chain that the
  # equivalent direct call (review_code = FALSE) would have produced
  # at the same outer seed.
  if (!is.null(x$seed)) {
    assign(".Random.seed", x$seed, envir = globalenv())
  }

  args <- x$proceed_args
  emit_fn <- switch(
    x$backend,
    "stan_via_brms" = emit_brms,
    stop(
      "proceed(): unrecognised backend \"",
      x$backend,
      "\" on review object.",
      call. = FALSE
    )
  )
  fit <- do.call(emit_fn, args)

  x$fit <- fit
  x$proceeded <- TRUE
  .fb_warn_poor_convergence(fit)
  fit
}
