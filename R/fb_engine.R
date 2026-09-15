# fb_engine() -- constructor noun for inference engines.
#
# One of the four flexyBayes constructor nouns. `fb_engine(name, opts)`
# names a concrete inference engine and carries its tuning options. It
# is consumed as the `backend` argument of the fitting verbs:
# `flexybayes(..., backend = fb_engine("brms", chains = 4L))`. The bare
# string shorthand (`backend = "brms"`) continues to work and resolves
# to the default engine via `.resolve_engine_string()` in dispatch.R.
#
# `name` and `paradigm` are closed vocabularies at this release: new
# engines extend them only through the design-record amendment process,
# the same posture the four registries take. `auto` is a routing
# directive rather than an engine, so it is not an `fb_engine` name --
# pass the bare string `backend = "auto"` for automatic routing.

# Closed engine vocabulary: each engine maps to its inference paradigm,
# the package that provides it, and (for toolchain status) that
# package's installability.
.fb_engine_paradigm <- c(
  inla = "laplace",
  brms = "mcmc"
)
.fb_engine_package <- c(
  inla = "INLA",
  brms = "brms"
)

# Tuning options an engine may carry. Restricted to the Monte-Carlo
# sampler controls; an unrecognised option is rejected at construction
# rather than silently dropped.
# `adapt_delta` joined the set at 0.11.0. The selection module's reference
# fit defaults to 0.95 and the convergence contract refuses on a single
# divergent transition, which is the contract working -- but without a lever
# a caller who hits one has nothing to reach for but a different seed. It is
# read by the routes that emit through Stan and ignored by the rest.
.fb_engine_opt_names <- c("n_samples", "warmup", "chains", "adapt_delta")


#' Construct an inference-engine specification
#'
#' Names a concrete inference engine and its tuning options. The result
#' is passed as the `backend` argument of the fitting verbs:
#' `flexybayes(..., backend = fb_engine("brms", chains = 4L))`. The bare
#' string form (`backend = "brms"`) remains valid and is equivalent to
#' the default `fb_engine()` for that engine.
#'
#' `name` and the derived `paradigm` are closed vocabularies. `"auto"`
#' is a routing directive, not an engine; pass `backend = "auto"`
#' directly for automatic routing.
#'
#' @section Recognised options:
#'
#' All four are sampler controls, read by the routes that emit through
#' Stan and ignored by the rest. An unrecognised name is an error at
#' construction rather than a silently dropped argument.
#'
#' \describe{
#'   \item{`n_samples`}{Post-warmup draws to keep, per chain.}
#'   \item{`warmup`}{Warmup iterations, per chain.}
#'   \item{`chains`}{Number of chains.}
#'   \item{`adapt_delta`}{Target acceptance probability for Stan's
#'   step-size adaptation, between 0 and 1. Raising it makes the sampler
#'   take smaller steps, which costs time and removes the divergent
#'   transitions a smaller step would have avoided. brms's own default is
#'   0.8; [fb_select()]'s reference fit starts from 0.95 instead, and a
#'   value given here replaces that.}
#' }
#'
#' `adapt_delta` is the lever to reach for when a selection run refuses
#' with `penalty_sampler_diagnostics_failed`. The adequacy gate refuses on
#' a single divergent transition, because a selection read off a posterior
#' the sampler did not explore is a statement about the sampler; a
#' heavy-tailed prior on a small trial reaches that threshold routinely.
#' A regularised horseshoe on a 200-row design ended two transitions of
#' six thousand in a divergence at the 0.95 default and completed cleanly
#' at `fb_engine("brms", adapt_delta = 0.999)`.
#'
#' @param name Character(1): the engine, one of `"inla"`, `"brms"`.
#' @param opts Named list of tuning options. Recognised names are
#'   `n_samples`, `warmup`, `chains` and `adapt_delta`, described under
#'   "Recognised options" below; an unrecognised name is an error.
#' @param ... Tuning options given individually, merged into `opts`
#'   (e.g. `fb_engine("brms", chains = 4L)`).
#' @return An `fb_engine` object: a classed list with elements `name`,
#'   `paradigm` (one of `mcmc`, `laplace`, `vb`, `map`),
#'   `toolchain_status` (one of `ready`, `requires_install`,
#'   `unavailable`), and `opts`.
#' @seealso [fb_approx()]
#' @examples
#' e <- fb_engine("brms", chains = 4L)
#' e$paradigm
#' e$toolchain_status
#' @export
fb_engine <- function(name, opts = list(), ...) {
  .check_string_scalar(
    name,
    "fb_engine(): `name` must be a non-empty single string."
  )
  if (!name %in% names(.fb_engine_paradigm)) {
    stop(
      "fb_engine(): unknown engine '",
      name,
      "'. Supported engines: ",
      paste(names(.fb_engine_paradigm), collapse = ", "),
      ". ('auto' is a routing directive -- pass backend = \"auto\" ",
      "directly.)",
      call. = FALSE
    )
  }

  if (!is.list(opts)) {
    stop("fb_engine(): `opts` must be a named list.", call. = FALSE)
  }
  opts <- utils::modifyList(opts, list(...))
  if (
    length(opts) &&
      (is.null(names(opts)) || any(!nzchar(names(opts))))
  ) {
    stop(
      "fb_engine(): every option in `opts` / `...` must be named.",
      call. = FALSE
    )
  }
  unknown <- setdiff(names(opts), .fb_engine_opt_names)
  if (length(unknown)) {
    stop(
      "fb_engine(): unrecognised option",
      if (length(unknown) != 1L) "s" else "",
      " ",
      paste(shQuote(unknown), collapse = ", "),
      ". Recognised options: ",
      paste(.fb_engine_opt_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  pkg <- .fb_engine_package[[name]]
  status <- if (requireNamespace(pkg, quietly = TRUE)) {
    "ready"
  } else {
    "requires_install"
  }

  structure(
    list(
      name = name,
      paradigm = unname(.fb_engine_paradigm[[name]]),
      toolchain_status = status,
      opts = opts
    ),
    class = c("fb_engine", "list")
  )
}

#' Test whether an object is an `fb_engine` specification
#'
#' @param x Any R object. The test is a class check, not a structural
#'   one.
#' @returns `TRUE` if `x` is an `fb_engine` object, `FALSE` otherwise.
#' @section Deprecated:
#' `is_fb_engine()` is deprecated in flexyBayes 0.13.0 and will be removed in
#' a later release. Use `inherits(x, "fb_engine")`. A call warns with a
#' condition of class `flexybayes_deprecated`.
#' @examples
#' inherits(fb_engine("inla"), "fb_engine")
#' inherits("inla", "fb_engine")
#' @export
is_fb_engine <- function(x) {
  .fb_warn_deprecated_export(
    "is_fb_engine",
    "Use `inherits(x, \"fb_engine\")`."
  )
  inherits(x, "fb_engine")
}

#' Print an `fb_engine` specification
#'
#' @param x An `fb_engine` object as returned by [fb_engine()].
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the one-line engine
#'   description it prints.
#' @export
print.fb_engine <- function(x, ...) {
  cat(
    "<fb_engine> ",
    x$name,
    " (",
    x$paradigm,
    ", ",
    x$toolchain_status,
    ")\n",
    sep = ""
  )
  if (length(x$opts)) {
    parts <- vapply(
      seq_along(x$opts),
      function(i) {
        paste0(names(x$opts)[i], " = ", format(x$opts[[i]]))
      },
      character(1L)
    )
    cat("  opts: ", paste(parts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
