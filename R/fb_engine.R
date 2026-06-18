# fb_engine() -- constructor noun for inference engines.
#
# One of the four flexyBayes constructor nouns. `fb_engine(name, opts)`
# names a concrete inference engine and carries its tuning options. It
# is consumed as the `backend` argument of the fitting verbs:
# `flexybayes(..., backend = fb_engine("greta", chains = 4L))`. The bare
# string shorthand (`backend = "greta"`) continues to work and resolves
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
  greta = "mcmc",
  inla = "laplace",
  brms = "mcmc"
)
.fb_engine_package <- c(
  greta = "greta",
  inla = "INLA",
  brms = "brms"
)

# Tuning options an engine may carry. Restricted to the Monte-Carlo
# sampler controls; an unrecognised option is rejected at construction
# rather than silently dropped.
.fb_engine_opt_names <- c("n_samples", "warmup", "chains")


#' Construct an inference-engine specification
#'
#' Names a concrete inference engine and its tuning options. The result
#' is passed as the `backend` argument of the fitting verbs:
#' `flexybayes(..., backend = fb_engine("greta", chains = 4L))`. The bare
#' string form (`backend = "greta"`) remains valid and is equivalent to
#' the default `fb_engine()` for that engine.
#'
#' `name` and the derived `paradigm` are closed vocabularies. `"auto"`
#' is a routing directive, not an engine; pass `backend = "auto"`
#' directly for automatic routing.
#'
#' @param name Character(1): the engine, one of `"greta"`, `"inla"`,
#'   `"brms"`.
#' @param opts Named list of tuning options. Recognised names are
#'   `n_samples`, `warmup`, `chains`; an unrecognised name is an error.
#' @param ... Tuning options given individually, merged into `opts`
#'   (e.g. `fb_engine("greta", chains = 4L)`).
#' @return An `fb_engine` object: a classed list with elements `name`,
#'   `paradigm` (one of `mcmc`, `laplace`, `vb`, `map`),
#'   `toolchain_status` (one of `ready`, `requires_install`,
#'   `unavailable`), and `opts`.
#' @seealso [fb_approx()]
#' @examples
#' e <- fb_engine("greta", chains = 4L)
#' e$paradigm
#' e$toolchain_status
#' @export
fb_engine <- function(name, opts = list(), ...) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop(
      "fb_engine(): `name` must be a non-empty single string.",
      call. = FALSE
    )
  }
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
#' @param x An object.
#' @return `TRUE` if `x` is an `fb_engine` object.
#' @export
is_fb_engine <- function(x) inherits(x, "fb_engine")

#' Print an `fb_engine` specification
#'
#' @param x An `fb_engine` object.
#' @param ... Unused.
#' @return `x`, invisibly.
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
