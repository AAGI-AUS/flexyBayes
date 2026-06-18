# fb_terms -- flexyBayes intermediate representation (IR)
#
# An S3 class that bridges asreml-format ingest (via parse_formula.R)
# and brms-formula ingest (via brms::brmsterms()) to
# the three v0.1 emit backends: greta (existing), brms / Stan (new),
# and INLA -- Integrated Nested Laplace Approximations (new).
#
# Internal for v0.1. Constructed inside `flexybayes()` (the asreml
# entry, retained) or `fb()` (the brms entry) -- not
# for direct user instantiation. The print and format S3 methods
# below are exported solely for dispatch under
# `devtools::load_all()` and `R CMD check --as-cran`.
#
# Scope of this file: class skeleton, constructor,
# validator, predicate, accessors, print, format. No ingest, no
# backend logic, no semantic interpretation of term descriptors.
# The ingest layers wire this into the existing asreml -> greta
# pathway without altering codegen.R.

# ---------------------------------------------------------------- #
# Constructor                                                      #
# ---------------------------------------------------------------- #

# Build an fb_terms object from validated fields.
#
# All inputs are expected to be pre-shaped by the caller
# (`flexybayes()` or `fb()`); the constructor performs structural
# validation only -- types, lengths, allowed values. Term-level
# semantic validation is the responsibility of `lgm_gate()`.
#
# @param response  length-1 non-empty character: response variable name.
# @param family    a base-R `family` object, or a length-1 character
#                  family name (e.g., "gaussian", "binomial").
# @param link      length-1 character or NULL -- link override.
# @param intercept logical(1).
# @param fixed_terms,random_terms,rcov_terms,addition_terms
#                  list of term descriptors. Each element must be a
#                  list with a non-empty character `type` field.
#                  parse_formula.R types (e.g., "factor", "simple",
#                  "fa_gxe", "ar1_spatial") are accepted; brms-derived
#                  types ("gr", "mm", "gp", "s", "cs") are added
#                  by the brms ingest layer.
# @param priors    an `fb_prior` object or NULL --
#                  defaults applied downstream by the emit backend.
# @param data_summary list with at least `n` (sample size) and
#                  factor-level metadata. Populated by the ingest
#                  layer.
# @param capabilities character vector of capability flags
#                  populated by `lgm_gate()`. Empty
#                  at construction.
# @param source    "asreml" or "brms" -- provenance, used by
#                  `triangulate()` for clear ingest-path labelling
#                  in cross-engine reports.
# @return an fb_terms object (S3, inherits from list).
new_fb_terms <- function(
  response,
  family,
  link = NULL,
  intercept = TRUE,
  fixed_terms = list(),
  random_terms = list(),
  rcov_terms = list(),
  addition_terms = list(),
  priors = NULL,
  data_summary = list(),
  capabilities = character(),
  source = c("asreml", "brms", "greta"),
  greta_meta = NULL
) {
  source <- match.arg(source)

  obj <- list(
    response = response,
    family = family,
    link = link,
    intercept = intercept,
    fixed_terms = fixed_terms,
    random_terms = random_terms,
    rcov_terms = rcov_terms,
    addition_terms = addition_terms,
    priors = priors,
    data_summary = data_summary,
    capabilities = capabilities,
    source = source,
    greta_meta = greta_meta
  )
  class(obj) <- c("fb_terms", "list")
  validate_fb_terms(obj)
  obj
}

# ---------------------------------------------------------------- #
# Validator                                                        #
# ---------------------------------------------------------------- #

# Invariant check. Called by the constructor and re-callable on a
# constructed object after mutation (e.g., after lgm_gate populates
# `capabilities`).
validate_fb_terms <- function(x) {
  stopifnot(inherits(x, "fb_terms"))

  if (
    !is.character(x$response) ||
      length(x$response) != 1L ||
      is.na(x$response) ||
      !nzchar(x$response)
  ) {
    stop(
      "`response` must be a non-empty length-1 character vector.",
      call. = FALSE
    )
  }

  if (
    !inherits(x$family, "family") &&
      !(is.character(x$family) && length(x$family) == 1L)
  ) {
    stop(
      "`family` must be a `family` object or a length-1 ",
      "character family name.",
      call. = FALSE
    )
  }

  if (
    !is.null(x$link) &&
      !(is.character(x$link) && length(x$link) == 1L)
  ) {
    stop("`link` must be NULL or a length-1 character vector.", call. = FALSE)
  }

  if (!is.logical(x$intercept) || length(x$intercept) != 1L) {
    stop("`intercept` must be TRUE or FALSE.", call. = FALSE)
  }
  # NA permitted only when source = "greta": greta-direct
  # models do not declare an intercept syntactically. Otherwise the
  # original "TRUE or FALSE" contract holds.
  if (is.na(x$intercept) && !identical(x$source, "greta")) {
    stop(
      "`intercept` must be TRUE or FALSE (NA only permitted on ",
      "the greta-direct entry, where the IR carries ",
      "`source = \"greta\"`).",
      call. = FALSE
    )
  }

  for (slot in c(
    "fixed_terms",
    "random_terms",
    "rcov_terms",
    "addition_terms"
  )) {
    if (!is.list(x[[slot]])) {
      stop("`", slot, "` must be a list.", call. = FALSE)
    }
    for (i in seq_along(x[[slot]])) {
      el <- x[[slot]][[i]]
      if (
        !is.list(el) ||
          is.null(el$type) ||
          !is.character(el$type) ||
          length(el$type) != 1L ||
          !nzchar(el$type)
      ) {
        stop(
          "`",
          slot,
          "[[",
          i,
          "]]` must be a list with a ",
          "non-empty character `type` field.",
          call. = FALSE
        )
      }
    }
  }

  if (!is.list(x$data_summary)) {
    stop("`data_summary` must be a list.", call. = FALSE)
  }

  if (!is.character(x$capabilities)) {
    stop("`capabilities` must be a character vector.", call. = FALSE)
  }

  if (!x$source %in% c("asreml", "brms", "greta")) {
    stop('`source` must be "asreml", "brms", or "greta".', call. = FALSE)
  }

  # greta_meta is required when source = "greta" and forbidden otherwise.
  # The slot holds the post-hoc IR built from a
  # user-supplied greta_model object.
  if (identical(x$source, "greta")) {
    if (!is.list(x$greta_meta)) {
      stop(
        "`greta_meta` must be a list when `source = \"greta\"`.",
        call. = FALSE
      )
    }
    for (req in c("arrays", "canonical_map", "model_dim", "n_data")) {
      if (is.null(x$greta_meta[[req]])) {
        stop(
          "`greta_meta$",
          req,
          "` is required when source = \"greta\".",
          call. = FALSE
        )
      }
    }
  } else {
    if (!is.null(x$greta_meta)) {
      stop(
        "`greta_meta` is only valid when `source = \"greta\"`.",
        call. = FALSE
      )
    }
  }

  invisible(x)
}

# ---------------------------------------------------------------- #
# Predicate                                                        #
# ---------------------------------------------------------------- #

is_fb_terms <- function(x) inherits(x, "fb_terms")

# ---------------------------------------------------------------- #
# Accessors                                                        #
# ---------------------------------------------------------------- #

# Lightweight field accessors. Internal -- used inside emit_*()
# helpers so they don't reach into list slots directly. A shared
# accessor surface insulates the rest of the package from any
# future representation change (e.g., S7 migration in v0.2+).
fb_response <- function(x) x$response
fb_family <- function(x) x$family
fb_link <- function(x) x$link
fb_intercept <- function(x) x$intercept
fb_fixed_terms <- function(x) x$fixed_terms
fb_random_terms <- function(x) x$random_terms
fb_rcov_terms <- function(x) x$rcov_terms
fb_addition_terms <- function(x) x$addition_terms
fb_priors <- function(x) x$priors
fb_data_summary <- function(x) x$data_summary
fb_capabilities <- function(x) x$capabilities
fb_source <- function(x) x$source

# ---------------------------------------------------------------- #
# S3 methods (exported for dispatch only)                          #
# ---------------------------------------------------------------- #

#' Print method for fb_terms (intermediate representation)
#'
#' Internal S3 method, registered for dispatch only. Used during
#' development and inside the `flexybayes()` / `fb()` flow to inspect
#' the parsed model object before backend dispatch.
#'
#' @param x   an `fb_terms` object.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @keywords internal
#' @export
print.fb_terms <- function(x, ...) {
  cat("<fb_terms> flexyBayes intermediate representation\n")
  cat("  source:    ", x$source, "\n", sep = "")
  cat("  response:  ", x$response, "\n", sep = "")

  fam_str <- if (inherits(x$family, "family")) {
    x$family$family
  } else {
    as.character(x$family)
  }
  link_str <- if (inherits(x$family, "family")) {
    x$family$link
  } else if (!is.null(x$link)) {
    x$link
  } else {
    "<default>"
  }
  cat("  family:    ", fam_str, " (", link_str, " link)\n", sep = "")
  cat("  intercept: ", x$intercept, "\n", sep = "")

  cat("  terms:\n")
  cat("    fixed:    ", length(x$fixed_terms), "\n", sep = "")
  cat("    random:   ", length(x$random_terms), "\n", sep = "")
  cat("    rcov:     ", length(x$rcov_terms), "\n", sep = "")
  cat("    addition: ", length(x$addition_terms), "\n", sep = "")

  cat(
    "  priors:    ",
    if (is.null(x$priors)) "<defaults>" else "<user-supplied>",
    "\n",
    sep = ""
  )
  cat(
    "  capabilities: ",
    if (length(x$capabilities)) {
      paste(x$capabilities, collapse = ", ")
    } else {
      "<not yet evaluated by lgm_gate>"
    },
    "\n",
    sep = ""
  )

  invisible(x)
}

#' Format method for fb_terms -- one-line summary
#'
#' Internal S3 method, registered for dispatch only. Used by R's
#' default print/format machinery for compact display in lists and
#' `data.frame` columns.
#'
#' @param x   an `fb_terms` object.
#' @param ... unused.
#' @return character(1).
#' @keywords internal
#' @export
format.fb_terms <- function(x, ...) {
  paste0(
    "<fb_terms: ",
    x$source,
    " ingest, ",
    x$response,
    " ~ ",
    length(x$fixed_terms),
    " fixed + ",
    length(x$random_terms),
    " random>"
  )
}
