# fb_cov() -- constructor noun for structured-covariance carriers.
#
# One of the four flexyBayes constructor nouns.
# `fb_cov(M, type, ...)` wraps a known covariance / precision / Cholesky
# / block-diagonal / low-rank carrier in a classed object the random-
# effect terms consume:
#
#   y ~ vm(geno, cov = fb_cov(G, type = "chol"))
#
# It is the v0.4.0 replacement for the v0.3.7 legacy kwarg
# forms (`vm(geno, chol = L)`, `vm(geno, precision = Q)`,
# `vm(geno, blocks = list(...))`, `vm(geno, low_rank_factor = U,
# low_rank_scheme = "pca")`), which deprecate at v0.4.0 and error at
# v0.5.0.
#
# Two faces, mirroring fb_approx():
#   (1) As a standalone runtime constructor it returns a classed list
#       carrying the matrix M, its declared `type`, optional `levels`,
#       a `representation_class` attribute (the locked-registry class
#       the carrier maps to), and a `validation_summary` attribute (a
#       light structural check run at construction). This is the
#       object str() / print() / introspection sees and the
#       architectural-invariant snapshot test exercises.
#   (2) Inside a formula the parser does NOT evaluate the call: it
#       deparse-extracts the matrix symbol and the `type` literal,
#       mapping them onto the existing `cov_representation` IR slot so
#       the matrix continues to resolve via `known_matrices` at fit
#       time exactly as the legacy kwarg forms did. The internal emit /
#       codegen / preflight path is unchanged.
#
# The five carrier types map onto the locked representation registry:
#
#   "dense"     -> dense_cov           V supplied directly.
#   "chol"      -> chol_cov            user supplies the Cholesky L.
#   "precision" -> sparse_precision    user supplies the precision Q.
#   "blocks"    -> block_diagonal      list of per-block V_k.
#   "low_rank"  -> low_rank            rank-K factor U + a registered
#                                      approximation scheme. Reserved
#                                      at v0.4.0 -- the fit path lands
#                                      with the approximate-covariance
#                                      route; the carrier vocabulary is
#                                      active so the surface is stable.

# Carrier type -> locked representation-registry class.
.FB_COV_TYPE_TO_REPRESENTATION <- c(
  dense = "dense_cov",
  chol = "chol_cov",
  precision = "sparse_precision",
  blocks = "block_diagonal",
  low_rank = "low_rank"
)

# Carrier type -> the cov_representation$format token the IR / emit
# path reads. The format vocabulary predates fb_cov(); the
# constructor is a typed front door onto the same tokens.
.FB_COV_TYPE_TO_FORMAT <- c(
  dense = "dense",
  chol = "chol",
  precision = "precision",
  blocks = "blocks",
  low_rank = "low_rank"
)

#' Construct a structured-covariance carrier
#'
#' Wraps a known covariance, precision, Cholesky factor, block-diagonal
#' list, or low-rank factor in a classed carrier object that the
#' structured random-effect terms consume, most directly
#' `vm(geno, cov = fb_cov(G, type = "chol"))`. It is the v0.4.0
#' constructor-noun replacement for the legacy `vm()` keyword forms
#' (`chol = `, `precision = `, `blocks = `, `low_rank_factor = `), which
#' deprecate at this release and are removed at v0.5.0.
#'
#' The carrier's `type` selects how the downstream emit path derives a
#' covariance square root. `"dense"` takes the lower Cholesky of the
#' supplied covariance; `"chol"` uses the supplied factor directly;
#' `"precision"` inverts via the Cholesky of the precision matrix;
#' `"blocks"` assembles a block-diagonal covariance from a list of
#' per-block matrices; `"low_rank"` pairs a rank-K factor with a
#' registered approximation scheme (reserved at v0.4.0 -- the carrier
#' vocabulary is active and validated, while the approximate-covariance
#' fit route activates in a subsequent release).
#'
#' Inside a model formula the carrier is written inline -- the matrix
#' argument names a matrix passed through `known_matrices`, exactly as
#' the legacy keyword forms did. The construction-time check is a light
#' structural probe (shape, lower-triangular pattern, symmetry, block
#' count); the full level-aware validation runs at fit time against the
#' grouping factor.
#'
#' @param M The covariance carrier. For `type = "blocks"` a base-R list
#'   of K square covariance matrices; otherwise a square numeric matrix
#'   (base-R or \pkg{Matrix}). Inside a formula this argument names a
#'   matrix resolved through `known_matrices`.
#' @param type Character(1): the carrier type, one of `"dense"`,
#'   `"chol"`, `"precision"`, `"blocks"`, `"low_rank"`. Defaults to
#'   `"dense"`.
#' @param levels Optional character vector of grouping-factor level
#'   labels the carrier's rows / columns align to. Carried as metadata;
#'   the fit-time validator checks alignment against the fitted factor.
#' @param scheme Character(1), required when `type = "low_rank"`: the
#'   name of a registered approximation scheme (see
#'   [validate_approximation()]). Ignored for the exact types.
#' @param ... Reserved for forward carrier-type options; currently
#'   carried verbatim on the object.
#' @return An `fb_cov` object: a classed list with the matrix `M`, the
#'   `type`, optional `levels`, and (for low-rank) `scheme` as elements,
#'   plus `representation_class` and `validation_summary` attributes.
#' @seealso [fb_approx()], [fb_engine()], [fb_prior()],
#'   [validate_approximation()]
#' @examples
#' G <- crossprod(matrix(rnorm(9L), 3L, 3L))
#' fb_cov(G, type = "dense")
#' L <- t(chol(G))
#' fb_cov(L, type = "chol")
#' @export
fb_cov <- function(M, type = "dense", levels = NULL, scheme = NULL, ...) {
  if (!is.character(type) || length(type) != 1L || !nzchar(type)) {
    stop("fb_cov(): `type` must be a non-empty single string.", call. = FALSE)
  }

  if (!type %in% names(.FB_COV_TYPE_TO_REPRESENTATION)) {
    stop(.fb_refusal_condition(
      reason_code = "fb_cov_type_unknown",
      message = paste0(
        "fb_cov(): '",
        type,
        "' is not a known carrier type. Supported ",
        "types: ",
        paste(names(.FB_COV_TYPE_TO_REPRESENTATION), collapse = ", "),
        "."
      ),
      family_class = "flexybayes_structured_cov_refusal",
      supplied = type
    ))
  }

  if (missing(M) || is.null(M)) {
    stop(.fb_refusal_condition(
      reason_code = "fb_cov_missing_matrix",
      message = paste0(
        "fb_cov(type = \"",
        type,
        "\"): the carrier matrix `M` is ",
        "required."
      ),
      family_class = "flexybayes_structured_cov_refusal"
    ))
  }

  if (type == "blocks" && (!is.list(M) || is.data.frame(M))) {
    stop(.fb_refusal_condition(
      reason_code = "blocks_not_a_list",
      message = paste0(
        "fb_cov(type = \"blocks\"): `M` must be a base-R list of K ",
        "covariance matrices; got class '",
        paste(class(M), collapse = "/"),
        "'."
      ),
      family_class = "flexybayes_structured_cov_refusal"
    ))
  }

  if (type == "low_rank") {
    if (
      is.null(scheme) ||
        !is.character(scheme) ||
        length(scheme) != 1L ||
        !nzchar(scheme)
    ) {
      stop(.fb_refusal_condition(
        reason_code = "low_rank_scheme_required",
        message = paste0(
          "fb_cov(type = \"low_rank\"): a `scheme` naming a registered ",
          "approximation is required. See ?validate_approximation for ",
          "the available schemes."
        ),
        family_class = "flexybayes_structured_cov_refusal"
      ))
    }
    entry <- tryCatch(.lookup_approximation(scheme), error = function(e) NULL)
    if (is.null(entry)) {
      known <- sort(ls(envir = .approximation_registry, all.names = FALSE))
      stop(.fb_refusal_condition(
        reason_code = "approximation_scheme_unknown",
        message = paste0(
          "fb_cov(type = \"low_rank\", scheme = \"",
          scheme,
          "\"): '",
          scheme,
          "' is not a registered approximation scheme. ",
          "Supported scheme",
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
  }

  out <- c(
    list(M = M, type = type, scheme = scheme, levels = levels),
    list(...)
  )
  structure(
    out,
    class = c("fb_cov", "list"),
    type = type,
    levels = levels,
    representation_class = .representation_class(
      .FB_COV_TYPE_TO_REPRESENTATION[[type]]
    ),
    validation_summary = .fb_cov_validation_summary(M, type)
  )
}

#' Test whether an object is an `fb_cov` carrier
#'
#' @param x An object.
#' @return `TRUE` if `x` is an `fb_cov` object.
#' @export
is_fb_cov <- function(x) inherits(x, "fb_cov")

# Light, level-agnostic structural probe run at construction time. The
# fit-time validators (.validate_chol_input, .validate_precision_input,
# .validate_blocks_input) carry the authoritative, group-level-aware
# checks; this summary is the introspection-surface companion. Returns
# a single human-readable string.
.fb_cov_validation_summary <- function(M, type) {
  if (type == "blocks") {
    if (!is.list(M)) {
      return("blocks: malformed carrier")
    }
    K <- length(M)
    sizes <- vapply(
      M,
      function(V) {
        d <- dim(V)
        if (is.null(d) || length(d) != 2L) NA_integer_ else as.integer(d[[1L]])
      },
      integer(1L)
    )
    return(paste0(
      K,
      " block",
      if (K != 1L) "s" else "",
      "; sizes ",
      paste(ifelse(is.na(sizes), "?", sizes), collapse = "+")
    ))
  }

  d <- dim(M)
  if (is.null(d) || length(d) != 2L) {
    return(paste0(type, ": non-matrix carrier (no fit-time shape check yet)"))
  }
  shape <- paste0(d[[1L]], "x", d[[2L]])
  if (d[[1L]] != d[[2L]]) {
    return(paste0(type, ": ", shape, " (non-square -- refused at fit time)"))
  }

  detail <- switch(
    type,
    "chol" = if (.is_lower_triangular(M)) {
      "lower-triangular"
    } else {
      "NOT lower-triangular (refused at fit time)"
    },
    "precision" = if (isSymmetric(unname(as.matrix(M)))) {
      "symmetric"
    } else {
      "NOT symmetric (refused at fit time)"
    },
    "low_rank" = "factor",
    ""
  )
  if (nzchar(detail)) {
    paste0(type, ": ", shape, " ", detail)
  } else {
    paste0(type, ": ", shape)
  }
}

#' Print an `fb_cov` carrier
#'
#' @param x An `fb_cov` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.fb_cov <- function(x, ...) {
  cat("<fb_cov> type = \"", x$type, "\"", sep = "")
  if (!is.null(x$scheme)) {
    cat(" | scheme = \"", x$scheme, "\"", sep = "")
  }
  cat("\n")
  cat("  representation: ", attr(x, "representation_class"), "\n", sep = "")
  vs <- attr(x, "validation_summary")
  if (!is.null(vs) && !is.na(vs)) {
    cat("  carrier: ", vs, "\n", sep = "")
  }
  lv <- attr(x, "levels")
  if (!is.null(lv)) {
    cat(
      "  levels: ",
      length(lv),
      " label",
      if (length(lv) != 1L) "s" else "",
      "\n",
      sep = ""
    )
  }
  invisible(x)
}
