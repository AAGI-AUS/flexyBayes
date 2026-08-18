# prior_summary() -- user-facing accessor for the resolved prior
# on a flexyBayes fit.
#
# The fb_terms IR carries the prior slot (`fit$extras$fb_terms$priors`)
# in one of three shapes:
#
#   (1) an `fb_prior` object (user-supplied OR the auto-default
#       uniform-on-SD constructed by `.default_uniform_prior()`);
#   (2) the legacy-bridge list `list(fixed_sd, vc_sd, legacy = TRUE)`
#       carried when neither a user `fb_prior()` nor the auto-default
#       fired (rare; default-fire is the standard path);
#   (3) NULL (pre-IR-recording fits; surfaced as `no_prior_recorded`).
#
# Closest precedent: `brms::prior_summary()` / `rstanarm::prior_summary()`.
# The flexyBayes implementation differs in two intentional ways:
#
#   - it surfaces whether the prior came from the auto-default
#     (so users know the bounded uniform on SD is the working
#     prior, not a placeholder);
#   - it flags `flexybayes_direct_greta` fits as declaration-only,
#     because `fb_greta()` accepts an `fb_prior()` only as a
#     declaration of what the user-built model graph encodes ---
#     `flexyBayes` never re-priors the user's model on that path.

#' Resolved-prior summary for a flexyBayes fit
#'
#' Returns the resolved priors used to fit the model -- either the
#' user-supplied `fb_prior()` object, the auto-default bounded
#' uniform on the standard-deviation scale (Gelman, 2006), or the
#' legacy scalar bridge (`prior_fixed_sd` + `prior_vc_sd`). The
#' return value is an S3 object with a `print()` method; the
#' underlying `fb_prior` (when applicable) is exposed under
#' `$fb_prior` for programmatic access.
#'
#' For `flexybayes_direct_greta` fits the priors are a *declaration*
#' of what the user-built model graph encodes; the summary flags
#' this with `declaration_only = TRUE`.
#'
#' @param object A `flexybayes`, `flexybayes_inla`, or
#'   `flexybayes_direct_greta` object.
#' @param ... Ignored by current methods (reserved for future
#'   per-component selection).
#'
#' @returns A `prior_summary_flexybayes` object, a list carrying the
#'   components below.
#'
#'   \describe{
#'     \item{`kind`}{One of `"fb_prior"`, `"legacy_scalar"`,
#'       `"no_prior_recorded"`.}
#'     \item{`backend`}{The backend the fit ran on:
#'       `"inla"`, `"brms"`, `"greta"`, or `"greta-direct"`.}
#'     \item{`fb_prior`}{The `fb_prior` object (when
#'       `kind == "fb_prior"`).}
#'     \item{`default_origin`}{`"auto"` when the prior was
#'       constructed by the bounded-uniform auto-default;
#'       `"user"` when supplied via the `prior` argument; `NA` for
#'       legacy / no-prior cases.}
#'     \item{`default_scale`, `default_basis`}{Attributes carried by
#'       the auto-default prior naming the response-scale upper
#'       bound and its basis (response-scale `sd(y)`, logit-scale
#'       constant, log-scale constant). `NULL` when the prior was
#'       user-supplied.}
#'     \item{`fixed_sd`, `vc_sd`}{Legacy scalar values (when
#'       `kind == "legacy_scalar"`).}
#'     \item{`declaration_only`}{`TRUE` on `flexybayes_direct_greta`
#'       fits -- the prior is a declaration of the user's
#'       greta-built model, not an enforcement.}
#'   }
#'
#' @examples
#' \dontrun{
#' # live brms (Stan) fit -- needs a working Stan toolchain
#' data(sleepstudy, package = "lme4")
#' fit <- fb_brms(Reaction ~ Days + (1 | Subject),
#'                data = sleepstudy,
#'                n_samples = 100, warmup = 100, chains = 1,
#'                verbose = FALSE, mcmc_verbose = FALSE)
#' prior_summary(fit)
#' }
#'
#' @export
prior_summary <- function(object, ...) UseMethod("prior_summary")

#' @rdname prior_summary
#' @export
prior_summary.flexybayes <- function(object, ...) {
  .prior_summary_impl(object, backend_label = "greta", declaration_only = FALSE)
}

#' @rdname prior_summary
#' @export
prior_summary.flexybayes_inla <- function(object, ...) {
  .prior_summary_impl(object, backend_label = "inla", declaration_only = FALSE)
}

#' @rdname prior_summary
#' @export
prior_summary.flexybayes_brms <- function(object, ...) {
  .prior_summary_impl(object, backend_label = "brms", declaration_only = FALSE)
}

#' @rdname prior_summary
#' @export
prior_summary.flexybayes_direct_greta <- function(object, ...) {
  .prior_summary_impl(
    object,
    backend_label = "greta-direct",
    declaration_only = TRUE
  )
}

#' @rdname prior_summary
#' @export
prior_summary.default <- function(object, ...) {
  stop(
    "`prior_summary()` does not know how to extract priors from ",
    "an object of class ",
    paste(class(object), collapse = "/"),
    ". Define a `prior_summary.<class>` method or pass a ",
    "flexyBayes fit.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------- #
# Implementation                                                   #
# ---------------------------------------------------------------- #

# .prior_summary_engine_default() --- the parameters this package did not
# prior, and why.
#
# The declared `fb_prior` covers the terms the default-prior walker
# reaches. Everything else runs under the engine's own hyperprior, and
# before 0.9.0 `prior_summary()` printed the declaration alone, so a
# combined GxE model reported the shared default on one variance
# component while the other silently carried brms's student_t(3, 0, 2.5).
# The fingerprint already records the gap for triangulate()'s
# matched-prior gate; this reads the same record.
.prior_summary_engine_default <- function(object) {
  fp <- object$extras$fingerprint
  if (!is.null(fp) && length(fp$engine_default_params) > 0L) {
    return(fp$engine_default_params)
  }
  fb_terms <- object$extras$fb_terms
  if (is.null(fb_terms)) {
    return(character(0))
  }
  rec <- tryCatch(.fb_prior_record(fb_terms), error = function(e) NULL)
  if (is.null(rec)) character(0) else rec$engine_default
}

# .prior_summary_residual_lowering() --- name the residual parameter the
# model actually has.
#
# A sectioned residual (`dsum(~ units | f)`, `at(f):units`) makes the
# residual a distributional predictor: the model has no scalar `sigma`,
# and the declared uniform on the SD scale is retargeted onto the
# log-sigma coefficients. Printing `sigma ~ uniform(0, U)` for such a fit
# names a parameter the model does not contain and a distribution Stan
# did not use, so the summary carries the retarget explicitly.
.prior_summary_residual_lowering <- function(object) {
  fb_terms <- object$extras$fb_terms
  if (is.null(fb_terms)) {
    return(NULL)
  }
  sectioned <- Filter(
    function(t) identical(t$type %||% "", "at_units"),
    fb_terms$residual_terms %||% list()
  )
  if (length(sectioned) == 0L) {
    return(NULL)
  }
  term <- sectioned[[1L]]
  term$var %||% term$outer %||% NA_character_
}

.prior_summary_impl <- function(object, backend_label, declaration_only) {
  fb_terms <- object$extras$fb_terms
  priors <- if (!is.null(fb_terms)) fb_terms$priors else NULL

  out <- if (is.null(priors)) {
    list(
      kind = "no_prior_recorded",
      backend = backend_label,
      fb_prior = NULL,
      default_origin = NA_character_,
      default_scale = NULL,
      default_basis = NULL,
      fixed_sd = NA_real_,
      vc_sd = NA_real_
    )
  } else if (isTRUE(priors$legacy)) {
    list(
      kind = "legacy_scalar",
      backend = backend_label,
      fb_prior = NULL,
      default_origin = NA_character_,
      default_scale = NULL,
      default_basis = NULL,
      fixed_sd = priors$fixed_sd,
      vc_sd = priors$vc_sd
    )
  } else if (inherits(priors, "fb_prior")) {
    default_basis <- attr(priors, "fb_prior_default_basis")
    default_scale <- attr(priors, "fb_prior_default_scale")
    origin <- if (!is.null(default_basis)) "auto" else "user"
    list(
      kind = "fb_prior",
      backend = backend_label,
      fb_prior = priors,
      default_origin = origin,
      default_scale = default_scale,
      default_basis = default_basis,
      fixed_sd = NA_real_,
      vc_sd = NA_real_
    )
  } else {
    # Unknown shape -- defensive surfacing rather than silent
    # mis-classification.
    list(
      kind = "unknown_shape",
      backend = backend_label,
      fb_prior = NULL,
      default_origin = NA_character_,
      default_scale = NULL,
      default_basis = NULL,
      fixed_sd = NA_real_,
      vc_sd = NA_real_,
      raw_class = class(priors)
    )
  }

  out$declaration_only <- isTRUE(declaration_only)

  # What the declaration does not cover. Both slots exist so a reader of
  # the printed summary sees the whole prior, not the half this package
  # chose.
  out$engine_default <- .prior_summary_engine_default(object)
  out$residual_lowered_to <- .prior_summary_residual_lowering(object)

  # brms is the authority on what reached Stan, so on a brms fit the
  # engine's own table is carried rather than reconstructed. This is what
  # closes the gap the declaration alone leaves: a retargeted residual
  # prior and an engine default both appear here under their real names.
  out$engine_prior_table <- if (
    identical(backend_label, "brms") && !is.null(object$brms)
  ) {
    tryCatch(brms::prior_summary(object$brms), error = function(e) NULL)
  } else {
    NULL
  }

  structure(out, class = c("prior_summary_flexybayes", "list"))
}


# ---------------------------------------------------------------- #
# Print method                                                     #
# ---------------------------------------------------------------- #

#' @export
print.prior_summary_flexybayes <- function(x, ...) {
  cat("<prior_summary>  backend = ", x$backend, sep = "")
  if (isTRUE(x$declaration_only)) {
    cat("  (declaration only)", sep = "")
  }
  cat("\n")

  if (isTRUE(x$declaration_only)) {
    cat(
      "  Note: fb_greta() does not modify the user's model graph; ",
      "the prior below records what the user declared the model ",
      "encodes.\n",
      sep = ""
    )
  }

  switch(
    x$kind,
    no_prior_recorded = {
      cat("  No prior record attached to the fit.\n")
    },
    legacy_scalar = {
      cat("  Source: legacy scalar bridge\n")
      cat(
        "    prior_fixed_sd = ",
        format(x$fixed_sd),
        "  ",
        "(beta ~ N(0, prior_fixed_sd))\n",
        sep = ""
      )
      cat(
        "    prior_vc_sd    = ",
        format(x$vc_sd),
        "  ",
        "(sigma ~ Lognormal(0, prior_vc_sd))\n",
        sep = ""
      )
    },
    fb_prior = {
      if (identical(x$default_origin, "auto")) {
        cat(
          "  Source: auto-default bounded uniform on SD ",
          "(weakly-informative; half-Cauchy advised for small J)\n",
          sep = ""
        )
        if (!is.null(x$default_scale)) {
          cat("    Upper bound U = ", format(x$default_scale), "\n", sep = "")
        }
        if (!is.null(x$default_basis)) {
          cat("    Scale basis    = ", x$default_basis, "\n", sep = "")
        }
      } else {
        cat("  Source: user-supplied fb_prior()\n")
      }
      cat("\n")
      print(x$fb_prior)
    },
    unknown_shape = {
      cat(
        "  Unknown prior shape (class = ",
        paste(x$raw_class, collapse = "/"),
        "). Inspect fit$extras$fb_terms$priors directly.\n",
        sep = ""
      )
    }
  )

  # ---- What the declaration above does not say ------------------------

  if (!is.null(x$residual_lowered_to) && !is.na(x$residual_lowered_to)) {
    cat(
      "\n  Residual: this model has no scalar `sigma`. The residual is a\n",
      "  distributional predictor with one log-sigma coefficient per level\n",
      "  of `", x$residual_lowered_to, "`, and a declared uniform on the SD\n",
      "  scale is retargeted onto those coefficients on the log scale --\n",
      "  it is not applied as written above.\n",
      sep = ""
    )
  }

  if (length(x$engine_default) > 0L) {
    cat(
      "\n  Carrying the engine's own default (not this package's):\n",
      sep = ""
    )
    for (nm in names(x$engine_default)) {
      cat("    ", nm, " -- ", x$engine_default[[nm]], "\n", sep = "")
    }
  }

  if (!is.null(x$engine_prior_table)) {
    cat("\n  As it reached Stan (brms::prior_summary):\n")
    print(x$engine_prior_table)
  }

  invisible(x)
}
