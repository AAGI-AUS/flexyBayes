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
#' @return A `prior_summary_flexybayes` object (list). Components:
#'   \describe{
#'     \item{`kind`}{One of `"fb_prior"`, `"legacy_scalar"`,
#'       `"no_prior_recorded"`.}
#'     \item{`backend`}{The backend the fit ran on:
#'       `"greta"`, `"inla"`, or `"greta-direct"`.}
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
  invisible(x)
}
