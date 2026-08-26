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
# The flexyBayes implementation differs in one intentional way: it
# surfaces whether the prior came from the auto-default (so users know
# the bounded uniform on SD is the working prior, not a placeholder).

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
#' @param object A `flexybayes_inla` or `flexybayes_brms` object.
#' @param ... Ignored by current methods (reserved for future
#'   per-component selection).
#'
#' @returns A `prior_summary_flexybayes` object, a list carrying the
#'   components below.
#'
#'   \describe{
#'     \item{`kind`}{One of `"fb_prior"`, `"legacy_scalar"`,
#'       `"no_prior_recorded"`.}
#'     \item{`backend`}{The backend the fit ran on: `"inla"` or
#'       `"brms"`.}
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
#'     \item{`scalars_supplied`}{Named logical: whether the caller
#'       supplied `prior_fixed_sd` and `prior_vc_sd`. An argument left
#'       unsupplied leaves the engine's own default in force, which
#'       `fixed_sd_engine_default` names.}
#'     \item{`legacy_vc_applied`}{Named character: the density the
#'       legacy scalar bridge put on each variance component it reached,
#'       empty when `prior_vc_sd` was not supplied.}
#'     \item{`not_applied`}{Named character: a prior this package
#'       declared that the engine has no parameter for -- a residual
#'       prior on a family whose dispersion is a function of the mean.}
#'     \item{`declaration_only`}{Reserved; always `FALSE` for every
#'       current backend.}
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
  # An empty record is an answer, not a missing one: it says every
  # variance component carries a prior this package chose. Treating
  # length 0 as "not recorded" and recomputing sent the summary down a
  # path that had no access to the scalar the fit ran under, so a fit
  # whose components were all covered still printed them as engine
  # defaults, three lines under the header naming the prior they carry.
  if (!is.null(fp) && !is.null(fp$engine_default_params)) {
    return(fp$engine_default_params)
  }
  fb_terms <- object$extras$fb_terms
  if (is.null(fb_terms)) {
    return(character(0))
  }
  rec <- tryCatch(
    .fb_prior_record(
      fb_terms,
      prior_vc_sd = .fb_prior_scalar_value(fb_terms, "vc_sd")
    ),
    error = function(e) NULL
  )
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

# .prior_summary_fixed_engine_default() --- what carries the fixed
# effects when `prior_fixed_sd` was not supplied. Named rather than left
# blank: "this package set no fixed-effect prior" and "the fixed effects
# have no prior" are different statements, and only the first is true.
.prior_summary_fixed_engine_default <- function(backend_label) {
  switch(
    backend_label,
    "brms" = paste0(
      "brms's own defaults -- flat on the population-level coefficients, ",
      "student_t on the intercept, centred on the response"
    ),
    "inla" = paste0(
      "INLA's own control.fixed defaults -- prec = 0.001 on the slopes ",
      "(a standard deviation near 32), and prec.intercept = 0, which is ",
      "flat"
    ),
    "the engine's own fixed-effect default"
  )
}

# .prior_summary_not_applied() --- a declared prior the engine has no
# parameter for.
#
# A residual-scale prior is declared for every fit the default-prior
# walker touches, including families whose brms counterpart carries no
# `sigma`: Gamma parameterises `shape`, Beta parameterises `phi`, and
# neither is on the standard-deviation scale the DSL lives on. The emit
# drops the row -- it has to, or brms refuses the fit -- and the drop is
# reported here rather than left for a reader to infer from its absence.
.prior_summary_not_applied <- function(priors, fb_terms, backend_label) {
  if (!identical(backend_label, "brms") || is.null(fb_terms)) {
    return(character(0))
  }
  fam <- tolower(as.character(fb_terms$family %||% "gaussian")[[1L]])
  if (.fb_family_has_brms_sigma(fam)) {
    return(character(0))
  }
  declares_sigma <- if (inherits(priors, "fb_prior")) {
    any(vapply(
      priors$specs,
      function(s) identical(s$target$type %||% "", "sigma"),
      logical(1)
    ))
  } else {
    isTRUE(priors$legacy)
  }
  if (!isTRUE(declares_sigma)) {
    return(character(0))
  }
  c(
    sigma = paste0(
      "family \"",
      fam,
      "\" has no residual scale parameter in brms -- its ",
      "dispersion is a function of the mean -- so the declared sigma prior ",
      "is not applied, and that dispersion parameter keeps brms's own default"
    )
  )
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

  # Which scalar prior arguments the caller wrote, and therefore which of
  # them this fit actually carries. Before 0.9.2 the printed summary
  # asserted `prior_fixed_sd` and `prior_vc_sd` unconditionally while
  # neither reached the engine on some routes, so the accessor built to
  # answer "what prior did this fit use?" named priors the fit did not
  # use. Absent on a fit built by a direct emit call, which is the same
  # as not supplied.
  ir <- fb_terms %||% list()
  out$scalars_supplied <- c(
    fixed_sd = .fb_prior_scalar_supplied(ir, "fixed_sd"),
    vc_sd = .fb_prior_scalar_supplied(ir, "vc_sd")
  )
  out$fixed_sd_applied <- .fb_prior_scalar_value(ir, "fixed_sd")

  # Which variance components the legacy scalar bridge actually priored,
  # and with what density. Read off the fit's own prior record rather
  # than rebuilt from the scalar, and empty unless the caller supplied
  # `prior_vc_sd`. summary()'s variance-component table projects this:
  # on brms the engine's own table answers first, and on INLA -- where
  # there is no engine prior table to read -- this is what lets the cell
  # name the prior the fit carries instead of the two words that were
  # true only while the bridge handed INLA nothing.
  out$legacy_vc_applied <- if (
    isTRUE(unname(out$scalars_supplied[["vc_sd"]]))
  ) {
    object$extras$fingerprint$priors %||% character(0)
  } else {
    character(0)
  }
  out$fixed_sd_engine_default <- .prior_summary_fixed_engine_default(
    backend_label
  )
  # Which fixed-effect coefficients an fb_prior() covers with a `b()`
  # row. The engine-default sentence below is true only of the
  # coefficients nothing names: since 0.9.2 a `b()` row reaches brms as a
  # coef-keyed prior row and INLA through a per-coefficient
  # `control.fixed` entry, so saying "the fixed effects carry the
  # engine's own defaults" without qualification would be the same
  # one-story breach in the other direction.
  out$fixed_b_named <- if (inherits(priors, "fb_prior")) {
    nm <- vapply(
      priors$specs,
      function(s) {
        if (identical(s$target$type, "b")) s$target$name else NA_character_
      },
      character(1)
    )
    unique(nm[!is.na(nm)])
  } else {
    character(0)
  }

  # A declared residual prior the engine cannot carry, because the family
  # has no residual scale parameter for it to apply to.
  out$not_applied <- .prior_summary_not_applied(
    priors,
    fb_terms,
    backend_label
  )

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

# One line for the fixed-effect prior: the scalar when the caller wrote
# it, the engine's own default by name when they did not.
.prior_summary_cat_fixed_line <- function(x) {
  supplied <- isTRUE(unname((x$scalars_supplied %||%
    c(fixed_sd = FALSE))[["fixed_sd"]]))
  if (supplied) {
    cat(
      "    prior_fixed_sd = ",
      format(x$fixed_sd_applied),
      "  ",
      "(every fixed effect, intercept included, ~ N(0, ",
      format(x$fixed_sd_applied),
      "))\n",
      sep = ""
    )
  } else {
    named <- x$fixed_b_named %||% character(0)
    cat(
      "    prior_fixed_sd:  not supplied -- ",
      if (length(named)) {
        paste0(
          "every fixed effect except ",
          paste(paste0("`", named, "`"), collapse = ", "),
          "\n      (given a prior by the fb_prior() row below) carries "
        )
      } else {
        "the fixed effects carry "
      },
      x$fixed_sd_engine_default,
      "\n",
      sep = ""
    )
  }
  invisible(NULL)
}

#' @export
print.prior_summary_flexybayes <- function(x, ...) {
  cat("<prior_summary>  backend = ", x$backend, sep = "")
  cat("\n")

  switch(
    x$kind,
    no_prior_recorded = {
      cat("  No prior record attached to the fit.\n")
    },
    legacy_scalar = {
      cat("  Source: legacy scalar bridge\n")
      .prior_summary_cat_fixed_line(x)
      if (isTRUE(unname(x$scalars_supplied[["vc_sd"]]))) {
        cat(
          "    prior_vc_sd    = ",
          format(x$vc_sd),
          "  ",
          "(sigma and every variance component ~ Lognormal(0, ",
          format(x$vc_sd),
          ") on the SD scale)\n",
          sep = ""
        )
      } else {
        cat(
          "    prior_vc_sd:     not supplied -- the variance components ",
          "carry each engine's own hyperprior\n",
          sep = ""
        )
      }
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
      .prior_summary_cat_fixed_line(x)
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

  if (length(x$not_applied) > 0L) {
    cat("\n  Declared but not applied on this engine:\n", sep = "")
    for (nm in names(x$not_applied)) {
      cat("    ", nm, " -- ", x$not_applied[[nm]], "\n", sep = "")
    }
  }

  if (!is.null(x$residual_lowered_to) && !is.na(x$residual_lowered_to)) {
    cat(
      "\n  Residual: this model has no scalar `sigma`. The residual is a\n",
      "  distributional predictor with one log-sigma coefficient per level\n",
      "  of `",
      x$residual_lowered_to,
      "`, and a declared uniform on the SD\n",
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
