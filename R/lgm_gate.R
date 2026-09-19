# lgm_gate -- LGM (latent Gaussian model) feasibility filter.
#
# Structural-first hybrid filter. Operates on an fb_terms
# object (the IR -- intermediate representation) and decides whether
# the model is feasible for INLA (Integrated Nested Laplace
# Approximations) dispatch.
#
# v0.2.5 implements 9 of the 10 checks (the structural-first set):
# the original six (family, predictor, distributional, re_prior,
# latent_class, hyperparam_budget) plus three INLA-emit-support
# checks (fixed_term_type_inla, random_term_type_inla,
# residual_term_type_inla) that fold emit_inla()
# structural refusals into the gate. Check 10 -- the post-fit
# numerical-confirm gate -- is stubbed; it lands alongside
# emit_inla()'s numerical-confirm pass where it can inspect
# mode.status, mlik finiteness, hyperpar boundary proximity, and
# posterior skewness.
#
# Twin objectives: low false-negative rate (don't refuse INLA-
# feasible models) and zero silent failure (don't accept models
# where INLA gives a wrong answer). Non-silent by construction --
# every refusal carries the rule id that fired plus a one-line
# human gloss; overrides are two-key-armed.
#
# Internal in v0.1 -- called from flexybayes() / fb() / emit_inla()
# during dispatch. Power users can reach it via `flexyBayes:::`.
# print.lgm_refusal is exported so dispatch works in installed
# builds.

# ---------------------------------------------------------------- #
# Public-style entry                                               #
# ---------------------------------------------------------------- #

# Run the LGM-OK filter on an fb_terms object.
#
# Returns one of:
# - the input fb_terms with $capabilities augmented to include
#   "lgm_compatible" plus any soft-warning rule ids (when all
#   structural checks pass);
# - the input fb_terms with $capabilities augmented to include
#   "lgm_force_overridden" plus the bypassed rule ids (when force =
#   "inla" + acknowledge_silent_bias_risk = TRUE + non-empty reason);
# - an lgm_refusal object (when one or more structural checks fail
#   and the override path is not taken).
#
# Callers (emit_inla; flexybayes(... backend = "inla"); fb(... backend
# = "auto")) dispatch on the return value: fb_terms -> proceed;
# lgm_refusal -> stop() the user-facing call with the formatted
# refusal, or attempt re-route per backend = "auto" semantics.
#
# @param fb an fb_terms object.
# @param force NULL or "inla". When "inla" and the structural
#   checks failed, override the refusal (subject to the two-key
#   guard).
# @param acknowledge_silent_bias_risk logical(1). Must be TRUE for
#   override to take effect. Default FALSE -- refusals are non-silent
#   by construction.
# @param reason character(1) non-empty. Logged into capabilities so
#   downstream tooling can audit forced-INLA fits.
# @param preflight Optional `<fb_preflight>` result. When supplied,
#   the gate runs an 11th rule
#   `.lgm_check_memory_feasibility_inla()` that refuses when the
#   INLA path's projected memory exceeds the active ceiling. When
#   NULL (the v0.3.5 calling convention), the gate behaves exactly
#   as today -- 10 structural rules, no memory awareness.
# @return fb_terms (pass / override) or lgm_refusal (default).
lgm_gate <- function(
  fb,
  force = NULL,
  acknowledge_silent_bias_risk = FALSE,
  reason = NULL,
  preflight = NULL
) {
  .check_fb_terms(
    fb,
    "`fb` must be an fb_terms object (see fb_from_asreml)."
  )

  checks <- list(
    .lgm_check_family(fb), # check 1
    # check 1b -- estimated correlation matrices go to brms. Ahead of the
    # structural checks so that the refusal a user reads first is the one
    # that names where the model is fitted.
    .lgm_check_correlated_covariance(fb),
    .lgm_check_predictor(fb), # check 2
    .lgm_check_distributional(fb), # check 3
    .lgm_check_re_prior(fb), # check 4
    .lgm_check_latent_class(fb), # check 5
    .lgm_check_hyperparam_budget(fb), # check 6
    .lgm_check_fixed_term_inla_support(fb), # check 7
    .lgm_check_random_term_inla_support(fb), # check 8
    .lgm_check_residual_term_inla_support(fb), # check 9
    # check 9b -- us() order, capped at five levels on INLA (see below).
    .lgm_check_us_gxe_order(fb),
    # check 11 -- INLA-
    # specific memory-feasibility rule. Trivially passes when the
    # caller did not supply a preflight result (v0.3.5 backward-
    # compatible default); otherwise estimates the INLA-path
    # design memory as a fixed multiple of the indexed estimate
    # and refuses if the INLA projection exceeds the active
    # ceiling. See .lgm_check_memory_feasibility_inla() below for
    # the multiplier rationale.
    .lgm_check_memory_feasibility_inla(fb, preflight),
    # check 12 -- joint shrinkage priors. INLA's latent field is Gaussian
    # given the hyperparameters, and a scale-mixture prior over the fixed
    # effects is not: a regularised horseshoe is a global scale times
    # per-coefficient local scales, and there is no hyperparameter
    # configuration under which the nested Laplace approximation integrates
    # it. The check exists because the alternative is a fit that runs, looks
    # normal, and reports a posterior under a prior nobody asked for.
    .lgm_check_prior_family_inla_support(fb)
    # post-fit numerical-confirm gate -- see emit_inla.
  )

  failures <- Filter(function(r) !isTRUE(r$pass), checks)
  warnings <- Filter(function(r) isTRUE(r$pass) && !is.null(r$warning), checks)

  # All structural checks passed.
  if (length(failures) == 0L) {
    caps <- c(fb$capabilities, "lgm_compatible")
    for (w in warnings) {
      caps <- c(caps, paste0("lgm_warning:", w$rule_id, ":", w$warning))
    }
    fb$capabilities <- caps
    return(fb)
  }

  # Override path -- two-key armed. The structural checks are
  # user-overridable here; the numerical-confirm check runs post-fit in
  # the INLA emitter (see emit_inla.R) and is not overridable.
  if (identical(force, "inla") && isTRUE(acknowledge_silent_bias_risk)) {
    if (
      is.null(reason) ||
        !is.character(reason) ||
        length(reason) != 1L ||
        !nzchar(reason)
    ) {
      stop(
        "`reason` (non-empty character(1)) is required when forcing ",
        "INLA dispatch with acknowledge_silent_bias_risk = TRUE.",
        call. = FALSE
      )
    }
    bypassed <- vapply(failures, function(f) f$rule_id, character(1))
    fb$capabilities <- c(
      fb$capabilities,
      "lgm_force_overridden",
      paste0("lgm_force_reason:", reason),
      paste0("lgm_force_bypassed:", paste(bypassed, collapse = ","))
    )
    return(fb)
  }

  # Default -- return a structured refusal.
  new_lgm_refusal(failures, warnings, fb)
}

# ---------------------------------------------------------------- #
# Refusal class                                                    #
# ---------------------------------------------------------------- #

new_lgm_refusal <- function(failures, warnings = list(), fb = NULL) {
  obj <- list(
    failures = failures,
    warnings = warnings,
    n_failures = length(failures),
    primary_rule = if (length(failures)) {
      failures[[1]]$rule_id
    } else {
      NA_character_
    },
    fb_source = if (!is.null(fb)) fb$source else NA_character_,
    fb_response = if (!is.null(fb)) fb$response else NA_character_,
    # Whether flexyBayes fits the same model on brms; NA when no IR was
    # given. The print method reads it to decide whether to point at brms.
    brms_fits = if (!is.null(fb)) .fb_brms_fits_model(fb) else NA
  )
  class(obj) <- c("lgm_refusal", "list")
  obj
}

is_lgm_refusal <- function(x) inherits(x, "lgm_refusal")

#' Print method for lgm_refusal -- LGM feasibility refusal
#'
#' Internal S3 method, registered for dispatch only. Prints why flexyBayes
#' cannot fit the model on INLA and what to do about it, and nothing else.
#'
#' The refusal used to print sixteen lines: a per-rule diagnostic string,
#' four lines explaining what `backend = "auto"` already does, an
#' unconditional sentence about a roadmap item for mixtures -- printed
#' whether or not the model contained one -- and a two-argument override
#' incantation. The R Journal review of 0.10.0 read one of these and said
#' so. A refusal is read by someone who wanted a fit, so it now carries the
#' reasons and one route onward; `rule_id` and `diagnostic` stay on the
#' object for anything that wants to branch on them, and `?lgm_gate`
#' documents the override.
#'
#' @param x   An `lgm_refusal` object, as returned by the feasibility
#'   gate when a model falls outside the latent Gaussian class.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the refusal it prints.
#' @seealso `vignette("flexyBayes-03-syntax")` lists every spelling the
#'   package covers and what each refusal means.
#' @keywords internal
#' @export
print.lgm_refusal <- function(x, ...) {
  cat("flexyBayes cannot fit this model on INLA.\n")
  for (f in x$failures) {
    cat("  - [", f$rule_id, "] ", f$reason, "\n", sep = "")
  }
  for (w in x$warnings) {
    cat("  ! [", w$rule_id, "] ", w$warning, "\n", sep = "")
  }
  # Point at brms only when flexyBayes fits the model there. The footer used
  # to say so unconditionally, including for models the brms route refuses.
  footer <- if (isFALSE(x$brms_fits)) {
    "See ?lgm_gate"
  } else {
    "Fit it with backend = \"brms\", or see ?lgm_gate"
  }
  cat(footer, " for the rule ids and the override.\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------- #
# Check helpers -- internal                                         #
# ---------------------------------------------------------------- #
#
# Each check returns a list:
#   list(rule_id     = character(1),
#        pass        = logical(1),
#        reason      = character(1) or NULL,
#        diagnostic  = character(1) or NULL,
#        warning     = character(1) or NULL,
#        reason_code = character(1) or NULL)
#
# `reason_code` names a registered refusal. When the primary failure
# carries one, an explicit backend = "inla" request raises that code rather
# than the generic inla_gate_refused.

.lgm_pass <- function(rule_id, warning = NULL) {
  list(
    rule_id = rule_id,
    pass = TRUE,
    reason = NULL,
    diagnostic = NULL,
    warning = warning
  )
}

.lgm_fail <- function(rule_id, reason, diagnostic = NULL, reason_code = NULL) {
  list(
    rule_id = rule_id,
    pass = FALSE,
    reason = reason,
    diagnostic = diagnostic,
    warning = NULL,
    reason_code = reason_code
  )
}

# Check 1 -- Family allowlist
#
# family must be in INLA's likelihood roster. With INLA loaded we
# query inla.models()$likelihood; otherwise fall back to a static
# allowlist refreshed at v0.1 from INLA documentation.
#
# NB: this allowlist mirrors INLA's own likelihood roster and is a
# routing-eligibility check, NOT flexyBayes's notion of supported
# families. The authoritative family gate is .resolve_family()
# (R/utils.R), which every user entry (asreml via fb.R, brms via
# fb_from_brms.R, flexybayes_stream() via R/stream_aggregate.R) passes
# through first, for a family name and a family object alike, and which
# admits only the families flexyBayes can emit. Families INLA recognises
# but flexyBayes cannot emit (e.g. survival / time-to-event) never reach
# this check -- .resolve_family() refuses them up front. When adding
# a family, register its emit path AND widen .resolve_family(); a family
# admitted for brms only that INLA's roster also carries is refused below.
.lgm_check_family <- function(fb) {
  fam_name <- if (inherits(fb$family, "family")) {
    fb$family$family
  } else {
    as.character(fb$family)
  }

  # Admitted at entry for brms only, although INLA's roster carries it: the
  # brms fit is validated (tests/testthat/test-brms-missing-response.R), and
  # no INLA fit is. hurdle_gamma, the other brms-only family, needs no entry
  # here, because INLA's roster has no counterpart and the roster check
  # refuses it.
  if (identical(tolower(fam_name), "lognormal")) {
    return(.lgm_fail(
      "family_allowlist",
      paste0(
        "family \"",
        fam_name,
        "\" is admitted on brms only: INLA's likelihood roster carries it, ",
        "but no INLA fit of it is validated in flexyBayes."
      ),
      diagnostic = paste0(
        "fb$family = \"",
        fam_name,
        "\"; a brms-only family"
      )
    ))
  }

  # Test the spelling INLA will be handed, not the one the user wrote.
  # The two vocabularies differ -- flexyBayes says negative_binomial,
  # negbinom and binary where INLA says nbinomial and binomial -- and the
  # reconciler that maps between them (.resolve_inla_family(), in
  # R/emit_inla.R) runs after this gate. Comparing the flexyBayes
  # spelling against INLA's own roster refused three of the eight
  # supported family spellings for a naming reason, on an engine that
  # carries every one of those likelihoods.
  resolved <- .resolve_inla_family(fb)

  allowed <- .inla_likelihood_allowlist()
  if (tolower(resolved) %in% tolower(allowed)) {
    # INLA takes each family's default link only. The INLA emit passes no
    # link to INLA::inla(), so INLA fits the likelihood's own default, the
    # first link on the family gate's list (.fb_family_links(), R/utils.R).
    # Any other link would be fitted as that default while the fit recorded
    # the link the call asked for. It is refused here, and `auto` routes the
    # model to brms, which receives the link as written.
    link_name <- if (inherits(fb$family, "family")) {
      fb$family$link
    } else {
      fb$link
    }
    default_link <- .fb_family_links()[[tolower(fam_name)]][1L]
    if (
      !is.null(link_name) &&
        !is.null(default_link) &&
        !identical(tolower(link_name), default_link)
    ) {
      return(.lgm_fail(
        "family_allowlist",
        paste0(
          "family \"",
          fam_name,
          "\" with link \"",
          link_name,
          "\" is fitted on brms only: flexyBayes's INLA route passes the ",
          "family's default link, \"",
          default_link,
          "\", and no other."
        ),
        diagnostic = paste0(
          "fb$family = \"",
          fam_name,
          "\"; fb$link = \"",
          link_name,
          "\"; INLA default link = \"",
          default_link,
          "\""
        )
      ))
    }
    return(.lgm_pass("family_allowlist"))
  }

  spelling_note <- if (!identical(tolower(resolved), tolower(fam_name))) {
    paste0(" (INLA spelling \"", resolved, "\")")
  } else {
    ""
  }
  .lgm_fail(
    "family_allowlist",
    paste0(
      "family \"",
      fam_name,
      "\"",
      spelling_note,
      " is not among INLA's likelihoods, so flexyBayes cannot fit it on ",
      "INLA."
    ),
    diagnostic = paste0(
      "fb$family = \"",
      fam_name,
      "\"; resolved INLA family = \"",
      resolved,
      "\""
    )
  )
}

.inla_likelihood_allowlist <- function() {
  # Prefer the live INLA roster when the package is installed
  # (Suggests + Additional_repositories declared in DESCRIPTION).
  # Fall back to the static list otherwise.
  if (requireNamespace("INLA", quietly = TRUE)) {
    dyn <- tryCatch(names(INLA::inla.models()$likelihood), error = function(e) {
      NULL
    })
    if (!is.null(dyn) && length(dyn)) return(dyn)
  }
  .inla_likelihood_allowlist_static()
}

# Static allowlist -- frozen at v0.1; refresh on each INLA major
# release if substantive changes are introduced. Covers the common
# agricultural / health / ecological likelihoods.
.inla_likelihood_allowlist_static <- function() {
  c(
    # Continuous
    "gaussian",
    "stdnormal",
    "lognormal",
    "logistic",
    "T",
    "exponential",
    "loggaussian",
    "loglogistic",
    "loggamma",
    "gamma",
    "beta",
    # Counts
    "poisson",
    "zeroinflatedpoisson0",
    "zeroinflatedpoisson1",
    "nbinomial",
    "nbinomial2",
    "zeroinflatednbinomial0",
    "zeroinflatednbinomial1",
    # Binary / binomial
    "binomial",
    "binary",
    "betabinomial",
    "cbinomial",
    "zeroinflatedbinomial0",
    "zeroinflatedbinomial1",
    # Survival
    "exponentialsurv",
    "weibullsurv",
    "loggaussiansurv",
    "lognormalsurv",
    "coxph",
    # Negative-binomial alias common in asreml workflows
    "negative_binomial",
    "negbinom"
  )
}

# Check 2 -- Predictor linearity
#
# asreml ingest is always linear (the asreml-format DSL has no
# non-linear syntax in v0.1). brms ingest (deliverable 1) sets a
# `nl = TRUE` flag or surfaces a non-linear term type -- refuse
# there.
.lgm_check_predictor <- function(fb) {
  has_nl <- any(vapply(
    fb$fixed_terms,
    function(t) {
      isTRUE(t$type == "non_linear") || isTRUE(t$nl)
    },
    logical(1)
  ))
  if (has_nl) {
    return(.lgm_fail(
      "predictor_linearity",
      paste0(
        "non-linear predictor detected (brms `nl = TRUE` or ",
        "equivalent); flexyBayes's INLA route takes a predictor that is ",
        "linear in its terms."
      ),
      diagnostic = "fixed_terms with type = \"non_linear\""
    ))
  }
  .lgm_pass("predictor_linearity")
}

# Check 3 -- Distributional regression
#
# brms-style: refuse if any auxiliary parameter (sigma, phi, nu,
# zi, hu) has a non-trivial RHS. asreml v0.1 has none -- addition
# slot may carry weights but never dpar regressions.
.lgm_check_distributional <- function(fb) {
  dpar_terms <- Filter(
    function(t) {
      isTRUE(grepl("^dpar_", t$type)) && !isTRUE(t$is_intercept_only)
    },
    fb$addition_terms
  )
  if (length(dpar_terms)) {
    return(.lgm_fail(
      "distributional_regression",
      paste0(
        "an auxiliary parameter has its own regression (",
        paste(
          vapply(dpar_terms, function(t) t$type, character(1)),
          collapse = ", "
        ),
        "); flexyBayes's INLA route takes one value for each auxiliary ",
        "parameter."
      ),
      diagnostic = "addition_terms with type matching ^dpar_"
    ))
  }
  .lgm_pass("distributional_regression")
}

# Check 4 -- Random-effect prior is Gaussian / MVN / GMRF
#
# Legacy asreml priors (lognormal scalar wrapped via
# fb$priors$legacy = TRUE) map to Gaussian-equivalent on sigma and pass
# trivially. fb_prior() declarations are inspected term-by-term below:
# an explicit non-Gaussian random-effect prior fails the check, since
# INLA's latent field is Gaussian by construction.
.lgm_check_re_prior <- function(fb) {
  if (isTRUE(fb$priors$legacy)) {
    return(.lgm_pass("re_gaussian_prior"))
  }
  if (is.null(fb$priors)) {
    return(.lgm_pass("re_gaussian_prior"))
  }

  # If priors is an fb_prior object (deliverable 3), look for
  # explicit non-Gaussian RE priors. Until that DSL ships, accept.
  bad <- character(0)
  if (inherits(fb$priors, "fb_prior") && !is.null(fb$priors$re)) {
    for (entry in fb$priors$re) {
      if (
        !is.null(entry$family) &&
          !tolower(entry$family) %in%
            c("normal", "gaussian", "mvn", "gmrf", "half_normal", "exponential")
      ) {
        bad <- c(bad, entry$family)
      }
    }
  }
  if (length(bad)) {
    return(.lgm_fail(
      "re_gaussian_prior",
      paste0(
        "non-Gaussian random-effect prior detected (",
        paste(unique(bad), collapse = ", "),
        "); flexyBayes fits random effects on INLA as Gaussian fields only."
      ),
      diagnostic = "fb$priors$re entries with non-Gaussian family"
    ))
  }
  .lgm_pass("re_gaussian_prior")
}

# Check 5 -- Latent-class detection
#
# Refuse models that introduce a discrete latent variable: finite
# mixtures, hidden Markov models, multistate transitions. Also
# catches family = "mixture(...)" from brms.
.lgm_check_latent_class <- function(fb) {
  bad_types <- c("mixture", "hmm", "multistate")
  found <- character(0)

  for (slot_name in c("fixed_terms", "random_terms", "residual_terms")) {
    slot <- fb[[slot_name]]
    for (term in slot) {
      if (is.character(term$type) && term$type %in% bad_types) {
        found <- c(found, term$type)
      }
    }
  }
  fam_name <- if (inherits(fb$family, "family")) {
    fb$family$family
  } else {
    as.character(fb$family)
  }
  if (
    is.character(fam_name) &&
      grepl("^mixture", tolower(fam_name))
  ) {
    found <- c(found, paste0("family:", fam_name))
  }

  if (length(found)) {
    return(.lgm_fail(
      "latent_class",
      paste0(
        "latent-class structure detected (",
        paste(unique(found), collapse = ", "),
        "); flexyBayes fits no finite mixture or hidden Markov model on ",
        "INLA."
      ),
      diagnostic = "fb$family or fb$<slot>_terms"
    ))
  }
  .lgm_pass("latent_class")
}

# Check 6 -- Hyperparameter budget
#
# Count hyperparameters: 1 likelihood (depending on family) + per
# random-effect contribution (1 for simple, 2 for fa, k+ for fa_gxe,
# n*(n+1)/2 for us_gxe, etc.) + per heterogeneous-residual term. Refuse
# if total > 15 (CCD/grid intractable). Warn at >10 -- emit_inla
# (deliverable 5+) will set int.strategy = "ccd" automatically.
.lgm_check_hyperparam_budget <- function(
  fb,
  hard_limit = 15L,
  warn_limit = 10L
) {
  n_hyp <- .lgm_count_hyperparams(fb)

  if (n_hyp > hard_limit) {
    return(.lgm_fail(
      "hyperparam_budget",
      paste0(
        "hyperparameter count = ",
        n_hyp,
        " is above the ",
        hard_limit,
        " flexyBayes allows on INLA, which integrates over the ",
        "hyperparameters numerically."
      ),
      diagnostic = paste0(
        "count = likelihood (",
        .lgm_family_hypers(fb),
        ") + random + residual contributions"
      )
    ))
  }

  if (n_hyp > warn_limit) {
    return(.lgm_pass(
      "hyperparam_budget",
      warning = paste0(
        "hyperparameter count = ",
        n_hyp,
        " exceeds soft limit ",
        warn_limit,
        "; emit_inla will set int.strategy = \"ccd\"."
      )
    ))
  }

  .lgm_pass("hyperparam_budget")
}

# Likelihood-side hyperparameter count, by family.
.lgm_family_hypers <- function(fb) {
  fam_name <- if (inherits(fb$family, "family")) {
    fb$family$family
  } else {
    as.character(fb$family)
  }
  switch(
    tolower(as.character(fam_name)),
    "gaussian" = 1L,
    "stdnormal" = 0L,
    "binomial" = 0L,
    "binary" = 0L,
    "poisson" = 0L,
    "negative_binomial" = 1L,
    "negbinom" = 1L,
    "nbinomial" = 1L,
    "gamma" = 1L,
    "beta" = 1L,
    "t" = 2L,
    "lognormal" = 1L,
    "exponential" = 0L,
    1L
  )
}

.lgm_count_hyperparams <- function(fb) {
  fam_hypers <- .lgm_family_hypers(fb)

  re_hypers <- 0L
  for (term in fb$random_terms) {
    contrib <- switch(
      term$type,
      "simple" = 1L,
      "ide" = 1L,
      "id" = 1L,
      "vm" = 1L,
      "ped" = 1L,
      "fa" = 2L,
      "fa_gxe" = if (!is.null(term$k)) {
        as.integer(2L * term$k)
      } else {
        2L
      },
      "us_gxe" = {
        n <- if (!is.null(term$n_outer)) term$n_outer else 4L
        as.integer(n * (n + 1L) / 2L)
      },
      "at_simple" = if (!is.null(term$n_outer)) {
        as.integer(term$n_outer)
      } else {
        2L
      },
      "at_units" = if (!is.null(term$n_var)) {
        as.integer(term$n_var)
      } else {
        2L
      },
      # The autoregressive latent fields carry more than one hyperparameter
      # each, and the budget is what decides whether INLA's numerical
      # integration is tractable, so an under-count here buys a fit that
      # should have been refused. INLA's f(model = "ar1") has a marginal
      # precision and a lag-1 correlation (2); the grouped form adds the
      # group correlation (3). The observation-level nugget that completes
      # the four-parameter separable model is the likelihood's own
      # precision, already counted by .lgm_family_hypers().
      "ar1" = 2L,
      "ar1_spatial" = 3L,
      "spline" = 1L,
      "polynomial" = 1L,
      1L
    )
    re_hypers <- re_hypers + contrib
  }

  residual_hypers <- 0L
  for (term in fb$residual_terms) {
    if (term$type != "units") {
      residual_hypers <- residual_hypers + 1L
    }
  }

  fam_hypers + re_hypers + residual_hypers
}

# ---------------------------------------------------------------- #
# Checks 7-9 -- INLA emit-support allowlists                       #
# ---------------------------------------------------------------- #
#
# Gate truth = emit truth. Each check
# mirrors the corresponding switch() default in
# `.build_inla_formula()` (R/emit_inla.R) and refuses at the gate
# what emit_inla would refuse at emit time. With these checks in
# place, `lgm_gate(fb) == accept` is a sufficient condition for
# `.build_inla_formula(gated)` to succeed structurally; the
# emit-side `stop()`s become internal contract-violation
# assertions.
#
# Allowlists are maintained as flat character vectors next to the
# checks so the gate's source-of-truth is co-located with the
# refusal text. Adding a new IR term type in the codegen layer
# requires updating both the emit switch() and the corresponding
# allowlist here.

# Allowlist: fixed-term types that emit_inla()'s fixed-effect
# switch() accepts. Anything else is refused at the gate with the
# offending type named. The v0.2.6 release adds
# `factor_numeric_interaction` -- treatment-coded indexed slopes for
# factor x continuous interactions -- to the allowlist. INLA's
# native `f:x` notation handles this term shape directly for the
# gaussian-identity case; non-gaussian cases have no active-engine
# route to the engine that represents the term.
.inla_fixed_term_type_allowlist <- function() {
  c(
    "factor",
    "continuous",
    "interaction",
    "factor_interaction",
    "factor_numeric_interaction",
    "expression"
  )
}

# Allowlist: random-term types that emit_inla()'s random-effect
# switch() accepts. The refused set (vm/ped/at/us/fa/ar1 +
# variants) has no INLA representation; the refusal
# text names the structured-covariance class.
.inla_random_term_type_allowlist <- function() {
  # "simple_slope_uncor" is an ordinary member: emit_inla() writes it as
  # f(<g>, model = "iid") + f(<g>_for_slope, <x>, model = "iid"). It used
  # to be conditional on a runtime verification artefact, retired at
  # 0.10.0.9000 -- the agreement it asserted is checked in the suite
  # against lme4, which is where every other capability claim here lives.
  # "simple_slope_cor" is not a member: the INLA emit still builds its
  # iid2d field, and .lgm_check_correlated_covariance() refuses it.
  c("simple", "ide", "id", "spline", "simple_slope_uncor")
}

# Check 1b -- estimated correlation matrices are fitted on brms.
#
# INLA carries an estimated correlation matrix as an iid2d or iidKd field
# whose hyperparameters take one joint Wishart prior. That prior couples
# the variances and the correlations whatever its scale, and fitted on INLA
# neither structure reproduced REML's covariance nor stayed invariant to
# rescaling the response: us(env):gen on agridat::besag.met returned a
# genetic correlation of 0.142 against ASReml's 0.930, and (Days | Subject)
# on lme4::sleepstudy a slope variance a third above REML's. brms fits both
# under separate standard-deviation and correlation priors, so the gate
# refuses them here and backend = "auto" routes them there. The INLA emit
# branches remain; readmitting either structure is its own decision.
.FB_INLA_CORRELATED_COV_TYPES <- c("simple_slope_cor", "us_gxe")

.lgm_check_correlated_covariance <- function(fb) {
  rule_id <- "correlated_covariance_requires_brms"
  hits <- Filter(
    function(t) (t$type %||% "") %in% .FB_INLA_CORRELATED_COV_TYPES,
    fb$random_terms %||% list()
  )
  if (length(hits) == 0L) {
    return(.lgm_pass(rule_id))
  }
  spellings <- vapply(hits, .fb_correlated_cov_spelling, character(1L))
  .lgm_fail(
    rule_id,
    paste0(
      paste(spellings, collapse = " and "),
      " estimates a correlation matrix, and flexyBayes fits estimated ",
      "correlation matrices on brms only. On INLA the term is one iid2d or ",
      "iidKd field under a joint Wishart prior that couples the variances ",
      "and the correlations, and on that route the fitted covariance did ",
      "not match REML's and changed when the response was rescaled. Pass ",
      "backend = \"brms\", or leave backend = \"auto\", which routes there."
    ),
    diagnostic = paste0(
      "fb$random_terms type(s): ",
      paste(vapply(hits, function(t) t$type, character(1L)), collapse = ", ")
    ),
    reason_code = rule_id
  )
}

# The user's spelling of a correlated term, for the refusal message.
.fb_correlated_cov_spelling <- function(term) {
  if (identical(term$type, "us_gxe")) {
    return(paste0("`us(", term$outer, "):", term$inner, "`"))
  }
  paste0("`(", term$slope_var, " | ", term$var, ")`")
}

# Allowlist: residual-term types that emit_inla()'s residual-term guard
# accepts. INLA folds residual variance into the likelihood and
# does not represent structured-residual forms.
.inla_residual_term_type_allowlist <- function() {
  # 0.9.0: the separable AR1 field moved to the random side. INLA emits it
  # as a latent autoregressive field plus the observation-level Gaussian
  # nugget -- four parameters against the three of ASReml's nugget-free
  # residual -- so the residual spelling is refused by name at dispatch
  # (.refuse_unrepresentable_structures(), R/dispatch.R) and only the
  # homogeneous residual reaches the INLA emit.
  c("units", "at_units")
}

# A sectioned residual reaches INLA through the multiple-likelihood
# interface, which carries one observation precision per section -- exactly
# `dsum(~ units | f)`. The `level` argument conditions the residual on
# selected levels instead, which is a different model, and stays refused.
.fb_inla_sectioned_residual <- function(fb) {
  for (term in fb$residual_terms %||% list()) {
    if (identical(term$type %||% "", "at_units") && is.null(term$level)) {
      return(term)
    }
  }
  NULL
}

# Check 9b -- the order of an unstructured GxE covariance on INLA.
#
# INLA carries an unstructured k-dimensional covariance as iid2d..iid5d and,
# beyond k = 5, as iidkd, whose parameterisation and hyperparameter naming
# differ. A k-level unstructured covariance costs k(k+1)/2 parameters -- 210
# at k = 20. flexyBayes caps it at five levels on INLA and refuses by name
# above that. (Check 1b refuses every us() term on INLA in any case.)
.FB_INLA_US_MAX_ORDER <- 5L

.lgm_check_us_gxe_order <- function(fb) {
  rule_id <- "us_gxe_order_inla"
  for (term in fb$random_terms %||% list()) {
    if (!identical(term$type %||% "", "us_gxe")) {
      next
    }
    k <- term$n_outer %||% length(term$outer_levels %||% character(0))
    if (is.numeric(k) && k > .FB_INLA_US_MAX_ORDER) {
      return(.lgm_fail(
        rule_id,
        paste0(
          "us(",
          term$outer,
          "):",
          term$inner,
          " asks for an unstructured covariance across ",
          k,
          " levels, which is ",
          k * (k + 1L) / 2L,
          " parameters, and flexyBayes caps it at ",
          .FB_INLA_US_MAX_ORDER,
          " levels on INLA. Use backend = \"brms\", which fits the full ",
          "structure, or diag(",
          term$outer,
          "):",
          term$inner,
          " for one variance per level."
        ),
        diagnostic = paste0("n_outer = ", k)
      ))
    }
  }
  .lgm_pass(rule_id)
}

# What the OTHER active engine does with this model, asked rather than
# assumed.
#
# The gate's refusals used to close with "and no other active backend
# represents it either", which was true when they were written and became
# false the moment brms learned to emit a term. A gate refusal for us(f):g
# or a nested random effect would then tell a user no engine could fit a
# model that `backend = "auto"` fits on the next line. The predicate is the
# authority on what brms represents, so the sentence asks it.
.lgm_other_engine_note <- function(fb) {
  if (isTRUE(.fb_brms_fits_model(fb))) {
    return(paste0(
      "flexyBayes fits this model on brms: pass backend = \"brms\", or ",
      "leave backend = \"auto\", which routes there when INLA refuses."
    ))
  }
  "flexyBayes cannot fit it on brms either."
}

# Whether flexyBayes fits the model on brms: the brms capability predicate
# passes and the brms formula reconstruction lowers every term. The
# predicate alone is not enough -- it passed corh(), an mgcv s() term and a
# level-conditioned at() while the reconstruction refused each by name, so
# a gate refusal told the user to fit on brms a model brms would refuse. An
# error that is not a flexyBayes refusal (brms not installed, say) says
# nothing about the lowering and leaves the predicate's answer standing.
.fb_brms_fits_model <- function(fb) {
  if (!isTRUE(.capability_brms(fb))) {
    return(FALSE)
  }
  lowered <- tryCatch(
    .fb_to_brms_formula(fb),
    flexybayes_refusal = function(e) e,
    error = function(e) NULL
  )
  !inherits(lowered, "flexybayes_refusal")
}

# Human-readable label for a structured-covariance random-term
# type, used in the refusal diagnostic so the user sees the
# semantic class (e.g., "variance-matrix / GBLUP") not just the
# bare type string.
.inla_random_term_class_label <- function(term_type) {
  switch(
    term_type,
    "vm" = "variance-matrix (e.g., GBLUP / kinship-driven random effect)",
    "ped" = "pedigree-based (numerator relationship matrix)",
    "at" = "heterogeneous variance by factor level",
    "at_simple" = "heterogeneous variance by simple factor",
    "us" = "unstructured covariance",
    "fa" = "factor-analytic",
    "fa_gxe" = "factor-analytic GxE",
    "us_gxe" = "unstructured GxE",
    "ar1" = "autoregressive lag-1",
    "ar1_spatial" = "AR1 spatial",
    "corh_gxe" = "equicorrelated GxE",
    "nested" = "nested / interaction random effect",
    "combo" = "multi-way interaction random effect",
    "polynomial" = "polynomial random effect",
    paste0("random term type \"", term_type, "\"")
  )
}

# Human-readable label for a structured-residual term type.
.inla_residual_term_class_label <- function(term_type) {
  switch(
    term_type,
    "at_units" = "heterogeneous residual by factor",
    "dsum" = "diagonal heterogeneous residual",
    "ar1_units" = "AR1 residual",
    paste0("residual term type \"", term_type, "\"")
  )
}

# Check 7 -- Fixed-term type allowlist for the INLA emit path.
#
# emit_inla()'s fixed-effect switch() recognises factor /
# continuous / interaction / factor_interaction / expression.
# Anything else falls through to a stop(); we
# catch that at the gate instead. In v0.2 the codegen does not
# produce out-of-allowlist fixed-term types; this check is
# defensive against codegen drift.
.lgm_check_fixed_term_inla_support <- function(fb) {
  allowed <- .inla_fixed_term_type_allowlist()
  bad <- character(0)
  for (term in fb$fixed_terms) {
    if (!is.character(term$type) || !(term$type %in% allowed)) {
      bad <- c(bad, as.character(term$type))
    }
  }
  if (length(bad)) {
    return(.lgm_fail(
      "fixed_term_type_inla",
      paste0(
        "fixed term type \"",
        bad[[1L]],
        "\" is outside the INLA emit allowlist. ",
        "Re-route via backend = \"brms\" (Stan passthrough)."
      ),
      diagnostic = paste0(
        "fb$fixed_terms type(s): ",
        paste(unique(bad), collapse = ", "),
        "; allowed: ",
        paste(allowed, collapse = ", ")
      )
    ))
  }
  .lgm_pass("fixed_term_type_inla")
}

# Check 8 -- Random-term type / format allowlist for the INLA
# emit path.
#
# emit_inla()'s random-effect switch() recognises simple / ide /
# id (mapped to INLA "iid"), spline (mapped to INLA "rw2"), and
# simple_slope_uncor (two independent iid fields). vm / ped enter the allowlist *conditionally*: the
# precision and pedigree_sparse_precision formats route through
# INLA's f(model = "generic0", Cmatrix = Q) interface. The blocks format
# (v0.3.10) is refused here. The emitter writes K independent
# f(<var>_id_block_<k>, model = "generic0", Cmatrix = solve(V_k)) calls,
# each with its own precision, which for two or more blocks is not the one
# block-diagonal covariance the carrier represents; the one-block fit, which
# equals the precision carrier of the same matrix, did not meet its brms
# oracle within the tolerance fixed for it (U-VC-1). The emitter branch
# stays for a forced override. The dense and chol formats stay refused on INLA (no
# native interface for a user-supplied dense covariance or
# Cholesky factor) and the user is pointed to the actionable
# workaround (re-express as precision via solve(V), or route via
# backend = "brms").
.lgm_check_random_term_inla_support <- function(fb) {
  bad_terms <- list()
  for (term in fb$random_terms) {
    # Refused by check 1b with its own reason; listing it here as well
    # would add a second message calling it an emit gap, which it is not.
    if ((term$type %||% "") %in% .FB_INLA_CORRELATED_COV_TYPES) {
      next
    }
    if (!.is_supported_random_term_for_inla(term)) {
      bad_terms[[length(bad_terms) + 1L]] <- term
    }
  }
  if (length(bad_terms)) {
    first <- bad_terms[[1L]]
    first_type <- as.character(first$type)
    first_label <- .inla_random_term_class_label(first_type)
    fmt <- first$cov_representation$format
    detail <- if (
      first_type %in% c("vm", "ped") && identical(fmt, "blocks")
    ) {
      paste0(
        "A blocks carrier is refused on INLA. flexyBayes's INLA route would ",
        "fit each block as its own field with its own standard deviation, ",
        "which for two or more ",
        "blocks is not the one block-diagonal covariance the carrier ",
        "represents, and a carrier of one block is the precision carrier ",
        "of its matrix. The precision carrier states the covariance with ",
        "one standard deviation: ",
        first_type,
        "(",
        first$var,
        ", cov = fb_cov(Q, type = \"precision\")) with ",
        "Q = Matrix::bdiag(lapply(blocks, solve)) in known_matrices, ",
        "where `blocks` is the list of block matrices. brms fits the same ",
        "model through the dense carrier of the block-diagonal matrix."
      )
    } else if (
      first_type %in%
        c("vm", "ped") &&
        !is.null(fmt) &&
        fmt %in% c("dense", "chol")
    ) {
      paste0(
        "flexyBayes's INLA route takes vm/ped only as a sparse precision ",
        "or a pedigree sparse precision ",
        "(format = \"precision\" or \"pedigree_sparse_precision\"). ",
        "The supplied format is \"",
        fmt,
        "\". Re-express as vm(",
        first$var,
        ", precision = solve(V)) to use INLA's generic0 ",
        "sparse-precision path, or route via backend = \"brms\", ",
        "which accepts the dense / chol format directly."
      )
    } else {
      paste0(
        "flexyBayes has no INLA lowering for this term. ",
        .lgm_other_engine_note(fb)
      )
    }
    return(.lgm_fail(
      "random_term_type_inla",
      paste0(
        "random term type \"",
        first_type,
        "\" (",
        first_label,
        ") is outside the INLA emit ",
        "allowlist. ",
        detail
      ),
      diagnostic = paste0(
        "fb$random_terms type(s): ",
        paste(
          vapply(
            bad_terms,
            function(t) {
              as.character(t$type)
            },
            character(1L)
          ),
          collapse = ", "
        )
      )
    ))
  }
  .lgm_pass("random_term_type_inla")
}

# Predicate: does a single random-effect term land on the INLA
# emit path? Gate truth = emit truth: every term that
# this returns TRUE for must have a matching branch in
# .build_inla_formula()'s switch, and vice versa. The predicate
# accepts vm / ped only when the
# cov_representation slot carries the precision /
# pedigree_sparse_precision format.
.is_supported_random_term_for_inla <- function(term) {
  if (!is.character(term$type)) {
    return(FALSE)
  }
  if (
    term$type %in%
      c("simple", "ide", "id", "spline", "simple_slope_uncor")
  ) {
    return(TRUE)
  }
  # Nested / interaction random effects. INLA fits these as an iid field
  # over the interaction index; .build_inla_formula() has the branch and
  # emit_inla() builds the index. Verified against agridat::john.alpha,
  # the alpha-design that reported this as a bug.
  if (term$type %in% c("nested", "combo")) {
    return(TRUE)
  }
  # diag(f):g / idh(f):g / at(f):g -- k independent iid fields, one per level
  # of the outer factor. The `level` argument conditions on selected levels,
  # which is a different model, and stays refused.
  if (identical(term$type, "at_simple") && is.null(term$level)) {
    return(TRUE)
  }
  # us(f):g is not accepted here: the INLA emit still builds its iidKd
  # field, and .lgm_check_correlated_covariance() refuses it.
  # fa(f, k):g -- a factor-analytic covariance. INLA writes each common
  # factor once as an iid field and again per level as a SCALED COPY, and
  # the copy's scaling parameter IS the loading. Identifiability (k <
  # n_outer) is enforced at parse time.
  if (identical(term$type, "fa_gxe")) {
    return(TRUE)
  }
  # WP16: AR1 (1D) and separable AR1xAR1 spatial latent fields on INLA.
  # emit_inla() validates the one-observation-per-node grid condition and
  # refuses an incomplete / replicated grid there (gate accepts the class;
  # emit validates the instance -- the simple_slope_uncor pattern).
  if (term$type %in% c("ar1", "ar1_spatial")) {
    return(TRUE)
  }
  if (term$type %in% c("vm", "ped")) {
    cov <- term$cov_representation
    if (is.null(cov)) {
      return(FALSE)
    }
    return(
      cov$format %in% c("precision", "pedigree_sparse_precision")
    )
  }
  FALSE
}

# Check 9 -- Residual-term type allowlist for the INLA emit path.
#
# The INLA emit currently produces only `units`, the homogeneous residual.
# That is a flexyBayes gap rather than an INLA one: INLA carries a residual
# variance per section through its multiple-likelihood interface (a matrix
# response with family = rep("gaussian", k)), verified 2026-09-09 against
# a known truth. Until the emit branch lands, only `units` passes. A heterogeneous residual -- dsum(~ units | f), at(f):units
# -- is emitted by brms as distributional regression on sigma; a separable
# autoregressive residual is refused by name before the gate sees it, since
# the field that would fit it belongs on the random side.
.lgm_check_residual_term_inla_support <- function(fb) {
  allowed <- .inla_residual_term_type_allowlist()
  bad <- character(0)
  labels <- character(0)
  for (term in fb$residual_terms) {
    ok <- is.character(term$type) && term$type %in% allowed
    # at(f, level):units conditions the residual on selected levels, which
    # is a different model from one variance per level, and has no
    # multiple-likelihood lowering.
    if (ok && identical(term$type, "at_units") && !is.null(term$level)) {
      ok <- FALSE
    }
    if (!ok) {
      bad <- c(bad, as.character(term$type))
      labels <- c(labels, .inla_residual_term_class_label(term$type))
    }
  }
  if (length(bad)) {
    first_type <- bad[[1L]]
    first_label <- labels[[1L]]
    return(.lgm_fail(
      "residual_term_type_inla",
      paste0(
        "residual term type \"",
        first_type,
        "\" (",
        first_label,
        ") is outside the INLA emit ",
        "allowlist: flexyBayes has no INLA lowering for this residual ",
        "term. ",
        .lgm_other_engine_note(fb)
      ),
      diagnostic = paste0(
        "fb$residual_terms type(s): ",
        paste(unique(bad), collapse = ", "),
        "; allowed: ",
        paste(allowed, collapse = ", ")
      )
    ))
  }
  .lgm_pass("residual_term_type_inla")
}



# ---------------------------------------------------------------- #
# Check 11 -- memory-feasibility rule                              #
# ---------------------------------------------------------------- #
#
# INLA-specific design-memory feasibility check. The hard ceiling
# probe lives in `.maybe_preflight()` (R/dispatch.R) and refuses
# dispatch upstream when the indexed total exceeds the active
# ceiling. This rule, by contrast, fires on the SOFTER ceiling
# where the indexed total fits but the INLA path's projected
# memory does not -- the INLA dense intermediates inflate the
# indexed estimate by a fixed multiplier per the
# 12-class representation taxonomy (dense_factor /
# dense_continuous / dense_factor_interaction representations
# allocate 3-4x the indexed_* variants when INLA materialises
# the full design matrix).
#
# The rule is trivially-pass when `preflight` is NULL (the v0.3.5
# calling convention; the gate is structural-only). When
# `preflight` is non-NULL, the rule reads
# `preflight$total_estimate_bytes` + `preflight$ceiling_bytes` and
# refuses if the inflated INLA projection exceeds the ceiling.
#
# Decision rationale: the structural-vs-
# memory distinction is the central new piece of dispatch
# information. memory_infeasibility_inla is the dedicated reason
# code for memory refusals; the nine structural reason codes
# carry the structural refusals. The two are
# distinguishable from `rejected_routes` alone, without re-running
# the gate.
#
# Multiplier sourcing: the v0.3.6 ship is a single-multiplier
# approximation. Per-term INLA memory estimators (the more
# accurate model) land at v0.3.7 alongside the routing-trace
# vignette refresh. The
# multiplier value (3) reflects the typical INLA dense-design
# inflation observed on the v0.3.x test corpus; it can be
# overridden per-session via getOption("flexyBayes.inla_memory_multiplier",
# default 3).
.lgm_check_memory_feasibility_inla <- function(fb, preflight) {
  rule_id <- "memory_feasibility_inla"
  if (is.null(preflight)) {
    return(.lgm_pass(rule_id))
  }
  if (!inherits(preflight, "fb_preflight")) {
    return(.lgm_pass(rule_id))
  }
  if (
    is.null(preflight$ceiling_bytes) ||
      is.na(preflight$ceiling_bytes)
  ) {
    return(.lgm_pass(rule_id))
  }

  # v0.3.10: when the per-term
  # INLA memory estimator landed on `preflight$memory_estimate`, the
  # check decomposes the estimate by representation class and reports
  # the breakdown on refusal. Falls back to the v0.3.6 single-
  # multiplier path when memory_estimate is absent (legacy callers).
  mem <- preflight$memory_estimate
  if (
    inherits(mem, "fb_memory_estimate") &&
      !is.null(mem$total) &&
      !is.na(mem$total)
  ) {
    inla_estimate <- as.numeric(mem$total)
    if (inla_estimate <= preflight$ceiling_bytes) {
      return(.lgm_pass(rule_id))
    }

    inla_gb <- inla_estimate / 1024^3
    ceiling_gb <- preflight$ceiling_bytes / 1024^3
    bd <- mem$breakdown
    top_share <- if (nrow(bd) > 0L) {
      bd[order(-bd$bytes), ][1L, ]
    } else {
      NULL
    }
    diag_text <- if (!is.null(top_share)) {
      sprintf(
        paste0(
          "per-term total = %.2f GB (overhead %.1fx); ",
          "largest contributor: %s (%s, %.0f%%)"
        ),
        inla_gb,
        mem$overhead_factor %||% 2,
        top_share$term_label,
        top_share$representation,
        (top_share$share %||% 0) * 100
      )
    } else {
      sprintf(
        "per-term total = %.2f GB; ceiling = %.2f GB",
        inla_gb,
        ceiling_gb
      )
    }
    return(.lgm_fail(
      "memory_feasibility_inla_per_term",
      sprintf(
        paste0(
          "INLA path's per-term memory estimate (%.2f GB) ",
          "exceeds the active ceiling (%.2f GB). Re-route via ",
          "backend = \"brms\", or leave backend = \"auto\", which ",
          "routes there when INLA refuses; or reduce structured-cov ",
          "cardinality."
        ),
        inla_gb,
        ceiling_gb
      ),
      diagnostic = diag_text
    ))
  }

  # Legacy single-multiplier path (v0.3.6 default; retained for
  # backward-compat when memory_estimate is unavailable).
  if (
    is.null(preflight$total_estimate_bytes) ||
      is.na(preflight$total_estimate_bytes)
  ) {
    return(.lgm_pass(rule_id))
  }
  multiplier <- getOption("flexyBayes.inla_memory_multiplier", 3)
  if (!is.numeric(multiplier) || length(multiplier) != 1L || multiplier <= 0) {
    multiplier <- 3
  }
  inla_estimate <- preflight$total_estimate_bytes * multiplier
  if (inla_estimate <= preflight$ceiling_bytes) {
    return(.lgm_pass(rule_id))
  }

  inla_gb <- inla_estimate / 1024^3
  ceiling_gb <- preflight$ceiling_bytes / 1024^3
  indexed_gb <- preflight$total_estimate_bytes / 1024^3
  .lgm_fail(
    rule_id,
    sprintf(
      paste0(
        "INLA path's projected design memory (%.2f GB; ",
        "~%.0fx the indexed estimate) exceeds the active ",
        "ceiling (%.2f GB). The indexed representation that would ",
        "avoid this ceiling has no active-engine emit path (see ",
        "NEWS.md). Reduce the model or data ",
        "scale, or raise the memory ceiling if the projection is ",
        "conservative for your hardware."
      ),
      inla_gb,
      multiplier,
      ceiling_gb
    ),
    diagnostic = sprintf(
      "indexed estimate = %.2f GB; INLA multiplier = %.0fx; ceiling = %.2f GB",
      indexed_gb,
      multiplier,
      ceiling_gb
    )
  )
}


# Check 12 -- joint shrinkage prior families
#
# `fb$penalty_blocks` would hold the resolved `fb_prior_blocks` list, one
# entry per declared shrinkage block over the candidate fixed effects. No
# public call reaches this check: nothing in R/ writes that field, and
# fb_select(), the only verb that takes a shrinkage prior, re-emits its
# reference model on brms whatever engine fitted the base (R/fb_select.R).
# The check is kept as an internal guard, so that a joint block handed to
# the INLA route by future code refuses at the gate instead of reaching the
# emitter. A model that carries none passes trivially.
#
# A block whose family is `none` is an explicit opt-OUT of shrinkage for the
# coordinates it names and does not refuse. A Laplace, regularised horseshoe
# or R2D2 prior is Gaussian given one local scale per coefficient, so the
# model can be written as a latent Gaussian model, but it then carries a
# hyperparameter for every coefficient and INLA's integration is built for a
# few. Ridge and elastic-net blocks have no INLA implementation.
.lgm_check_prior_family_inla_support <- function(fb) {
  blocks <- fb$penalty_blocks %||% list()
  if (length(blocks) == 0L) {
    return(.lgm_pass("prior_family_unsupported_inla"))
  }
  families <- vapply(
    blocks,
    function(b) as.character(b$family %||% "none")[[1L]],
    character(1L)
  )
  shrinkage <- setdiff(unique(families), "none")
  if (length(shrinkage) == 0L) {
    return(.lgm_pass("prior_family_unsupported_inla"))
  }
  .lgm_fail(
    "prior_family_unsupported_inla",
    reason = paste0(
      "the joint shrinkage prior(s) ",
      paste(shrinkage, collapse = ", "),
      " are fitted on brms only: flexyBayes has no INLA route for a joint ",
      "shrinkage block. Fit the reference model on brms, as fb_select() ",
      "does, or select with penalty = NULL, which keeps every candidate ",
      "protected."
    ),
    diagnostic = paste0(
      "penalty block families: ",
      paste(families, collapse = ", ")
    )
  )
}
