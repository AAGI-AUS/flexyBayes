# emit_brms -- flexyBayes Stan-passthrough emit via brms.
#
# Stan-passthrough emit on fb_brms(backend = "brms"). The contract:
# flexyBayes builds the IR; this emit reconstructs the brms formula,
# translates the prior via .priors_to_brms(), calls brms::brm(), and
# wraps the result as c("flexybayes_brms", "flexybayes", "list").
# brms then authors the Stan code, compiles, and samples. flexyBayes
# does NOT author Stan code -- the wrap-don't-rewrite pattern.
#
# Three contracts the wrap honours:
#
#   (1) Result-shape compatibility. The c("flexybayes_brms",
#       "flexybayes", "list") dispatch order lets the parent S3
#       methods (coef, vcov, fitted, residuals, formula, family,
#       nobs, model.matrix, summary) work on a populated $glm shim.
#       confint and logLik are overridden in this file to read brms
#       posterior draws via posterior::as_draws_matrix() instead of
#       the bare parent method (which no active engine reaches).
#
#   (2) backend_decision() uniformity. $extras$backend_decision is
#       populated by .dispatch_backend() after this emit returns;
#       the slot's shape mirrors the INLA path.
#
#   (3) triangulate() peer status. The new class is the dispatch
#       key for fb_as_draws_simple.flexybayes_brms (defined in
#       R/triangulate.R) and canonical_names.flexybayes_brms
#       (defined in R/canonical_names.R).
#
# Stan compile latency. brms's first call in a session compiles
# Stan; cmdstanr backend can shave seconds via cache reuse. We do
# not promote a compile-cache option to the v0.2 surface; users
# wanting it pass `.brms_args = list(backend = "cmdstanr")` through
# the dots argument once emit_brms() accepts it (v0.2.5+).
#
# brms compile + sample times are quality-degrading if review_code
# is the actual user intent. emit_brms supports return_code = TRUE
# via brms::make_stancode(), which authors the Stan code without
# compiling or sampling. This is also how fb_brms(backend = "brms",
# review_code = TRUE) populates the <flexybayes_review> code slot.

# The column name .fb_to_brms_formula()'s sigma-distributional offset
# term (C6) references, and emit_brms() materialises on a copy of
# `data`. flexyBayes lowers `weights =` for the Gaussian family as a
# known per-observation PRECISION multiplier, Var(y_i) = sigma^2 / w_i
# -- the ASReml / lme4 / glm(weights=) / INLA `scale=` convention (see
# R/emit_inla.R). brms's own `weights()` addition term is a different
# quantity: per brms's own documentation (brmsformula.Rd, "Additional
# response information" -- "weighted regression ... is implemented by
# multiplying the log-posterior values of each observation by their
# corresponding weights") and confirmed via the Stan code
# brms::make_stancode() generates (`target += weights[n] *
# normal_lpdf(...)`), it is a LIKELIHOOD-POWER weighting: rescaling
# every weight by a common constant leaves the log-posterior optimum
# for sigma UNCHANGED (a normalised weighted average), not scaled by
# that constant -- confirmed empirically against
# lme4::lmer(weights=) on identical simulated data, where the
# addition-term route's sigma diverged from lme4's by about 25%
# (not Monte Carlo noise; stable across repeat fits). Precision
# weighting is instead lowered on brms's sigma distributional
# parameter (log link) as a known OFFSET,
# log(sigma_i) = log(sigma_base) - 0.5 * log(w_i), which reproduces
# lme4::lmer(weights=)'s sigma to within Monte Carlo error (grounding
# spike: lme4 sigma 0.9819 vs brms-offset 0.983 on the same data,
# against 0.7342 for the addition-term route -- WP-C report, item C6,
# carries the full trail). A fixed internal name rather than threading
# a caller-chosen name through both functions; namespaced against
# colliding with a real column the way the rest of this emit's
# synthesised names are. NOT a leading-dot / double-underscore /
# trailing-underscore name -- verified live: brms::brm() refuses
# "Variable names may not contain double underscores or underscores
# at the end."
.FB_BRMS_WEIGHTS_OFFSET_COL <- "flexybayes_obs_weight_offset"

# ---------------------------------------------------------------- #
# Public-style entry -- internal in v0.2                            #
# ---------------------------------------------------------------- #

emit_brms <- function(
  fb,
  data,
  known_matrices = list(),
  weights = NULL,
  n_samples = 1000L,
  warmup = 500L,
  chains = 4L,
  seed = NULL,
  control = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  the_call = NULL,
  fixed = NULL,
  random = NULL,
  residual = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_,
  na_action = NULL,
  requested_backend = NULL,
  requested_aggregate = NULL,
  ...
) {
  .check_fb_terms(fb, "`fb` must be an fb_terms object.")

  .check_installed(
    "brms",
    "Package 'brms' is required for backend = \"brms\". ",
    "Install via:\n  install.packages('brms')\n",
    "A working C++ toolchain (rstan or cmdstanr) is required ",
    "for Stan compilation."
  )
  .check_installed(
    "posterior",
    "Package 'posterior' is required for the brms posterior ",
    "extraction shim. Install via install.packages('posterior')."
  )

  # -- formula + family + prior ------------------------------------
  # A missing response has to be handled explicitly here: brms drops those
  # rows otherwise, silently turning augmentation into deletion.
  mi_mode <- .fb_brms_missing_response_mode(fb, data)
  if (identical(mi_mode, "mi")) {
    attr(fb, "brms_mi_response") <- TRUE
  }
  # C6: materialise the sigma-offset column .fb_to_brms_formula()'s
  # distributional sigma sub-formula references -- see the constant's
  # definition above for why this is an offset on log(sigma), not the
  # raw weights, and not brms's own weights() addition term. `[[<-` on
  # a plain data.frame is copy-on-modify (the caller's own frame is
  # untouched); this follows the same convention
  # .inla_legalise_factor_levels() (R/emit_inla.R, C4) already uses
  # for its own `data[[v]] <- ...`. Only when the weights are actually
  # non-constant -- a constant vector (any value, not only 1) is the
  # unweighted model under a different spelling, and a non-zero offset
  # built from it would silently rescale the reported sigma rather than
  # leaving it unweighted (and, for a family with no sigma at all, the
  # distributional sub-formula this triggers below would error outright
  # on a request that was never really asking to be weighted).
  ir_weights <- .fb_ir_weights(fb) %||% weights
  if (.fb_weights_nonconstant(ir_weights)) {
    data[[.FB_BRMS_WEIGHTS_OFFSET_COL]] <- -0.5 * log(ir_weights)
  }
  brms_form <- .fb_to_brms_formula(fb)
  # Materialised relationship covariances for any vm() / ped() term,
  # threaded into brms via data2 (empty list when there are none).
  brms_data2 <- .fb_brms_data2(fb, known_matrices, data)
  brms_family <- .fb_family_to_brms(fb$family, fb$link)
  brms_prior <- .priors_to_brms(
    fb$priors %||%
      list(legacy = TRUE, fixed_sd = prior_fixed_sd, vc_sd = prior_vc_sd),
    fb,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd
  )
  # A prior naming a term the model does not have is caught here, against
  # brms's own parameter list, rather than left to brms's parser. The
  # parser reports a synthesised Stan name (`b_nonexistent_term`) the
  # user never wrote, untyped (field-sweep FS-20).
  .check_brms_prior_rows_reachable(
    brms_prior,
    formula = brms_form,
    data = data,
    data2 = brms_data2,
    family = brms_family
  )

  # -- return_code: skip the fit; emit the Stan code only -----------
  if (isTRUE(return_code)) {
    stancode <- tryCatch(
      brms::make_stancode(
        formula = brms_form,
        data = data,
        data2 = brms_data2,
        family = brms_family,
        prior = brms_prior,
        silent = 2L
      ),
      error = function(e) {
        stop(
          "brms::make_stancode() failed: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    return(invisible(as.character(stancode)))
  }

  if (verbose) {
    cat(
      "\n-- flexyBayes: brms (Stan passthrough) fit ",
      paste(rep("-", 18), collapse = ""),
      "\n",
      sep = ""
    )
    cat("  formula: ", deparse(brms_form), "\n", sep = "")
    cat("  family : ", .deparse_brms_family(brms_family), "\n", sep = "")
    if (!is.null(brms_prior)) {
      cat("  priors : ", nrow(brms_prior), " row(s)\n", sep = "")
    }
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
  }

  # -- fit ----------------------------------------------------------
  iter <- as.integer(warmup) + as.integer(n_samples)

  t0 <- proc.time()
  # brms warns "It appears as if you have specified an
  # upper bounded prior on a parameter that has no natural upper
  # bound. ... please specify argument 'ub' of 'set_prior'" when a
  # coef-keyed sd row carries a `uniform(0, X)` prior. brms refuses
  # the same row with an explicit `ub` slot ("Prior argument 'coef'
  # may not be specified when using boundaries"), so the warning is
  # structurally unavoidable for the slope-variance row.
  # The bounded prior expression `uniform(0, X)` is still applied
  # correctly on the Stan side; the warning is purely informational
  # and would otherwise leak into every Stan-passthrough fit on
  # (x || g). Suppress only that exact warning text via
  # withCallingHandlers; other brms warnings pass through.
  brmsfit <- tryCatch(
    withCallingHandlers(
      brms::brm(
        formula = brms_form,
        data = data,
        data2 = brms_data2,
        family = brms_family,
        prior = brms_prior,
        chains = as.integer(chains),
        iter = iter,
        warmup = as.integer(warmup),
        silent = if (isTRUE(verbose)) 0L else 2L,
        refresh = if (isTRUE(mcmc_verbose)) 200L else 0L,
        # Reproducibility and sampler tuning reach Stan from the public
        # entry point. brms's own defaults are seed = NA (draw a seed)
        # and control = NULL (Stan's own adapt_delta and max_treedepth),
        # so a NULL from flexybayes() maps onto exactly those.
        seed = if (is.null(seed)) NA else as.integer(seed),
        control = control,
        ...
      ),
      warning = function(w) {
        if (
          grepl(
            paste0(
              "upper bounded prior on a parameter that has ",
              "no natural upper bound"
            ),
            conditionMessage(w),
            fixed = TRUE
          )
        ) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) {
      stop("brms::brm() failed: ", conditionMessage(e), call. = FALSE)
    }
  )
  elapsed <- unname((proc.time() - t0)["elapsed"])

  # -- post-fit summaries ------------------------------------------
  draws_arr <- tryCatch(
    posterior::as_draws_matrix(brmsfit),
    error = function(e) {
      stop(
        "posterior::as_draws_matrix() failed on the brms fit: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  draws_mat <- as.matrix(draws_arr)

  # Posterior summary in a shape compatible with the parent print /
  # summary methods. Single-chain treatment matches the way the
  # generic flexybayes summary expects mean / sd / quantiles.
  post_summary <- tryCatch(
    summary(brmsfit),
    error = function(e) NULL
  )

  # Convergence diagnostics: brms / posterior expose ess_bulk and
  # rhat per parameter. We translate to the parent print method's
  # convention (`gelman$psrf[, "Point est."]` for Rhat; numeric
  # vector for n_eff).
  conv_diag <- .brms_convergence(brmsfit, draws_mat)

  # GLM shim -- minimal scaffolding that lets coef / vcov / fitted /
  # residuals / family / formula / model.matrix dispatch through the
  # parent .flexybayes methods. Fixed-effect coef names follow brms's
  # b_<term> convention with the b_ prefix stripped (brms convention
  # for canonical names; see canonical_names.flexybayes_brms in
  # R/canonical_names.R for the matching identity mapper).
  glm_obj <- .brms_glm_shim(
    brmsfit,
    draws_mat,
    data,
    family = brms_family,
    formula_used = brms_form
  )

  # Variance components: extract sd_<group>__Intercept rows from the
  # brms posterior + the residual sigma. Same shape as the existing
  # variance_comps table that print.flexybayes / summary expect.
  vc_table <- .brms_variance_comps(brmsfit, draws_mat)

  extras <- structure(
    list(
      summary = post_summary,
      convergence = conv_diag,
      variance_comps = vc_table,
      blups = NULL, # brms ranef() lives on $brms slot
      predictions = NULL, # brms posterior_predict on demand
      code = NULL, # populated only on review_code path
      param_names = colnames(draws_mat),
      parse_info = list(
        fixed = list(
          response = fb$response,
          intercept = isTRUE(fb$intercept),
          terms = fb$fixed_terms
        ),
        random = fb$random_terms,
        residual = fb$residual_terms,
        family = list(family = fb$family, link = fb$link),
        # C6: the observation-weight vector lowered onto brms's sigma
        # offset (see .FB_BRMS_WEIGHTS_OFFSET_COL), NULL when the model
        # carries none. Read back at summary() time by
        # .brms_weights_sigma_row() so the single "sigma" varcomp row
        # can be recovered from the b_sigma_Intercept draws exactly as
        # .brms_residual_by_level_table() already does for `sigma ~
        # 0 + f`.
        weights = ir_weights
      ),
      call_info = list(
        fixed = fixed,
        random = random,
        residual = residual,
        data_name = data_name,
        family = family,
        link = link,
        known_matrices = known_matrices,
        weights = weights,
        n_samples = n_samples,
        warmup = warmup,
        chains = chains,
        seed = seed,
        control = control,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd,
        # The missing-response policy the call requested, as the native
        # word. An `omit` fit that came back as an `augment` re-fit would
        # be a different model under the same name, so the policy belongs
        # in the call record and not only in the na_action summary.
        na_action = na_action,
        # The engine, the representation and the reporting the call
        # ASKED for -- not the ones it resolved to. update() re-issues
        # each of them, and with all three absent an identity update()
        # of a Stan fit came back as an aggregated INLA fit: a different
        # inference engine and a different object shape, under the same
        # name and with nothing said.
        #
        # `backend` is recorded as the request ("auto") rather than the
        # engine the request landed on, because the recorded value is a
        # policy: a re-fit whose model has changed must be free to route
        # again, while a fit that named its engine must come back on
        # that engine.
        backend = requested_backend,
        aggregate = requested_aggregate,
        verbose = verbose
      ),
      run_time = elapsed,
      # What has to match before this fit can be compared with another
      # engine's. Built here rather than at comparison time because the
      # analysis data are what the fit actually saw, and by the time
      # triangulate() runs the caller may hold a different object of the
      # same name. See R/model_fingerprint.R.
      fingerprint = .fb_model_fingerprint(
        fb,
        data,
        prior_vc_sd = prior_vc_sd
      ),
      model_info = list(
        n_obs = nrow(data),
        # brms_form may carry `(x || g)` which
        # stats::model.matrix() rejects (`||` requires logical(1)
        # operands). Compute n_fixed off a fixed-only formula
        # reconstructed from the IR's fixed_terms slot so the
        # random-effects block never reaches model.matrix.
        n_fixed = ncol(stats::model.matrix(
          .fb_to_fixed_only_formula(fb),
          data = data
        )),
        n_random = length(fb$random_terms),
        n_params = ncol(draws_mat),
        family = fb$family,
        link = fb$link
      ),
      the_call = the_call,
      formula = brms_form,
      brms_prior = brms_prior
    ),
    class = "flexybayes_extras"
  )

  out <- structure(
    list(
      glm = glm_obj,
      brms = brmsfit,
      extras = extras
    ),
    class = c("flexybayes_brms", "flexybayes", "list")
  )

  # Read off the same canonical table the INLA path warns from, so a
  # component at the boundary is reported identically on both engines.
  .fb_warn_boundary_collapse(out)

  out
}


# ---------------------------------------------------------------- #
# Prior-row reachability (field-sweep FS-20)                        #
# ---------------------------------------------------------------- #

# Refuse a prior row that names no parameter of the model, before
# brms's parser does.
#
# The oracle is brms's own `default_prior()`, read from the same
# formula / data / family the fit will use: it enumerates every
# (class, coef, group) triple the model carries. A flexyBayes prior row
# whose triple is absent from that enumeration is the user naming a term
# that does not exist. flexyBayes holds the model's term list at this
# point and can say so in the user's vocabulary; before 0.9.2 it left
# the check to brms, which answered with a synthesised Stan parameter
# name (`b_nonexistent_term`) and an untyped `simpleError`.
#
# A row matches when brms lists its exact (class, coef, group) triple.
# The match has to be exact rather than falling back to the class-wide
# row: brms lists a blank-`coef` row for every class the model carries,
# and treating that as a wildcard would let `b(\"nonexistent_term\")`
# match the class-wide `b` row -- which is how the mismatch reached
# brms's parser in the first place. `default_prior()` compiles nothing,
# so the check costs a formula walk.
.check_brms_prior_rows_reachable <- function(
  prior,
  formula,
  data,
  data2,
  family
) {
  if (is.null(prior) || !nrow(prior)) {
    return(invisible(TRUE))
  }
  available <- tryCatch(
    as.data.frame(brms::default_prior(
      formula,
      data = data,
      data2 = data2,
      family = family
    )),
    error = function(e) NULL
  )
  # No oracle, no verdict: an unreadable default_prior() must not turn
  # into a refusal of a model that would otherwise fit.
  if (is.null(available) || !nrow(available)) {
    return(invisible(TRUE))
  }
  norm <- function(x) {
    x <- as.character(x %||% "")
    x[is.na(x)] <- ""
    x
  }
  available$coef <- norm(available$coef)
  available$group <- norm(available$group)
  for (i in seq_len(nrow(prior))) {
    cls <- prior$class[[i]]
    cf <- norm(prior$coef[[i]])
    grp <- norm(prior$group[[i]])
    same_class <- available$class == cls
    if (any(same_class & available$coef == cf & available$group == grp)) {
      next
    }
    if (!nzchar(grp) && !nzchar(cf)) {
      # A class-wide row on a class the model does not carry at all.
      .fb_stop_prior_target_absent(
        target_label = paste0("class \"", cls, "\""),
        kind = "parameter class",
        available = unique(available$class),
        engine = "brms"
      )
    }
    if (nzchar(grp)) {
      .fb_stop_prior_target_absent(
        target_label = paste0("sd(group = \"", grp, "\")"),
        kind = "grouping factor",
        available = unique(available$group[same_class]),
        engine = "brms"
      )
    }
    .fb_stop_prior_target_absent(
      target_label = paste0("b(\"", cf, "\")"),
      kind = "fixed-effect coefficient",
      available = unique(available$coef[same_class]),
      engine = "brms"
    )
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------- #
# IR -> brms formula reconstruction                                 #
# ---------------------------------------------------------------- #

# Walk the IR back to a brms-shaped formula. The corpus we accept on
# this path is the v0.2 five-shape brms corpus: fixed-only Gaussian;
# fixed + one random intercept; fixed + crossed random intercepts;
# binomial RI (single-column 0/1); Poisson RI. Anything outside the
# corpus was already refused at fb_from_brms() ingest -- so the IR
# we see here only carries supported shapes.
.fb_to_brms_formula <- function(fb) {
  response <- fb$response
  if (is.null(response) || !nzchar(response)) {
    stop(
      "IR is missing the response slot; cannot reconstruct the ",
      "brms formula.",
      call. = FALSE
    )
  }

  fixed_labels <- vapply(
    fb$fixed_terms,
    function(t) {
      if (!is.null(t$label)) {
        t$label
      } else if (!is.null(t$var)) {
        t$var
      } else if (!is.null(t$expr)) {
        deparse(t$expr)
      } else {
        NA_character_
      }
    },
    character(1)
  )
  fixed_labels <- fixed_labels[!is.na(fixed_labels)]

  fixed_rhs <- if (length(fixed_labels) > 0L) {
    paste(fixed_labels, collapse = " + ")
  } else {
    "1"
  }
  if (!isTRUE(fb$intercept)) {
    fixed_rhs <- paste0(fixed_rhs, " - 1")
  }

  re_terms <- character(0)
  for (term in fb$random_terms) {
    ttype <- term$type %||% "<unknown>"

    # Genomic / pedigree relationship random effect. Emit brms's native
    # known-covariance group term, (1 | gr(<var>, cov = <K>)). brms
    # Cholesky-factors the supplied covariance internally -- the
    # K = L L' decorrelation that turns a structured genetic random
    # effect into one Stan fits cheaply -- so vm() / ped() reach the
    # brms backend and GBLUP becomes three-engine triangulatable. The
    # covariance matrix itself is threaded into brms's data2 by
    # .fb_brms_data2().
    if (ttype %in% c("vm", "ped")) {
      grp <- term$var
      if (is.null(grp) || !nzchar(grp)) {
        stop("vm() / ped() term missing the group variable.", call. = FALSE)
      }
      covname <- .fb_brms_covname(term)
      re_terms <- c(
        re_terms,
        paste0("(1 | gr(", grp, ", cov = ", covname, "))")
      )
      next
    }

    # Interaction random intercepts (nested A:B / combo A:B:C): brms forms
    # the interaction grouping natively as (1 | A:B). This is the faithful
    # full-HMC path for multi-stratum designed experiments -- validated on
    # besag.met, where brms recovers every variance component while INLA
    # collapses the finest strata. auto therefore keeps INLA's
    # refusal for these models and routes them here.
    if (ttype %in% c("nested", "combo")) {
      int_vars <- if (identical(ttype, "nested")) {
        c(term$outer, term$inner)
      } else {
        as.character(term$vars)
      }
      re_terms <- c(
        re_terms,
        paste0("(1 | ", paste(int_vars, collapse = ":"), ")")
      )
      next
    }

    # Heterogeneous variance over an outer factor: ASReml's diag(f):g, and
    # its synonyms idh(f):g and at(f):g. A separate variance for every level
    # of the outer factor, with no covariance between them.
    #
    # brms writes this as an UNCORRELATED random slope, (0 + f || g): for each
    # level of the grouping factor a vector of level-specific effects, whose
    # standard deviations are free and whose correlations are suppressed.
    #
    # The mapping is validated against ASReml rather than assumed --
    # design-preserving-missingness/oracle_heterogeneous.R fits
    # `diag(site):gen` in ASReml, the explicit dummy() expansion in lme4 and
    # this emit in brms on the same data, and checks the PARAMETER COUNT
    # before the values: a diagonal over k levels is exactly k variances and 0
    # covariances. The two REML arms agree to 3.5e-05 and all three fit the
    # diagonal structure.
    #
    # The count matters because values alone cannot distinguish diag() from
    # us() when the true cross-level correlations are zero. It was checking
    # the count that caught lme4's `||` quietly fitting a correlated block for
    # a factor slope, which a value comparison had passed.
    if (identical(ttype, "at_simple")) {
      if (!is.null(term$level)) {
        stop(.fb_refusal_condition(
          reason_code = "at_level_conditioning_unsupported",
          message = paste0(
            "at(",
            term$outer,
            ", ",
            paste(term$level, collapse = ", "),
            "):",
            term$inner,
            " conditions the random effect on selected ",
            "levels of ",
            term$outer,
            ", which is a different model from a ",
            "heterogeneous variance across all of them. Drop the level to ",
            "fit diag(",
            term$outer,
            "):",
            term$inner,
            "."
          )
        ))
      }
      re_terms <- c(
        re_terms,
        paste0("(0 + ", term$outer, " || ", term$inner, ")")
      )
      next
    }

    # corh(f):g asks for heterogeneous variances with ONE correlation shared
    # by every pair of levels -- k + 1 parameters, sitting between diag()'s k
    # and us()'s k(k+1)/2. brms offers no such structure: a group-level term
    # is either uncorrelated (`||`) or fully unstructured (`|`), with prior
    # classes `sd` and `cor` and nothing between them. (brms does have a
    # compound-symmetry structure, cosy(), but it applies to the RESIDUAL
    # autocorrelation, not to group-level effects.)
    #
    # Approximating it with us() would fit k(k+1)/2 - 1 parameters nobody
    # asked for and report correlations the model was meant to constrain
    # equal; approximating it with diag() would drop the correlation
    # entirely. Both are a different model wearing this one's name.
    if (identical(ttype, "corh_gxe")) {
      stop(.fb_refusal_condition(
        reason_code = "corh_no_equicorrelation_representation",
        message = paste0(
          "corh(",
          term$outer,
          "):",
          term$inner,
          " requests heterogeneous ",
          "variances with a single shared correlation. No active backend ",
          "represents that structure: brms group-level effects are either ",
          "uncorrelated or fully unstructured, with nothing between. Use ",
          "diag(",
          term$outer,
          "):",
          term$inner,
          " for heterogeneous ",
          "variances with no correlation, or us(",
          term$outer,
          "):",
          term$inner,
          " to estimate every pairwise correlation freely."
        )
      ))
    }

    # Unstructured covariance over an outer factor: ASReml's us(f):g. The
    # same construction with correlations RETAINED -- a single character
    # apart from the diagonal case above, and the whole difference between
    # k parameters and k(k+1)/2 of them.
    if (identical(ttype, "us_gxe")) {
      re_terms <- c(
        re_terms,
        paste0("(0 + ", term$outer, " | ", term$inner, ")")
      )
      next
    }

    if (!ttype %in% c("simple", "ide", "id", "simple_slope_uncor")) {
      stop(.fb_refusal_condition(
        reason_code = "brms_cannot_represent_term",
        message = .brms_term_refusal_message(term),
        term_type = ttype
      ))
    }
    grp <- term$var
    if (is.null(grp) || !nzchar(grp)) {
      stop("Random-effect term missing the group variable.", call. = FALSE)
    }
    # Uncorrelated random slope. Round-trip to the brms
    # double-pipe shape -- (x || g) for the intercept + slope form
    # (the lme4 / brms-default semantics for the bare slope),
    # (0 + x || g) for the slope-only form (with explicit
    # intercept suppression).
    if (identical(term$type, "simple_slope_uncor")) {
      sv <- term$slope_var
      if (is.null(sv) || !nzchar(sv)) {
        stop("Uncorrelated random-slope term missing slope_var.", call. = FALSE)
      }
      re_terms <- c(
        re_terms,
        paste0(
          "(",
          if (isTRUE(term$with_intercept)) "" else "0 + ",
          sv,
          " || ",
          grp,
          ")"
        )
      )
      next
    }
    re_terms <- c(re_terms, paste0("(1 | ", grp, ")"))
  }

  rhs <- if (length(re_terms) > 0L) {
    paste(c(fixed_rhs, re_terms), collapse = " + ")
  } else {
    fixed_rhs
  }

  # `mi` marks the response as partially missing so brms samples the absent
  # values as parameters instead of dropping their rows. Added only when the
  # caller asks for it -- see .fb_brms_missing_response_mode(). Observation
  # weights (C6) are NOT an addition term here -- they are lowered on the
  # sigma distributional sub-formula below as a known offset, for the
  # reasons the .FB_BRMS_WEIGHTS_OFFSET_COL definition above documents.
  addition <- character(0)
  if (isTRUE(attr(fb, "brms_mi_response"))) {
    addition <- c(addition, "mi()")
  }
  lhs <- if (length(addition)) {
    paste(response, "|", paste(addition, collapse = " + "))
  } else {
    response
  }

  main <- stats::as.formula(paste(lhs, "~", rhs))

  # Heterogeneous residual variance: ASReml's dsum(~ units | site), and the
  # at(site):units spelling that parses to the same node. A separate residual
  # variance per level of the sectioning factor.
  #
  # brms expresses it as distributional regression on sigma. Two details make
  # or break the mapping, and both are read off a live fit rather than assumed
  # (design-preserving-missingness, the residual-heterogeneity probe):
  #
  #   * the predictor is on LOG sigma, so a coefficient is not a variance --
  #     the variance for a level is exp(2 * b). Comparing brms's `b_sigma_*`
  #     against ASReml's `site!R` directly would compare a log-scale contrast
  #     against a variance.
  #   * `0 + f` gives one coefficient PER LEVEL rather than contrasts against
  #     a reference, which is what lines up with ASReml's per-section
  #     variances.
  #
  # Against ASReml's dsum on the same simulated data, fitted through
  # flexybayes() and its default priors, this returns per-site
  # posterior-mean variances of 0.1146 / 1.1516 / 4.6981 where ASReml
  # gives 0.1093 / 1.1248 / 4.6079 -- within 4.8%, largest on the
  # smallest variance, which is where a posterior mean and a REML point
  # estimate are least alike. Probe 2 of
  # design-preserving-missingness/oracle_heterogeneous.R, run
  # 2026-08-15. The entry point exposes no sampler seed, so those digits
  # are one run's realisation at a bulk effective size above 8,000; the
  # per-site percentages are stable, the last digit is not.
  het <- Filter(
    function(t) identical(t$type %||% "", "at_units"),
    fb$residual_terms %||% list()
  )
  # C6: observation weights ALSO put flexyBayes on brms's sigma
  # distributional sub-formula (a known offset, not a heterogeneous-
  # residual factor) -- so the early return below is reachable only
  # when NEITHER a heterogeneous residual NOR weights are in play.
  has_weight_offset <- .fb_weights_nonconstant(.fb_ir_weights(fb))
  if (length(het) == 0L && !has_weight_offset) {
    return(main)
  }
  if (length(het) > 1L) {
    stop(.fb_refusal_condition(
      reason_code = "heterogeneous_residual_multiple_factors",
      message = paste0(
        "More than one heterogeneous-residual term was supplied (",
        paste(
          vapply(het, function(t) t$var %||% "?", character(1)),
          collapse = ", "
        ),
        "). A residual is sectioned by one factor; nesting several would ",
        "need their interaction, which must be stated explicitly."
      )
    ))
  }
  # sigma_rhs accumulates the two independent reasons the sigma
  # sub-formula can exist: a heterogeneous-residual factor (0 + factor,
  # one coefficient per level) and/or the C6 weight offset. Either can
  # be present alone; both can be present together (a per-section
  # residual variance on data that also carries known observation
  # weights) -- the two compose by `+` on the linear predictor for
  # log(sigma).
  sigma_rhs <- NULL
  if (length(het) == 1L) {
    term <- het[[1L]]
    if (!is.null(term$level)) {
      stop(.fb_refusal_condition(
        reason_code = "at_level_conditioning_unsupported",
        message = paste0(
          "at(",
          term$var,
          ", ",
          paste(term$level, collapse = ", "),
          "):units conditions the residual on selected levels of ",
          term$var,
          ", which is a different model from a heterogeneous residual ",
          "across all of them. Drop the level to fit dsum(~ units | ",
          term$var,
          ")."
        )
      ))
    }

    # A residual variance only exists for a family that HAS one. Poisson
    # and binomial carry no sigma, so `sigma ~ f` would be a predictor on
    # a parameter the model does not have -- brms would error, or worse,
    # accept it for a family where sigma means something else.
    fam <- tolower(as.character(fb$family %||% "gaussian"))
    if (!.fb_family_has_brms_sigma(fam)) {
      stop(.fb_refusal_condition(
        reason_code = "heterogeneous_residual_family_has_no_sigma",
        message = paste0(
          "A heterogeneous residual variance was requested for family '",
          fam,
          "', which has no residual scale parameter to vary -- its ",
          "dispersion is a function of the mean. Model the dispersion ",
          "directly for that family, or fit a Gaussian on a transformed ",
          "response."
        )
      ))
    }
    sigma_rhs <- paste0("0 + ", term$var)
  } else {
    sigma_rhs <- "1"
  }
  if (has_weight_offset) {
    sigma_rhs <- paste0(
      sigma_rhs, " + offset(", .FB_BRMS_WEIGHTS_OFFSET_COL, ")"
    )
  }

  brms::bf(main, stats::as.formula(paste("sigma ~", sigma_rhs)))
}

# .brms_term_refusal_message() --- what to tell a user whose random term has
# no branch in the reconstruction above.
#
# The reconstruction rebuilds the fixed and random blocks from the IR, so a
# term it does not recognise would simply not appear in the emitted formula:
# a model missing a component, returned without complaint. The refusal is
# therefore mandatory, and it names the term type and the route that does
# fit rather than a bare "not supported".
.brms_term_refusal_message <- function(term) {
  ttype <- term$type %||% "<unknown>"
  var_name <- term$var %||% "x"
  route <- switch(
    ttype,
    "spline" = paste0(
      "A penalised spline is emitted by INLA as a second-order random ",
      "walk, so fit spl(",
      var_name,
      ") with backend = \"inla\"."
    ),
    "smooth_mgcv" = paste0(
      "An mgcv smooth basis has no brms lowering in the ASReml grammar, ",
      "and no active engine fits s(",
      var_name,
      ") at all: INLA's rw2 route is spl(",
      var_name,
      "), a different term type built on an ASReml-style basis, not a ",
      "lowering of the mgcv smooth. Use spl(",
      var_name,
      ") with backend = \"inla\" instead of s(",
      var_name,
      ")."
    ),
    "polynomial" = paste0(
      "Fit pol(",
      var_name,
      ") with backend = \"inla\", or write the ",
      "polynomial terms in the fixed part of the formula."
    ),
    "continuous" = paste0(
      "A continuous variable in the random part is a random regression. ",
      "Write the slope as (",
      var_name,
      " || <grouping factor>) in the ",
      "brms-style grammar."
    ),
    "at" = paste0(
      "Write the heterogeneous variance as diag(",
      var_name,
      "):<inner factor>, which brms emits as (0 + ",
      var_name,
      " || <inner factor>)."
    ),
    "us" = paste0(
      "Write the unstructured covariance as us(",
      var_name,
      "):<inner factor>, which brms emits as (0 + ",
      var_name,
      " | <inner factor>)."
    ),
    "ar1" = ,
    "ar1_spatial" = paste0(
      "An autoregressive latent field is emitted by INLA only. Fit it ",
      "with backend = \"inla\" on a complete grid with one observation ",
      "per node."
    ),
    "vm_gxe" = paste0(
      "A known-covariance term crossed with a second factor has no brms ",
      "lowering. Fit vm(",
      term$inner %||% var_name,
      ", K) on its own, or use ASReml for the crossed form."
    ),
    paste0(
      "No active engine emits this structure. fb_terms() will describe ",
      "what was parsed, and the capability table in the README lists what ",
      "each engine fits."
    )
  )
  paste0(
    "brms cannot represent the random term type \"",
    ttype,
    "\": the brms formula reconstruction has no lowering for it, and ",
    "emitting the model without the term would answer a different ",
    "question. ",
    route
  )
}

# .fb_brms_missing_response_mode() --- decide how brms should treat a
# missing response, and refuse where it cannot treat it faithfully.
#
# brms drops rows whose response is NA. Silently: a 48-row dataset with six
# missing responses reaches Stan as N = 42, with a warning that is easy to
# miss and no effect on the returned object. That makes
# `na_action = "augment"` on brms mean complete-case deletion -- the
# argument promises the design is preserved and delivers the opposite.
#
# The repair depends on the family, because brms's `mi()` addition term is
# Gaussian-only. `bf(y | mi() ~ ...)` with family = poisson raises
# "Argument 'mi' is not supported for family 'poisson(log)'", and the same
# for bernoulli. So:
#
#   gaussian     -> `y | mi() ~ ...`, and brms samples the missing responses.
#   non-Gaussian -> refuse, and point at INLA, which carries a missing
#                   response as a native prediction target for every family.
#
# Refusing is not a limitation being papered over. Under ignorability the
# parameter posterior is the same whether the cell is imputed or omitted, so
# a user who wants the non-Gaussian fit can pass na_action = "omit" -- what
# they cannot have is deletion wearing the name of augmentation. The
# identity is a statement about the posterior, not about what an optimiser
# returns from a ragged design; see the note in R/na_action.R.
.fb_brms_missing_response_mode <- function(fb, data) {
  resp <- fb$response
  has_missing <- !is.null(resp) &&
    nzchar(resp) &&
    resp %in% names(data) &&
    anyNA(data[[resp]])
  if (!has_missing) {
    return("none")
  }
  fam <- tolower(as.character(fb$family %||% "gaussian"))
  if (fam %in% c("gaussian", "lognormal", "student", "t")) {
    return("mi")
  }
  entry <- .lookup_refusal("brms_cannot_augment_nongaussian")
  stop(.fb_refusal_condition(
    reason_code = "brms_cannot_augment_nongaussian",
    message = if (!is.null(entry)) {
      entry$message_template
    } else {
      paste0(
        "backend = \"brms\" cannot carry a missing response for family \"",
        fam,
        "\"."
      )
    },
    family_class = "flexybayes_brms_cannot_augment_nongaussian",
    backend = "brms"
  ))
}

# Resolve the covariance-matrix symbol a vm() / ped() term references on
# the brms route. brms reaches GBLUP / pedigree BLUP only through an
# exact dense-able carrier (dense / chol / precision); the block-
# diagonal and low-rank carriers are INLA-only and are refused loudly
# rather than mis-rendered.
.fb_brms_covname <- function(term) {
  cov <- term$cov_representation
  fmt <- (cov$format %||% "dense")
  if (fmt %in% c("blocks", "low_rank")) {
    stop(
      "emit_brms() supports vm() / ped() only with an exact dense-able ",
      "covariance carrier (dense / chol / precision); the \"",
      fmt,
      "\" carrier is INLA-only. Re-route via backend = \"inla\".",
      call. = FALSE
    )
  }
  sym <- cov$data %||% term$mat
  if (is.null(sym) || is.na(sym) || !nzchar(sym)) {
    stop(
      "vm() / ped() term carries no resolvable covariance-matrix symbol ",
      "for the brms route.",
      call. = FALSE
    )
  }
  sym
}

# Build brms's data2 list: the materialised dense covariance for every
# vm() / ped() relationship term, keyed by the carrier symbol the
# gr(<var>, cov = <symbol>) formula references. brms wants a covariance
# matrix (it Cholesky-factors it and scales by the estimated group sd),
# so each carrier is densified to K -- dense as-is, chol -> L L',
# precision -> solve(Q). When the supplied matrix has no dimnames they
# are set from the grouping factor's levels so brms aligns the matrix
# positionally (the same level-alignment contract the INLA
# path enforces); existing dimnames are preserved so brms aligns by
# name.
.fb_brms_data2 <- function(fb, known_matrices, data) {
  d2 <- list()
  for (term in fb$random_terms) {
    if (!(term$type %||% "") %in% c("vm", "ped")) {
      next
    }
    covname <- .fb_brms_covname(term)
    cov <- term$cov_representation
    fmt <- (cov$format %||% "dense")
    sym <- cov$data %||% term$mat
    M <- known_matrices[[sym]]
    if (is.null(M)) {
      stop(
        "emit_brms(): the covariance matrix '",
        sym,
        "' for vm(",
        term$var,
        ") is not in known_matrices. Pass it via known_matrices = list(",
        sym,
        " = <your relationship matrix>).",
        call. = FALSE
      )
    }
    K <- switch(
      fmt,
      "dense" = as.matrix(M),
      "chol" = tcrossprod(as.matrix(M)),
      "precision" = solve(as.matrix(M)),
      "pedigree_sparse_precision" = as.matrix(solve(M)),
      stop(
        "emit_brms(): unsupported covariance carrier '",
        fmt,
        "' for the brms vm() / ped() route.",
        call. = FALSE
      )
    )
    if (is.null(rownames(K)) || is.null(colnames(K))) {
      lev <- levels(as.factor(data[[term$var]]))
      if (nrow(K) == length(lev)) {
        dimnames(K) <- list(lev, lev)
      }
    }
    d2[[covname]] <- K
  }
  d2
}

# Build a fixed-only formula from the IR's fixed_terms
# slot. Used by emit_brms()'s model_info$n_fixed calculation --
# stats::model.matrix() rejects the lme4 double-pipe `(x || g)`
# random-effects shorthand (`||` requires logical(1) operands), so
# we cannot feed the full brms_form into model.matrix(). The
# fixed-only formula mirrors .fb_to_brms_formula() minus the
# random-effects block.
.fb_to_fixed_only_formula <- function(fb) {
  response <- fb$response
  fixed_labels <- vapply(
    fb$fixed_terms,
    function(t) {
      if (!is.null(t$label)) {
        t$label
      } else if (!is.null(t$var)) {
        t$var
      } else if (!is.null(t$expr)) {
        deparse(t$expr)
      } else {
        NA_character_
      }
    },
    character(1)
  )
  fixed_labels <- fixed_labels[!is.na(fixed_labels)]

  fixed_rhs <- if (length(fixed_labels) > 0L) {
    paste(fixed_labels, collapse = " + ")
  } else {
    "1"
  }
  if (!isTRUE(fb$intercept)) {
    fixed_rhs <- paste0(fixed_rhs, " - 1")
  }
  stats::as.formula(paste(response, "~", fixed_rhs))
}

# Translate fb's family + link into a brms family. Conventions:
# Gaussian / Poisson / Gamma use the base stats::family() call (brms
# does not re-export these). Binomial single-column maps to
# brms::bernoulli() (the corpus refuses `y | trials(n)`, so binomial
# is single-trial). Lognormal / Beta / negbinomial are brms-defined
# families pulled by name from the brms namespace.
.fb_family_to_brms <- function(fam, link) {
  fam <- tolower(as.character(fam %||% "gaussian"))
  link_clean <- if (!is.null(link)) tolower(link) else NULL

  ns <- asNamespace("brms")

  if (fam == "gaussian") {
    return(
      if (is.null(link_clean)) {
        stats::gaussian()
      } else {
        stats::gaussian(link = link_clean)
      }
    )
  }
  if (fam %in% c("binomial", "binary", "bernoulli")) {
    # Single-column 0/1 corpus shape -> Bernoulli logit.
    return(
      if (is.null(link_clean)) {
        get("bernoulli", envir = ns)()
      } else {
        get("bernoulli", envir = ns)(link = link_clean)
      }
    )
  }
  if (fam == "poisson") {
    return(
      if (is.null(link_clean)) {
        stats::poisson()
      } else {
        stats::poisson(link = link_clean)
      }
    )
  }
  if (fam %in% c("negbin", "negbinom", "negative_binomial", "nbinomial")) {
    return(
      if (is.null(link_clean)) {
        get("negbinomial", envir = ns)()
      } else {
        get("negbinomial", envir = ns)(link = link_clean)
      }
    )
  }
  if (fam == "gamma") {
    return(
      if (is.null(link_clean)) {
        stats::Gamma(link = "log")
      } else {
        stats::Gamma(link = link_clean)
      }
    )
  }
  if (fam == "lognormal") {
    return(get("lognormal", envir = ns)())
  }
  if (fam == "beta") {
    return(get("Beta", envir = ns)())
  }
  if (fam == "hurdle_gamma") {
    # brms-native: `brms::brmsfamily("hurdle_gamma")` declares
    # dpars mu, shape, hu on brms 2.23.0. The mean model carries the
    # link; `hu` (the zero-mass probability) keeps brms's own logit link
    # and its own prior, which is what `.fb_default_prior_targets()`
    # records as an engine default rather than a flexyBayes choice.
    return(
      if (is.null(link_clean)) {
        get("hurdle_gamma", envir = ns)()
      } else {
        get("hurdle_gamma", envir = ns)(link = link_clean)
      }
    )
  }
  stop(
    "emit_brms() does not yet translate family = \"",
    fam,
    "\". Supported families: gaussian, binomial (single-column ",
    "Bernoulli), poisson, negative_binomial, gamma, lognormal, ",
    "beta, hurdle_gamma. Other families are deferred to a future ",
    "release.",
    call. = FALSE
  )
}

.deparse_brms_family <- function(fam) {
  if (inherits(fam, "family") || inherits(fam, "brmsfamily")) {
    nm <- fam$family %||% "<unknown>"
    lk <- fam$link %||% "default"
    return(paste0(nm, " (", lk, " link)"))
  }
  as.character(fam)
}


# ---------------------------------------------------------------- #
# GLM shim + convergence + variance components                       #
# ---------------------------------------------------------------- #

# .brms_dpar_names() --- the distributional parameters this fit models
# with their own linear predictor.
#
# A sectioned residual emits `sigma ~ 0 + f`, which brms records under
# `$formula$pforms`. Empty on every ordinary fit.
#
# @noRd
# @keywords internal
.brms_dpar_names <- function(brmsfit) {
  nm <- names(brmsfit$formula$pforms %||% list())
  if (is.null(nm)) {
    return(character(0L))
  }
  nm[nzchar(nm)]
}

# .brms_mean_b_cols() --- the population-level columns of the MEAN model.
#
# brms writes every linear-predictor coefficient with a `b_` prefix, and
# a distributional parameter's coefficients as `b_<dpar>_<term>`. Sweeping
# `^b_` therefore collects the log-sigma coefficients of a sectioned
# residual alongside the mean effects, which is how `coef()` came to
# report `sigma_EnvE1` as if it were a fixed effect on the response, and
# how the fixed-effect design matrix (6 columns) came to be reconciled
# against a 12-element coefficient basis -- the failure a caller met as an
# estimability error from `predict(classify = )` on any sectioned-residual
# fit, including for a plain fixed factor.
#
# @noRd
# @keywords internal
.brms_mean_b_cols <- function(cn, dpars) {
  b_cols <- cn[grepl("^b_", cn)]
  if (length(dpars) == 0L || length(b_cols) == 0L) {
    return(b_cols)
  }
  prefixes <- paste0("b_", dpars, "_")
  keep <- !Reduce(
    `|`,
    lapply(prefixes, function(p) startsWith(b_cols, p))
  )
  b_cols[keep]
}

# Build a $glm shim that satisfies coef / vcov / fitted / residuals /
# family / formula / model.matrix on the brms path. Coefficient
# names follow the canonical (brms-stripped) convention: b_<term>
# rows from the posterior are exposed as <term>; the Intercept row
# becomes "(Intercept)".
.brms_glm_shim <- function(brmsfit, draws_mat, data, family, formula_used) {
  cn <- colnames(draws_mat)

  b_cols <- .brms_mean_b_cols(cn, .brms_dpar_names(brmsfit))
  canon_names <- vapply(
    b_cols,
    function(nm) {
      bare <- sub("^b_", "", nm)
      if (identical(bare, "Intercept")) "(Intercept)" else bare
    },
    character(1),
    USE.NAMES = FALSE
  )

  if (length(b_cols) > 0L) {
    fixed_draws <- draws_mat[, b_cols, drop = FALSE]
    coefs <- colMeans(fixed_draws)
    vcov_mx <- stats::cov(fixed_draws)
    rownames(vcov_mx) <- colnames(vcov_mx) <- canon_names
    names(coefs) <- canon_names
    # The draws travel with the summary statistics so summary() and
    # confint() can report a posterior quantile interval and a
    # probability of direction counted from the draws, rather than
    # normal approximations to both presented under the exact
    # quantities' names.
    colnames(fixed_draws) <- canon_names
  } else {
    coefs <- numeric(0)
    vcov_mx <- matrix(numeric(0), 0L, 0L)
    fixed_draws <- NULL
  }

  # fitted values: dispatch through the stats::fitted generic to the
  # brms-registered method (brms registers fitted.brmsfit via
  # NAMESPACE; the S3 dispatch path is portable across brms versions
  # whereas brms::fitted is not always re-exported). Wrap in tryCatch
  # so a brms version skew (e.g., signature change) does not break
  # construction; fall back to NA placeholders that the parent
  # methods will surface explicitly.
  fitted_vals <- tryCatch(
    as.numeric(stats::fitted(brmsfit, summary = TRUE)[, "Estimate"]),
    error = function(e) rep(NA_real_, nrow(data))
  )
  linear_pred <- tryCatch(
    as.numeric(stats::fitted(brmsfit, scale = "linear", summary = TRUE)[,
      "Estimate"
    ]),
    error = function(e) fitted_vals
  )
  response_col <- all.vars(.brms_main_formula(formula_used)[[2L]])[1L]
  y_vec <- if (!is.null(response_col) && response_col %in% names(data)) {
    as.numeric(data[[response_col]])
  } else {
    rep(NA_real_, nrow(data))
  }
  residual_vec <- y_vec - fitted_vals

  fixed_formula <- .brms_fixed_only_formula(formula_used)

  glm_obj <- structure(
    list(
      coefficients = coefs,
      fitted.values = fitted_vals,
      linear.predictors = linear_pred,
      residuals = residual_vec,
      family = family,
      formula = fixed_formula,
      data = data,
      y = y_vec
    ),
    class = c("flexybayes_glm", "list"),
    posterior_vcov = vcov_mx,
    posterior_draws = fixed_draws
  )
  glm_obj
}

# .brms_main_formula() --- the mean-model formula, whatever wrapper it
# arrived in.
#
# A sectioned residual makes the emitted model distributional, and the
# object handed down the chain is then a `brmsformula` rather than a
# plain one. Indexing a `brmsformula` by position does not reach the
# formula: `form[[2L]]` and `form[[3L]]` return its `pforms` and `pfix`
# slots, both of which are lists. Every positional reader below therefore
# has to unwrap first, or it silently works on the wrong object -- which
# is what left `glm$formula` carrying its random-effect bars, `glm$y` and
# `glm$residuals` all-NA, and `predict(classify = )` dying inside the
# estimability seam with "Perhaps a 'data' or 'params' argument is
# needed" on every sectioned-residual fit.
#
# @noRd
# @keywords internal
.brms_main_formula <- function(form) {
  if (inherits(form, "brmsformula") && !is.null(form$formula)) {
    return(form$formula)
  }
  form
}

# Drop the random-effect terms from a brms formula so the parent
# model.matrix / predict path can use a plain fixed-effect formula.
.brms_fixed_only_formula <- function(form) {
  form <- .brms_main_formula(form)
  rhs <- form[[3L]]
  no_re <- .strip_re_terms(rhs)
  if (is.null(no_re)) {
    no_re <- quote(1)
  }
  out <- form
  out[[3L]] <- no_re
  out
}

.strip_re_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("+"))) {
    l <- .strip_re_terms(expr[[2L]])
    r <- .strip_re_terms(expr[[3L]])
    if (is.null(l)) {
      return(r)
    }
    if (is.null(r)) {
      return(l)
    }
    return(call("+", l, r))
  }
  if (is.call(expr) && identical(expr[[1L]], as.name("|"))) {
    return(NULL)
  }
  if (is.call(expr) && identical(expr[[1L]], as.name("||"))) {
    return(NULL)
  }
  if (is.call(expr) && identical(expr[[1L]], as.name("("))) {
    inner <- .strip_re_terms(expr[[2L]])
    if (is.null(inner)) {
      return(NULL)
    }
    return(call("(", inner))
  }
  expr
}

# Build the convergence-diag shape print.flexybayes expects.
# gelman$psrf is a matrix with rownames = parameter, columns
# including "Point est."; we synthesise it from posterior::rhat.
# n_eff is a numeric vector of ess_bulk values keyed by parameter.
.brms_convergence <- function(brmsfit, draws_mat) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(list(gelman = NULL, n_eff = NULL))
  }

  # Tail ESS is recorded alongside bulk because they answer different
  # questions: bulk covers the centre of the marginal, tail the quantiles
  # a credible interval is read off. A fit can mix well in the middle and
  # badly at the 2.5% bound, which is precisely the number a user quotes.
  summ <- tryCatch(
    posterior::summarise_draws(
      posterior::as_draws_array(brmsfit),
      "rhat",
      "ess_bulk",
      "ess_tail"
    ),
    error = function(e) NULL
  )
  if (is.null(summ)) {
    return(list(gelman = NULL, n_eff = NULL))
  }

  psrf <- matrix(
    c(summ$rhat, rep(NA_real_, nrow(summ))),
    nrow = nrow(summ),
    ncol = 2L
  )
  colnames(psrf) <- c("Point est.", "Upper C.I.")
  rownames(psrf) <- summ$variable

  # Divergent transitions are a separate signal from R-hat and ESS: a
  # chain can mix well over a region it never reached. Counted here so the
  # fit carries it, since triangulate()'s diagnostics gate reads the fit
  # rather than re-deriving the sampler state. NA when the sampler does
  # not report NUTS parameters.
  n_div <- tryCatch(
    {
      np <- brms::nuts_params(brmsfit)
      sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
    },
    error = function(e) NA_integer_
  )

  list(
    gelman = list(psrf = psrf),
    n_eff = stats::setNames(summ$ess_bulk, summ$variable),
    n_eff_tail = stats::setNames(summ$ess_tail, summ$variable),
    n_divergent = n_div
  )
}

# Build the variance_comps table -- same shape as the INLA path
# (component, estimate, sd, q2.5, q97.5).
.brms_variance_comps <- function(brmsfit, draws_mat) {
  cn <- colnames(draws_mat)
  vc_cols <- cn[grepl("^sd_", cn) | cn == "sigma"]
  if (length(vc_cols) == 0L) {
    return(data.frame(
      component = character(0),
      estimate = numeric(0),
      sd = numeric(0),
      q2.5 = numeric(0),
      q97.5 = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(vc_cols, function(nm) {
    x <- draws_mat[, nm]
    qs <- stats::quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    list(
      component = .brms_vc_canonical_name(nm),
      estimate = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      q2.5 = qs[1L],
      q97.5 = qs[2L]
    )
  })
  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))

  # Posterior medians travel as an attribute rather than a sixth column:
  # the five column names are a broom contract (R/tidiers.R). The
  # boundary-collapse display flag in summary() reads them, and computes
  # them from the same draws the row itself was summarised from.
  attr(out, "posterior_median") <- stats::setNames(
    vapply(
      vc_cols,
      function(nm) stats::median(draws_mat[, nm], na.rm = TRUE),
      numeric(1L)
    ),
    vapply(vc_cols, .brms_vc_canonical_name, character(1L))
  )
  out
}

# brms VC names: sd_<group>__Intercept -> sd_<group>; sigma stays.
.brms_vc_canonical_name <- function(nm) {
  if (identical(nm, "sigma")) {
    return("sigma")
  }
  bare <- sub("^sd_", "", nm)
  bare <- sub("__Intercept$", "", bare)
  paste0("sd_", bare)
}


# ---------------------------------------------------------------- #
# Per-level residual variances on a sectioned residual              #
# ---------------------------------------------------------------- #
#
# A sectioned residual leaves no scalar `sigma` behind, so the variance-
# component table above -- which reads `sd_*` rows and `sigma` -- reports
# nothing at all for the residual. What the model carries instead is a
# vector of coefficients on LOG sigma, one per level of the sectioning
# factor. A reader who takes `b_sigma_envE2` for an ASReml `env_2!R` is
# off by a logarithm and a square.
#
# The block below closes that gap on the surface a user actually reads.
# Both quantities are summarised FROM THE DRAWS of exp(b) and exp(2 * b),
# not by transforming a posterior mean: exp() is convex, so
# exp(2 * mean(b)) is a different number from mean(exp(2 * b)) and the
# gap grows with the coefficient's posterior spread. The point summary is
# the median because it is the one summary the monotone transform carries
# through exactly, and the interval is the transformed quantile pair for
# the same reason.

# The sectioned-residual terms carried by a fit's IR, if any.
.brms_at_units_terms <- function(object) {
  Filter(
    function(t) identical(t$type %||% "", "at_units"),
    object$extras$parse_info$residual %||% list()
  )
}

# One row per level of the sectioning factor: the residual variance and
# the residual standard deviation, each as a posterior median and a 95%
# credible interval on the natural scale.
#
# Returns NULL when the fit carries no sectioned residual, when the draws
# are unavailable, or when the emitted coefficients cannot be matched to
# the factor -- a partial table would be worse than none, because the
# level labels are the whole point.
.brms_residual_by_level_table <- function(object, probs = c(0.025, 0.975)) {
  terms <- .brms_at_units_terms(object)
  if (length(terms) != 1L || is.null(object$brms)) {
    return(NULL)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(NULL)
  }
  var_name <- terms[[1L]]$var %||% NA_character_
  if (is.na(var_name)) {
    return(NULL)
  }

  draws <- tryCatch(
    as.matrix(posterior::as_draws_matrix(object$brms)),
    error = function(e) NULL
  )
  if (is.null(draws)) {
    return(NULL)
  }
  cols <- grep("^b_sigma_", colnames(draws), value = TRUE)
  if (length(cols) == 0L) {
    return(NULL)
  }

  # `sigma ~ 0 + f` names each coefficient <factor><level>, so the level
  # is what remains once the factor's own name is removed. Ordering
  # follows the data's factor levels where they are recoverable, so the
  # printed block reads in the same order as every other per-level
  # summary the user has seen.
  labels <- sub(paste0("^b_sigma_", var_name), "", cols)
  declared <- tryCatch(
    levels(object$brms$data[[var_name]]),
    error = function(e) NULL
  )
  if (!is.null(declared) && setequal(labels, declared)) {
    ord <- match(declared, labels)
    cols <- cols[ord]
    labels <- labels[ord]
  }

  summarise <- function(x) {
    q <- stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)
    c(stats::median(x, na.rm = TRUE), q[1L], q[2L])
  }
  rows <- lapply(seq_along(cols), function(i) {
    b <- draws[, cols[[i]]]
    v <- summarise(exp(2 * b))
    s <- summarise(exp(b))
    data.frame(
      level = labels[[i]],
      variance = v[[1L]],
      variance_lower = v[[2L]],
      variance_upper = v[[3L]],
      sd = s[[1L]],
      sd_lower = s[[2L]],
      sd_upper = s[[3L]],
      # The posterior spread of the SD itself, so the same rows can be
      # carried into summary()$varcomp, whose `std.error` column every
      # other variance component fills. The printed block does not show
      # it -- it selects its columns by name.
      sd_se = stats::sd(exp(b), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  attr(out, "factor") <- var_name
  attr(out, "probs") <- probs
  out
}

# .brms_weights_sigma_row() --- the single "sigma" row when weights (C6)
# put sigma on brms's distributional sub-formula without a heterogeneous
# residual alongside it.
#
# `sigma ~ 1 + offset(...)` (R/emit_brms.R's .fb_to_brms_formula()) has
# no scalar `sigma` draw either, for the same reason .brms_at_units_terms()
# fixes above -- it names the intercept coefficient `b_sigma_Intercept`
# instead, on the LOG scale, and OFFSET at 0 (i.e. weight = 1). Exp'd,
# that is the model's residual SD at unit weight -- the same "sigma at
# w = 1" convention INLA's summary.hyperpar and lme4::lmer(weights =)
# report (see .FB_BRMS_WEIGHTS_OFFSET_COL's banner for the full
# grounding trail against both). Median point summary and
# quantile-transformed interval for the same convexity reason as the
# per-level rows above.
#
# NULL whenever the fit carries no weights, carries a heterogeneous
# residual (that case's per-level rows already include the weight
# offset net of level, since the two compose additively on log sigma --
# .fb_varcomp_residual_by_level() covers it), or the Intercept
# coefficient cannot be found.
.brms_weights_sigma_row <- function(object, probs = c(0.025, 0.975)) {
  het <- .brms_at_units_terms(object)
  w <- object$extras$parse_info$weights
  if (length(het) > 0L || !.fb_weights_nonconstant(w) || is.null(object$brms)) {
    return(NULL)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(NULL)
  }
  draws <- tryCatch(
    as.matrix(posterior::as_draws_matrix(object$brms)),
    error = function(e) NULL
  )
  if (is.null(draws) || !("b_sigma_Intercept" %in% colnames(draws))) {
    return(NULL)
  }
  b <- draws[, "b_sigma_Intercept"]
  s <- exp(b)
  q <- stats::quantile(s, probs = probs, names = FALSE, na.rm = TRUE)
  data.frame(
    level = "sigma",
    variance = NA_real_,
    variance_lower = NA_real_,
    variance_upper = NA_real_,
    sd = stats::median(s, na.rm = TRUE),
    sd_lower = q[1L],
    sd_upper = q[2L],
    sd_se = stats::sd(s, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# Render the block. Shared by print() and summary() so the two cannot
# drift apart; a no-op on every fit without a sectioned residual.
.print_brms_residual_by_level <- function(object) {
  tab <- .brms_residual_by_level_table(object)
  if (is.null(tab)) {
    return(invisible(NULL))
  }
  probs <- attr(tab, "probs")
  cat(
    "\n-- Residual by level of `",
    attr(tab, "factor"),
    "` (posterior median, ",
    round(100 * diff(probs)),
    "% interval) ",
    strrep("-", 8),
    "\n",
    sep = ""
  )
  body <- data.frame(
    level = tab$level,
    variance = round(tab$variance, 4),
    lower = round(tab$variance_lower, 4),
    upper = round(tab$variance_upper, 4),
    SD = round(tab$sd, 4),
    SD.lower = round(tab$sd_lower, 4),
    SD.upper = round(tab$sd_upper, 4),
    stringsAsFactors = FALSE
  )
  print(body, row.names = FALSE)
  cat(
    "  Summarised from the draws of exp(2 * b_sigma_<level>); the fitted\n",
    "  coefficients are on log sigma, so neither they nor their posterior\n",
    "  means are variances. There is no scalar residual variance in this\n",
    "  model -- the sectioning replaces it.\n",
    sep = ""
  )
  invisible(tab)
}


# ---------------------------------------------------------------- #
# Subclass S3 overrides on the brms path                            #
# ---------------------------------------------------------------- #

#' Print method for the brms-passthrough flexybayes subclass
#'
#' Opens with the header every engine's print shares, then adds the
#' sampler diagnostics and a brms-specific footer (the live `brmsfit`
#' lives at `$brms`; the GLM shim at `$glm`; `$extras` carries the same
#' diagnostics as the INLA path).
#'
#' @param x A `flexybayes_brms` object.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the description it
#'   prints.
#' @export
print.flexybayes_brms <- function(x, ...) {
  mi <- x$extras$model_info

  .fb_print_header(x, "Bayesian mixed model", "-")

  cat(
    "  Params   :",
    mi$n_params,
    "monitored;",
    mi$n_fixed,
    "fixed,",
    mi$n_random,
    "random terms\n"
  )

  if (!is.null(x$extras$convergence$gelman)) {
    rhat <- x$extras$convergence$gelman$psrf[, "Point est."]
    if (length(rhat)) {
      mx <- max(rhat, na.rm = TRUE)
      flag <- if (mx < 1.05) {
        " [OK]"
      } else if (mx < 1.10) {
        " [borderline]"
      } else {
        " [!]"
      }
      cat("  Max Rhat:", round(mx, 3), flag, "\n")
    }
  }
  if (!is.null(x$extras$convergence$n_eff)) {
    mn <- min(x$extras$convergence$n_eff, na.rm = TRUE)
    if (is.finite(mn)) cat("  Min ESS:", round(mn, 0), "\n")
  }

  .print_brms_residual_by_level(x)

  cat(strrep("-", 55), "\n")
  cat("  $glm    -- GLM-compatible shim (coef, vcov, fitted, etc.)\n")
  cat("  $brms   -- live brmsfit (loo, posterior_predict, summary)\n")
  cat("  $extras -- diagnostics, variance components, call info\n")

  invisible(x)
}

#' Credible intervals on the brms path
#'
#' Uses the brms posterior draws directly (the parent
#' `confint.flexybayes` refuses unconditionally with
#' `fit_lacks_posterior_draws`, since no active engine reaches it
#' without its own override). Returns quantile-based credible
#' bounds over the `b_<term>` rows; row names are stripped of the
#' brms `b_` prefix to align with `coef()`.
#'
#' @param object A `flexybayes_brms` object.
#' @param parm Subset of fixed-effect names to return (NULL = all).
#' @param level Credible level (default 0.95).
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns A numeric matrix with one row per fixed-effect term and two
#'   columns holding the lower and upper credible bounds at `level`. Row
#'   names are the term names with the brms `b_` prefix removed.
#' @export
confint.flexybayes_brms <- function(object, parm = NULL, level = 0.95, ...) {
  .check_installed(
    "posterior",
    "Package 'posterior' is required for confint() on brms ",
    "fits."
  )
  draws <- as.matrix(posterior::as_draws_matrix(object$brms))
  cn <- colnames(draws)
  # Mean-model coefficients only, so this table and coef() / vcov()
  # describe the same parameter set on a distributional fit.
  b_cols <- .brms_mean_b_cols(cn, .brms_dpar_names(object$brms))
  if (length(b_cols) == 0L) {
    return(matrix(numeric(0), 0L, 2L))
  }

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  ci_mat <- t(apply(
    draws[, b_cols, drop = FALSE],
    2L,
    stats::quantile,
    probs = probs,
    names = FALSE,
    na.rm = TRUE
  ))
  rownames(ci_mat) <- vapply(
    b_cols,
    function(nm) {
      bare <- sub("^b_", "", nm)
      if (identical(bare, "Intercept")) "(Intercept)" else bare
    },
    character(1),
    USE.NAMES = FALSE
  )
  colnames(ci_mat) <- paste0(round(probs * 100, 1), "%")

  if (!is.null(parm)) {
    ci_mat <- ci_mat[parm, , drop = FALSE]
  }
  ci_mat
}

#' Predict from a brms-passthrough flexybayes fit
#'
#' Delegates to brms's `posterior_epred()` (response-scale posterior
#' mean) or `posterior_linpred()` (linear-predictor scale) on the
#' live `brmsfit` carried at `$brms`. The parent
#' `predict.flexybayes()` path uses a `$glm$linear.predictors`
#' point estimate that handles only the original-data case; this
#' subclass override accepts `newdata` and returns the posterior-mean
#' prediction (per-row mean over draws), or the full posterior
#' matrix when `summary = FALSE`.
#'
#' Population-level vs. group-level prediction follows brms's
#' `re_formula` convention: the default `re_formula = NULL` includes
#' all random effects; pass `re_formula = NA` for population-level
#' predictions only.
#'
#' @param object A `flexybayes_brms` object.
#' @param newdata Optional data.frame at which to predict. When
#'   omitted, returns the in-sample posterior summary.
#' @param type `"response"` (default; posterior_epred) or `"link"`
#'   (posterior_linpred).
#' @param re_formula Forwarded to brms; `NULL` (default) includes
#'   all random effects, `NA` excludes them (population-level).
#' @param se.fit Logical: if `TRUE`, returns a list with `fit`
#'   (posterior mean) and `se.fit` (posterior SD).
#' @param summary Logical: if `TRUE` (default), summarise across
#'   draws to a numeric vector; if `FALSE`, return the
#'   `draws x rows` posterior matrix.
#' @param classify The factors to break a marginal-means table down by:
#'   a character value (`"Variety"`, `"Variety:env"`) or a one-sided
#'   formula (`~ Variety`). `NULL` (the default) is the historical
#'   behaviour. See [predict.flexybayes_inla()] for how the two active
#'   engines' prediction paths differ.
#' @param level Credible level for the classify table's interval, as a
#'   proportion. Default `0.95`.
#' @param ... Forwarded to `brms::posterior_epred()` /
#'   `brms::posterior_linpred()`.
#' @returns One of four shapes, by argument. With `classify` set, the
#'   marginal-means table for those factors. Otherwise with
#'   `summary = FALSE`, the `draws x rows` posterior matrix; with
#'   `se.fit = TRUE`, a list of `fit` (posterior mean per row) and
#'   `se.fit` (posterior SD per row); and by default a numeric vector of
#'   posterior means, one per row of `newdata` or of the fitted data.
#' @export
predict.flexybayes_brms <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  re_formula = NULL,
  se.fit = FALSE,
  summary = TRUE,
  classify = NULL,
  level = 0.95,
  ...
) {
  type <- match.arg(type)
  if (!is.null(classify)) {
    .fb_classify_newdata_note(newdata)
    return(.fb_predict_classify(object, classify, level))
  }
  if (is.null(object$brms)) {
    stop(
      "Cannot predict from a flexybayes_brms object with an ",
      "empty $brms slot.",
      call. = FALSE
    )
  }

  # brms exports `posterior_epred` (response-scale) and
  # `posterior_linpred` (linear-predictor scale) as generics; the
  # `.brmsfit` methods dispatch automatically via the generic.
  # Pull the exported generic directly from the brms namespace so
  # the call is portable whether or not the user has loaded brms.
  .check_installed(
    "brms",
    "Package 'brms' is required for predict() on a ",
    "flexybayes_brms object."
  )
  brms_fn <- if (type == "link") {
    get("posterior_linpred", envir = asNamespace("brms"))
  } else {
    get("posterior_epred", envir = asNamespace("brms"))
  }

  post_mat <- if (is.null(newdata)) {
    brms_fn(object$brms, re_formula = re_formula, ...)
  } else {
    brms_fn(object$brms, newdata = newdata, re_formula = re_formula, ...)
  }

  if (!isTRUE(summary)) {
    return(post_mat)
  }

  fit_vec <- as.numeric(colMeans(post_mat))
  if (isTRUE(se.fit)) {
    se_vec <- as.numeric(apply(post_mat, 2L, stats::sd))
    return(list(fit = fit_vec, se.fit = se_vec))
  }
  fit_vec
}

#' Log-likelihood on the brms path
#'
#' Delegates to `brms::log_lik()` then sums pointwise log-likelihood
#' across observations and averages across draws. The `df` attribute
#' carries the parameter count from `$extras$model_info`; `nobs`
#' carries the observation count.
#'
#' @param object A `flexybayes_brms` object.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns An object of class `logLik`: the pointwise log-likelihood
#'   summed over observations and averaged across draws, carrying `df`
#'   (parameter count) and `nobs` (observation count) attributes.
#' @export
logLik.flexybayes_brms <- function(object, ...) {
  # Dispatch to brms's `log_lik` exported generic. brms ships
  # log_lik.brmsfit as an S3 method on the exported generic;
  # pulling the generic directly from the brms namespace works
  # whether or not the user has loaded brms.
  if (!requireNamespace("brms", quietly = TRUE)) {
    ll_val <- NA_real_
  } else {
    ll_fn <- get("log_lik", envir = asNamespace("brms"))
    ll_mat <- tryCatch(ll_fn(object$brms), error = function(e) NULL)
    ll_val <- if (is.null(ll_mat)) {
      NA_real_
    } else {
      mean(rowSums(ll_mat, na.rm = TRUE))
    }
  }

  structure(
    ll_val,
    df = object$extras$model_info$n_params,
    nobs = object$extras$model_info$n_obs,
    class = "logLik"
  )
}
