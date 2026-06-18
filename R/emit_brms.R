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
#       Methods that read $greta$draws (confint, logLik) are
#       overridden in this file to read brms posterior draws via
#       posterior::as_draws_matrix() instead.
#
#   (2) backend_decision() uniformity. $extras$backend_decision is
#       populated by .dispatch_backend() after this emit returns;
#       the slot's shape mirrors the greta / INLA path.
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
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  the_call = NULL,
  fixed = NULL,
  random = NULL,
  rcov = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_,
  ...
) {
  if (!is_fb_terms(fb)) {
    stop("`fb` must be an fb_terms object.", call. = FALSE)
  }

  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for backend = \"brms\". ",
      "Install via:\n  install.packages('brms')\n",
      "A working C++ toolchain (rstan or cmdstanr) is required ",
      "for Stan compilation.",
      call. = FALSE
    )
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package 'posterior' is required for the brms posterior ",
      "extraction shim. Install via install.packages('posterior').",
      call. = FALSE
    )
  }

  # -- formula + family + prior ------------------------------------
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
        rcov = fb$rcov_terms,
        family = list(family = fb$family, link = fb$link)
      ),
      call_info = list(
        fixed = fixed,
        random = random,
        rcov = rcov,
        data_name = data_name,
        family = family,
        link = link,
        known_matrices = known_matrices,
        weights = weights,
        n_samples = n_samples,
        warmup = warmup,
        chains = chains,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd
      ),
      run_time = elapsed,
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

  structure(
    list(
      glm = glm_obj,
      brms = brmsfit,
      extras = extras
    ),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
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

    if (
      !ttype %in% c("simple", "ide", "id", "simple_slope_uncor")
    ) {
      stop(
        "emit_brms() does not yet support random term type \"",
        ttype,
        "\". Re-route via backend = \"greta\".",
        call. = FALSE
      )
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

  stats::as.formula(paste(response, "~", rhs))
}

# Resolve the covariance-matrix symbol a vm() / ped() term references on
# the brms route. brms reaches GBLUP / pedigree BLUP only through an
# exact dense-able carrier (dense / chol / precision); the block-
# diagonal and low-rank carriers are greta / INLA-only and are refused
# loudly rather than mis-rendered.
.fb_brms_covname <- function(term) {
  cov <- term$cov_representation
  fmt <- (cov$format %||% "dense")
  if (fmt %in% c("blocks", "low_rank")) {
    stop(
      "emit_brms() supports vm() / ped() only with an exact dense-able ",
      "covariance carrier (dense / chol / precision); the \"", fmt,
      "\" carrier is greta / INLA-only. Re-route via backend = \"greta\".",
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
# positionally (the same level-alignment contract the greta / INLA
# paths enforce); existing dimnames are preserved so brms aligns by
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
        "emit_brms(): the covariance matrix '", sym, "' for vm(", term$var,
        ") is not in known_matrices. Pass it via known_matrices = list(",
        sym, " = <your relationship matrix>).",
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
        "emit_brms(): unsupported covariance carrier '", fmt,
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
  stop(
    "emit_brms() does not yet translate family = \"",
    fam,
    "\". Supported families: gaussian, binomial (single-column ",
    "Bernoulli), poisson, negative_binomial, gamma, lognormal, ",
    "beta. Other families are deferred to a future release.",
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

# Build a $glm shim that satisfies coef / vcov / fitted / residuals /
# family / formula / model.matrix on the brms path. Coefficient
# names follow the canonical (brms-stripped) convention: b_<term>
# rows from the posterior are exposed as <term>; the Intercept row
# becomes "(Intercept)".
.brms_glm_shim <- function(brmsfit, draws_mat, data, family, formula_used) {
  cn <- colnames(draws_mat)

  b_cols <- cn[grepl("^b_", cn)]
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
    coefs <- colMeans(draws_mat[, b_cols, drop = FALSE])
    vcov_mx <- stats::cov(draws_mat[, b_cols, drop = FALSE])
    rownames(vcov_mx) <- colnames(vcov_mx) <- canon_names
    names(coefs) <- canon_names
  } else {
    coefs <- numeric(0)
    vcov_mx <- matrix(numeric(0), 0L, 0L)
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
  response_col <- all.vars(formula_used[[2L]])[1L]
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
    posterior_vcov = vcov_mx
  )
  glm_obj
}

# Drop the random-effect terms from a brms formula so the parent
# model.matrix / predict path can use a plain fixed-effect formula.
.brms_fixed_only_formula <- function(form) {
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

  summ <- tryCatch(
    posterior::summarise_draws(
      posterior::as_draws_array(brmsfit),
      "rhat",
      "ess_bulk"
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

  list(
    gelman = list(psrf = psrf),
    n_eff = stats::setNames(summ$ess_bulk, summ$variable)
  )
}

# Build the variance_comps table -- same shape as the greta path
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
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
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
# Subclass S3 overrides on the brms path                            #
# ---------------------------------------------------------------- #

#' Print method for the brms-passthrough flexybayes subclass
#'
#' Mirrors `print.flexybayes` (call info + run time + diagnostics)
#' with a brms-specific footer (the live `brmsfit` lives at
#' `$brms`; the GLM shim at `$glm`; `$extras` carries the same
#' diagnostics as the greta path).
#'
#' @param x A `flexybayes_brms` object.
#' @param ... Ignored.
#' @export
print.flexybayes_brms <- function(x, ...) {
  ci <- x$extras$call_info
  mi <- x$extras$model_info

  cat("Bayesian mixed model  [flexyBayes / brms passthrough]\n")
  cat(strrep("-", 55), "\n")
  cat("  Fixed  :", deparse(ci$fixed), "\n")
  if (!is.null(ci$random)) {
    cat("  Random :", deparse(ci$random), "\n")
  }
  cat("  Family :", mi$family, "(", mi$link, "link )\n")

  nch <- ci$chains
  ns <- ci$n_samples
  cat(
    "  MCMC   :",
    nch,
    "chain(s) x",
    ns,
    "samples",
    "(warmup =",
    ci$warmup,
    ") --",
    round(x$extras$run_time, 1),
    "sec (Stan via brms)\n"
  )

  cat(
    "  Params :",
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

  cat(strrep("-", 55), "\n")
  cat("  $glm    -- GLM-compatible shim (coef, vcov, fitted, etc.)\n")
  cat("  $brms   -- live brmsfit (loo, posterior_predict, summary)\n")
  cat("  $extras -- diagnostics, variance components, call info\n")

  invisible(x)
}

#' Credible intervals on the brms path
#'
#' Uses the brms posterior draws directly (the parent
#' `confint.flexybayes` reads `$greta$draws`, which is `NULL` on
#' the brms-passthrough path). Returns quantile-based credible
#' bounds over the b_<term> rows; row names are stripped of the
#' brms `b_` prefix to align with `coef()`.
#'
#' @param object A `flexybayes_brms` object.
#' @param parm Subset of fixed-effect names to return (NULL = all).
#' @param level Credible level (default 0.95).
#' @param ... Ignored.
#' @export
confint.flexybayes_brms <- function(object, parm = NULL, level = 0.95, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package 'posterior' is required for confint() on brms ",
      "fits.",
      call. = FALSE
    )
  }
  draws <- as.matrix(posterior::as_draws_matrix(object$brms))
  cn <- colnames(draws)
  b_cols <- cn[grepl("^b_", cn)]
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
#' @param ... Forwarded to `brms::posterior_epred()` /
#'   `brms::posterior_linpred()`.
#' @export
predict.flexybayes_brms <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  re_formula = NULL,
  se.fit = FALSE,
  summary = TRUE,
  ...
) {
  type <- match.arg(type)
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
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for predict() on a ",
      "flexybayes_brms object.",
      call. = FALSE
    )
  }
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
#' @param ... Ignored.
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
