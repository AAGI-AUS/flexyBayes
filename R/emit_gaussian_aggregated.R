# emit_gaussian_aggregated -- aggregated backend wiring (v0.3.2).
#
# Aggregated emit path: gaussian-identity LMMs in scope route through the
# `<fb_aggregated>` sufficient-statistics object instead of the per-row
# likelihood. The aggregated form is algebraically identical to the
# per-row form (not an approximation).
#
# INLA path (the only active aggregated backend). Cell-mean weighted
# gaussian with a custom-prior
#     correction on the log-precision hyperparameter that absorbs the
#     sigma-dependent within-cell SS term:
#
#       prior_expr <- "expression:
#                       a = (WSS / 2) + 5e-5;
#                       b = (N - K) / 2 + 1;
#                       log_dens = -a * exp(theta) + b * theta;
#                       return(log_dens);"
#
#     The constants `a` and `b` collapse the default `loggamma(1, 5e-5)`
#     precision prior with the `Delta(theta)` correction. The combined
#     prior is fed via
#     `control.family$hyper$prec$prior` and INLA's deterministic
#     numerical engine recovers the per-row posterior on beta + sigma
#     to bit-exact precision (spike: differences <= 1e-4 on a 1000-row
#     example).
#
# Heterogeneous residual: when the IR's residual_terms include an
# `at_units` term (per-level residual variance, e.g. `at(env):units`),
# the cell-mean machinery extends naturally PROVIDED the residual-
# grouping factor is itself a cell key (so each cell has a single
# residual sigma). When the residual factor is NOT a cell key, the
# aggregation closure breaks and the emit refuses with reason code
# `heterogeneous_residual_factor_not_in_cell_key`.
#
# Internal -- not exported. Called from .dispatch_backend() when the
# `<fb_aggregation_plan>` declares eligibility.

emit_gaussian_aggregated <- function(
  fb,
  fb_aggregated,
  data,
  backend = "inla",
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
  residual = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_,
  known_matrices = NULL,
  weights = NULL,
  seed = NULL,
  control = NULL,
  na_action = NULL,
  requested_backend = NULL,
  requested_aggregate = NULL
) {
  backend <- match.arg(backend)

  .check_fb_terms(fb, "emit_gaussian_aggregated() requires an <fb_terms> IR.")
  .check_fb_aggregated(
    fb_aggregated,
    "emit_gaussian_aggregated() requires an <fb_aggregated> ",
    "from .fb_aggregate_gaussian()."
  )

  # Defensive scope re-check (the caller -- .dispatch_backend() --
  # already gates on <fb_aggregation_plan>$eligible, but the emit
  # path also documents its own invariants so a future direct caller
  # cannot misuse it).
  if (!identical(fb$family, "gaussian")) {
    stop(
      "emit_gaussian_aggregated(): family must be gaussian; got '",
      fb$family,
      "'.",
      call. = FALSE
    )
  }
  if (!(is.null(fb$link) || identical(fb$link, "identity"))) {
    stop(
      "emit_gaussian_aggregated(): link must be identity; got '",
      fb$link,
      "'.",
      call. = FALSE
    )
  }

  # Sufficient-statistics sanity gates.
  agg <- fb_aggregated$sufficient_stats
  if (any(agg$n_k <= 0L)) {
    stop(
      "emit_gaussian_aggregated(): non-positive n_k in cell ",
      which(agg$n_k <= 0L)[1L],
      ".",
      call. = FALSE
    )
  }
  wss_per_cell <- agg$S2_k - agg$S1_k^2 / agg$n_k
  # Allow tiny negative drift from floating-point summation.
  if (any(wss_per_cell < -1e-8 * abs(agg$S2_k))) {
    stop(
      "emit_gaussian_aggregated(): negative within-cell SS in cell ",
      which.min(wss_per_cell),
      " (likely numerical instability ",
      "upstream).",
      call. = FALSE
    )
  }
  wss_per_cell <- pmax(wss_per_cell, 0)
  WSS_total <- sum(wss_per_cell)

  # Heterogeneous-residual handling: detect residual_terms shape and
  # validate that any at_units residual factor lies in the cell key.
  residual_plan <- .agg_emit_residual_plan(fb, fb_aggregated)

  # Build per-RI level indices for the cell key (one per random_cols
  # entry; foundation already places each RI grouping factor into the
  # cell key as a column on the sufficient_stats data.table).
  ri_plan <- .agg_emit_ri_plan(fb, fb_aggregated)

  if (return_code) {
    # The aggregated emit at v0.3.2 does NOT yet expose review_code =
    # TRUE for the aggregated path; the per-row review path remains the
    # documented v0.2 surface. Refused by name:
    stop(
      "emit_gaussian_aggregated(): return_code = TRUE is not yet ",
      "supported on the aggregated path (deferred to a follow-on ",
      "release). Pass aggregate = FALSE if you need the per-row ",
      "Stan source.",
      call. = FALSE
    )
  }

  t0 <- Sys.time()
  fit_engine_out <- .emit_gaussian_aggregated_inla(
    fb = fb,
    fb_aggregated = fb_aggregated,
    data = data,
    wss_per_cell = wss_per_cell,
    WSS_total = WSS_total,
    ri_plan = ri_plan,
    residual_plan = residual_plan,
    verbose = verbose
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  .agg_emit_build_fit(
    fb = fb,
    fb_aggregated = fb_aggregated,
    data = data,
    backend = backend,
    engine_out = fit_engine_out,
    ri_plan = ri_plan,
    residual_plan = residual_plan,
    elapsed = elapsed,
    the_call = the_call,
    fixed = fixed,
    random = random,
    residual = residual,
    family = family,
    link = link,
    data_name = data_name,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    known_matrices = known_matrices,
    weights = weights,
    seed = seed,
    control = control,
    na_action = na_action,
    verbose = verbose,
    requested_backend = requested_backend,
    requested_aggregate = requested_aggregate
  )
}


# ---------------------------------------------------------------- #
# Residual-structure planning                                       #
# ---------------------------------------------------------------- #
# Returns a list:
#   $kind        "homogeneous" | "at_units"
#   $factor      NULL (homogeneous) | <character> column name (at_units)
#   $level_col   NULL | <integer vector length K> mapping cell -> level
#                of $factor. Used to attach per-level scales in INLA.
#   $n_levels    NULL | integer count of distinct residual levels.
#   $levels      NULL | character vector of level labels.
#
# Refuses (structured error) when an at_units term references a factor
# that is NOT in the cell-key set: the cell-constant mu property still
# holds, but the cell-constant sigma property does not, so the per-row
# vs per-cell algebraic identity breaks. Reason code
# `heterogeneous_residual_factor_not_in_cell_key`.
.agg_emit_residual_plan <- function(fb, fb_aggregated) {
  residual <- fb$residual_terms
  if (length(residual) == 0L) {
    return(list(
      kind = "homogeneous",
      factor = NULL,
      level_col = NULL,
      n_levels = NULL,
      levels = NULL
    ))
  }

  homogeneous_types <- c("units", "id", "ide", "simple")
  hetero_types <- c("at_units")

  for (term in residual) {
    ttype <- term$type %||% "units"
    if (ttype %in% homogeneous_types) {
      next
    }
    if (ttype %in% hetero_types) {
      f <- as.character(term$var)
      if (!(f %in% fb_aggregated$cell_key_cols)) {
        stop(.fb_refusal_condition(
          reason_code = "heterogeneous_residual_factor_not_in_cell_key",
          message = paste0(
            "emit_gaussian_aggregated(): heterogeneous residual ",
            "at(",
            f,
            "):units refused -- '",
            f,
            "' is not in the ",
            "cell key {",
            paste(fb_aggregated$cell_key_cols, collapse = ", "),
            "}. The cell-constant sigma property does not hold, ",
            "so the per-row / per-cell algebraic identity breaks. ",
            "Pass aggregate = FALSE for the per-row path."
          ),
          family_class = "flexybayes_aggregate_emit_refusal",
          factor = f,
          cell_key = fb_aggregated$cell_key_cols
        ))
      }
      col_data <- fb_aggregated$sufficient_stats[[f]]
      if (!is.factor(col_data)) {
        col_data <- factor(col_data)
      }
      return(list(
        kind = "at_units",
        factor = f,
        level_col = as.integer(col_data),
        n_levels = nlevels(col_data),
        levels = levels(col_data)
      ))
    }
    # Anything else (us_units, fa_units, ar1_units, ...) is out of scope.
    stop(.fb_refusal_condition(
      reason_code = "residual_type_unsupported_for_aggregation",
      message = paste0(
        "emit_gaussian_aggregated(): residual term type '",
        ttype,
        "' is not supported by the aggregated path. Only homogeneous ",
        "(units / id) and at_units heterogeneous residual are ",
        "supported. Pass aggregate = FALSE for the per-row path."
      ),
      family_class = "flexybayes_aggregate_emit_refusal",
      residual_type = ttype
    ))
  }
  list(
    kind = "homogeneous",
    factor = NULL,
    level_col = NULL,
    n_levels = NULL,
    levels = NULL
  )
}


# Random-intercept planning. Returns a list of per-RI-term records:
#   $col          character: column name in sufficient_stats
#   $level_idx    integer vector length K: per-cell level index
#   $n_levels     integer: factor's full level count
.agg_emit_ri_plan <- function(fb, fb_aggregated) {
  cols <- fb_aggregated$random_cols
  if (!length(cols)) {
    return(list())
  }
  out <- vector("list", length(cols))
  for (i in seq_along(cols)) {
    g <- cols[i]
    col_data <- fb_aggregated$sufficient_stats[[g]]
    if (!is.factor(col_data)) {
      col_data <- factor(col_data)
    }
    out[[i]] <- list(
      col = g,
      level_idx = as.integer(col_data),
      n_levels = nlevels(col_data)
    )
  }
  out
}



# ---------------------------------------------------------------- #
# INLA path                                                         #
# ---------------------------------------------------------------- #
.emit_gaussian_aggregated_inla <- function(
  fb,
  fb_aggregated,
  data,
  wss_per_cell,
  WSS_total,
  ri_plan,
  residual_plan,
  verbose
) {
  .check_installed(
    "INLA",
    "Package 'INLA' is required for the INLA-backed aggregated ",
    "emit. Install from https://inla.r-inla-download.org/."
  )

  if (identical(residual_plan$kind, "at_units")) {
    stop(
      "emit_gaussian_aggregated() INLA: heterogeneous residual ",
      "at_units on INLA requires the multi-likelihood INLA stack ",
      "API (deferred), and no other active backend has an aggregated ",
      "heterogeneous-residual path, so pass aggregate = FALSE for the ",
      "per-row path -- on INLA, or on brms, which represents ",
      "the sectioned residual as a distributional predictor.",
      call. = FALSE
    )
  }

  agg <- as.data.frame(fb_aggregated$sufficient_stats)
  K <- fb_aggregated$K
  N <- fb_aggregated$N
  agg$ybar_k <- agg$S1_k / agg$n_k

  # Pull per-RI hyperpriors from the same priors_to_inla() translation
  # the per-row emit_inla() uses, so the aggregated + per-row INLA fits
  # share the random-intercept prior shape. The residual prior is NOT
  # plumbed through at v0.3.2: the custom-prior expression below
  # composes Delta with INLA's default loggamma(1, 5e-5) base. A user-
  # supplied residual prior on `sigma` is therefore ignored on the
  # aggregated INLA path -- documented limitation; matched-prior
  # composition is a v0.3.3 follow-on.
  hyper_ctrl <- if (inherits(fb$priors, "fb_prior")) {
    priors_to_inla(fb$priors)
  } else if (.fb_prior_scalar_supplied(fb, "vc_sd")) {
    # The legacy scalar route, keyed exactly as priors_to_inla() keys an
    # fb_prior, so the aggregated and per-row INLA fits carry the same
    # random-intercept prior. They have to: their agreement is an
    # algebraic identity that only holds when the two share one prior.
    .priors_legacy_to_inla(fb, .fb_prior_scalar_value(fb, "vc_sd"))
  } else {
    list()
  }

  fixed_form <- .agg_emit_inla_formula(fb, fb_aggregated, ri_plan, hyper_ctrl)

  # Custom-prior expression on the gaussian-precision hyperparameter:
  # combines INLA's default loggamma(1, 5e-5) prior with the Delta
  # correction that absorbs the within-cell SS into the prior so the
  # cell-mean weighted likelihood + custom prior recovers the per-row
  # posterior.
  a <- WSS_total / 2 + 5e-5
  b <- (N - K) / 2 + 1
  # Recenter the log-density to ~0 at its mode theta* = log(b / a). The
  # additive constant `c` does not change the prior (it is absorbed in the
  # normalising constant), but it keeps the values INLA's expression
  # evaluator sees O(1) near the mode. Without it the log-density is ~ -b
  # at the mode (e.g. ~ -2.5e9 at N = 5e9), at which magnitude INLA's
  # hyperparameter integration becomes numerically unstable and can
  # corrupt the fixed-effect posterior at extreme aggregated N.
  c0 <- b * (log(b / a) - 1)
  custom_prior_expr <- sprintf(
    paste0(
      "expression:",
      " a = %.16e;",
      " b = %.16e;",
      " c = %.16e;",
      " log_dens = -a * exp(theta) + b * theta - c;",
      " return(log_dens);"
    ),
    a,
    b,
    c0
  )

  control_family <- list(
    hyper = list(prec = list(prior = custom_prior_expr))
  )

  inla_call <- list(
    formula = fixed_form,
    data = agg,
    family = "gaussian",
    scale = agg$n_k,
    control.family = control_family,
    control.compute = list(
      config = TRUE,
      return.marginals = TRUE,
      dic = FALSE,
      waic = FALSE
    )
  )

  # Spliced in only when it carries something, so a call that supplies no
  # `prior_fixed_sd` reaches INLA exactly as it did before the argument
  # was wired.
  # `coef_names` is the aggregated design's own column vocabulary: the
  # aggregated INLA formula names model-matrix columns, so a `b()` prior
  # row has to be checked against those rather than against the per-row
  # term labels, or a row naming a factor level would arrive as a
  # `control.fixed` entry INLA silently ignores.
  control_fixed <- .build_inla_control_fixed(
    fb,
    .fb_prior_scalar_value(fb, "fixed_sd"),
    coef_names = setdiff(fb_aggregated$fixed_cols, "(Intercept)")
  )
  if (length(control_fixed)) {
    inla_call$control.fixed <- control_fixed
  }

  inla_fit <- do.call(INLA::inla, inla_call)

  list(
    backend = "inla",
    inla = inla_fit,
    fixed_form = fixed_form,
    custom_prior = custom_prior_expr,
    K = K,
    N = N
  )
}


# Build an INLA-side formula on the cell-level data.table. The
# `<fb_aggregated>$sufficient_stats` data.table carries the model-
# matrix contrast columns (`(Intercept)`, `f1b`, `f1c`, ...) rather
# than the original factor variables, so the formula references
# those contrast columns directly. `(Intercept)` is dropped because
# R's formula machinery auto-adds the intercept (unless `fb$intercept
# == FALSE`, in which case we suppress with `- 1`). RI grouping
# columns enter via INLA's `f(<col>, model = "iid")` syntax. Backticks
# guard against syntactically-awkward contrast names.
.agg_emit_inla_formula <- function(
  fb,
  fb_aggregated,
  ri_plan,
  hyper_ctrl = list()
) {
  fixed_cols <- setdiff(fb_aggregated$fixed_cols, "(Intercept)")
  rhs_fixed <- if (length(fixed_cols)) {
    paste(paste0("`", fixed_cols, "`"), collapse = " + ")
  } else {
    "1"
  }

  # When hyper_ctrl carries a translated user prior for this RI, splice
  # it into the f(...) call. In the default case the synthesised
  # uniform-on-SD prior (.default_uniform_prior() in flexybayes()) keys
  # every random group, so hyper_ctrl carries the same faithful uniform
  # expression prior the per-row emit_inla() path uses -- both paths
  # share one prior (scale fixed once from the per-row data), which is
  # what makes the per-row and aggregated INLA fits algebraically
  # equivalent. The former per-row/aggregated divergence under defaults
  # is resolved by that unification (priors_to_inla() now emits the exact
  # uniform on both paths instead of a per-row PC approximation).
  #
  # The fallback below fires only when no fb_prior is in play (e.g. the
  # legacy `prior_vc_sd` path) or a user supplies a partial prior that
  # omits this group. It pins the scale-invariant loggamma(1, 5e-5) so
  # the aggregated fit matches the per-row path, which inherits the same
  # INLA default when hyper_ctrl is empty.
  default_iid_hyper <- list(prior = "loggamma", param = c(1, 5e-05))
  ri_pieces <- vapply(
    ri_plan,
    function(r) {
      entry <- hyper_ctrl[[r$col]] %||% default_iid_hyper
      hyper_arg <- .inla_hyper_arg(entry)
      sprintf(
        "f(`%s`, model = \"iid\"%s)",
        r$col,
        if (nzchar(hyper_arg)) paste0(", ", hyper_arg) else ""
      )
    },
    character(1L)
  )
  rhs_ri <- if (length(ri_pieces)) {
    paste(ri_pieces, collapse = " + ")
  } else {
    ""
  }

  rhs <- if (nzchar(rhs_ri)) {
    paste(rhs_fixed, rhs_ri, sep = " + ")
  } else {
    rhs_fixed
  }
  if (!isTRUE(fb$intercept)) {
    rhs <- paste(rhs, "- 1")
  }

  stats::as.formula(paste("ybar_k ~", rhs))
}


# ---------------------------------------------------------------- #
# Shared fit-object constructor                                     #
# ---------------------------------------------------------------- #
# Builds the `<flexybayes>` fit shape that downstream methods
# (print/summary/predict/triangulate) consume. Mirrors emit_inla()
# construction so the new aggregated fits are
# byte-interchangeable for the existing method surface.
.agg_emit_build_fit <- function(
  fb,
  fb_aggregated,
  data,
  backend,
  engine_out,
  ri_plan,
  residual_plan,
  elapsed,
  the_call,
  fixed,
  random,
  residual,
  family,
  link,
  data_name,
  n_samples,
  warmup,
  chains,
  prior_fixed_sd,
  prior_vc_sd,
  known_matrices = NULL,
  weights = NULL,
  seed = NULL,
  control = NULL,
  na_action = NULL,
  verbose = TRUE,
  requested_backend = NULL,
  requested_aggregate = NULL
) {
  # Posterior summary: pull the engine-side marginals into a uniform
  # shape (mean/sd/q025/q975 per parameter) from summary.fixed +
  # summary.hyperpar.
  posterior_summary <- .agg_inla_summarise(
    engine_out, fb_aggregated, ri_plan, residual_plan
  )

  # Per-row reconstructed fitted values for the $glm shim.
  fitted_row <- .agg_reconstruct_fitted_row(
    fb,
    fb_aggregated,
    data,
    posterior_summary,
    ri_plan
  )
  y_row <- data[[fb$response]]
  resid_row <- as.numeric(y_row) - fitted_row

  glm_obj <- structure(
    list(
      coefficients = posterior_summary$beta_means,
      vcov = posterior_summary$beta_vcov,
      fitted.values = fitted_row,
      linear.predictors = fitted_row,
      residuals = resid_row,
      y = y_row,
      family = list(family = fb$family, link = fb$link %||% "identity"),
      formula = the_call,
      data = data
    ),
    class = c("flexybayes_glm_shim", "lm", "list")
  )

  extras <- structure(
    list(
      summary = posterior_summary,
      convergence = posterior_summary$convergence,
      variance_comps = posterior_summary$variance_comps,
      run_time = elapsed,
      # Built from the PRE-aggregation data on purpose. Aggregation is a
      # computational route to the same likelihood, not a different model,
      # so an aggregated fit and a per-row fit of the same model must read
      # as comparable to triangulate(). See R/model_fingerprint.R.
      fingerprint = .fb_model_fingerprint(fb, data, prior_vc_sd = prior_vc_sd),
      parse_info = list(
        # Mirror emit_inla()'s fixed-info shape so confint /
        # summary / methods.R downstream consumers see the same
        # `list(response, intercept, terms)` structure.
        fixed = list(
          response = fb$response,
          intercept = fb$intercept,
          terms = fb$fixed_terms
        ),
        random = fb$random_terms,
        residual = fb$residual_terms,
        family = list(family = fb$family, link = fb$link %||% "identity"),
        smooths = list() # the aggregated path excludes smooths.
      ),
      # The same nineteen fields, under the same names and in the same
      # order, as the per-row emits record. update() re-issues every
      # recorded field to flexybayes(), so a field this route left out
      # was a default silently substituted for what the user asked for --
      # which is why update() refused an aggregated fit outright rather
      # than re-fitting a different model under the same name.
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
        na_action = na_action,
        # The engine and the representation the call ASKED for -- the
        # requested `aggregate` above all, since this route exists only
        # because the request permitted it. Recording the resolved
        # engine instead would pin a re-fit of a changed model to a
        # representation the changed model may no longer be eligible
        # for.
        backend = requested_backend,
        aggregate = requested_aggregate,
        verbose = verbose
      ),
      model_info = list(
        n_obs = fb_aggregated$N,
        n_cells = fb_aggregated$K,
        n_fixed = length(fb$fixed_terms) + isTRUE(fb$intercept),
        n_random = length(fb$random_terms),
        family = fb$family,
        link = fb$link %||% "identity"
      ),
      aggregation_meta = list(
        N = fb_aggregated$N,
        K = fb_aggregated$K,
        compression = fb_aggregated$compression,
        residual = residual_plan$kind,
        # Which prior the posterior the user sees was taken under. The
        # aggregated likelihood plus the legacy scalar bridge's precision
        # prior recover the per-row posterior to numerical precision, so
        # that path is tagged "per_row_equivalent". An explicit prior
        # object breaks that matched-prior equivalence against a
        # default-prior per-row fit (and on the INLA aggregated path the
        # residual prior is moreover not plumbed through -- see the
        # "Matched priors" note on triangulate()), so it is tagged
        # "custom". The auto-default is a third case: it is an fb_prior
        # object, so the class test alone called it "custom" and the
        # commonest fit in the package announced an explicit prior nobody
        # supplied. It gets its own tag, which claims neither. Surfaced
        # via canonical_names() and the aggregated print / summary.
        prior_parametrization = .agg_prior_parametrization(fb$priors)
      ),
      the_call = the_call,
      formula = the_call
    ),
    class = "flexybayes_extras"
  )

  # The INLA slot holds the raw inla fit object directly, matching
  # fb_as_draws_simple.flexybayes_inla()'s `fit$inla` contract.
  raw_slot <- engine_out$inla
  raw_name <- "inla"

  # Backend identity in the class vector so generics that have no
  # `.flexybayes_aggregated` method (e.g., fb_as_draws_simple()) fall
  # through to the correct backend-specific method via S3 dispatch.
  # `flexybayes_inla` in the vector routes aggregated fits to the same
  # methods as per-row INLA fits.
  fit <- structure(
    list(
      glm = glm_obj,
      extras = extras
    ),
    class = c(
      "flexybayes_aggregated",
      "flexybayes_inla",
      "flexybayes"
    )
  )
  fit[[raw_name]] <- raw_slot
  fit
}


.agg_inla_summarise <- function(
  engine_out,
  fb_aggregated,
  ri_plan,
  residual_plan
) {
  inla_fit <- engine_out$inla
  fixed_summary <- inla_fit$summary.fixed
  beta_means <- fixed_summary$mean
  names(beta_means) <- rownames(fixed_summary)

  # `nrow =` guards the intercept-only case: diag() of a length-one
  # vector is otherwise read as diag(n) and returns an n x n identity
  # (a 0 x 0 matrix here), not a 1 x 1 covariance. Mirrors the count path.
  beta_vcov <- diag(fixed_summary$sd^2, nrow = length(beta_means))
  dimnames(beta_vcov) <- list(names(beta_means), names(beta_means))

  hyper <- inla_fit$summary.hyperpar
  prec_y_mean <- hyper$mean[grepl(
    "^Precision for the Gaussian observations",
    rownames(hyper)
  )]
  sigma_means <- if (length(prec_y_mean)) {
    sqrt(1 / prec_y_mean)
  } else {
    NA_real_
  }

  tau_means <- if (length(ri_plan)) {
    out <- numeric(length(ri_plan))
    for (i in seq_along(ri_plan)) {
      pat <- paste0("^Precision for ", ri_plan[[i]]$col, "$")
      prec_i <- hyper$mean[grepl(pat, rownames(hyper))]
      out[i] <- if (length(prec_i)) sqrt(1 / prec_i) else NA_real_
    }
    out
  } else {
    numeric(0)
  }

  list(
    backend = "inla",
    beta_means = beta_means,
    beta_vcov = beta_vcov,
    sigma_means = sigma_means,
    tau_means = tau_means,
    convergence = list(gelman = list(psrf = NULL)),
    variance_comps = list(sigma = sigma_means, tau = tau_means)
  )
}


# Reconstruct per-row fitted values from the cell-level posterior +
# original data. Pure linear-predictor expansion at the fixed-effect
# level; random-intercept BLUP shrinkage is a follow-on extension
# (would join `predict.flexybayes_aggregated` if/when added). The
# formula uses the ORIGINAL IR variable names (e.g. `f1 + f2`) so the
# per-row data.frame's factor columns resolve cleanly through
# `model.matrix()`; the resulting contrast columns must match the
# cell_design's contrast columns by construction (same R formula
# semantics applied to the same factors).
.agg_reconstruct_fitted_row <- function(
  fb,
  fb_aggregated,
  data,
  posterior_summary,
  ri_plan
) {
  fixed_form <- .fb_aggregate_fixed_formula(fb)
  if (is.null(fixed_form)) {
    X_row <- matrix(
      1,
      nrow = nrow(data),
      ncol = 1L,
      dimnames = list(NULL, "(Intercept)")
    )
  } else {
    X_row <- stats::model.matrix(fixed_form, data = data)
  }
  # Defensive: align column order between X_row and the beta vector
  # (same factor order through model.matrix, but be explicit).
  beta_aligned <- posterior_summary$beta_means[colnames(X_row)]
  as.numeric(X_row %*% beta_aligned)
}
