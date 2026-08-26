# emit_count_aggregated.R -- Exact aggregated emit for the binomial and
# poisson families.
#
# Count families aggregate exactly with no nuisance scale parameter. A
# cell whose linear predictor is constant has, summed over its rows:
#
#   binomial -- sum of successes ~ Binomial(sum of trials, p_cell), since
#               independent Binomials with a shared probability add. The
#               binomial coefficient that distinguishes the cell-total
#               likelihood from the per-row product is free of the model
#               parameters, so the posterior is unchanged.
#   poisson  -- sum of counts ~ Poisson(lambda_cell * sum of exposure),
#               since independent Poissons with a shared rate add. The
#               multinomial coefficient is parameter-free, so the
#               posterior is unchanged.
#
# Both are therefore exact compression, not approximation: INLA is given
# the cell totals with `Ntrials` (binomial) or `E` (poisson). Unlike
# the gaussian path (emit_gaussian_aggregated.R) there is no within-cell
# sum-of-squares term and no custom-prior correction -- the only
# hyperparameters are the random-intercept precisions, identical to the
# per-row model, so matched priors give a matched posterior.
#
# Internal -- not exported. Called from .fb_stream_emit() for the count
# families.

emit_count_aggregated <- function(
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
  fam <- fb$family

  .check_fb_terms(fb, "emit_count_aggregated() requires an <fb_terms> IR.")
  .check_fb_aggregated(
    fb_aggregated,
    "emit_count_aggregated() requires an <fb_aggregated> object."
  )
  if (!fam %in% c("binomial", "poisson")) {
    stop(
      "emit_count_aggregated(): family must be binomial or poisson; ",
      "got '",
      fam,
      "'.",
      call. = FALSE
    )
  }

  ri_plan <- .agg_emit_ri_plan(fb, fb_aggregated)

  t0 <- Sys.time()
  engine_out <- .emit_count_aggregated_inla(fb, fb_aggregated, ri_plan, verbose)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  .agg_count_build_fit(
    fb = fb,
    fb_aggregated = fb_aggregated,
    data = data,
    backend = backend,
    engine_out = engine_out,
    ri_plan = ri_plan,
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
# INLA path                                                         #
# ---------------------------------------------------------------- #
.emit_count_aggregated_inla <- function(fb, fb_aggregated, ri_plan, verbose) {
  .check_installed(
    "INLA",
    "Package 'INLA' is required for the INLA-backed aggregated ",
    "count emit. Install from https://inla.r-inla-download.org/."
  )

  fam <- fb$family
  agg <- as.data.frame(fb_aggregated$sufficient_stats)

  # Response + family-specific exposure column.
  if (identical(fam, "binomial")) {
    agg$.resp <- agg$succ_k
    n_trials <- agg$trials_k
    e_arg <- NULL
  } else {
    agg$.resp <- agg$count_k
    n_trials <- NULL
    e_arg <- agg$expo_k
  }

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
  inla_form <- .agg_count_inla_formula(
    fb,
    fb_aggregated,
    ri_plan,
    hyper_ctrl,
    lhs = ".resp"
  )

  inla_args <- list(
    formula = inla_form,
    data = agg,
    family = fam,
    control.compute = list(
      config = TRUE,
      return.marginals = TRUE,
      dic = FALSE,
      waic = FALSE
    )
  )
  if (!is.null(n_trials)) {
    inla_args$Ntrials <- n_trials
  }
  if (!is.null(e_arg)) {
    inla_args$E <- e_arg
  }

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
    inla_args$control.fixed <- control_fixed
  }

  inla_fit <- do.call(INLA::inla, inla_args)

  list(
    backend = "inla",
    inla = inla_fit,
    fixed_form = inla_form,
    K = fb_aggregated$K,
    N = fb_aggregated$N
  )
}

# Build the cell-level INLA formula for a count family. Mirrors
# .agg_emit_inla_formula() but parameterises the response and omits the
# gaussian scale handling.
.agg_count_inla_formula <- function(
  fb,
  fb_aggregated,
  ri_plan,
  hyper_ctrl = list(),
  lhs
) {
  fixed_cols <- setdiff(fb_aggregated$fixed_cols, "(Intercept)")
  rhs_fixed <- if (length(fixed_cols)) {
    paste(paste0("`", fixed_cols, "`"), collapse = " + ")
  } else {
    "1"
  }

  # In the default case the synthesised uniform-on-SD prior keys every
  # random group, so hyper_ctrl carries the same faithful uniform
  # expression prior the per-row emit_inla() path uses (both paths share
  # one prior). The fallback below fires only when no fb_prior is in play
  # (e.g. the legacy `prior_vc_sd` path) or a user supplies a partial
  # prior that omits this group; it pins the scale-invariant
  # loggamma(1, 5e-5) so the aggregated fit matches the per-row path,
  # which inherits the same INLA default when hyper_ctrl is empty.
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

  stats::as.formula(paste(lhs, "~", rhs))
}


# ---------------------------------------------------------------- #
# Posterior summary                                                 #
# ---------------------------------------------------------------- #
.agg_count_inla_summarise <- function(engine_out, fb_aggregated, ri_plan) {
  inla_fit <- engine_out$inla
  fixed_summary <- inla_fit$summary.fixed
  beta_means <- fixed_summary$mean
  names(beta_means) <- rownames(fixed_summary)
  beta_vcov <- diag(fixed_summary$sd^2, nrow = length(beta_means))
  dimnames(beta_vcov) <- list(names(beta_means), names(beta_means))

  hyper <- inla_fit$summary.hyperpar
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
    sigma_means = numeric(0),
    tau_means = tau_means,
    convergence = list(gelman = list(psrf = NULL)),
    variance_comps = list(sigma = numeric(0), tau = tau_means)
  )
}

# ---------------------------------------------------------------- #
# Fit-object constructor                                            #
# ---------------------------------------------------------------- #
# Builds the <flexybayes_aggregated> shape downstream methods consume,
# mirroring .agg_emit_build_fit() but for the count families (no
# residual scale; fitted values on the response scale via the inverse
# link).
.agg_count_build_fit <- function(
  fb,
  fb_aggregated,
  data,
  backend,
  engine_out,
  ri_plan,
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
  posterior_summary <- .agg_count_inla_summarise(
    engine_out, fb_aggregated, ri_plan
  )

  eta_row <- .agg_reconstruct_fitted_row(
    fb,
    fb_aggregated,
    data,
    posterior_summary,
    ri_plan
  )
  fitted_row <- .agg_count_inverse_link(eta_row, fb$family)
  y_row <- data[[fb$response]]
  resid_row <- as.numeric(y_row) - fitted_row

  glm_obj <- structure(
    list(
      coefficients = posterior_summary$beta_means,
      vcov = posterior_summary$beta_vcov,
      fitted.values = fitted_row,
      linear.predictors = eta_row,
      residuals = resid_row,
      y = y_row,
      family = list(
        family = fb$family,
        link = fb$link %||% .agg_count_link(fb$family)
      ),
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
      # Built from the PRE-aggregation data on purpose: cell-level
      # sufficient statistics are a computational route to the same
      # likelihood, not a different model, so an aggregated fit and a
      # per-row fit of the same model read as comparable to
      # triangulate(). See R/model_fingerprint.R.
      fingerprint = .fb_model_fingerprint(fb, data, prior_vc_sd = prior_vc_sd),
      parse_info = list(
        fixed = list(
          response = fb$response,
          intercept = fb$intercept,
          terms = fb$fixed_terms
        ),
        random = fb$random_terms,
        residual = fb$residual_terms,
        family = list(
          family = fb$family,
          link = fb$link %||% .agg_count_link(fb$family)
        ),
        smooths = list()
      ),
      # The same nineteen fields, under the same names and in the same
      # order, as the per-row emits and the aggregated Gaussian emit
      # record. update() re-issues every recorded field to flexybayes(),
      # so a field this route left out was a default silently substituted
      # for what the user asked for -- which is why an aggregated count
      # fit refused update() outright rather than re-fitting a different
      # model under the same name.
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
        # The engine and the representation the call ASKED for, and the
        # reporting it ran under. See the aggregated Gaussian emit for
        # why the request rather than the resolved engine is recorded.
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
        link = fb$link %||% .agg_count_link(fb$family)
      ),
      aggregation_meta = list(
        N = fb_aggregated$N,
        K = fb_aggregated$K,
        compression = fb_aggregated$compression,
        residual = "none",
        # Three cases, not two -- see .agg_prior_parametrization(). The
        # auto-default is an fb_prior object, so a class test alone
        # labelled an ordinary call "custom (explicit prior supplied)".
        prior_parametrization = .agg_prior_parametrization(fb$priors),
        streamed = isTRUE(fb_aggregated$streamed)
      ),
      the_call = the_call,
      formula = the_call
    ),
    class = "flexybayes_extras"
  )

  raw_slot <- engine_out$inla
  raw_name <- "inla"

  fit <- structure(
    list(glm = glm_obj, extras = extras),
    class = c(
      "flexybayes_aggregated",
      "flexybayes_inla",
      "flexybayes"
    )
  )
  fit[[raw_name]] <- raw_slot
  fit
}

#' @noRd
#' @keywords internal
.agg_count_link <- function(family) {
  switch(family, binomial = "logit", poisson = "log", "identity")
}

#' @noRd
#' @keywords internal
.agg_count_inverse_link <- function(eta, family) {
  switch(family, binomial = stats::plogis(eta), poisson = exp(eta), eta)
}
