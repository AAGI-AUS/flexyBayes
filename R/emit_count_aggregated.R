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
# the cell totals with `Ntrials` (binomial) or `E` (poisson), and greta
# attaches the matching count distribution to the cell-total data. Unlike
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
  backend = c("inla", "greta"),
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
  rcov = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_
) {
  backend <- match.arg(backend)
  fam <- fb$family

  if (!inherits(fb, "fb_terms")) {
    stop("emit_count_aggregated() requires an <fb_terms> IR.", call. = FALSE)
  }
  if (!inherits(fb_aggregated, "fb_aggregated")) {
    stop(
      "emit_count_aggregated() requires an <fb_aggregated> object.",
      call. = FALSE
    )
  }
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
  engine_out <- if (identical(backend, "inla")) {
    .emit_count_aggregated_inla(fb, fb_aggregated, ri_plan, verbose)
  } else {
    .emit_count_aggregated_greta(
      fb,
      fb_aggregated,
      ri_plan,
      n_samples,
      warmup,
      chains,
      prior_fixed_sd,
      prior_vc_sd,
      verbose,
      mcmc_verbose
    )
  }
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
    rcov = rcov,
    family = family,
    link = link,
    data_name = data_name,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd
  )
}


# ---------------------------------------------------------------- #
# INLA path                                                         #
# ---------------------------------------------------------------- #
.emit_count_aggregated_inla <- function(fb, fb_aggregated, ri_plan, verbose) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required for the INLA-backed aggregated ",
      "count emit. Install from https://inla.r-inla-download.org/.",
      call. = FALSE
    )
  }

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
# greta path                                                        #
# ---------------------------------------------------------------- #
.emit_count_aggregated_greta <- function(
  fb,
  fb_aggregated,
  ri_plan,
  n_samples,
  warmup,
  chains,
  prior_fixed_sd,
  prior_vc_sd,
  verbose,
  mcmc_verbose
) {
  if (!requireNamespace("greta", quietly = TRUE)) {
    stop(
      "Package 'greta' is required for the greta-backed aggregated ",
      "count emit.",
      call. = FALSE
    )
  }

  fam <- fb$family
  cell_design <- fb_aggregated$cell_design
  agg <- fb_aggregated$sufficient_stats
  K <- fb_aggregated$K
  p <- ncol(cell_design)

  # Cell-level linear predictor by a single matrix multiply (the same
  # S3-dispatch-safe construction the gaussian greta path uses).
  Z_per_term <- lapply(ri_plan, function(r) {
    Z <- matrix(0, K, r$n_levels)
    Z[cbind(seq_len(K), r$level_idx)] <- 1
    Z
  })
  M_combined <- if (length(Z_per_term)) {
    do.call(cbind, c(list(cell_design), Z_per_term))
  } else {
    cell_design
  }

  greta_ns <- asNamespace("greta")
  g <- greta_ns
  M_g <- g$as_data(M_combined)
  beta_g <- g$normal(0, prior_fixed_sd, dim = p)

  tau_per_term <- vector("list", length(ri_plan))
  u_per_term <- vector("list", length(ri_plan))
  for (i in seq_along(ri_plan)) {
    nL <- ri_plan[[i]]$n_levels
    tau <- g$normal(0, prior_vc_sd, truncation = c(0, Inf))
    u_ <- g$normal(0, tau, dim = nL)
    tau_per_term[[i]] <- tau
    u_per_term[[i]] <- u_
  }

  c_method <- get("c.greta_array", envir = greta_ns)
  matmul_method <- get("%*%.greta_array", envir = greta_ns)
  theta <- beta_g
  for (u_ in u_per_term) {
    theta <- c_method(theta, u_)
  }
  eta_cell <- matmul_method(M_g, theta)

  if (identical(fam, "binomial")) {
    prob_cell <- g$ilogit(eta_cell)
    succ_g <- g$as_data(as.integer(agg$succ_k))
    trials_g <- g$as_data(as.integer(agg$trials_k))
    g$`distribution<-`(succ_g, g$binomial(trials_g, prob_cell))
  } else {
    # `exp()` dispatches through greta's Math group-generic on the
    # greta array (greta exports no `exp` namespace function); `*`
    # dispatches through the Ops group-generic, as in the gaussian path.
    expo_g <- g$as_data(agg$expo_k)
    rate <- expo_g * exp(eta_cell)
    count_g <- g$as_data(as.integer(agg$count_k))
    g$`distribution<-`(count_g, g$poisson(rate))
  }

  model_env <- new.env(parent = greta_ns)
  model_env$beta <- beta_g
  tau_names <- character(length(tau_per_term))
  for (i in seq_along(tau_per_term)) {
    nm <- paste0("tau_", i)
    model_env[[nm]] <- tau_per_term[[i]]
    tau_names[[i]] <- nm
  }
  call_args <- c(list(quote(beta)), lapply(tau_names, as.name))
  model_call <- as.call(c(list(quote(model)), call_args))
  m <- base::eval(model_call, envir = model_env)
  draws <- g$mcmc(
    m,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    verbose = isTRUE(mcmc_verbose)
  )

  list(
    backend = "greta",
    draws = draws,
    model = m,
    beta = beta_g,
    tau_per_term = tau_per_term,
    cell_design = cell_design,
    K = K,
    N = fb_aggregated$N,
    p = p
  )
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

.agg_count_greta_summarise <- function(engine_out, fb_aggregated, ri_plan) {
  draws <- engine_out$draws
  m <- do.call(rbind, draws)

  beta_idx <- seq_len(engine_out$p)
  tau_idx <- if (length(ri_plan)) {
    (engine_out$p + 1L):(engine_out$p + length(ri_plan))
  } else {
    integer(0)
  }

  beta_means <- colMeans(m[, beta_idx, drop = FALSE])
  names(beta_means) <- colnames(fb_aggregated$cell_design)
  beta_vcov <- stats::cov(m[, beta_idx, drop = FALSE])
  dimnames(beta_vcov) <- list(names(beta_means), names(beta_means))
  tau_means <- if (length(tau_idx)) {
    colMeans(m[, tau_idx, drop = FALSE])
  } else {
    numeric(0)
  }

  psrf_mat <- tryCatch(
    {
      coda_ns <- asNamespace("coda")
      gd <- coda_ns$gelman.diag(draws, multivariate = FALSE, autoburnin = FALSE)
      gd$psrf
    },
    error = function(e) NULL
  )

  list(
    backend = "greta",
    beta_means = beta_means,
    beta_vcov = beta_vcov,
    sigma_means = numeric(0),
    tau_means = tau_means,
    convergence = list(gelman = list(psrf = psrf_mat)),
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
  rcov,
  family,
  link,
  data_name,
  n_samples,
  warmup,
  chains,
  prior_fixed_sd,
  prior_vc_sd
) {
  posterior_summary <- if (identical(backend, "inla")) {
    .agg_count_inla_summarise(engine_out, fb_aggregated, ri_plan)
  } else {
    .agg_count_greta_summarise(engine_out, fb_aggregated, ri_plan)
  }

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
      parse_info = list(
        fixed = list(
          response = fb$response,
          intercept = fb$intercept,
          terms = fb$fixed_terms
        ),
        random = fb$random_terms,
        rcov = fb$rcov_terms,
        family = list(
          family = fb$family,
          link = fb$link %||% .agg_count_link(fb$family)
        ),
        smooths = list()
      ),
      call_info = list(
        fixed = fixed,
        random = random,
        rcov = rcov,
        data_name = data_name,
        family = family,
        link = link,
        n_samples = n_samples,
        warmup = warmup,
        chains = chains,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd
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
        prior_parametrization = if (inherits(fb$priors, "fb_prior")) {
          "custom"
        } else {
          "per_row_equivalent"
        },
        streamed = isTRUE(fb_aggregated$streamed)
      ),
      the_call = the_call,
      formula = the_call
    ),
    class = "flexybayes_extras"
  )

  raw_slot <- if (identical(backend, "greta")) {
    list(draws = engine_out$draws)
  } else {
    engine_out$inla
  }
  raw_name <- if (identical(backend, "greta")) "greta" else "inla"

  fit <- structure(
    list(glm = glm_obj, extras = extras),
    class = c(
      "flexybayes_aggregated",
      if (identical(backend, "inla")) "flexybayes_inla",
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
