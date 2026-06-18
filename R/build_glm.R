# Build GLM-compatible output from posterior draws
# Not exported.

# Reconstruct a fixed-effects-only formula from a parsed fixed_info.
# Used by .build_glm() so that model.frame() / model.matrix() / terms()
# never see brms-style random-effects bars (e.g. `(1|g)`), which would
# otherwise be evaluated as bitwise-OR and warn on factor groups
# ("Ops.factor(1, g) : '|' not meaningful for factors"). For the asreml
# ingest path, fixed_info$terms already excludes random terms, so the
# reconstructed formula is equivalent to the caller-supplied `fixed`.
.fixed_only_formula <- function(fixed_info, env = parent.frame()) {
  resp <- fixed_info$response
  rhs_parts <- vapply(fixed_info$terms, function(t) t$label, character(1))
  rhs <- if (length(rhs_parts)) paste(rhs_parts, collapse = " + ") else "1"
  if (!isTRUE(fixed_info$intercept)) {
    rhs <- paste(rhs, "- 1")
  }
  stats::as.formula(paste(resp, "~", rhs), env = env)
}

# Construct a GLM-compatible object from flexyBayes posterior results
#
# @param fixed_info Parsed fixed formula
# @param data Original data.frame
# @param draws mcmc.list posterior draws
# @param ev Evaluation environment
# @param fam_link Family/link list
# @param the_call Original match.call()
# @param fixed Original fixed formula (may include brms-style random
#   bars; we strip them via .fixed_only_formula() before any
#   model.frame/model.matrix/terms call)
# @return Object of class c("flexybayes_glm", "glm", "lm")
.build_glm <- function(fixed_info, data, draws, ev, fam_link, the_call, fixed) {
  N <- nrow(data)
  resp <- fixed_info$response
  y <- data[[resp]]

  # Extract fixed effect coefficients from posterior
  coef_info <- .extract_fixed_coefs(fixed_info, draws)
  beta_hat <- coef_info$coefficients
  beta_vcov <- coef_info$vcov

  # Reconstruct a fixed-only formula so brms-style `(1|g)` terms
  # don't reach model.frame() / model.matrix() / terms().
  fixed_only <- .fixed_only_formula(
    fixed_info,
    env = if (inherits(fixed, "formula")) environment(fixed) else parent.frame()
  )

  # Build model frame and model matrix
  mf <- tryCatch(
    model.frame(fixed_only, data = data, na.action = na.omit),
    error = function(e) data
  )

  mm <- tryCatch(
    {
      X <- model.matrix(fixed_only, data = data)
      X
    },
    error = function(e) {
      matrix(1, nrow = N, ncol = 1, dimnames = list(NULL, "(Intercept)"))
    }
  )

  # Compute fitted values on response scale
  fitted_vals <- .compute_fitted(draws, ev, N, fam_link)
  resid_vals <- y - fitted_vals

  # Linear predictor (on link scale)
  linpred <- .compute_linear_predictor(draws, ev, N)

  # Family object
  fam_obj <- .get_stats_family(fam_link)

  # QR decomposition of model matrix
  qr_mm <- qr(mm)

  # Deviance approximation
  dev <- sum(resid_vals^2)

  # Number of fixed parameters
  n_fixed <- length(beta_hat)

  # Construct the GLM object
  glm_obj <- list(
    coefficients = beta_hat,
    residuals = resid_vals,
    fitted.values = fitted_vals,
    effects = NULL,
    R = NULL,
    rank = n_fixed,
    qr = qr_mm,
    family = fam_obj,
    linear.predictors = linpred,
    deviance = dev,
    aic = NA_real_,
    null.deviance = sum((y - mean(y))^2),
    iter = NA_integer_,
    weights = rep(1, N),
    prior.weights = rep(1, N),
    df.residual = N - n_fixed,
    df.null = N - 1L,
    y = y,
    converged = TRUE,
    boundary = FALSE,
    model = mf,
    call = the_call,
    formula = fixed_only,
    terms = tryCatch(terms(fixed_only, data = data), error = function(e) {
      terms(fixed_only)
    }),
    data = data,
    offset = NULL,
    control = list(),
    method = "flexyBayes",
    contrasts = attr(mm, "contrasts"),
    xlevels = .getXlevels(
      tryCatch(terms(fixed_only, data = data), error = function(e) {
        terms(fixed_only)
      }),
      mf
    ),
    na.action = NULL
  )

  # Store posterior vcov matrix
  attr(glm_obj, "posterior_vcov") <- beta_vcov

  class(glm_obj) <- c("flexybayes_glm", "glm", "lm")
  glm_obj
}

# Extract fixed effect posterior means and vcov from draws
.extract_fixed_coefs <- function(fixed_info, draws) {
  # Combine all chains into a single matrix
  all_draws <- do.call(rbind, lapply(draws, as.matrix))
  col_names <- colnames(all_draws)

  # Identify fixed effect parameters
  fixed_params <- character(0)
  coef_names <- character(0)

  if (fixed_info$intercept) {
    if ("mu_atg" %in% col_names) {
      fixed_params <- c(fixed_params, "mu_atg")
      coef_names <- c(coef_names, "(Intercept)")
    }
  }

  for (term in fixed_info$terms) {
    if (term$type == "factor") {
      tag <- term$var
      prefix <- if (fixed_info$intercept) {
        paste0("tau_", tag)
      } else {
        paste0("alpha_", tag)
      }
      # Find matching columns (may have [1], [2], etc.)
      matches <- grep(paste0("^", prefix), col_names, value = TRUE)
      if (length(matches) > 0) {
        fixed_params <- c(fixed_params, matches)
        if (length(matches) == 1 && !grepl("\\[", matches)) {
          coef_names <- c(coef_names, tag)
        } else {
          lvls <- if (!is.null(term$levels)) {
            term$levels
          } else {
            paste0(seq_along(matches))
          }
          coef_names <- c(coef_names, paste0(tag, lvls))
        }
      }
    } else if (term$type == "continuous") {
      param_nm <- paste0("beta_", term$var)
      matches <- grep(paste0("^", param_nm, "$"), col_names, value = TRUE)
      if (length(matches) > 0) {
        fixed_params <- c(fixed_params, matches)
        coef_names <- c(coef_names, term$var)
      }
    } else if (term$type == "factor_interaction") {
      tag <- paste(term$vars, collapse = "_x_")
      prefix <- if (fixed_info$intercept) {
        paste0("tau_", tag)
      } else {
        paste0("alpha_", tag)
      }
      matches <- grep(paste0("^", prefix), col_names, value = TRUE)
      if (length(matches) > 0) {
        fixed_params <- c(fixed_params, matches)
        coef_names <- c(coef_names, paste0(tag, "[", seq_along(matches), "]"))
      }
    } else if (term$type %in% c("interaction", "expression")) {
      tag <- if (term$type == "interaction") {
        paste0(term$vars[1], "_x_", term$vars[2])
      } else {
        gsub("[^A-Za-z0-9_]", "_", term$label)
      }
      param_nm <- paste0("beta_", tag)
      matches <- grep(paste0("^", param_nm, "$"), col_names, value = TRUE)
      if (length(matches) > 0) {
        fixed_params <- c(fixed_params, matches)
        coef_names <- c(coef_names, term$label)
      }
    }
  }

  if (length(fixed_params) == 0) {
    return(list(
      coefficients = setNames(numeric(0), character(0)),
      vcov = matrix(nrow = 0, ncol = 0)
    ))
  }

  # Subset draws to fixed params only
  # Handle potential missing columns gracefully
  avail <- fixed_params[fixed_params %in% col_names]
  if (length(avail) == 0) {
    return(list(
      coefficients = setNames(numeric(0), character(0)),
      vcov = matrix(nrow = 0, ncol = 0)
    ))
  }

  draws_fixed <- all_draws[, avail, drop = FALSE]
  beta_hat <- colMeans(draws_fixed)
  beta_vcov <- cov(draws_fixed)

  # Rename to interpretable names
  names(beta_hat) <- coef_names[seq_along(avail)]
  rownames(beta_vcov) <- colnames(beta_vcov) <- coef_names[seq_along(avail)]

  list(coefficients = beta_hat, vcov = beta_vcov)
}

# Compute fitted values from posterior mean of mu_i_atg
.compute_fitted <- function(draws, ev, N, fam_link) {
  tryCatch(
    {
      # Use posterior means of the linear predictor
      linpred <- .compute_linear_predictor(draws, ev, N)

      # Apply inverse link
      switch(
        fam_link$link,
        "identity" = linpred,
        "log" = exp(linpred),
        "logit" = 1 / (1 + exp(-linpred)),
        "probit" = pnorm(linpred),
        linpred # fallback
      )
    },
    error = function(e) {
      rep(NA_real_, N)
    }
  )
}

# Compute linear predictor from posterior means
.compute_linear_predictor <- function(draws, ev, N) {
  tryCatch(
    {
      # Reconstruct linear predictor from posterior mean parameter values
      all_draws <- do.call(rbind, lapply(draws, as.matrix))
      col_names <- colnames(all_draws)
      post_means <- colMeans(all_draws)

      # Start with zeros
      linpred <- rep(0, N)

      # Add intercept if present
      if ("mu_atg" %in% col_names) {
        linpred <- linpred + post_means["mu_atg"]
      }

      # Add factor effects (look for tau_* and alpha_* with _id vectors)
      env_names <- ls(ev)
      id_vars <- grep("_id$", env_names, value = TRUE)
      for (id_var in id_vars) {
        base_name <- sub("_id$", "", id_var)
        id_vec <- get(id_var, envir = ev, inherits = FALSE)

        # Check for tau_ or alpha_ parameters
        for (prefix in c("tau_", "alpha_", "u_")) {
          param_pattern <- paste0("^", prefix, base_name)
          matching <- grep(param_pattern, col_names, value = TRUE)
          if (length(matching) > 0) {
            vals <- post_means[matching]
            if (length(vals) == length(unique(id_vec))) {
              linpred <- linpred + vals[id_vec]
            }
          }
        }
      }

      # Add continuous effects
      for (nm in col_names) {
        if (grepl("^beta_", nm)) {
          var_name <- sub("^beta_", "", nm)
          if (var_name %in% env_names) {
            x_vals <- get(var_name, envir = ev, inherits = FALSE)
            if (is.numeric(x_vals) && length(x_vals) == N) {
              linpred <- linpred + post_means[nm] * x_vals
            }
          }
        }
      }

      linpred
    },
    error = function(e) {
      rep(NA_real_, N)
    }
  )
}

# Build variance components table from posterior summary
.build_variance_comps <- function(post_summary, params) {
  if (is.null(post_summary)) {
    return(NULL)
  }

  stats <- post_summary$statistics
  quants <- post_summary$quantiles
  if (is.null(stats) || is.null(quants)) {
    return(NULL)
  }

  # Find sigma_ and sd_ parameters (variance components)
  vc_params <- grep("^(sigma_|sd_|psi_|Lambda_)", rownames(stats), value = TRUE)
  if (length(vc_params) == 0) {
    return(NULL)
  }

  data.frame(
    component = vc_params,
    estimate = stats[vc_params, "Mean"],
    sd = stats[vc_params, "SD"],
    q2.5 = quants[vc_params, "2.5%"],
    q50 = quants[vc_params, "50%"],
    q97.5 = quants[vc_params, "97.5%"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Build BLUPs from posterior summary
.build_blups <- function(post_summary, random_terms) {
  if (is.null(post_summary) || length(random_terms) == 0) {
    return(list())
  }

  stats <- post_summary$statistics
  if (is.null(stats)) {
    return(list())
  }

  blups <- list()
  for (term in random_terms) {
    prefix <- switch(
      term$type,
      "simple" = ,
      "ide" = ,
      "id" = paste0("u_", term$var),
      "vm" = paste0("u_", term$var),
      "ped" = paste0("u_", term$var),
      "nested" = paste0("u_", term$inner, "_in_", term$outer),
      NULL
    )
    if (!is.null(prefix)) {
      matches <- grep(
        paste0("^", prefix, "(\\[|$)"),
        rownames(stats),
        value = TRUE
      )
      if (length(matches) > 0) {
        blups[[prefix]] <- stats[matches, "Mean"]
      }
    }
  }
  blups
}

# Build predictions data frame
.build_predictions <- function(data, fixed_info, post_summary, ev, fam_link) {
  N <- nrow(data)
  y <- data[[fixed_info$response]]

  tryCatch(
    {
      all_draws <- if (!is.null(ev$atg_draws)) {
        do.call(rbind, lapply(ev$atg_draws, as.matrix))
      } else {
        NULL
      }

      if (is.null(all_draws)) {
        return(NULL)
      }

      fitted <- .compute_fitted(ev$atg_draws, ev, N, fam_link)
      resid <- y - fitted

      data.frame(
        obs = seq_len(N),
        observed = y,
        fitted = fitted,
        residual = resid,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) NULL
  )
}
