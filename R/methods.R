# S3 methods for the flexybayes class

#' Print a compact description of a flexyBayes fit
#'
#' Reports the call the fit came from, the engine that ran it, the model
#' dimensions, and the headline posterior quantities. Use [summary()] for
#' the full coefficient and variance-component tables.
#'
#' @param x A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the description it
#'   prints.
#' @export
print.flexybayes <- function(x, ...) {
  ci <- x$extras$call_info
  mi <- x$extras$model_info

  cat("Bayesian mixed model  [flexyBayes]\n")
  cat(strrep("-", 55), "\n")

  cat("  Fixed  :", deparse(ci$fixed), "\n")
  if (!is.null(ci$random)) {
    cat("  Random :", deparse(ci$random), "\n")
  }
  if (!is.null(ci$residual) && !identical(deparse(ci$residual), "~units")) {
    cat("  Residual :", deparse(ci$residual), "\n")
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
    "sec\n"
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

  # The truth-display surface (v0.3.8): "Representation:" and
  # "Engine:" as two adjacent lines. Replaces the single-line
  # "Exact.: <exactness>" rendering that conflated the representation
  # regime (exact vs aggregated_exact) with the inference engine
  # (greta MCMC vs INLA Laplace etc.). Older fits without $exactness
  # fall through silently.
  if (!is.null(x$exactness)) {
    bd <- x$extras$backend_decision
    cat("  Representation: ", .repr_label_for_fit(x, bd), "\n", sep = "")
    cat("  Engine:         ", .engine_label_for_fit(x, bd), "\n", sep = "")
  }

  # Quick convergence
  if (!is.null(x$extras$convergence$gelman)) {
    rhat <- x$extras$convergence$gelman$psrf[, "Point est."]
    max_rhat <- max(rhat, na.rm = TRUE)
    flag <- if (max_rhat < 1.05) {
      " [OK]"
    } else if (max_rhat < 1.10) {
      " [borderline]"
    } else {
      " [!]"
    }
    cat("  Max Rhat:", round(max_rhat, 3), flag, "\n")
  }
  min_eff <- if (!is.null(x$extras$convergence$n_eff)) {
    min(x$extras$convergence$n_eff, na.rm = TRUE)
  } else {
    NA
  }
  if (!is.na(min_eff)) {
    cat("  Min ESS:", round(min_eff, 0), "\n")
  }

  # The component list is read off the object instead of being declared,
  # so a fit only ever advertises a slot it actually carries. The former
  # fixed block named `$greta` on every object printed through this
  # method, including brms and INLA fits that have no such slot.
  cat(strrep("-", 55), "\n")
  .print_fit_components(
    x,
    c(
      glm = "GLM-compatible (summary, emmeans, etc.)",
      brms = "live brmsfit (loo, posterior_predict, summary)",
      inla = "raw INLA fit (use INLA's summary, plot, etc.)",
      greta = "native greta (draws, model, calculate)",
      extras = "diagnostics, BLUPs, variance components"
    ),
    width = 8L
  )

  invisible(x)
}


#' Print the component list a fitted object actually carries
#'
#' Prints one line per slot the object holds, taking the wording and the
#' print order from the caller's label table. A slot named in the table
#' but absent from the object is skipped, so a print method cannot
#' advertise a component the backend never built.
#'
#' @param x A fitted object stored as a list, whose names are the slot
#'   names to be matched against `labels`.
#' @param labels A named character vector of one-line slot descriptions,
#'   named by slot and ordered as the lines should print.
#' @param width Integer field width for the `$slot` column. Callers pass
#'   the width their surrounding output already uses so the `--`
#'   separators stay aligned.
#' @returns Invisibly, a character vector of the slot names printed, in
#'   print order. Called for the lines it writes to the console.
#'
#' @noRd
#' @keywords internal
.print_fit_components <- function(x, labels, width = 8L) {
  present <- intersect(names(labels), names(x))
  for (slot in present) {
    cat(
      "  ",
      formatC(paste0("$", slot), width = -width),
      "-- ",
      labels[[slot]],
      "\n",
      sep = ""
    )
  }
  invisible(present)
}

#' Summarise a flexyBayes fit
#'
#' Prints the fixed-effect posterior summaries, the variance components,
#' and the sampler or approximation diagnostics the backend recorded. On a
#' brms fit with a sectioned residual it additionally prints the residual
#' variance and standard deviation for each level of the sectioning
#' factor, computed from the draws rather than by transforming a
#' posterior mean.
#'
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Invisibly, a list holding the posterior summary tables the
#'   method printed, so a caller can use the numbers without re-parsing
#'   the output.
#' @export
summary.flexybayes <- function(object, ...) {
  ci <- object$extras$call_info
  mi <- object$extras$model_info

  cat("Bayesian mixed model summary  [flexyBayes]\n")
  cat(strrep("=", 60), "\n")
  cat("  Fixed  :", deparse(ci$fixed), "\n")
  if (!is.null(ci$random)) {
    cat("  Random :", deparse(ci$random), "\n")
  }
  cat("  Family :", mi$family, "/", mi$link, "\n")
  cat(
    "  N =",
    mi$n_obs,
    ", chains =",
    ci$chains,
    ", samples =",
    ci$n_samples,
    "\n"
  )

  # Representation:/Engine: two-line truth display (v0.3.8),
  # replacing the single-line "Exact.:". The conditional
  # compression-ratio line avoids misleading claims on small or
  # non-aggregable data by only showing the ratio when N/K >= 2.
  if (!is.null(object$exactness)) {
    bd <- object$extras$backend_decision
    cat("  Representation: ", .repr_label_for_fit(object, bd), "\n", sep = "")
    cat("  Engine:         ", .engine_label_for_fit(object, bd), "\n", sep = "")
  }
  am <- object$extras$aggregation_meta
  if (
    identical(object$exactness, "aggregated_exact") &&
      !is.null(am) &&
      am$N / am$K >= 2
  ) {
    cat(sprintf(
      "  Agg.   : N = %s rows -> K = %s cells (ratio %.0f:1)\n",
      format(am$N, big.mark = " ", scientific = FALSE),
      format(am$K, big.mark = " ", scientific = FALSE),
      am$N / am$K
    ))
  }
  cat("\n")

  # Fixed effects
  cat("-- Fixed effects (posterior) ", strrep("-", 33), "\n")
  beta <- coef(object)
  if (length(beta) > 0) {
    ci_mat <- confint(object)
    fx_df <- data.frame(
      Estimate = beta,
      Post.SD = sqrt(diag(vcov(object))),
      `2.5%` = ci_mat[, 1],
      `97.5%` = ci_mat[, 2],
      check.names = FALSE
    )
    print(round(fx_df, 4))
  } else {
    cat("  (none)\n")
  }

  # Variance components
  cat("\n-- Variance components ", strrep("-", 38), "\n")
  vc <- object$extras$variance_comps
  if (!is.null(vc) && nrow(vc) > 0) {
    print(data.frame(
      Component = vc$component,
      Estimate = round(vc$estimate, 4),
      SD = round(vc$sd, 4),
      `2.5%` = round(vc$q2.5, 4),
      `97.5%` = round(vc$q97.5, 4),
      check.names = FALSE,
      row.names = NULL
    ))
  } else {
    cat("  (none available)\n")
  }

  # A sectioned residual has no scalar sigma for the table above to
  # report, so the per-level variances are printed in their own block
  # (a no-op on every other fit). See R/emit_brms.R.
  .print_brms_residual_by_level(object)

  # Convergence
  cat("\n-- Convergence ", strrep("-", 45), "\n")
  conv <- object$extras$convergence
  if (!is.null(conv$gelman)) {
    rhat <- conv$gelman$psrf[, "Point est."]
    cat(
      "  Rhat range:",
      round(min(rhat, na.rm = TRUE), 3),
      "-",
      round(max(rhat, na.rm = TRUE), 3),
      "\n"
    )
  }
  if (!is.null(conv$n_eff)) {
    cat(
      "  ESS  range:",
      round(min(conv$n_eff, na.rm = TRUE), 0),
      "-",
      round(max(conv$n_eff, na.rm = TRUE), 0),
      "\n"
    )
  }
  cat("  Run time  :", round(object$extras$run_time, 1), "sec\n")

  invisible(object$extras$summary)
}

# Human-readable label for the aggregated-fit prior parametrization,
# shared by the aggregated print + summary methods. The
# "per_row_equivalent" case reassures the user that the matched-prior
# guarantee holds (aggregated posterior == per-row posterior under the
# default prior); the "custom" case flags that an explicit prior was
# supplied and points at prior_summary().
.agg_fit_label <- function(mi) {
  # The aggregated header hard-coded "aggregated-gaussian" on every
  # aggregated fit, including the binomial and Poisson emits the count
  # aggregator produces -- a wrong family name on the surface the
  # package uses to signal exactness. Read it off the fit.
  fam <- mi$family %||% NA_character_
  if (!is.character(fam) || length(fam) != 1L || is.na(fam) || !nzchar(fam)) {
    return("aggregated")
  }
  paste0("aggregated-", fam)
}

.agg_prior_label <- function(pp) {
  switch(
    pp,
    per_row_equivalent = paste0(
      "per-row-equivalent (default prior; ",
      "aggregated posterior matches per-row)"
    ),
    custom = "custom (explicit prior supplied; see prior_summary())",
    pp
  )
}

#' Print a flexybayes_aggregated object
#'
#' Brief one-screen summary of a fit produced by
#' `flexybayes(..., aggregate = "auto"/TRUE)` or `fb_brms(..., aggregate
#' = ...)`. The header names the fit's own family
#' (`aggregated-gaussian`, `aggregated-binomial`, `aggregated-poisson`).
#' Includes the `exactness` field and the cell compression ratio (when
#' N/K >= 2).
#'
#' @param x   A `<flexybayes_aggregated>` object, as returned by a fit
#'   run on the aggregated representation.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the one-screen summary
#'   it prints.
#' @export
print.flexybayes_aggregated <- function(x, ...) {
  mi <- x$extras$model_info
  am <- x$extras$aggregation_meta
  bd <- x$extras$backend_decision

  cat(sprintf("Bayesian mixed model  [flexyBayes / %s]\n", .agg_fit_label(mi)))
  cat(strrep("-", 60), "\n")
  cat("  family:    ", mi$family, "(", mi$link, "link)\n")
  cat("  N obs:     ", mi$n_obs, "\n")
  cat("  K cells:   ", mi$n_cells, "\n")
  cat("  fixed:     ", mi$n_fixed, "\n")
  cat("  random:    ", mi$n_random, "\n")
  cat("  backend:   ", bd$backend, "(path =", bd$path, ")\n")
  cat("  runtime:   ", round(x$extras$run_time, 2), "sec\n")
  # Truth display as adjacent lines (v0.3.8).
  cat("  Representation: ", .repr_label_for_fit(x, bd), "\n", sep = "")
  cat("  Engine:         ", .engine_label_for_fit(x, bd), "\n", sep = "")
  if (!is.null(am) && am$N / am$K >= 2) {
    cat(sprintf(
      "  aggregation: N = %s rows -> K = %s cells (ratio %.0f:1)\n",
      format(am$N, big.mark = " ", scientific = FALSE),
      format(am$K, big.mark = " ", scientific = FALSE),
      am$N / am$K
    ))
  }
  if (!is.null(am$prior_parametrization)) {
    cat("  priors:    ", .agg_prior_label(am$prior_parametrization), "\n")
  }
  # An aggregated fit carries `$inla` or `$greta` depending on the engine
  # that ran it, so the component list is read off the object rather than
  # assuming the INLA path.
  cat(strrep("-", 60), "\n")
  .print_fit_components(
    x,
    c(
      glm = "per-row reconstructed fitted values + coef shim",
      inla = "raw aggregated INLA fit (use INLA's summary etc.)",
      greta = "raw aggregated greta draws (coda mcmc.list)",
      extras = "summary, aggregation_meta, backend_decision"
    ),
    width = 10L
  )
  invisible(x)
}


#' Summarise a flexybayes_aggregated object
#'
#' Posterior summary read off the aggregated INLA fit's
#' `summary.fixed` + `summary.hyperpar` slots. The header names the
#' fit's own family. Shows the compression line when N/K >= 2.
#'
#' @param object A `<flexybayes_aggregated>` object, as returned by a fit
#'   run on the aggregated representation.
#' @param ...    Ignored. Present for compatibility with the generic.
#' @returns Invisibly, a list holding the posterior summary tables the
#'   method printed, so a caller can use the numbers without re-parsing
#'   the output.
#' @export
summary.flexybayes_aggregated <- function(object, ...) {
  mi <- object$extras$model_info
  am <- object$extras$aggregation_meta
  bd <- object$extras$backend_decision
  ps <- object$extras$summary

  cat(sprintf(
    "Bayesian mixed model summary  [flexyBayes / %s]\n",
    .agg_fit_label(mi)
  ))
  cat(strrep("=", 65), "\n")
  cat("  family:    ", mi$family, "/", mi$link, "\n")
  cat("  N =", mi$n_obs, ", K =", mi$n_cells, "\n")
  cat("  backend:   ", bd$backend, "\n")
  cat("  exactness: ", object$exactness, "\n")
  if (!is.null(am$prior_parametrization)) {
    cat("  priors:    ", .agg_prior_label(am$prior_parametrization), "\n")
  }
  if (!is.null(am) && am$N / am$K >= 2) {
    cat(sprintf(
      "  aggregation: N = %s rows -> K = %s cells (ratio %.0f:1)\n",
      format(am$N, big.mark = " ", scientific = FALSE),
      format(am$K, big.mark = " ", scientific = FALSE),
      am$N / am$K
    ))
  }
  cat("\n")

  cat("-- Fixed effects (posterior) ", strrep("-", 35), "\n")
  beta <- ps$beta_means
  if (length(beta) > 0L) {
    sds <- sqrt(diag(ps$beta_vcov))
    fx_df <- data.frame(
      Estimate = round(beta, 4),
      Post.SD = round(sds, 4),
      row.names = names(beta)
    )
    print(fx_df)
  } else {
    cat("  (none)\n")
  }

  cat("\n-- Variance components ", strrep("-", 40), "\n")
  vc <- list(sigma = ps$sigma_means, tau = ps$tau_means)
  if (length(vc$sigma)) {
    cat(sprintf("  sigma (residual SD): %.4f\n", vc$sigma))
  }
  if (length(vc$tau)) {
    for (i in seq_along(vc$tau)) {
      cat(sprintf("  tau_%d  (random SD): %.4f\n", i, vc$tau[i]))
    }
  }

  cat("\n  Run time:", round(object$extras$run_time, 2), "sec\n")
  invisible(ps)
}


#' Extract fixed effect coefficients
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A named numeric vector of the fixed effects' posterior means,
#'   on the treatment-contrast basis, named for the design-matrix
#'   columns. Empty when the model carries no fixed effects.
#' @export
coef.flexybayes <- function(object, ...) {
  object$glm$coefficients
}

#' Extract variance-covariance matrix of fixed effects
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns The posterior covariance matrix of the fixed effects, square
#'   with one row and column per coefficient and dimnames taken from
#'   [coef()]. This is a posterior covariance, not a sampling-theory
#'   variance estimate, though the downstream packages that consume it
#'   treat it as one.
#' @export
vcov.flexybayes <- function(object, ...) {
  attr(object$glm, "posterior_vcov")
}

#' Credible intervals for the fixed effects of a flexyBayes fit
#'
#' Returns posterior quantile-based credible intervals, not frequentist
#' confidence intervals: the bounds are empirical quantiles of the
#' fixed-effect posterior draws the fit carries.
#'
#' This method reads the draws slot written by the greta-shaped emits
#' (`fit$greta$draws`). A fit from an engine that stores its posterior in
#' another shape reaches its own method -- [confint.flexybayes_inla()] for
#' INLA, [confint.flexybayes_brms()] for brms. An object that carries
#' neither is refused by name rather than returned empty, because an empty
#' interval matrix reads as "no fixed effects" and not as "this fit cannot
#' answer the question".
#'
#' @param object A `flexybayes` fit carrying a `$greta$draws` posterior
#'   slot.
#' @param parm Character vector of parameter names to return, or `NULL`
#'   (the default) for every fixed effect.
#' @param level Credible level for the interval, as a proportion. The
#'   default `0.95` returns the 2.5th and 97.5th posterior percentiles.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A matrix with one row per fixed effect and two columns holding
#'   the lower and upper credible bounds, named for the percentiles used.
#' @export
confint.flexybayes <- function(object, parm = NULL, level = 0.95, ...) {
  draws <- object$greta$draws
  if (is.null(draws)) {
    stop(.fb_refusal_condition(
      reason_code = "fit_lacks_posterior_draws",
      message = paste0(
        "confint() cannot form credible intervals for this fit: it carries ",
        "no posterior-draw slot (`$greta$draws`). The fit's class is <",
        paste(class(object), collapse = ", "), ">. Fixed-effect intervals ",
        "are available from summary() on every active engine, and an INLA ",
        "fit answers confint() from its own marginals."
      )
    ))
  }
  all_draws <- do.call(rbind, lapply(draws, as.matrix))

  # Get fixed effect column names
  beta_names <- names(coef(object))
  if (length(beta_names) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }

  # Use the raw draw columns
  col_names <- colnames(all_draws)

  # Find the draws columns for fixed effects
  fixed_info <- object$extras$parse_info$fixed
  fixed_draw_cols <- .get_fixed_draw_columns(fixed_info, col_names)

  if (length(fixed_draw_cols) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  avail <- fixed_draw_cols[fixed_draw_cols %in% col_names]
  ci_mat <- t(apply(
    all_draws[, avail, drop = FALSE],
    2,
    quantile,
    probs = probs
  ))
  rownames(ci_mat) <- beta_names[seq_along(avail)]
  colnames(ci_mat) <- paste0(round(probs * 100, 1), "%")

  if (!is.null(parm)) {
    ci_mat <- ci_mat[parm, , drop = FALSE]
  }

  ci_mat
}

# Helper: get draw column names for fixed effects
.get_fixed_draw_columns <- function(fixed_info, col_names) {
  result <- character(0)
  if (fixed_info$intercept && "mu_atg" %in% col_names) {
    result <- c(result, "mu_atg")
  }
  for (term in fixed_info$terms) {
    if (term$type == "factor") {
      prefix <- if (fixed_info$intercept) {
        paste0("tau_", term$var)
      } else {
        paste0("alpha_", term$var)
      }
      matches <- grep(paste0("^", prefix), col_names, value = TRUE)
      result <- c(result, matches)
    } else if (term$type == "continuous") {
      nm <- paste0("beta_", term$var)
      if (nm %in% col_names) result <- c(result, nm)
    } else if (term$type == "factor_interaction") {
      tag <- paste(term$vars, collapse = "_x_")
      prefix <- if (fixed_info$intercept) {
        paste0("tau_", tag)
      } else {
        paste0("alpha_", tag)
      }
      matches <- grep(paste0("^", prefix), col_names, value = TRUE)
      result <- c(result, matches)
    } else if (term$type %in% c("interaction", "expression")) {
      tag <- if (term$type == "interaction") {
        paste0(term$vars[1], "_x_", term$vars[2])
      } else {
        gsub("[^A-Za-z0-9_]", "_", term$label)
      }
      nm <- paste0("beta_", tag)
      if (nm %in% col_names) result <- c(result, nm)
    }
  }
  result
}

#' Extract fitted values
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A numeric vector of in-sample fitted values on the response
#'   scale, one per row of the fitted data, each the posterior mean of
#'   that observation's conditional expectation. Rows carried as latent
#'   under `na_action = "augment"` receive a fitted value like any other.
#' @export
fitted.flexybayes <- function(object, ...) {
  object$glm$fitted.values
}

#' Extract residuals
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A numeric vector of response residuals, the observed value
#'   minus the posterior-mean fitted value, one per row of the fitted
#'   data. A row whose response was missing and carried as latent has no
#'   observed value and returns `NA`.
#' @export
residuals.flexybayes <- function(object, ...) {
  object$glm$residuals
}

#' Predict from a flexybayes model
#'
#' @param object A fitted `flexybayes` object of any backend.
#' @param newdata Optional new data frame for prediction. If NULL, returns
#'   fitted values from the original data.
#' @param type `"link"` for linear predictor, `"response"` for response scale.
#' @param se.fit A single logical. `TRUE` returns standard errors
#'   alongside the predictions.
#' @param ... Ignored. Present for compatibility with the generic.
#' @param chunk_size Optional integer. When supplied, and when
#'   `newdata` has more rows than `chunk_size`, prediction iterates
#'   over chunks of this size and concatenates the results. Default
#'   `NULL` preserves the single-pass behaviour. Per-chunk prediction
#'   uses the same factor-level dictionary so `chunk_size` does not
#'   change the numerical output -- only the wall-time and peak
#'   memory profile. A typical setting at large `newdata` is
#'   `chunk_size = 10000L`.
#' @param allow_new_levels One of `"population"` (default), `"sample"`,
#'   or `"refuse"`. Policy for handling factor levels in `newdata`
#'   that are not in the fit-time dictionary. `"population"` sets
#'   unknown-level rows to NA on the
#'   affected column and emits a warning naming the count; downstream
#'   prediction returns NA for these rows. `"refuse"` raises a
#'   structured stop on the first unknown level. `"sample"` (active
#'   since v0.3.5) layers a fresh `Normal(0, tau_<group>)`
#'   random-effect realisation onto each unknown row per posterior
#'   draw -- caller's `set.seed()` controls reproducibility. The
#'   in-memory return reports the posterior-mean prediction
#'   (sampled-RE contribution averages toward zero across draws);
#'   the file-backed return (`output_file = ...`) captures the
#'   proper per-row posterior interval reflecting sampled-RE
#'   uncertainty. Only consulted when the fit carries an
#'   `extras$fb_dataset` slot (fits produced under v0.3.4+); legacy
#'   fits skip dictionary resolution entirely.
#' @param output_file Optional character path. When supplied,
#'   prediction is written to disk under the format resolved by
#'   `format`. The file is a
#'   tabular structure with columns `point`, `lower`, `upper`
#'   (95\% posterior credible bounds), optional `se.fit`, plus the
#'   columns of `newdata`. Refuses if `output_file` already exists
#'   (no silent overwrite). Returns the path invisibly. Requires
#'   a fit with posterior draws on `$greta$draws` (greta-backend
#'   fits including `fb_brms(..., backend = "greta")`); INLA fits
#'   do not currently expose per-draw access and route to a
#'   structured refusal. Default `NULL` (in-memory return).
#' @param format One of `"auto"` (default), `"csv"`, `"rds"`, `"fst"`.
#'   Only consulted when `output_file` is supplied. `"auto"`
#'   resolves to: `"csv"` when `interop = TRUE`; `"fst"` when
#'   `nrow(newdata) >= 1e6` and `fst` is installed; `"rds"`
#'   otherwise. `"fst"` requested without `fst` installed raises a
#'   structured refusal naming the install command. The fst path
#'   is 30--40x faster than rds at N >= 1e6 per the
#'   Stage-3B-shape benchmark (`benchmark_results/fst_stage3b_2026-05-23`).
#' @param interop Logical. When `TRUE`, the format-resolution rule
#'   under `format = "auto"` prefers `"csv"` (universally readable)
#'   over rds / fst (R-only / fst-only). Useful for handing the
#'   prediction grid to a non-R consumer. Default `FALSE`. Only
#'   consulted when `output_file` is supplied and `format = "auto"`.
#' @return If `output_file` is supplied: invisibly returns the path
#'   that was written to. Otherwise: if `se.fit = FALSE`, a numeric
#'   vector. If `se.fit = TRUE`, a list with `fit` and `se.fit`.
#' @export
predict.flexybayes <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  se.fit = FALSE,
  chunk_size = NULL,
  allow_new_levels = c("population", "sample", "refuse"),
  output_file = NULL,
  format = c("auto", "csv", "rds", "fst"),
  interop = FALSE,
  ...
) {
  type <- match.arg(type)
  allow_new_levels <- match.arg(allow_new_levels)
  format <- match.arg(format)

  if (is.null(newdata)) {
    if (!is.null(output_file)) {
      stop(
        "predict.flexybayes(): `output_file` is only supported ",
        "with `newdata`. For fitted-value persistence, write the ",
        "fit's $glm$fitted.values to disk directly.",
        call. = FALSE
      )
    }
    if (type == "link") {
      pred <- object$glm$linear.predictors
    } else {
      pred <- object$glm$fitted.values
    }
    if (se.fit) {
      return(list(fit = pred, se.fit = rep(NA_real_, length(pred))))
    }
    return(pred)
  }

  # Dictionary-backed factor handling. When the fit carries a
  # persisted <fb_dataset> descriptor (v0.3.4+ fits), resolve factor
  # levels in newdata against the fit-time codes per allow_new_levels
  # policy. Legacy fits (no extras$fb_dataset slot) skip this step and
  # fall through to the v0.3.3 behaviour.
  # "sample" was reserved at v0.3.4 (deferred-stop); v0.3.5 activates
  # the branch -- the resolver attaches `attr(newdata,
  # "sample_re_summary")` so downstream paths can layer a sampled
  # RE realisation onto each unknown row per draw.
  fb_dataset_meta <- object$extras$fb_dataset
  if (!is.null(fb_dataset_meta)) {
    newdata <- .predict_resolve_factors(
      newdata,
      fb_dataset_meta,
      allow_new_levels
    )
  }

  # File-backed output path. Routes before the chunked-iteration
  # branch because .predict_to_file()
  # handles its own chunking with a per-chunk per-draw posterior
  # interval computation, accumulates point/lower/upper vectors
  # across chunks, and writes the result via .predict_write_file()
  # under the format resolved from `format` / `interop` /
  # `nrow(newdata)` / fst availability. Returns invisible(path).
  if (!is.null(output_file)) {
    return(.predict_to_file(
      object = object,
      newdata = newdata,
      output_file = output_file,
      format = format,
      interop = interop,
      type = type,
      se.fit = se.fit,
      chunk_size = chunk_size,
      allow_new_levels = allow_new_levels
    ))
  }

  # In-memory "sample" path: when the resolver tagged unknown rows
  # via `sample_re_summary`, route through .predict_in_memory_sample()
  # which uses per-draw computation to layer sampled-RE
  # contributions onto unknown rows. Single-pass; sample mode does
  # not chunk in the in-memory return path (sample + chunked
  # iteration would re-sample per chunk and break bitwise
  # equivalence; chunked sample is only supported via
  # `output_file`).
  if (
    identical(allow_new_levels, "sample") &&
      !is.null(attr(newdata, "sample_re_summary"))
  ) {
    return(.predict_in_memory_sample(
      object = object,
      newdata = newdata,
      type = type,
      se.fit = se.fit
    ))
  }

  # Chunked iteration. Routes through the chunk iterator when
  # chunk_size is supplied AND the row count exceeds
  # it. Chunked path recurses into predict.flexybayes(chunk_size =
  # NULL) per chunk so the chunked branch only fires once at the
  # outer call. Dictionary resolution above has already happened on
  # the full newdata; per-chunk recursion passes through unchanged
  # factor columns. (The recursive call would re-resolve harmlessly
  # since known-level factor columns round-trip as identity, but the
  # outer-level resolution gives a single warning + single attribute
  # rather than n_chunks worth.)
  if (
    !is.null(chunk_size) &&
      is.numeric(chunk_size) &&
      length(chunk_size) == 1L &&
      chunk_size > 0L &&
      nrow(newdata) > chunk_size
  ) {
    return(.predict_chunked_iterate(
      object,
      newdata,
      chunk_size = as.integer(chunk_size),
      type = type,
      se.fit = se.fit,
      allow_new_levels = allow_new_levels,
      ...
    ))
  }

  # Smooth-aware prediction on newdata. The smooth terms live on the
  # IR's random_terms slot (asreml-style fits put s() in `random`, not
  # in the fixed formula). The truth source is the IR's
  # parse_info$random list plus the parse_info$smooths slot populated
  # at codegen time. If the IR carries smooth_mgcv random terms but
  # the smooths slot is empty/missing (legacy fits saved before smooth
  # support shipped), refuse cleanly. If the smooths slot is
  # populated, use mgcv::PredictMat() per smooth and column-bind to
  # the linear-only design matrix.
  fixed_formula <- object$glm$formula
  random_terms <- object$extras$parse_info$random
  smooths_slot <- object$extras$parse_info$smooths
  has_smooth_in_ir <- !is.null(random_terms) &&
    any(vapply(
      random_terms,
      function(t) !is.null(t$type) && t$type == "smooth_mgcv",
      logical(1)
    ))

  if (
    has_smooth_in_ir &&
      (is.null(smooths_slot) || length(smooths_slot) == 0L)
  ) {
    stop(
      "predict.flexybayes(): the fit was produced before smooth-",
      "basis retention shipped (or via fb_greta() where smooth ",
      "basis is user-managed). Refusing to use stats::model.matrix() ",
      "on smooth terms -- it would produce silently wrong ",
      "predictions on newdata. Re-fit via flexybayes() / fb_brms() ",
      "on the current package version to enable predict() with ",
      "newdata.",
      call. = FALSE
    )
  }

  tryCatch(
    {
      if (has_smooth_in_ir && length(smooths_slot) > 0L) {
        # Route the smooth in-memory path through the shared kernel
        # so it matches the file-backed (.predict_to_file) and
        # sampled-RE
        # (.predict_in_memory_sample) paths byte-for-byte. The kernel
        # builds the linear-only design matrix via .linear_formula(),
        # multiplies the fixed-effect per-draw beta block, layers the
        # per-smooth Bnew %*% (raw * sigma) contribution per draw, and
        # returns the n_rows x n_draws_eff linear predictor matrix.
        # Posterior aggregation (mean point + Monte Carlo se.fit) is
        # applied here so the kernel stays contract-neutral.
        per_draw_lin <- .predict_linear_draws(
          object = object,
          newdata = newdata,
          n_draws_cap = 4000L,
          sampled_re = NULL,
          include = c("fixed", "smooth")
        )
        fam_link <- object$extras$parse_info$family
        per_draw <- if (identical(type, "response")) {
          .apply_link_inverse(per_draw_lin, fam_link$link)
        } else {
          per_draw_lin
        }
        pred <- rowMeans(per_draw)
        if (se.fit) {
          se <- apply(per_draw, 1L, stats::sd)
          return(list(fit = pred, se.fit = se))
        }
        return(pred)
      }

      # No-smooth path. Preserves the coef()-based posterior-mean
      # computation so INLA fits (which do not currently expose
      # per-draw access on $greta$draws) continue to predict() on
      # newdata. The file-backed and sampled-RE paths remain
      # greta-only; the kernel refuses cleanly there.
      mm <- stats::model.matrix(
        stats::delete.response(stats::terms(fixed_formula)),
        data = newdata
      )
      beta <- coef(object)

      common <- intersect(colnames(mm), names(beta))
      if (length(common) == 0) {
        warning("No matching coefficients for newdata prediction.")
        pred <- rep(NA_real_, nrow(newdata))
      } else {
        pred <- as.numeric(mm[, common, drop = FALSE] %*% beta[common])
      }

      if (type == "response") {
        fam_link <- object$extras$parse_info$family
        pred <- switch(
          fam_link$link,
          "identity" = pred,
          "log" = exp(pred),
          "logit" = 1 / (1 + exp(-pred)),
          "probit" = pnorm(pred),
          pred
        )
      }

      if (se.fit) {
        V <- vcov(object)
        common <- intersect(colnames(mm), names(beta))
        if (length(common) > 0 && nrow(V) > 0) {
          se <- sqrt(rowSums(
            (mm[, common, drop = FALSE] %*% V[common, common]) *
              mm[, common, drop = FALSE]
          ))
        } else {
          se <- rep(NA_real_, nrow(newdata))
        }
        return(list(fit = pred, se.fit = se))
      }

      pred
    },
    error = function(e) {
      # The legacy-fit correctness refusal above is a deliberate
      # stop() with a structured message; re-raise it rather than
      # swallow into a warning + NA vector.
      msg <- conditionMessage(e)
      if (grepl("smooth-basis retention", msg, fixed = TRUE)) {
        stop(msg, call. = FALSE)
      }
      warning("Prediction with newdata failed: ", msg)
      rep(NA_real_, nrow(newdata))
    }
  )
}

# Internal helper: TRUE when the rhs of `f` contains at least one
# s() smooth specification. Used by predict.flexybayes() to drive
# the smooth-aware branch + the legacy-fit refusal.
.formula_has_smooth <- function(f) {
  if (!inherits(f, "formula")) {
    return(FALSE)
  }
  rhs <- f[[length(f)]]
  any(grepl("\\bs\\(", deparse(rhs)))
}

# Internal helper: strip s() terms from the rhs of `f` and
# return the linear-only formula. Preserves the lhs, the intercept
# convention, and any non-s() terms (factor / continuous /
# interaction / I() expressions). Used by predict.flexybayes() to
# build the linear-only design matrix before column-binding the
# per-smooth basis from mgcv::PredictMat().
.linear_formula <- function(f) {
  if (!inherits(f, "formula")) {
    return(f)
  }
  tt <- stats::terms(f)
  labs <- attr(tt, "term.labels")
  keep <- !grepl("^s\\(", labs)
  # Preserve intercept presence (attr(tt, "intercept") is 0 or 1).
  intercept_str <- if (isTRUE(attr(tt, "intercept") == 1)) "1" else "0"
  rhs_terms <- c(intercept_str, labs[keep])
  rhs <- paste(rhs_terms, collapse = " + ")
  lhs <- if (length(f) == 3L) deparse(f[[2]]) else NULL
  stats::as.formula(
    if (!is.null(lhs)) paste(lhs, "~", rhs) else paste("~", rhs),
    env = environment(f)
  )
}

#' Plug-in conditional log-likelihood of a flexyBayes fit
#'
#' Evaluates the conditional log-likelihood at the posterior-mean fitted
#' values -- a plug-in quantity, not a posterior summary, and not the
#' marginal likelihood. It exists so `AIC()`-style comparisons of two fits
#' from the same engine have something to read.
#'
#' The method computes what the object carries and refuses by name when it
#' carries too little. Three requirements are checked before any arithmetic
#' runs: the response vector in `$glm$y`, the fitted values from
#' [fitted()], and a recorded family. A Gaussian fit additionally sharpens
#' its residual scale from the posterior draws when they are present, and
#' otherwise falls back to the residual root-mean-square. Families outside
#' the Gaussian, binomial and Poisson set have no plug-in form here and are
#' refused rather than reported as `NA`, because a silent `NA` propagates
#' into `anova()` and `AIC()` as though the comparison had been made.
#'
#' INLA fits reach [logLik.flexybayes_inla()] instead, which states that
#' INLA reports a marginal log-likelihood and returns `NA` with that
#' message.
#'
#' @param object A `flexybayes` fit carrying a response vector, fitted
#'   values, and a recorded response family.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A `logLik` object: the scalar log-likelihood with `df` and
#'   `nobs` attributes taken from the fit's recorded model information.
#' @export
logLik.flexybayes <- function(object, ...) {
  y <- object$glm$y
  fam_link <- object$extras$parse_info$family
  if (is.null(y) || is.null(fam_link$family)) {
    stop(.fb_refusal_condition(
      reason_code = "conditional_loglik_not_available",
      message = paste0(
        "logLik() cannot be computed for this fit: it carries ",
        if (is.null(y)) "no response vector (`$glm$y`)" else
          "no recorded response family",
        ". The fit's class is <", paste(class(object), collapse = ", "),
        ">. Returning NA here would let AIC() and anova() report a ",
        "comparison that was never made."
      )
    ))
  }

  fitted_values <- stats::fitted(object)

  # ---- Plug-in evaluation, one branch per supported family -------------
  if (identical(fam_link$family, "gaussian")) {
    # The posterior draws sharpen the residual scale where they exist. A
    # fit without them (an engine storing its posterior in another shape)
    # still has a well-defined plug-in scale from the residuals, so this
    # is a fallback rather than a refusal.
    sigma_e <- NA_real_
    draws <- object$greta$draws
    if (!is.null(draws)) {
      all_draws <- do.call(rbind, lapply(draws, as.matrix))
      sigma_cols <- grep("^sigma_e_atg", colnames(all_draws), value = TRUE)
      if (length(sigma_cols) > 0L) {
        sigma_e <- mean(all_draws[, sigma_cols[1L]])
      }
    }
    if (!is.finite(sigma_e)) {
      sigma_e <- sqrt(mean((y - fitted_values)^2))
    }
    ll <- sum(stats::dnorm(y, mean = fitted_values, sd = sigma_e, log = TRUE))
  } else if (fam_link$family %in% c("binomial", "binary")) {
    p <- pmax(pmin(fitted_values, 1 - 1e-10), 1e-10)
    ll <- sum(stats::dbinom(y, size = 1L, prob = p, log = TRUE))
  } else if (identical(fam_link$family, "poisson")) {
    ll <- sum(stats::dpois(y, lambda = pmax(fitted_values, 1e-10), log = TRUE))
  } else {
    stop(.fb_refusal_condition(
      reason_code = "conditional_loglik_not_available",
      message = paste0(
        "logLik() has no plug-in conditional log-likelihood for the '",
        fam_link$family, "' family. The implemented families are ",
        "gaussian, binomial and poisson. Use the engine's own model ",
        "comparison instead -- loo::loo() on brms draws, or the DIC and ",
        "WAIC that summary() reports for an INLA fit."
      )
    ))
  }

  structure(
    ll,
    df = object$extras$model_info$n_params,
    nobs = object$extras$model_info$n_obs,
    class = "logLik"
  )
}

#' Number of observations a flexyBayes fit was fitted to
#'
#' Reads the observation count recorded at emit time, so the answer is the
#' number of rows the engine actually saw. On a fit whose missing responses
#' were augmented rather than dropped, that count includes the augmented
#' design cells.
#'
#' @param object A `flexybayes` fit of any engine.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A single integer, the number of observations in the fitted
#'   data.
#' @export
nobs.flexybayes <- function(object, ...) {
  object$extras$model_info$n_obs
}

#' Extract model family
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns The response family the model was fitted under, as the
#'   `family` object the fit recorded at emit time -- so it names the
#'   family and link the engine actually used, not the string the caller
#'   passed.
#' @export
family.flexybayes <- function(object, ...) {
  object$glm$family
}

#' Extract model formula
#' @param x A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns The fixed-effect (population-level) formula as a `formula`
#'   object. Random-effect and residual-structure terms are not part of
#'   it: read them from the fit's `fb_terms` intermediate representation,
#'   or from [fb_plan()] before fitting.
#' @export
formula.flexybayes <- function(x, ...) {
  x$glm$formula
}

#' Fixed-effect model matrix of a flexyBayes fit
#'
#' Rebuilds the population-level design matrix from the formula and data
#' the fit carries. Random-effect and residual-structure terms are not part
#' of it: this is the basis the fixed-effect coefficients are expressed in.
#'
#' The formula and data are resolved from whichever slots the object holds.
#' The greta-shaped and brms-shaped emits keep both under `$glm`; an INLA
#' fit keeps its data at `$data` and recovers its fixed-effect formula
#' through [formula.flexybayes_inla()]. An object supplying neither is
#' refused by name, since `model.matrix()` on a `NULL` formula silently
#' returns the intercept-only matrix of the calling frame's data.
#'
#' @param object A `flexybayes` fit carrying a fixed-effect formula and the
#'   data it was fitted to.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A numeric matrix with one row per observation and one column
#'   per fixed-effect basis column, carrying the usual `assign` attribute.
#' @export
model.matrix.flexybayes <- function(object, ...) {
  form <- object$glm$formula
  if (is.null(form)) {
    form <- tryCatch(stats::formula(object), error = function(e) NULL)
  }
  dat <- object$glm$data %||% object$data

  if (is.null(form) || is.null(dat)) {
    stop(.fb_refusal_condition(
      reason_code = "model_matrix_not_recoverable",
      message = paste0(
        "model.matrix() cannot rebuild the fixed-effect design for this ",
        "fit: it carries ",
        if (is.null(form)) "no fixed-effect formula" else "no fitted data",
        ". The fit's class is <", paste(class(object), collapse = ", "),
        ">. Rebuild the design from the original call instead."
      )
    ))
  }

  stats::model.matrix(form, data = dat)
}

#' Re-fit a flexyBayes model with modified arguments
#'
#' Rebuilds the original [flexybayes()] call from the argument record the
#' fit carries, applies the overrides supplied in `...`, and re-fits.
#'
#' The record has to be complete for this to be safe. An argument the fit
#' did not record would be re-supplied as its default, so a re-fit could
#' quietly drop a relationship matrix or a prior scale and return a
#' different model under the same name. When any of the thirteen recorded
#' arguments is absent, the method refuses and names what is missing rather
#' than re-fitting a model the user did not ask for. INLA fits currently
#' record six of the thirteen and therefore refuse.
#'
#' @param object A `flexybayes` fit carrying a complete argument record in
#'   `$extras$call_info`.
#' @param ... Named arguments overriding the recorded call, for example
#'   `n_samples = 2000L` or a replacement `data` frame.
#' @returns A new fitted `flexybayes` object of the same engine as the
#'   original, unless an override routes it elsewhere.
#' @export
update.flexybayes <- function(object, ...) {
  cl <- object$extras$call_info
  required <- c(
    "fixed", "random", "residual", "family", "link", "known_matrices",
    "weights", "n_samples", "warmup", "chains", "prior_fixed_sd",
    "prior_vc_sd"
  )
  missing_fields <- setdiff(required, names(cl))
  dat <- object$glm$data %||% object$data

  if (length(missing_fields) > 0L || is.null(dat)) {
    stop(.fb_refusal_condition(
      reason_code = "update_call_not_reconstructable",
      message = paste0(
        "update() cannot re-fit this model: its recorded call is ",
        "incomplete, so a re-fit would silently substitute defaults for ",
        if (length(missing_fields) > 0L) {
          paste0("`", paste(missing_fields, collapse = "`, `"), "`")
        } else {
          "the fitted data"
        },
        ". The fit's class is <", paste(class(object), collapse = ", "),
        ">. Re-issue the original flexybayes() call with the arguments ",
        "you want changed."
      )
    ))
  }

  dots <- list(...)
  args <- list(
    fixed = cl$fixed,
    random = cl$random,
    residual = cl$residual,
    data = dat,
    family = cl$family,
    link = cl$link,
    known_matrices = cl$known_matrices,
    weights = cl$weights,
    n_samples = cl$n_samples,
    warmup = cl$warmup,
    chains = cl$chains,
    # Carried so a re-fit repeats the sampler settings the original used.
    # Absent from a pre-0.9.0 record, where NULL is the old behaviour.
    seed = cl$seed,
    control = cl$control,
    prior_fixed_sd = cl$prior_fixed_sd,
    prior_vc_sd = cl$prior_vc_sd
  )

  # Override with user-supplied arguments
  for (nm in names(dots)) {
    args[[nm]] <- dots[[nm]]
  }

  do.call(flexybayes, args)
}

#' Compare flexyBayes models on a plug-in information criterion
#'
#' Ranks two or more fits by a DIC-shaped criterion built from the plug-in
#' conditional log-likelihood of [logLik.flexybayes()] and the recorded
#' parameter count. The criterion is approximate and the method says so on
#' every print: it penalises the nominal parameter count rather than an
#' effective one, so it is a coarse ordering, not a model-selection
#' procedure.
#'
#' Every fit compared must supply both ingredients. A fit whose
#' log-likelihood is `NA` -- an INLA fit reports a marginal likelihood, not
#' a conditional one -- or which records no parameter count is refused by
#' name, because ranking on an `NA` produces an ordering that looks
#' computed and is not.
#'
#' @param object A `flexybayes` fit supplying a conditional log-likelihood
#'   and a recorded parameter count.
#' @param ... Further `flexybayes` fits to compare against `object`.
#' @returns Invisibly, a data frame with one row per model carrying the
#'   log-likelihood, parameter count, criterion value, and the difference
#'   from the best model. Printed as a side effect.
#' @export
anova.flexybayes <- function(object, ...) {
  models <- c(list(object), list(...))
  n_models <- length(models)

  # ---- Both ingredients must exist for every fit ------------------------
  # logLik() on an INLA fit returns NA by design (marginal, not
  # conditional), and an INLA fit records no n_params. Either alone makes
  # the criterion unformable, so both are checked before any arithmetic.
  lls <- vapply(
    models,
    function(m) as.numeric(stats::logLik(m)),
    numeric(1)
  )
  n_params <- vapply(
    models,
    function(m) as.integer(m$extras$model_info$n_params %||% NA_integer_),
    integer(1)
  )
  unusable <- which(!is.finite(lls) | is.na(n_params))
  if (length(unusable) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "conditional_loglik_not_available",
      message = paste0(
        "anova() cannot rank these fits: model",
        if (length(unusable) == 1L) " " else "s ",
        paste(unusable, collapse = ", "),
        " supplied no finite conditional log-likelihood or no parameter ",
        "count. INLA fits are the usual case -- INLA reports a marginal ",
        "log-likelihood, and its DIC and WAIC are already in summary(). ",
        "For brms fits, loo::loo_compare() is the supported comparison."
      )
    ))
  }

  # Approximate DIC-like comparison
  dic <- -2 * lls + 2 * n_params

  model_names <- paste0("Model ", seq_len(n_models))

  result <- data.frame(
    Model = model_names,
    logLik = round(lls, 2),
    npar = n_params,
    DIC = round(dic, 2),
    delta = round(dic - min(dic), 2),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  cat("Bayesian model comparison\n")
  cat(strrep("-", 50), "\n")
  print(result)
  cat(
    "\nNote: this criterion is approximate -- it penalises the nominal\n",
    "parameter count, not an effective one. For a defensible comparison\n",
    "of brms fits use loo::loo_compare(); an INLA fit already reports\n",
    "DIC and WAIC in summary().\n"
  )

  invisible(result)
}


# ---------------------------------------------------------------- #
# flexybayes_direct_greta -- subclass method overrides              #
# ---------------------------------------------------------------- #
# fb_greta() returns objects of class
#   c("flexybayes_direct_greta", "flexybayes", "list")
# (subclass first so S3 dispatch finds the direct-greta overrides
# before the parent flexybayes methods). Three methods dispatch to
# greta-direct-aware overrides; everything else (summary, vcov,
# confint, prior_summary, posterior_samples) inherits verbatim from
# the parent "flexybayes" class.

#' Print a flexybayes object built via fb_greta() (direct greta entry)
#'
#' Distinguishes the direct-greta entry from the asreml-style /
#' brms-style formula entries by listing target parameters and the
#' inferred likelihood family.
#'
#' @param x A `flexybayes_direct_greta` object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`.
#' @export
print.flexybayes_direct_greta <- function(x, ...) {
  ci <- x$extras$call_info
  mi <- x$extras$model_info
  cm <- mi$canonical_map

  cat("Bayesian model  [flexyBayes / direct-greta entry]\n")
  cat(strrep("-", 55), "\n")
  cat("  Entry  : fb_greta() (user-built greta model)\n")
  cat("  Family :", mi$family, "(likelihood:", mi$likelihood, ")\n")
  cat(
    "  MCMC   :",
    ci$chains,
    "chain(s) x",
    ci$n_samples,
    "samples",
    "(warmup =",
    ci$warmup,
    ") --",
    round(x$extras$run_time, 1),
    "sec\n"
  )
  cat(
    "  Params :",
    mi$n_params,
    "monitored",
    "(",
    length(cm),
    "target greta_arrays)\n"
  )

  if (length(cm) > 0L) {
    cat("  Targets:", paste(names(cm), collapse = ", "), "\n")
    has_renames <- any(names(cm) != unname(cm))
    if (has_renames) {
      cat("  Canonical map:\n")
      for (i in seq_along(cm)) {
        if (names(cm)[i] != cm[i]) {
          cat("    ", names(cm)[i], " -> ", cm[i], "\n", sep = "")
        }
      }
    }
  }

  if (!is.null(x$extras$convergence$gelman)) {
    rhat <- x$extras$convergence$gelman$psrf[, "Point est."]
    max_rhat <- max(rhat, na.rm = TRUE)
    flag <- if (max_rhat < 1.05) {
      " [OK]"
    } else if (max_rhat < 1.10) {
      " [borderline]"
    } else {
      " [!]"
    }
    cat("  Max Rhat:", round(max_rhat, 3), flag, "\n")
  }
  min_eff <- if (!is.null(x$extras$convergence$n_eff)) {
    min(x$extras$convergence$n_eff, na.rm = TRUE)
  } else {
    NA
  }
  if (!is.na(min_eff)) {
    cat("  Min ESS:", round(min_eff, 0), "\n")
  }

  cat(strrep("-", 55), "\n")
  cat(
    "  $glm    -- minimal GLM-compatible shim ",
    "(posterior means of target params)\n",
    sep = ""
  )
  cat("  $greta  -- native greta (model, draws, greta_arrays)\n")
  cat("  $extras -- diagnostics, canonical_map, backend_decision\n")

  invisible(x)
}

#' Extract canonical-named posterior means from an fb_greta() fit
#'
#' Returns the posterior means of the target greta_arrays under
#' their canonical names. Differs from `coef.flexybayes()` only in
#' the source of the names (canonical map rather than fixed-effect
#' formula terms).
#'
#' @param object A `flexybayes_direct_greta` object.
#' @param ... Additional arguments (ignored).
#' @return Named numeric vector of posterior mean target parameters.
#' @export
coef.flexybayes_direct_greta <- function(object, ...) {
  # The parent class stores posterior means on object$glm$coefficients
  # already keyed by canonical names (renamed at fit time in
  # fb_greta()). Return that verbatim.
  object$glm$coefficients
}

#' Predict from a flexybayes_direct_greta fit
#'
#' Direct-greta fits do not encode a Wilkinson-Rogers formula on the
#' IR, so `predict()` cannot mechanically apply the model matrix to
#' `newdata`. The user must supply an explicit predictor function
#' `f(theta, newdata)` mapping a length-`p` named parameter vector
#' and a data frame to a length-`nrow(newdata)` vector of predicted
#' values. Posterior-mean prediction is computed by applying `f()`
#' to each posterior draw and averaging.
#'
#' @param object A `flexybayes_direct_greta` object.
#' @param newdata A data frame.
#' @param predictor A function `function(theta, newdata) -> numeric`
#'   where `theta` is a named vector of canonical-named target
#'   parameters and `newdata` is the data frame passed in.
#' @param type Character; `"response"` (default; on the response
#'   scale) or `"link"` (on the linear-predictor scale). Currently
#'   identity link only; non-Gaussian families queued for v0.3.
#' @param n_draws Integer; number of posterior draws to use (default
#'   500). The predictor function is applied to each draw and the
#'   mean is returned.
#' @param ... Additional arguments (ignored).
#' @return Numeric vector of length `nrow(newdata)`.
#' @export
predict.flexybayes_direct_greta <- function(
  object,
  newdata,
  predictor,
  type = c("response", "link"),
  n_draws = 500L,
  ...
) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    stop(
      "`newdata` must be supplied to predict() on an fb_greta() ",
      "fit (the IR does not encode a Wilkinson-Rogers formula ",
      "that predict() can mechanically apply).",
      call. = FALSE
    )
  }
  if (missing(predictor) || !is.function(predictor)) {
    stop(
      "`predictor` must be a function `function(theta, newdata) ",
      "-> numeric` mapping a named parameter vector and newdata to ",
      "a vector of predictions. See ?predict.flexybayes_direct_greta ",
      "for the contract.",
      call. = FALSE
    )
  }

  draws <- object$greta$draws
  all_draws <- do.call(rbind, lapply(draws, as.matrix))

  # Apply the canonical-name rename to the draw matrix so the
  # predictor sees the canonical names (consistent with coef()).
  cm <- object$extras$model_info$canonical_map
  greta_names <- names(cm)
  canon <- unname(cm)
  for (i in seq_along(greta_names)) {
    if (identical(greta_names[i], canon[i])) {
      next
    }
    pat <- paste0("^", greta_names[i], "(\\[|$)")
    repl <- paste0(canon[i], "\\1")
    colnames(all_draws) <- sub(pat, repl, colnames(all_draws))
  }

  n_total <- nrow(all_draws)
  if (n_draws >= n_total) {
    idx <- seq_len(n_total)
  } else {
    idx <- sort(sample.int(n_total, n_draws, replace = FALSE))
  }

  preds_mat <- vapply(
    idx,
    function(k) {
      theta <- setNames(as.numeric(all_draws[k, ]), colnames(all_draws))
      out <- predictor(theta, newdata)
      if (!is.numeric(out) || length(out) != nrow(newdata)) {
        stop(
          "`predictor(theta, newdata)` must return a numeric vector ",
          "of length nrow(newdata) = ",
          nrow(newdata),
          "; got length ",
          length(out),
          ".",
          call. = FALSE
        )
      }
      out
    },
    numeric(nrow(newdata))
  )

  if (is.null(dim(preds_mat))) {
    preds_mat <- matrix(preds_mat, nrow = 1L)
  }
  rowMeans(preds_mat)
}
