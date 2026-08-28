# S3 methods for the flexybayes class

#' Print a compact description of a flexyBayes fit
#'
#' Reports the call the fit came from, the model the representation
#' describes, the engine that ran it, the design and observation counts,
#' and the headline posterior quantities. Use [summary()] for the full
#' coefficient and variance-component tables.
#'
#' The header is shared with [print.flexybayes_inla()] and the brms
#' print, so the three cannot disagree about what the fit is. Sampler
#' lines print only on an engine that sampled: a nested Laplace
#' approximation has no chains and no warmup, and reporting them over one
#' tells a reader the fit is stochastic when it is not.
#'
#' @param x A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the description it
#'   prints.
#' @export
print.flexybayes <- function(x, ...) {
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
  # so a fit only ever advertises a slot it actually carries.
  cat(strrep("-", 55), "\n")
  .print_fit_components(
    x,
    c(
      glm = "GLM-compatible (summary, emmeans, etc.)",
      brms = "live brmsfit (loo, posterior_predict, summary)",
      inla = "raw INLA fit (use INLA's summary, plot, etc.)",
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
#' Prints the fixed-effect posterior summaries, the variance components
#' with the prior each one carried, the engine's own convergence
#' diagnostics, and whatever engine-native panels the fit earns. On a
#' brms fit with a sectioned residual it additionally prints the residual
#' variance and standard deviation for each level of the sectioning
#' factor, computed from the draws rather than by transforming a
#' posterior mean.
#'
#' The returned object is the same eleven-slot
#' `summary.flexybayes` object on every active engine, so
#' `summary(fit)$varcomp` answers whichever engine ran the fit. Before
#' 0.9.1 the two engines returned two incomparable objects -- INLA's
#' four-slot list and brms's `brmssummary` -- and neither carried a
#' variance-component table at all.
#'
#' A fit whose model carries an autoregressive latent field gains one
#' further slot, `spatial_field`, which is the field's own parameters on
#' the correlation and standard-deviation scales. It is the one
#' engine-native slot the object carries, it is present only where the
#' model has a field, and the eleven above it are on every fit.
#'
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Invisibly, an object of class
#'   `c("summary.flexybayes", "list")` carrying the slots below. Printed
#'   as a side effect.
#'
#'   \describe{
#'     \item{`fixed`}{Data frame: `term`, `estimate`, `std.error`,
#'       `conf.low`, `conf.high`. Posterior mean, posterior standard
#'       deviation and a credible interval, not a sampling-theory
#'       estimate and its confidence interval.}
#'     \item{`varcomp`}{Data frame: `component`, `estimate`,
#'       `std.error`, `conf.low`, `conf.high`, `prior`, `note`. One row
#'       per variance component (and per correlation on a fit carrying a
#'       latent field), on the standard-deviation scale, under the
#'       canonical component name. `prior` is a projection of
#'       [prior_summary()]. `note` is `"collapsed"` when the component's
#'       97.5% quantile sits below 1% of the posterior-median residual
#'       standard deviation -- a display heuristic for a posterior piled
#'       against zero, not a test.}
#'     \item{`random`}{Named list, one data frame per grouping factor,
#'       with the columns `group`, `level`, `estimate`, `std.error`,
#'       `conf.low`, `conf.high`. Empty when the model carries no random
#'       terms. The same object [ranef()] returns.}
#'     \item{`missing`}{Data frame of the unobserved design cells --
#'       ASReml's `mv` factor -- with the columns `row`, `estimate`,
#'       `std.error`, `conf.low`, `conf.high`, followed by any design
#'       index variables the fit recorded. `row` indexes the data the
#'       engine was handed. The posterior is the engine's own: INLA's
#'       fitted-value marginal at that row, brms's sampled `mi()`
#'       response. Zero rows -- never `NULL` -- when every response was
#'       observed. The same table [coef()] returns for
#'       `what = "missing"`.}
#'     \item{`converge`}{List, in the engine's own terms: R-hat,
#'       effective sample sizes and divergent transitions where a sampler
#'       ran; mode status, marginal likelihood and the largest
#'       Kullback-Leibler divergence where INLA's Laplace approximation
#'       did. No R-hat is invented for an approximation that has none.}
#'     \item{`n_design`}{Rows the engine was handed.}
#'     \item{`n_observed`}{How many of those carried an observed
#'       response.}
#'     \item{`na_action`}{The missing-response record the fit carries, or
#'       `NULL` on a fit assembled without that layer.}
#'     \item{`model`}{One-line description of the random and residual
#'       structure, derived from the model representation rather than
#'       from any engine's emitted formula.}
#'     \item{`engine`}{`"inla"` or `"brms"`.}
#'     \item{`call`}{The recorded call.}
#'     \item{`spatial_field`}{**Present only on a fit carrying an
#'       autoregressive latent field**, and absent -- not `NULL`-valued
#'       -- otherwise. A data frame of `parameter`, `median`, `lower`,
#'       `upper`: the field's correlation and standard-deviation
#'       parameters and the nugget standard deviation, on the scale a
#'       reader thinks in rather than INLA's precision scale. The point
#'       estimate is the posterior median, because the standard
#'       deviations are read off precision marginals and only a monotone
#'       summary survives the reciprocal square root exactly. The same
#'       table the printed field panel renders, unrounded.}
#'   }
#' @seealso [prior_summary()] for the full resolved prior, and
#'   [nobs()] for the design and observed counts on their own.
#' @export
summary.flexybayes <- function(object, ...) {
  out <- .fb_summary_object(object)
  print(out)
  invisible(out)
}

# .agg_prior_parametrization() --- which prior an aggregated fit ran
# under, as one token.
#
# The predicate used to be `inherits(fb$priors, "fb_prior")` alone, which
# was right while the only fb_prior on a fit was one the user wrote. Since
# the auto-default became an fb_prior object, that test says "custom" on
# an ordinary call nobody passed a prior to -- so the commonest fit in the
# package printed "explicit prior supplied" about a package default.
#
# The discriminator is the one prior_summary() already uses to report
# `default_origin`: the auto-default carries `fb_prior_default_basis`,
# recording which scale basis built it, and a hand-written fb_prior()
# does not. Reading the same attribute keeps the two surfaces from
# disagreeing about the same object.
#
# Three tokens, because there are three cases and only two of them used
# to be distinguishable:
#   "package_default"    the auto-default bounded-uniform-on-SD prior
#   "custom"             an explicit fb_prior() from the caller
#   "per_row_equivalent" the legacy scalar bridge (and no recorded prior)
#
# @noRd
# @keywords internal
.agg_prior_parametrization <- function(priors) {
  if (!inherits(priors, "fb_prior")) {
    return("per_row_equivalent")
  }
  if (!is.null(attr(priors, "fb_prior_default_basis"))) {
    return("package_default")
  }
  "custom"
}

# Human-readable label for the aggregated-fit prior parametrization,
# shared by the aggregated print + summary methods. The
# "per_row_equivalent" case reassures the user that the matched-prior
# guarantee holds (aggregated posterior == per-row posterior under the
# default prior); "package_default" names the prior the package injected
# without claiming that equivalence, which is not established for the
# bounded-uniform default on the aggregated route; the "custom" case
# flags that an explicit prior was supplied and points at
# prior_summary().
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
    package_default = paste0(
      "package default (bounded uniform on SD; see prior_summary())"
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
  # The aggregated path is INLA-only, so the component list carries the
  # `$inla` slot; still read off the object rather than assumed, in case a
  # future backend extends the aggregated path.
  cat(strrep("-", 60), "\n")
  .print_fit_components(
    x,
    c(
      glm = "per-row reconstructed fitted values + coef shim",
      inla = "raw aggregated INLA fit (use INLA's summary etc.)",
      extras = "summary, aggregation_meta, backend_decision"
    ),
    width = 10L
  )
  invisible(x)
}


#' Summarise a fit run on the aggregated representation
#'
#' Returns the same eleven-slot `summary.flexybayes` object every other
#' engine returns, so `summary(fit)$varcomp` answers on an aggregated fit
#' as it does on a per-row one. Aggregation is the default route for an
#' ordinary Gaussian call with a random term, so until 0.9.1 that slot
#' was `NULL` on the commonest fit the package produces: this method
#' returned its own list of the aggregated posterior's raw pieces and
#' carried no variance-component table at all.
#'
#' Two things are aggregated-specific and both appear in the printed
#' header rather than as extra slots. The banner names the representation
#' (`aggregated-gaussian`, `aggregated-binomial`, `aggregated-poisson`)
#' where a per-row fit names its engine, and the `aggregation` line gives
#' the row-to-cell compression the fit ran under. The `model` slot
#' carries the same compression in its text, because a reader comparing
#' the fixed-effect tables of an aggregated and a per-row fit otherwise
#' has nothing on the object saying the two were computed over different
#' numbers of rows.
#'
#' The variance components come from the engine's own posterior. On the
#' INLA route that is the hyperparameter marginal transformed to the
#' standard-deviation scale, the same construction the per-row route
#' uses, and not the reciprocal square root of a tabulated precision
#' mean. A route that recorded posterior means alone reports the means
#' with the three interval columns `NA` and the row's `note` cell
#' reading `"no interval recorded"`.
#'
#' The raw aggregated pieces this method used to return -- `beta_means`,
#' `beta_vcov`, `sigma_means`, `tau_means` -- are unchanged at
#' `fit$extras$summary`.
#'
#' @param object A `<flexybayes_aggregated>` object, as returned by a fit
#'   run on the aggregated representation.
#' @param ...    Ignored. Present for compatibility with the generic.
#' @returns Invisibly, an object of class
#'   `c("summary.flexybayes", "list")` with the slots documented at
#'   [summary.flexybayes()]. Printed as a side effect.
#' @seealso [summary.flexybayes()] for the slot-by-slot contract, and
#'   [print.flexybayes_aggregated()] for the one-screen fit description.
#' @export
summary.flexybayes_aggregated <- function(object, ...) {
  out <- .fb_summary_object(object)
  print(out)
  invisible(out)
}


#' Extract coefficients from a flexyBayes fit
#'
#' Returns the fixed effects by default, and the random-effect
#' predictions or the unobserved design cells on request. An ASReml user
#' reaches for `coef(fit)$random`; the equivalent here is
#' `coef(fit, what = "random")`, which is what [ranef()] calls.
#'
#' The default is the historical return -- a named numeric vector of the
#' fixed effects -- so every caller that treats `coef(fit)` as a numeric
#' vector keeps working, including the \pkg{emmeans} and
#' \pkg{marginaleffects} seams, [predict()], [vcov()] and the tidiers.
#'
#' @param object A fitted `flexybayes` object of any backend.
#' @param what Which part of the fit to return. `"fixed"` (the default)
#'   is the fixed-effect posterior means. `"random"` is one data frame
#'   per grouping factor, with the columns `group`, `level`, `estimate`,
#'   `std.error`, `conf.low` and `conf.high`; a grouping factor carrying
#'   more than one effect names the effect in `level`. `"missing"` is
#'   the table of unobserved design cells, in the columns `row`,
#'   `estimate`, `std.error`, `conf.low`, `conf.high` plus any design
#'   index variables the fit recorded. `"all"` is the named list of all
#'   three.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns For `what = "fixed"`, a named numeric vector of the fixed
#'   effects' posterior means on the treatment-contrast basis, named for
#'   the design-matrix columns, empty when the model carries no fixed
#'   effects. For `"random"`, a named list of data frames, empty when the
#'   model carries no random terms. For `"missing"`, a data frame with
#'   zero rows when every response was observed. For `"all"`, a list with
#'   the elements `fixed`, `random` and `missing`.
#' @seealso [ranef()] for the random-effect table on its own,
#'   [summary.flexybayes()] whose `random` and `missing` slots are the
#'   same two objects, and [nobs()] for the design and observed counts.
#' @export
coef.flexybayes <- function(
  object,
  what = c("fixed", "random", "missing", "all"),
  ...
) {
  .fb_coef_what(object, match.arg(what), object$glm$coefficients)
}

#' Random-effect predictions from a flexyBayes fit
#'
#' The random-effect table, one data frame per grouping factor, in the
#' same columns on every engine: `group`, `level`, `estimate`,
#' `std.error`, `conf.low`, `conf.high`. Identical to
#' `coef(object, what = "random")`, and present under this name because
#' it is the name the mixed-model ecosystem uses.
#'
#' The method is registered on both this package's own `ranef()` generic
#' and on \pkg{nlme}'s, which is the generic \pkg{lme4} and \pkg{brms}
#' re-export. Attaching any of those masks the local generic, and the
#' registration on nlme's is what keeps a bare `ranef(fit)` dispatching
#' there. In a session where the masking is ambiguous, call
#' `flexyBayes::ranef(fit)`.
#'
#' @param object A fitted `flexybayes` object of any backend.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A named list of data frames, one per grouping factor. Empty
#'   when the model carries no random terms.
#' @seealso [coef.flexybayes()], which this delegates to.
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' @rdname ranef
#' @export
#' @exportS3Method nlme::ranef
ranef.flexybayes <- function(object, ...) {
  coef(object, what = "random")
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
#' confidence intervals: the bounds would be empirical quantiles of the
#' fixed-effect posterior draws the fit carries.
#'
#' No active engine stores its posterior in the shape this bare fallback
#' method reads. A fit from an active engine reaches its own method --
#' [confint.flexybayes_inla()] for INLA, [confint.flexybayes_brms()] for
#' brms -- so this method always refuses by name rather than returning an
#' empty interval matrix, which would read as "no fixed effects" rather
#' than "this fit cannot answer the question".
#'
#' @param object A flexybayes fit.
#' @param parm Character vector of parameter names to return, or `NULL`
#'   (the default) for every fixed effect.
#' @param level Credible level for the interval, as a proportion. The
#'   default `0.95` returns the 2.5th and 97.5th posterior percentiles.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Does not return: raises the classed `fit_lacks_posterior_draws`
#'   refusal.
#' @export
confint.flexybayes <- function(object, parm = NULL, level = 0.95, ...) {
  stop(.fb_refusal_condition(
    reason_code = "fit_lacks_posterior_draws",
    message = paste0(
      "confint() cannot form credible intervals for this fit: it carries ",
      "no posterior-draw slot. The fit's class is <",
      paste(class(object), collapse = ", "), ">. Fixed-effect intervals ",
      "are available from summary() on every active engine, and an INLA ",
      "or brms fit answers confint() from its own method."
    )
  ))
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
#' [fitted()], and a recorded family. A Gaussian fit uses the residual
#' root-mean-square as its plug-in scale. Families outside the Gaussian,
#' binomial and Poisson set have no plug-in form here and are refused
#' rather than reported as `NA`, because a silent `NA` propagates into
#' `anova()` and `AIC()` as though the comparison had been made.
#'
#' brms fits reach [logLik.flexybayes_brms()] instead (brms's own
#' `log_lik()`), and INLA fits reach [logLik.flexybayes_inla()], which
#' states that INLA reports a marginal log-likelihood and returns `NA`
#' with that message.
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
    # The residual root-mean-square is the plug-in scale.
    sigma_e <- sqrt(mean((y - fitted_values)^2))
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
#' A fit whose missing responses were augmented rather than dropped was
#' handed more rows than it has observations: the design cell of a lost
#' plot is still a row, carried as a latent quantity so the index set a
#' structured covariance is built over survives. The two counts are
#' therefore different numbers and the argument says which one is wanted.
#'
#' `type = "design"` (the default, and the historical behaviour) is the
#' number of rows the engine saw. `type = "observed"` is how many of
#' those carried an observed response, read from the record the
#' missing-response layer left on the fit.
#'
#' @param object A `flexybayes` fit of any engine.
#' @param type Which count to return: `"design"` (default) or
#'   `"observed"`.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A single integer.
#' @seealso [summary.flexybayes()], whose `n_design` and `n_observed`
#'   slots are the same two numbers.
#' @export
nobs.flexybayes <- function(object, type = c("design", "observed"), ...) {
  type <- match.arg(type)
  n_design <- object$extras$model_info$n_obs
  if (identical(type, "design")) {
    return(n_design)
  }

  rec <- object$extras$na_action
  if (!is.null(rec$n_observed)) {
    return(rec$n_observed)
  }

  # No record: a fit built by calling an emit directly. Count the
  # observed responses from the data the fit carries rather than
  # returning the design count under the other name.
  dat <- object$glm$data %||% object$data
  resp <- .fb_fit_ir(object)$response %||%
    object$extras$parse_info$fixed$response
  if (!is.null(dat) && !is.null(resp) && resp %in% names(dat)) {
    return(sum(!is.na(dat[[resp]])))
  }
  warning(
    "nobs(type = \"observed\") cannot tell how many responses were ",
    "observed for this fit: it carries neither a missing-response record ",
    "nor recoverable data, so the count is NA rather than the design ",
    "count under another name.",
    call. = FALSE
  )
  NA_integer_
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
#' The brms-shaped emit keeps both under `$glm`; an INLA fit keeps its data
#' at `$data` and recovers its fixed-effect formula through
#' [formula.flexybayes_inla()]. An object supplying neither is refused by
#' name, since `model.matrix()` on a `NULL` formula silently returns the
#' intercept-only matrix of the calling frame's data.
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
#' quietly drop a relationship matrix, a prior scale, or a missing-response
#' policy, and return a different model under the same name. When any
#' recorded argument is absent, the method refuses and names what is
#' missing rather than re-fitting a model the user did not ask for. A fit
#' produced before its engine recorded the full set therefore refuses,
#' which is the intended behaviour: there is no pre-0.9.1 object to
#' special-case into a silent default.
#'
#' @section The prior a re-fit runs under:
#'
#' A re-fit uses the prior the original fit resolved to, taken from the
#' record on the fit itself rather than re-derived from the arguments.
#' This matters because the auto-default bounded uniform on the
#' standard-deviation scale fires only when the call supplies neither a
#' `prior` nor a `prior_vc_sd`, and a re-fit that re-issues the recorded
#' `prior_vc_sd` supplies one by construction. Before 0.9.1 the
#' auto-default therefore never fired on a re-fit: an identity `update()`
#' fell through to the engine's own hyperprior and could report a
#' variance component half the size of the one it re-fitted, with nothing
#' on the object to say so.
#'
#' Precedence, in order:
#'
#' \itemize{
#'   \item A `prior` in `...` is used as given. The recorded prior is not
#'     re-issued alongside it.
#'   \item A `prior_vc_sd` in `...` is used as given, and the recorded
#'     prior is *not* re-issued -- the two forms describe the same
#'     variance-component scale, and passing both would let the recorded
#'     prior silently override the scalar the caller just wrote.
#'   \item Otherwise the fit's own resolved prior decides, and *how* it
#'     decides depends on where that prior came from -- see below.
#'   \item A fit carrying the legacy scalar bridge, or no prior record at
#'     all, re-fits on the recorded scalars exactly as it did before.
#' }
#'
#' \strong{A policy re-fires; a bespoke prior carries.} The two things
#' `$extras$fb_terms$priors` can hold are not the same kind of thing, and
#' a re-fit treats them differently.
#'
#' \itemize{
#'   \item \strong{The auto-default} bounded uniform is the output of a
#'     policy -- one bounded uniform per variance component the model
#'     has, with the bound read off the data -- rather than a statement
#'     about particular terms. A re-fit passes neither `prior` nor either
#'     scalar, so the policy runs again over the *updated* model and
#'     data. On an identity `update()` the same data rebuild the same
#'     bound, so nothing moves. On `update(fit, random = ~ g + h)` the
#'     added term gets the same bounded uniform as its siblings.
#'   \item \strong{A user-supplied `fb_prior()`} is a statement about the
#'     terms it names, and is carried verbatim. A term the user did not
#'     name -- including one an override has just added -- keeps
#'     following the engine's own hyperprior, exactly as it did on the
#'     first fit. `prior_fixed_sd` and `prior_vc_sd` are re-issued
#'     alongside it, so a route that reads them for the terms the prior
#'     object does not name keeps the values the original ran under.
#' }
#'
#' The asymmetry is deliberate. Re-applying a policy's old output to a new
#' model leaves a mixed state -- the original terms priored, the added one
#' falling to the engine while its siblings keep the uniform -- which is
#' not what either the policy or the user asked for. Re-deciding a
#' *user's* explicit prior, on the other hand, would put words in their
#' mouth. Read [prior_summary()] on the re-fit to see what it ran under,
#' and [summary.flexybayes()]'s `varcomp` table to see it per component.
#'
#' @section The engine and representation a re-fit runs on:
#'
#' `backend` and `aggregate` are re-issued from the record alongside
#' everything else. Before 0.9.1 neither was recorded, so a re-fit took
#' the formal defaults: an identity `update()` of a Stan fit came back as
#' an aggregated INLA fit -- a different inference engine and a different
#' representation, under the same name and with nothing said.
#'
#' What is recorded is the value the *call* asked for, not the engine the
#' call resolved to. A fit made under `backend = "auto"` records `"auto"`,
#' and a re-fit routes again from there. A fit that named `"brms"` records
#' `"brms"` and comes back on Stan. The distinction is deliberate: the
#' recorded value is the user's policy, and a re-fit whose model has
#' changed -- `update(fit, random = ~ g + h)` -- has to be free to route
#' the changed model rather than be pinned to the engine that suited the
#' old one. `aggregate` behaves the same way. A recorded `"auto"` lets the
#' aggregation gate re-decide, a recorded `FALSE` keeps the re-fit
#' per-row.
#'
#' An override in `...` wins, as for every other recorded argument, so
#' `update(fit, backend = "brms")` moves the re-fit to Stan.
#'
#' What `update()` reproduces is the model and the policy behind it, not
#' the display settings of the session that first ran it. `verbose` is
#' recorded on the fit for completeness and is deliberately *not*
#' re-issued. A re-fit follows the current call's display default, and
#' `update(fit, verbose = FALSE)` quietens one on the spot. The
#' distinction is in what the argument can get wrong. A silently switched
#' engine is a different answer under the same name. A banner printed
#' where the first fit printed none is a banner.
#'
#' @param object A `flexybayes` fit carrying a complete argument record in
#'   `$extras$call_info`.
#' @param ... Named arguments overriding the recorded call, for example
#'   `n_samples = 2000L` or a replacement `data` frame.
#' @returns A new fitted `flexybayes` object of the same engine as the
#'   original, unless an override routes it elsewhere.
#' @seealso [prior_summary()] for the resolved prior a fit carries, and
#'   [summary.flexybayes()], whose `varcomp` table names the prior each
#'   component ran under.
#' @export
update.flexybayes <- function(object, ...) {
  cl <- object$extras$call_info
  required <- c(
    "fixed", "random", "residual", "family", "link", "known_matrices",
    "weights", "n_samples", "warmup", "chains", "prior_fixed_sd",
    "prior_vc_sd", "na_action", "backend", "aggregate"
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

  # D16 / WP-D (found by WP-G2's stress run). The override loop below
  # ("Override with user-supplied arguments") walks `names(dots)` and
  # assigns `args[[nm]] <- dots[[nm]]` for each. `stats::update()`'s
  # classic idiom is an UNNAMED replacement formula as the second
  # positional argument -- `update(fit, . ~ . + z)` -- and an unnamed
  # element has no entry in `names(dots)` to walk (a single unnamed
  # element makes `names(dots)` NULL outright; a mix of named and
  # unnamed makes it carry ""), so that loop silently visits zero
  # iterations for it. Confirmed live before this fix:
  # `update(fit, . ~ . + z)` returned a fit whose recorded `fixed`
  # formula, and whose coefficients, were IDENTICAL to the original --
  # a valid object answering a question nobody asked, not a re-fit.
  # This is not a documented spelling: every worked example in this
  # package (`README.md`'s accessor table, the roxygen above) uses a
  # named argument (`update(fit, random = ~ Block + Variety)`), so this
  # refuses the unsupported idiom by name rather than building a
  # formula-delta merger `stats::update.formula()`-style, which no
  # documented usage asks for.
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    has_unnamed <- is.null(dot_names) || any(!nzchar(dot_names))
    if (has_unnamed) {
      stop(.fb_refusal_condition(
        reason_code = "update_unnamed_argument_not_supported",
        message = paste0(
          "update() does not support the classic stats::update() idiom ",
          "of an unnamed replacement formula, for example ",
          "`update(fit, . ~ . + z)`. Every argument update() re-issues ",
          "is matched by name, so an unnamed argument has nothing to ",
          "match and was previously discarded silently, re-fitting the ",
          "unchanged model under the same call. Name the slot you want ",
          "changed instead: `update(fit, fixed = y ~ x + z)` for the ",
          "brms-style entry, or `random = ~ ...` / `residual = ~ ...` ",
          "on the ASReml grammar."
        ),
        family_class =
          "flexybayes_update_unnamed_argument_not_supported_refusal"
      ))
    }
  }

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
    prior_vc_sd = cl$prior_vc_sd,
    # An `omit` fit re-fits as `omit`. The default happens to be the same
    # word today, but a re-fit that quietly reverted to augmentation would
    # be computed on a different set of rows from the one it was named
    # after.
    na_action = cl$na_action,
    # The engine and the representation the original call asked for,
    # re-issued as the request rather than as the engine the request
    # resolved to -- see the Rd section.
    #
    # `verbose` is recorded next to these two and is deliberately NOT
    # re-issued. What update() reproduces is the model and the policy
    # behind it, not the display settings of the session that first ran
    # it: a re-fit on another engine, or in another representation, is a
    # different answer under the same name, while a re-fit that prints
    # its banner when the original did not is a re-fit that printed a
    # banner. The record keeps `verbose` because completeness there costs
    # nothing, and a re-fit follows the current call's display default.
    backend = cl$backend,
    aggregate = cl$aggregate
  )

  # ---- The prior the fit ran under, not the one its arguments imply ----
  #
  # flexybayes() decides whether the auto-default bounded uniform on the
  # SD scale fires from argument missingness:
  #
  #   default_prior_active <- is.null(prior) && missing(prior_vc_sd)
  #
  # and the block above re-issues `prior_vc_sd` from the record, so
  # `missing()` is FALSE on every re-fit and the gate never opens. The
  # re-fit then fell through to the legacy scalar bridge and let the
  # engine choose its own hyperprior -- an identity update() moved
  # sd_Block from 2.278 to 1.144 on a per-row INLA fit, from 1.558 to
  # 0.980 on the aggregated Gaussian route, and the same way on brms,
  # with $varcomp$prior changing from the resolved uniform to
  # "engine default" and nothing else saying so.
  #
  # Argument missingness is the wrong thing to re-derive the answer from,
  # because the fit already carries the answer: fb_terms$priors holds
  # whatever the original resolved to.
  #
  # WHAT IS CARRIED DEPENDS ON WHERE THE PRIOR CAME FROM, and the two
  # cases are not the same kind of thing.
  #
  # A resolved prior carries `fb_prior_default_basis` when the
  # auto-default built it -- the same attribute prior_summary() reads to
  # report `default_origin = "auto"`, and .fb_prior_record() to decide
  # `shared_default`. That object is not a choice the user made about
  # these terms. It is the output of a POLICY: a bounded uniform on the
  # SD scale, one per variance component the model has, with the bound
  # read off the data. Carrying the object re-applies yesterday's output
  # to today's model, and a model that gained a term gets a mixed state
  # -- the original terms priored, the new one falling to the engine
  # while its siblings keep the uniform. So the policy is re-fired
  # instead: neither `prior` nor either scalar is passed, `missing()`
  # holds on the gate above, and flexybayes() rebuilds the default over
  # the UPDATED model and data. On an identity re-fit the same data
  # rebuild the same bound, so this is indistinguishable from carrying
  # the object.
  #
  # A user-supplied fb_prior() is the opposite: a bespoke statement about
  # named terms, which is carried verbatim. A term the user did not name
  # keeps following the engine's own hyperprior, exactly as it did on the
  # first fit. The asymmetry is documented in the Rd.
  #
  # Dropping the scalars alongside the object is safe on both live
  # engines. INLA never reads `prior_fixed_sd`, and brms reads it only on
  # the legacy bridge, which by construction is not this branch.
  #
  # A caller who writes a prior of their own outranks the record. Both
  # spellings are honoured, and only one of them is passed: re-issuing
  # the recorded prior next to a caller's `prior_vc_sd` would let the
  # record win over the scalar just written, which is the same silent
  # substitution in the other direction.
  recorded_prior <- object$extras$fb_terms$priors
  caller_priced <- any(c("prior", "prior_vc_sd") %in% names(dots))
  if (inherits(recorded_prior, "fb_prior") && !caller_priced) {
    from_auto_default <- !is.null(
      attr(recorded_prior, "fb_prior_default_basis")
    )
    if (from_auto_default) {
      args$prior_fixed_sd <- NULL
      args$prior_vc_sd <- NULL
    } else {
      args$prior <- recorded_prior
    }
  }

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
