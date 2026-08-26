# Standard model interface for per-row INLA fits.
#
# A per-row INLA fit carries class
# `c("flexybayes_inla", "flexybayes", "list")`. Before 0.9.0 it stood
# outside the `flexybayes` parent, because its internal shape differs from
# the greta fit -- it has no `$glm` shim -- and that left it without the
# coef() / vcov() / predict() / formula() / family() interface the rest
# of the R modelling ecosystem (emmeans, marginaleffects) dispatches on.
# The methods below are the INLA-specific half of that interface and win
# dispatch over the parent; the parent methods with no sibling here
# resolve from the slots the object carries or refuse by name (see
# R/methods.R). This file supplies those methods, reading the
# fixed-effect posterior
# from the INLA fit's `summary.fixed` (means, treatment-contrast basis)
# and a Monte-Carlo joint covariance from `inla.posterior.sample()`
# (the fit is built with `control.compute = list(config = TRUE)`).
#
# Scope: the fixed-effect (population-level) surface. Random effects are
# held at their population mean (zero) for prediction on new data, which
# is the convention emmeans / marginaleffects assume for marginal means
# and average predictions.

# ---------------------------------------------------------------- #
# Internal: fixed-effect posterior on the treatment-contrast basis  #
# ---------------------------------------------------------------- #

# Monte-Carlo draws of the fixed-effect coefficients from an INLA fit,
# as a matrix (n_samples x p) with `summary.fixed` rownames as columns.
# INLA names the fixed effects in the latent field as "<name>:1".
.inla_fixef_draws <- function(object, n_samples = 2000L) {
  .check_installed(
    "INLA",
    "Package 'INLA' is required to sample from a ",
    "<flexybayes_inla> fit."
  )
  if (is.null(object$inla)) {
    stop("Cannot sample: the INLA fit object is missing.", call. = FALSE)
  }

  fixed_names <- rownames(object$inla$summary.fixed)
  # Force single-threaded sampling. inla.posterior.sample() otherwise
  # spawns one worker per requested sample, which (a) trips the core
  # limit R CMD check enforces via _R_CHECK_LIMIT_CORES_ and (b) is
  # wasteful for the modest sample counts the accessor draws. The fit
  # is built with control.compute = list(config = TRUE) (see
  # R/emit_inla.R), so config is not the failure mode here.
  samples <- tryCatch(
    INLA::inla.posterior.sample(
      as.integer(n_samples),
      object$inla,
      num.threads = "1:1"
    ),
    error = function(e) {
      stop(
        "INLA::inla.posterior.sample() failed: ",
        conditionMessage(e),
        ". Ensure the fit was built with ",
        "control.compute = list(config = TRUE).",
        call. = FALSE
      )
    }
  )

  latent_rows <- rownames(samples[[1L]]$latent)
  row_for <- vapply(
    fixed_names,
    function(nm) {
      hit <- which(latent_rows == paste0(nm, ":1"))
      if (length(hit) != 1L) {
        stop(
          "Could not locate the latent row for fixed effect '",
          nm,
          "' in the INLA posterior sample.",
          call. = FALSE
        )
      }
      latent_rows[hit]
    },
    character(1L)
  )

  mat <- vapply(
    samples,
    function(s) s$latent[row_for, 1L],
    numeric(length(fixed_names))
  )
  mat <- t(mat)
  # C4/FS-26: `fixed_names` above is deliberately the legalised name --
  # that is what the INLA latent field actually carries and what the
  # `paste0(nm, ":1")` lookup has to match. The column names on the
  # RETURNED matrix are the user's own, restored here at the very end,
  # so vcov.flexybayes_inla()'s dimnames agree with
  # coef.flexybayes_inla()'s (also restored) names -- the two are
  # indexed against each other (e.g. `V[colnames(X), colnames(X)]` in
  # the emmeans / marginaleffects seam) and must speak the same labels.
  colnames(mat) <- .inla_restore_term_labels(object$level_labels, fixed_names)
  mat
}


#' Coefficients of a per-row INLA fit
#'
#' Posterior means of the fixed effects, read from the INLA fit's
#' `summary.fixed` slot (treatment-contrast basis). These are the
#' coefficients consumed by [emmeans::emmeans()] and
#' [marginaleffects::predictions()] via the flexyBayes support methods.
#'
#' An INLA fit carries no `$glm` slot, so the fixed vector is read here
#' rather than inherited; every other value of `what` is resolved by the
#' shared body [coef.flexybayes()] uses, so the two engines answer the
#' same question the same way.
#'
#' A factor with a non-syntactic level (for example a level containing a
#' space, `"low N"`) is legalised internally (`make.names()`) before the
#' INLA emit, so the fit itself never dies on it; this method restores
#' the user's own, unlegalised label in the names it returns, as do
#' [ranef()], `summary()` and `predict(classify = )` on the same fit.
#'
#' @param object A `flexybayes_inla` fit.
#' @param what Which part of the fit to return: `"fixed"` (the default),
#'   `"random"`, `"missing"` or `"all"`. See [coef.flexybayes()] for the
#'   shape of each.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return For the default `what = "fixed"`, a named numeric vector of
#'   fixed-effect posterior means; otherwise as documented for
#'   [coef.flexybayes()].
#' @export
coef.flexybayes_inla <- function(
  object,
  what = c("fixed", "random", "missing", "all"),
  ...
) {
  sf <- object$inla$summary.fixed
  # C4/FS-26: `rownames(sf)` carries the legalised (make.names()) level
  # names the fit was actually built against when the model touches a
  # non-syntactic factor level; restore the user's own labels here so
  # coef() -- and ranef(), which delegates to coef(what = "random") via
  # R/methods.R -- print what the user wrote.
  fixed_names <- .inla_restore_term_labels(
    object$level_labels,
    rownames(sf)
  )
  .fb_coef_what(
    object,
    match.arg(what),
    stats::setNames(sf$mean, fixed_names)
  )
}

#' Posterior covariance of a per-row INLA fit's fixed effects
#'
#' Monte-Carlo estimate of the joint posterior covariance of the fixed
#' effects, computed from `inla.posterior.sample()`. The marginal
#' standard deviations match `summary.fixed$sd`; the off-diagonals carry
#' the joint dependence that contrast / marginal-mean standard errors
#' require. Because the estimate is sampling-based it varies slightly
#' between calls; raise `n_samples` for a tighter estimate.
#'
#' @param object A `flexybayes_inla` fit.
#' @param n_samples Posterior sample size for the covariance estimate
#'   (default 2000).
#' @param ... Ignored. Present for compatibility with the generic.
#' @return Posterior covariance matrix of the fixed effects, with
#'   `summary.fixed` rownames as dimnames.
#' @export
vcov.flexybayes_inla <- function(object, n_samples = 2000L, ...) {
  draws <- .inla_fixef_draws(object, n_samples = n_samples)
  v <- stats::cov(draws)
  dimnames(v) <- list(colnames(draws), colnames(draws))
  v
}

#' Credible intervals for the fixed effects of a per-row INLA fit
#'
#' Quantiles of INLA's own posterior marginals for the fixed effects, not
#' frequentist confidence intervals and not a normal approximation to
#' them. The bounds come from `INLA::inla.qmarginal()` applied to the
#' marginal densities the fit stores, so any credible level is available
#' rather than only the 0.95 that `summary.fixed` tabulates.
#'
#' Before 0.9.0 an INLA fit did not inherit the `flexybayes` parent class,
#' so `confint()` on one reached `stats::confint.default` and failed on a
#' missing `vcov` contract. This method is the INLA sibling of
#' [confint.flexybayes()].
#'
#' @param object A `flexybayes_inla` fit carrying fixed-effect marginals.
#' @param parm Character vector of coefficient names to return, or `NULL`
#'   (the default) for every fixed effect.
#' @param level Credible level for the interval, as a proportion. The
#'   default `0.95` returns the 2.5th and 97.5th posterior percentiles.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns A matrix with one row per fixed effect and two columns holding
#'   the lower and upper credible bounds, named for the percentiles used.
#' @export
confint.flexybayes_inla <- function(
  object,
  parm = NULL,
  level = 0.95,
  ...
) {
  marg <- object$inla$marginals.fixed
  if (is.null(marg) || length(marg) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "fit_lacks_posterior_draws",
      message = paste0(
        "confint() cannot form credible intervals for this INLA fit: it ",
        "carries no fixed-effect posterior marginals. Refit with ",
        "control.compute = list(return.marginals = TRUE), which is the ",
        "flexyBayes default, or read the 0.95 bounds that summary() ",
        "already reports."
      )
    ))
  }

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  bounds <- t(vapply(
    marg,
    function(m) vapply(probs, INLA::inla.qmarginal, numeric(1), marginal = m),
    numeric(2)
  ))
  # C4/FS-26: `names(marg)` is the legalised name; restore the user's
  # own here so confint()'s rownames stay consistent with coef()'s and
  # vcov()'s (also restored) names on the same fit.
  dimnames(bounds) <- list(
    .inla_restore_term_labels(object$level_labels, names(marg)),
    paste0(round(probs * 100, 1), "%")
  )

  if (!is.null(parm)) {
    bounds <- bounds[parm, , drop = FALSE]
  }
  bounds
}

#' Fixed-effect model formula of a per-row INLA fit
#'
#' The fixed-effect (population-level) formula, recovered from the
#' captured call. Random-effect and residual-structure terms are not
#' part of this formula.
#'
#' @param x A `flexybayes_inla` fit.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return A `formula`.
#' @export
formula.flexybayes_inla <- function(x, ...) {
  f <- x$extras$call_info$fixed
  if (is.null(f)) {
    stop("The INLA fit does not carry a fixed-effect formula.", call. = FALSE)
  }
  f <- stats::as.formula(f)
  # Enforce the documented contract: this is the *fixed-effect* formula.
  # A brms-grammar fit (`y ~ x + (1 | g)`) is stored whole in
  # call_info$fixed, so strip any random-effect bar terms here. The
  # downstream model-matrix reconstruction (predict / emmeans /
  # marginaleffects) must see only the fixed-effect basis; a leftover
  # `(1 | g)` evaluates as `1 | g` (a logical op on the factor) and
  # injects a spurious column that fails coefficient reconciliation.
  # Idempotent on an already-fixed-only formula.
  if (length(f) == 3L) {
    f <- .brms_fixed_only_formula(f)
  }
  f
}

#' Response family of a per-row INLA fit
#'
#' @param object A `flexybayes_inla` fit.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return A list with `family` and `link` entries.
#' @export
family.flexybayes_inla <- function(object, ...) {
  mi <- object$extras$model_info
  list(family = mi$family %||% "gaussian", link = mi$link %||% "identity")
}

#' In-sample fitted values from a per-row INLA fit
#'
#' Returns INLA's posterior-mean fitted values (response scale) for the
#' observed rows, taken from `summary.fitted.values`. Without this method
#' `fitted()` dispatched to `stats::fitted.default`, which silently returned
#' `NULL` for an INLA fit because the object does not populate a `$glm` slot.
#'
#' @param object A `flexybayes_inla` fit.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return A numeric vector of posterior-mean fitted values, one per
#'   observation.
#' @export
fitted.flexybayes_inla <- function(object, ...) {
  sf <- object$inla$summary.fitted.values
  if (is.null(sf) || is.null(sf$mean)) {
    stop(
      "fitted.flexybayes_inla(): this INLA fit does not carry fitted ",
      "values (no `$inla$summary.fitted.values`).",
      call. = FALSE
    )
  }
  out <- as.numeric(sf$mean)
  n_obs <- tryCatch(nrow(object$data), error = function(e) length(out))
  if (is.finite(n_obs) && length(out) > n_obs) {
    out <- out[seq_len(n_obs)]
  }
  out
}

#' Response residuals from a per-row INLA fit
#'
#' Observed response minus the posterior-mean fitted value (on the response
#' scale). The response is recovered from the fit's fixed-effect formula
#' evaluated against the stored data, so a transformed response
#' (`log(y) ~ ...`) residualises on the modelled scale.
#'
#' @param object A `flexybayes_inla` fit.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return A numeric vector of response residuals, one per observation.
#' @export
residuals.flexybayes_inla <- function(object, ...) {
  fit_vals <- fitted(object)
  f <- formula(object)
  y <- if (length(f) == 3L) {
    tryCatch(as.numeric(eval(f[[2L]], envir = object$data)),
             error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(y) || length(y) != length(fit_vals)) {
    stop(
      "residuals.flexybayes_inla(): could not recover the response on the ",
      "fitted scale from the stored data.",
      call. = FALSE
    )
  }
  y - fit_vals
}

#' Log-likelihood of a per-row INLA fit (not computed)
#'
#' INLA reports a *marginal* log-likelihood (the model evidence, available
#' through [summary()]), not the *conditional* model log-likelihood that the
#' `logLik()` generic denotes and that the greta / brms backends expose.
#' Returning the marginal quantity under the `logLik` name would conflate two
#' different things, so this method returns `NA` (with the degrees of
#' freedom and observation count filled in) and a one-line note. This also
#' lets downstream summaries (for example [glance()]) degrade gracefully
#' instead of erroring with "no applicable method".
#'
#' @param object A `flexybayes_inla` fit.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return A `logLik` object whose value is `NA_real_`, carrying `df` and
#'   `nobs` attributes.
#' @export
logLik.flexybayes_inla <- function(object, ...) {
  message(
    "logLik() is not defined for INLA fits: INLA reports a *marginal* ",
    "log-likelihood (model evidence; see summary()), not the conditional ",
    "log-likelihood logLik() denotes. Returning NA."
  )
  np <- tryCatch(nrow(object$inla$summary.fixed),
                 error = function(e) NA_integer_)
  n_obs <- tryCatch(nrow(object$data), error = function(e) NA_integer_)
  val <- NA_real_
  attr(val, "df") <- np
  attr(val, "nobs") <- n_obs
  class(val) <- "logLik"
  val
}

#' Population-level predictions from a per-row INLA fit
#'
#' Fixed-effect (population-level) predictions: the linear predictor is
#' `X beta` with random effects held at their population mean (zero).
#' On the identity link the response- and link-scale predictions
#' coincide. This is the prediction surface \pkg{marginaleffects} uses
#' for average predictions and slopes.
#'
#' The `classify` path is the same one [predict.flexybayes()] takes: a
#' marginal-means table built through the emmeans seam, whose interval on
#' this engine comes from the Gaussian approximation of the joint
#' fixed-effect posterior rather than from INLA's own marginals. The
#' printed table names that.
#'
#' @param object A `flexybayes_inla` fit.
#' @param newdata Optional data frame; defaults to the fit data. Ignored,
#'   with a warning, when `classify` is supplied.
#' @param type `"response"` or `"link"`.
#' @param se.fit Logical: also return delta-method standard errors from
#'   the fixed-effect covariance.
#' @param classify The factors to break a marginal-means table down by:
#'   a character value (`"Variety"`, `"Variety:env"`) or a one-sided
#'   formula (`~ Variety`). `NULL` (the default) is the historical
#'   behaviour.
#' @param level Credible level for the classify table's interval, as a
#'   proportion. Default `0.95`.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return With `classify`, a data frame of class `fb_predict_classify`.
#'   Otherwise a numeric vector of predictions, or a list `fit` /
#'   `se.fit` when `se.fit = TRUE`.
#' @export
predict.flexybayes_inla <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  se.fit = FALSE,
  classify = NULL,
  level = 0.95,
  ...
) {
  type <- match.arg(type)
  if (!is.null(classify)) {
    .fb_classify_newdata_note(newdata)
    return(.fb_predict_classify(object, classify, level))
  }
  fam <- family(object)
  if (!identical(fam$link, "identity")) {
    stop(
      "predict.flexybayes_inla() currently supports the identity ",
      "link only; this fit uses link '",
      fam$link,
      "'. Work with ",
      "the posterior draws via fb_as_draws_simple() for non-identity ",
      "links.",
      call. = FALSE
    )
  }

  if (is.null(newdata)) {
    # C4/FS-26: .fb_fit_data(), not the raw object$data -- the latter is
    # the level-legalised copy the fit was built against, and
    # coef(object) below now carries the user's own (restored) names.
    # Using the raw slot here would reintroduce the exact mismatch
    # .fb_fixef_model_matrix()'s own reconciliation check exists to
    # catch (the same one predict(classify = ) hit during development;
    # see R/ecosystem_support.R's .fb_fit_data()).
    newdata <- .fb_fit_data(object)
  }
  trms <- stats::delete.response(stats::terms(formula(object)))
  bhat <- coef(object)
  X <- .fb_fixef_model_matrix(trms, newdata, names(bhat), .fb_fit_data(object))
  eta <- as.numeric(X %*% bhat)

  if (!se.fit) {
    return(eta)
  }
  V <- vcov(object)
  se <- sqrt(rowSums((X %*% V) * X))
  list(fit = eta, se.fit = se)
}
