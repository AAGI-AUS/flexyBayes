# Standard model interface for per-row INLA fits.
#
# A per-row INLA fit carries class
# `c("flexybayes_inla", "flexybayes", "list")`. Before 0.9.0 it stood
# outside the `flexybayes` parent, because its internal shape differs from
# the brms fit (it has no `$glm` shim) and that left it without the
# coef() / vcov() / predict() / formula() / family() interface the rest
# of the R modelling ecosystem (emmeans, marginaleffects) dispatches on.
# The methods below are the INLA-specific half of that interface and win
# dispatch over the parent. The parent methods with no sibling here
# resolve from the slots the object carries or refuse by name (see
# R/methods.R). This file supplies those methods, reading the
# fixed-effect posterior
# from the INLA fit's `summary.fixed` (means, treatment-contrast basis)
# and the joint covariance from the Gaussian approximations INLA stores
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
    .fb_inla_posterior_sample(
      object,
      as.integer(n_samples),
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
  # `fixed_names` above is the legalised name: the INLA latent field
  # carries it and the `paste0(nm, ":1")` lookup has to match it. The
  # returned matrix gets the user's own labels, so its dimnames agree with
  # coef.flexybayes_inla(). The two are indexed against each other (e.g.
  # `V[colnames(X), colnames(X)]` in the emmeans / marginaleffects seam).
  colnames(mat) <- .inla_restore_term_labels(object$level_labels, fixed_names)
  # A column fixed_alias = "drop" fixed at zero draws within 1e-5 of it.
  aliased <- object$extras$fixed_alias$dropped$column
  mat[, colnames(mat) %in% aliased] <- 0
  mat
}

# Posterior weights of INLA's stored hyperparameter configurations.
.fb_inla_config_weights <- function(cs) {
  log_post <- vapply(cs$config, function(k) k$log.posterior, numeric(1L))
  w <- exp(log_post - max(log_post))
  w / sum(w)
}

# The Gaussian approximation of the latent field at configuration k: the
# Cholesky factor of its precision Q (the stored Qinv holds only the entries
# on the factor's sparsity pattern), its mean, and the kriging gain that
# conditions a draw on the linear constraints A x = e.
#
# `reuse`, from an earlier configuration, refactorises a Q of the same pattern
# on its ordering. The factor is the one Cholesky() would return.
#
# An intrinsic term (rw1, rw2, besag) leaves Q singular in a direction its
# sum-to-zero constraint fixes, such as a constant shared by the smooth and
# the intercept, and rounding can then make Q fail to factorise. On the
# constraint set A x = e, Q + c A'A gives the same exponent up to a constant
# for every c > 0, so conditioning a draw from it on A x = e gives the same
# distribution INLA::inla.posterior.sample() draws from Q.
.fb_inla_config_gaussian <- function(cs, k, reuse = NULL) {
  config <- cs$config[[k]]
  q <- Matrix::forceSymmetric(config$Q, uplo = "U")
  a_mat <- if (!is.null(cs$constr)) matrix(cs$constr$A, nrow = cs$constr$nc)
  factored <- .fb_inla_config_factor(q, reuse)
  if (is.null(factored) && !is.null(a_mat)) {
    # c puts c A'A on the scale of an average diagonal entry of Q.
    c_scale <- mean(Matrix::diag(q)) / max(rowSums(a_mat^2))
    q <- Matrix::forceSymmetric(
      q + Matrix::Matrix(c_scale * crossprod(a_mat), sparse = TRUE),
      uplo = "U"
    )
    factored <- .fb_inla_config_factor(q, reuse)
  }
  if (is.null(factored)) {
    stop(.fb_refusal_condition(
      reason_code = "inla_config_precision_not_positive_definite",
      message = sprintf(
        paste0(
          "The precision of the latent field INLA stored for ",
          "hyperparameter configuration %d is not positive definite%s. ",
          "The latent field's posterior is improper in a direction no ",
          "constraint of the model fixes, so vcov(), predict() standard ",
          "errors and posterior draws are not computed."
        ),
        k,
        if (is.null(a_mat)) "" else ", even on its linear constraints"
      )
    ))
  }
  chol_q <- factored$factor
  reuse <- factored$reuse
  out <- list(
    chol_q = chol_q,
    mu = config$improved.mean,
    gain = NULL,
    reuse = reuse
  )
  if (!is.null(a_mat)) {
    sigma_at <- as.matrix(Matrix::solve(chol_q, t(a_mat)))
    out$a_mat <- a_mat
    out$sigma_at <- sigma_at
    out$gain <- sigma_at %*% solve(a_mat %*% sigma_at)
  }
  out
}

# The Cholesky factor of q, or NULL when q is not positive definite.
.fb_inla_config_factor <- function(q, reuse) {
  same_pattern <- !is.null(reuse) &&
    identical(reuse$p, q@p) &&
    identical(reuse$i, q@i)
  not_pd <- function(cnd) grepl("not positive", conditionMessage(cnd))
  chol_q <- withCallingHandlers(
    tryCatch(
      if (same_pattern) {
        Matrix::update(reuse$factor, q)
      } else {
        Matrix::Cholesky(q, LDL = FALSE)
      },
      error = function(e) {
        if (not_pd(e)) NULL else stop(e)
      }
    ),
    # CHOLMOD warns before it fails on a matrix that is not positive definite.
    warning = function(w) {
      if (not_pd(w)) invokeRestart("muffleWarning")
    }
  )
  if (is.null(chol_q)) {
    return(NULL)
  }
  if (!same_pattern) {
    reuse <- list(factor = chol_q, p = q@p, i = q@i)
  }
  list(factor = chol_q, reuse = reuse)
}

# Joint posterior covariance of the fixed effects, computed exactly from the
# Gaussian approximations INLA stores for each hyperparameter configuration:
# the covariance of the mixture that inla.posterior.sample() draws from.
.inla_fixef_vcov <- function(object) {
  cs <- object$inla$misc$configs
  if (is.null(object$inla) || is.null(cs)) {
    stop(
      "vcov() needs an INLA fit built with ",
      "control.compute = list(config = TRUE).",
      call. = FALSE
    )
  }
  fixed_names <- rownames(object$inla$summary.fixed)
  p <- length(fixed_names)
  # The stored precision matrices omit the linear-predictor block that the
  # contents table lists first.
  offset <- sum(cs$contents$length) - cs$n
  idx <- cs$contents$start[match(fixed_names, cs$contents$tag)] - offset
  if (anyNA(idx)) {
    stop(
      "Could not locate the fixed effects in the INLA configurations.",
      call. = FALSE
    )
  }
  unit <- Matrix::sparseMatrix(
    i = idx,
    j = seq_len(p),
    x = 1,
    dims = c(cs$n, p)
  )
  w <- .fb_inla_config_weights(cs)

  second <- matrix(0, p, p)
  first <- numeric(p)
  reuse <- NULL
  for (k in seq_along(cs$config)) {
    g <- .fb_inla_config_gaussian(cs, k, reuse)
    reuse <- g$reuse
    cov_f <- as.matrix(Matrix::solve(g$chol_q, unit))[idx, , drop = FALSE]
    mu <- g$mu
    if (!is.null(g$gain)) {
      mu <- mu - drop(g$gain %*% (g$a_mat %*% mu - cs$constr$e))
      cov_f <- cov_f -
        g$gain[idx, , drop = FALSE] %*%
          t(g$sigma_at[idx, , drop = FALSE])
    }
    mu <- mu[idx]
    second <- second + w[k] * (cov_f + tcrossprod(mu))
    first <- first + w[k] * mu
  }
  v <- second - tcrossprod(first)
  v <- (v + t(v)) / 2

  # Back to the original scale, by the factors the posterior draws use.
  b <- .fb_inla_draw_map(object, fixed_names)$b
  v <- v * tcrossprod(b)

  labels <- .inla_restore_term_labels(object$level_labels, fixed_names)
  dimnames(v) <- list(labels, labels)
  aliased <- labels %in% object$extras$fixed_alias$dropped$column
  v[aliased, ] <- 0
  v[, aliased] <- 0
  v
}


#' Coefficients of a per-row INLA fit
#'
#' Posterior means of the fixed effects, read from the INLA fit's
#' `summary.fixed` slot (treatment-contrast basis). These are the
#' coefficients consumed by [emmeans::emmeans()] and
#' [marginaleffects::predictions()] via the flexyBayes support methods.
#'
#' An INLA fit carries no `$glm` slot, so the fixed vector is read here
#' and not inherited. Every other value of `what` is resolved by the
#' shared body [coef.flexybayes()] uses, so the two engines answer the
#' same question the same way.
#'
#' A factor with a non-syntactic level (for example a level containing a
#' space, `"low N"`) is legalised internally (`make.names()`) before the
#' INLA emit, so the fit itself never dies on it. This method restores
#' the user's own, unlegalised label in the names it returns, as do
#' [ranef()], `summary()` and `predict(classify = )` on the same fit.
#'
#' @param object A `flexybayes_inla` fit.
#' @param what Which part of the fit to return: `"fixed"` (the default),
#'   `"random"`, `"missing"` or `"all"`. See [coef.flexybayes()] for the
#'   shape of each.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return For the default `what = "fixed"`, a named numeric vector of
#'   fixed-effect posterior means. Otherwise as documented for
#'   [coef.flexybayes()].
#' @export
coef.flexybayes_inla <- function(
  object,
  what = c("fixed", "random", "missing", "all"),
  ...
) {
  sf <- object$inla$summary.fixed
  # `rownames(sf)` carries the legalised (make.names()) level names the
  # fit was built against when the model uses a non-syntactic factor
  # level. Restore the user's own labels, so coef() and ranef() (which
  # delegates to coef(what = "random")) print what the user wrote.
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

#' Joint posterior draws of the variance components of an INLA fit
#'
#' @description
#' Draws the hyperparameters of an INLA fit jointly and returns them on the
#' scale and under the names of the fit's variance-component table
#' (`summary(fit)$varcomp`): a standard deviation for the residual (`sigma`)
#' and for each random term (`sd_<term>`), and a correlation as it is.
#'
#' @details
#' The draws come from `INLA::inla.hyperpar.sample()` and are mapped to the
#' response's original scale. A hyperparameter held fixed, such as the
#' precision of a known variance, is not drawn. With `seed`, the draws from a
#' given fit are repeatable on one INLA thread
#' (`options(flexyBayes.inla_threads = "1:1")`). A refit of the same model in
#' another R session can reach a different mode (see [flexybayes()], section
#' *INLA threads*).
#'
#' @param object A `flexybayes_inla` fit.
#' @param n The number of draws.
#' @param seed An integer seed, or `NULL` (the default) to use the session's
#'   random number stream.
#'
#' @returns A data frame with `n` rows and one column per variance
#'   component.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(site = factor(rep(1:6, each = 5)))
#'   d$y <- rnorm(6)[d$site] + rnorm(30)
#'   fit <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
#'                     verbose = FALSE)
#'   head(fb_variance_draws(fit, 100L, seed = 1L))
#' }
#' }
#' @export
fb_variance_draws <- function(object, n, seed = NULL) {
  if (!inherits(object, "flexybayes_inla") || is.null(object$inla)) {
    stop("`object` must be a flexyBayes fit on INLA.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("`n` must be a single positive number of draws.", call. = FALSE)
  }
  draw <- function() .fb_inla_hyperpar_sample(object, as.integer(n))
  hyper <- if (is.null(seed)) draw() else .fb_with_seed(seed, draw())
  hyper <- as.matrix(hyper)
  map <- .mapper_inla(object, object$fb)$map
  out <- lapply(colnames(hyper), function(nm) {
    if (startsWith(nm, "Precision for ")) 1 / sqrt(hyper[, nm]) else hyper[, nm]
  })
  names(out) <- vapply(
    colnames(hyper),
    function(nm) unname(map[nm] %||% nm),
    character(1L)
  )
  as.data.frame(out, check.names = FALSE)
}

#' Posterior covariance of a per-row INLA fit's fixed effects
#'
#' The joint posterior covariance of the fixed effects, the matrix that
#' contrast and marginal-mean standard errors need. It is computed exactly
#' from the Gaussian approximations INLA stores for the fit, averaged over
#' INLA's hyperparameter configurations, so it involves no random draws:
#' every call on the same fit returns the same matrix and the session's
#' random number stream is left untouched. It is the covariance of the
#' draws `INLA::inla.posterior.sample()` would give, without their Monte
#' Carlo error. The standard deviations on its diagonal agree closely with
#' `summary.fixed$sd`, which comes from INLA's marginal densities.
#'
#' @param object A `flexybayes_inla` fit.
#' @param n_samples Not used, as the covariance is exact. Accepted so that
#'   calls written for earlier versions, which estimated the matrix from
#'   `n_samples` posterior draws, still run.
#' @param ... Ignored. Present for compatibility with the generic.
#' @return Posterior covariance matrix of the fixed effects, with the
#'   coefficient names of [coef()] as dimnames.
#' @export
vcov.flexybayes_inla <- function(object, n_samples = 2000L, ...) {
  .inla_fixef_vcov(object)
}

#' Credible intervals for the fixed effects of a per-row INLA fit
#'
#' Quantiles of INLA's own posterior marginals for the fixed effects, not
#' frequentist confidence intervals and not a normal approximation to
#' them. The bounds come from `INLA::inla.qmarginal()` applied to the
#' marginal densities the fit stores, so any credible level is available
#' and not only the 0.95 that `summary.fixed` tabulates.
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
  # `names(marg)` is the legalised name. Restore the user's own, so
  # confint() rownames match those of coef() and vcov() on the same fit.
  dimnames(bounds) <- list(
    .inla_restore_term_labels(object$level_labels, names(marg)),
    paste0(round(probs * 100, 1), "%")
  )
  aliased <- rownames(bounds) %in% object$extras$fixed_alias$dropped$column
  bounds[aliased, ] <- NA_real_

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
  # marginaleffects) must see only the fixed-effect basis. A leftover
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
  # A known variance per row is a row effect in INLA's linear predictor and
  # part of the likelihood on brms. The fitted value is the mean without it
  # (R/known_variance.R).
  known <- .fb_inla_known_variance_effect(object)
  if (!is.null(known)) {
    out <- out - known
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
    tryCatch(
      as.numeric(eval(f[[2L]], envir = object$data)),
      error = function(e) NULL
    )
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
#' `logLik()` generic denotes and that the brms backend exposes.
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
    "log-likelihood (model evidence, see summary()), not the conditional ",
    "log-likelihood logLik() denotes. Returning NA."
  )
  np <- tryCatch(nrow(object$inla$summary.fixed), error = function(e) {
    NA_integer_
  })
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
#' With `gen` random, `classify = "gen"` returns each genotype's fixed part
#' plus its random effect, and on a multi-environment fit with `gen:env`
#' random it also adds that genotype's `gen:env` effects averaged over
#' environments. This matches ASReml's default `predict()`: equal-weight
#' averaging over the fixed-term variables not in `classify`, a covariate
#' held at its mean, variables that appear only in random terms ignored, and
#' every random term whose variables are all in `classify` or averaged over
#' included. The intervals integrate over INLA's hyperparameter
#' configurations, so they include the uncertainty in the variance
#' components. A table with no random term is computed by emmeans, with
#' intervals from the Gaussian approximation of the fixed-effect
#' posterior. The arguments of
#' ASReml's `predict()` that change these rules (`present`, `average`,
#' `ignore`, `use`, `except`, `only`, `associate`, `levels`, `vcov`, `sed`,
#' `parallel`, `aliased`, `design.points`) are refused with
#' `classify_argument_not_supported`.
#'
#' @param object A `flexybayes_inla` fit.
#' @param newdata Optional data frame. Defaults to the fit data. Ignored,
#'   with a warning, when `classify` is supplied. On a model with an
#'   `offset()` the prediction includes the offset evaluated on `newdata`,
#'   so `newdata` needs every column the offset reads. An absent one
#'   refuses with `newdata_offset_missing` instead of being taken as zero.
#' @param type `"response"` or `"link"`.
#' @param se.fit Logical: also return delta-method standard errors from
#'   the fixed-effect covariance.
#' @param classify The factors to break a predicted-value table down by:
#'   a character value (`"Variety"`, `"Variety:env"`) or a one-sided
#'   formula (`~ Variety`). The factors may be fixed or random. `NULL`
#'   (the default) is the historical behaviour.
#' @param level Credible level for the classify table's interval, or for the
#'   predictive interval with `allow_new_levels = TRUE`, as a proportion.
#'   Default `0.95`.
#' @param allow_new_levels Logical. `FALSE` (the default) is the
#'   population-level prediction above. `TRUE` returns the posterior
#'   predictive of the rows of `newdata` (see Details), in which every random
#'   effect enters: a level the fit saw at its posterior draw, and a level it
#'   did not see drawn from its fitted variance component.
#' @param mean_by With `allow_new_levels = TRUE`, the columns of `newdata`
#'   whose combinations define the sets of rows to average: a character
#'   vector or a one-sided formula. `mean_by = ~ TRIAL` predicts each trial
#'   mean, `~ TRIAL + VARIETY + NRATE` each cell mean. `NULL` (the default)
#'   predicts each row.
#' @param n_draws With `allow_new_levels = TRUE`, the number of posterior
#'   draws. Default `1000`.
#' @param re_formula `NA` (the default): the prediction is population-level,
#'   the fixed part of the model, as `re_formula = NA` gives on a brms fit.
#'   Any other value refuses with `method_argument_not_supported`. [fitted()]
#'   returns the predictions at the fitted levels of the random effects.
#' @param summary `TRUE` (the default): the prediction is the posterior
#'   mean. `FALSE` refuses with `method_argument_not_supported`.
#'   [fb_as_draws_simple()] returns the posterior draws.
#' @param ... Ignored, except that an argument of ASReml's `predict()`
#'   given with `classify` is refused by name.
#'
#' @details
#' With `allow_new_levels = TRUE`, the prediction for a set of rows is the
#' mean response over them, observation noise included, as an observed trial
#' or cell mean at a new locality would be. Given one posterior draw it is
#' normal, with mean `x'beta` averaged over the rows plus the mean of the
#' effects at levels the fit saw, and variance the sum over random terms of
#' `sigma^2 * sum_l (n_l / n)^2` over the unseen levels `l` (`n_l` of the `n`
#' rows at level `l`), plus `sigma_field^2` times the mean AR1 x AR1
#' correlation over pairs of rows in a trial the fit did not see, plus the
#' residual variance over `n`. The prediction is the mixture of these normals
#' over the draws: `estimate` and `sd` are its exact mean and standard
#' deviation, and `lower` and `upper` its quantiles. The fixed effects and the
#' seen effects are drawn jointly from INLA's approximation, and the variances
#' from their continuous joint posterior (as [fb_variance_draws()] draws
#' them). Every random number comes from R's generator, so `set.seed()` before
#' the call makes the result repeatable.
#'
#' The new-level draw is built for random intercepts, interactions
#' (`TRIAL:REP`), `diag(f):g` and `id(trial):ar1(row):ar1(col)`, over an
#' independent residual or `dsum(~ units | f)`, for a Gaussian response on the
#' identity link. Any other structure is refused with
#' `predict_new_levels_structure_not_supported`, a new level with no fitted
#' variance (a new level of `f` in `diag(f):g`, a new residual section, a row
#' or column off the fitted lattice) with `predict_new_levels_no_variance`,
#' and another family with `predict_new_levels_family_not_supported`.
#'
#' On a fit with `fixed_alias = "drop"`, a row that is not a linear
#' combination of the rows the model was fitted to has no estimable
#' prediction: it is returned as `NA`, and the `status` attribute (a
#' `status` element with `se.fit = TRUE`, a `status` column with
#' `classify`) reads `"Aliased"` there and `"Estimable"` elsewhere, as in
#' asreml's `predict()`.
#'
#' @return With `classify`, a data frame of class `fb_predict_classify`.
#'   With `allow_new_levels = TRUE`, a data frame of class `fb_predict_new`:
#'   one row per set, holding the `mean_by` columns (or `row`), `n_rows`,
#'   `estimate`, `sd`, `lower` and `upper` (and `status` on a fit with
#'   `fixed_alias = "drop"`, where a set holding a row the fit does not
#'   estimate is `NA`), with the joint predictive draws, draws by sets, in
#'   `attr(, "draws")`. Otherwise a numeric vector of predictions, or a list
#'   `fit` / `se.fit` when `se.fit = TRUE`.
#' @export
predict.flexybayes_inla <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  se.fit = FALSE,
  classify = NULL,
  level = 0.95,
  allow_new_levels = FALSE,
  mean_by = NULL,
  n_draws = 1000L,
  re_formula = NA,
  summary = TRUE,
  ...
) {
  type <- match.arg(type)
  re_na <- is.atomic(re_formula) &&
    length(re_formula) == 1L &&
    is.na(re_formula)
  if (!re_na) {
    .fb_refuse_method_argument(
      "predict() on an INLA fit",
      "re_formula",
      paste0(
        "it returns population-level predictions, the fixed part of the ",
        "model, which is re_formula = NA. fitted() returns the predictions ",
        "at the fitted levels of the random effects."
      )
    )
  }
  if (!isTRUE(summary)) {
    .fb_refuse_method_argument(
      "predict() on an INLA fit",
      "summary",
      paste0(
        "it returns posterior means, which is summary = TRUE. ",
        "fb_as_draws_simple() returns the posterior draws."
      )
    )
  }
  .fb_predict_new_route(allow_new_levels, mean_by, classify, type, se.fit)
  if (isTRUE(allow_new_levels)) {
    return(.fb_predict_new_inla(object, newdata, mean_by, n_draws, level))
  }
  if (!is.null(classify)) {
    .fb_classify_newdata_note(newdata)
    return(.fb_predict_classify(object, classify, level, ...))
  }
  fam <- family(object)
  if (!identical(fam$link, "identity")) {
    stop(.fb_refusal_condition(
      reason_code = "predict_link_not_supported",
      message = paste0(
        "predict.flexybayes_inla() currently supports the identity ",
        "link only, and this fit uses link '",
        fam$link,
        "'. Work with ",
        "the posterior draws via fb_as_draws_simple() for non-identity ",
        "links."
      )
    ))
  }

  supplied <- !is.null(newdata)
  if (is.null(newdata)) {
    # .fb_fit_data(), not the raw object$data: the raw slot is the
    # level-legalised copy the fit was built against, while coef(object)
    # below carries the user's own names. The raw slot would trigger the
    # mismatch that .fb_fixef_model_matrix() checks for.
    newdata <- .fb_fit_data(object)
  }
  # The offset on the rows predicted at, read before the design matrix so an
  # absent offset column refuses by name, not inside model.frame().
  off <- .fb_offset_values(.fb_fit_ir(object), newdata, newdata = supplied)
  trms <- stats::delete.response(stats::terms(formula(object)))
  bhat <- coef(object)
  X <- .fb_fixef_model_matrix(trms, newdata, names(bhat), .fb_fit_data(object))
  eta <- as.numeric(X %*% bhat)
  if (!is.null(off)) {
    eta <- eta + off
  }
  status <- .fb_fixed_alias_status(object, X)
  if (!is.null(status)) {
    eta[status == "Aliased"] <- NA_real_
    attr(eta, "status") <- status
  }

  if (!se.fit) {
    return(eta)
  }
  V <- vcov(object)
  se <- sqrt(rowSums((X %*% V) * X))
  if (is.null(status)) {
    return(list(fit = eta, se.fit = se))
  }
  se[status == "Aliased"] <- NA_real_
  list(fit = eta, se.fit = se, status = status)
}


# .fb_with_seed(): evaluate under a named seed and put the caller's random
# stream back exactly as it was found.
.fb_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      suppressWarnings(rm(
        list = ".Random.seed",
        envir = globalenv()
      )),
      add = TRUE
    )
  }
  set.seed(seed)
  force(expr)
}
