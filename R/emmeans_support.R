# emmeans integration for flexyBayes
#
# Provides recover_data() and emm_basis() methods so that
# emmeans::emmeans(fit, ~ factor) works on greta- and INLA-backed fits.
# Both methods are registered for the foreign emmeans generics via
# @exportS3Method (delayed S3 registration), so emmeans stays in
# Suggests.
#
# The fixed-effect coefficients and covariance come from the fit's own
# coef() / vcov() methods; the design matrix is reconciled to the
# coefficient basis by .fb_fixef_model_matrix() (see ecosystem_support.R).
# Posterior summaries are reported on z-intervals (df = Inf), the
# Bayesian convention emmeans uses for sampled posteriors.

# Shared recover_data core for either backend.
.fb_recover_data <- function(object, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for this method.", call. = FALSE)
  }
  trms <- .fb_fixef_terms(object)
  data <- .fb_fit_data(object)
  if (is.null(data)) {
    stop(
      "emmeans support: the fit does not carry its model data.",
      call. = FALSE
    )
  }
  # A synthetic model call lets emmeans::recover_data.call() extract the
  # predictor variables from `data` using the fixed-effect terms. `...`
  # is deliberately not forwarded: emmeans threads its own `data` /
  # `params` arguments through it, which would clash with the explicit
  # `data` here ("formal argument 'data' matched by multiple actual
  # arguments").
  fcall <- call("lm", formula = stats::formula(object), data = quote(data))
  emmeans::recover_data(fcall, trms, na.action = NULL, data = data)
}

# Shared emm_basis core. Builds X to match the fitted coefficient names,
# carries the posterior covariance as V, and supplies a non-estimability
# basis so emmeans rejects non-estimable combinations on the over-
# parameterised (greta) basis. df = Inf -> z-based intervals.
.fb_emm_basis <- function(object, trms, xlev, grid, ...) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for this method.", call. = FALSE)
  }
  bhat <- coef(object)
  V <- vcov(object)
  data <- .fb_fit_data(object)
  X <- .fb_fixef_model_matrix(trms, grid, names(bhat), data)
  bhat <- bhat[colnames(X)]
  V <- V[colnames(X), colnames(X), drop = FALSE]
  list(
    X = X,
    bhat = bhat,
    nbasis = estimability::nonest.basis(X),
    V = V,
    dffun = function(k, dfargs) Inf,
    dfargs = list()
  )
}


#' emmeans support: recover model data (greta backend)
#'
#' @param object A `flexybayes` fit.
#' @param ... Passed to [emmeans::recover_data()].
#' @return A data frame of predictors for the reference grid.
#' @exportS3Method emmeans::recover_data
recover_data.flexybayes <- function(object, ...) {
  .fb_recover_data(object, ...)
}

#' @rdname recover_data.flexybayes
#' @param object A `flexybayes_inla` fit.
#' @exportS3Method emmeans::recover_data
recover_data.flexybayes_inla <- function(object, ...) {
  .fb_recover_data(object, ...)
}

#' emmeans support: estimation basis (greta backend)
#'
#' @param object A `flexybayes` fit.
#' @param trms Fixed-effect terms supplied by emmeans.
#' @param xlev Factor levels supplied by emmeans.
#' @param grid Reference grid supplied by emmeans.
#' @param ... Ignored.
#' @return A list with `X`, `bhat`, `nbasis`, `V`, `dffun`, `dfargs`.
#' @exportS3Method emmeans::emm_basis
emm_basis.flexybayes <- function(object, trms, xlev, grid, ...) {
  .fb_emm_basis(object, trms, xlev, grid, ...)
}

#' @rdname emm_basis.flexybayes
#' @param object A `flexybayes_inla` fit.
#' @exportS3Method emmeans::emm_basis
emm_basis.flexybayes_inla <- function(object, trms, xlev, grid, ...) {
  .fb_emm_basis(object, trms, xlev, grid, ...)
}
