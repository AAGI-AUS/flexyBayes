# marginaleffects integration for flexyBayes
#
# marginaleffects discovers a model through a small set of S3 generics:
# get_coef(), set_coef(), get_vcov(), get_predict() and get_data().
# We register all five for both backends via @exportS3Method (delayed
# registration -- marginaleffects stays in Suggests).
#
# get_predict() computes population-level (fixed-effect) predictions
# from the SAME design-matrix / coefficient basis that get_coef() and
# set_coef() expose, so the delta-method standard errors marginaleffects
# builds (perturb coef via set_coef -> re-predict) are internally
# consistent. Random effects are held at their population mean; the
# identity link is assumed (the Gaussian fixed-effect surface).

# Coefficient override stash used by set_coef() -> get_predict() so the
# marginaleffects delta-method Jacobian sees perturbed coefficients.
.fb_coef_override_attr <- "fb_coef_override"

.fb_get_coef <- function(model, ...) coef(model)

.fb_set_coef <- function(model, coefs, ...) {
  attr(model, .fb_coef_override_attr) <- coefs
  model
}

.fb_get_vcov <- function(model, ...) {
  v <- vcov(model)
  if (is.null(v)) {
    stop(
      "marginaleffects support: the fit has no fixed-effect ",
      "covariance.",
      call. = FALSE
    )
  }
  v
}

.fb_get_data <- function(model, ...) .fb_fit_data(model)

.fb_get_predict <- function(model, newdata = NULL, type = "response", ...) {
  fam <- tryCatch(family(model), error = function(e) list(link = "identity"))
  if (!identical(fam$link %||% "identity", "identity")) {
    stop(
      "marginaleffects support currently covers the identity link ",
      "only; this fit uses link '",
      fam$link,
      "'. Work with the ",
      "posterior draws via fb_as_draws_simple() instead.",
      call. = FALSE
    )
  }

  if (is.null(newdata)) {
    newdata <- .fb_fit_data(model)
  }
  bhat <- attr(model, .fb_coef_override_attr) %||% coef(model)
  trms <- .fb_fixef_terms(model)
  data <- .fb_fit_data(model)
  X <- .fb_fixef_model_matrix(trms, newdata, names(bhat), data)
  eta <- as.numeric(X %*% bhat[colnames(X)])
  data.frame(rowid = seq_len(nrow(newdata)), estimate = eta)
}


#' marginaleffects support: fixed-effect coefficients (greta backend)
#' @param model A `flexybayes` fit.
#' @param ... Ignored.
#' @return Named numeric vector of coefficients.
#' @exportS3Method marginaleffects::get_coef
get_coef.flexybayes <- function(model, ...) .fb_get_coef(model, ...)

#' @rdname get_coef.flexybayes
#' @param model A `flexybayes_inla` fit.
#' @exportS3Method marginaleffects::get_coef
get_coef.flexybayes_inla <- function(model, ...) .fb_get_coef(model, ...)

#' marginaleffects support: set coefficients (greta backend)
#' @param model A `flexybayes` fit.
#' @param coefs Replacement coefficient vector.
#' @param ... Ignored.
#' @return The fit with a coefficient override attached.
#' @exportS3Method marginaleffects::set_coef
set_coef.flexybayes <- function(model, coefs, ...) {
  .fb_set_coef(model, coefs, ...)
}

#' @rdname set_coef.flexybayes
#' @param model A `flexybayes_inla` fit.
#' @exportS3Method marginaleffects::set_coef
set_coef.flexybayes_inla <- function(model, coefs, ...) {
  .fb_set_coef(model, coefs, ...)
}

#' marginaleffects support: covariance (greta backend)
#' @param model A `flexybayes` fit.
#' @param ... Ignored.
#' @return Fixed-effect covariance matrix.
#' @exportS3Method marginaleffects::get_vcov
get_vcov.flexybayes <- function(model, ...) .fb_get_vcov(model, ...)

#' @rdname get_vcov.flexybayes
#' @param model A `flexybayes_inla` fit.
#' @exportS3Method marginaleffects::get_vcov
get_vcov.flexybayes_inla <- function(model, ...) .fb_get_vcov(model, ...)

#' Model data accessor (greta backend)
#'
#' Registered for [insight::get_data()] so marginaleffects (which
#' discovers a model's data through insight) can build reference grids
#' and average predictions without an explicit `newdata`.
#'
#' @param x A `flexybayes` fit.
#' @param ... Ignored.
#' @return The model data frame.
#' @exportS3Method insight::get_data
get_data.flexybayes <- function(x, ...) .fb_get_data(x, ...)

#' @rdname get_data.flexybayes
#' @param x A `flexybayes_inla` fit.
#' @exportS3Method insight::get_data
get_data.flexybayes_inla <- function(x, ...) .fb_get_data(x, ...)

#' marginaleffects support: population-level predictions (greta backend)
#' @param model A `flexybayes` fit.
#' @param newdata Data frame to predict on (default: fit data).
#' @param type Prediction scale (identity link only).
#' @param ... Ignored.
#' @return A data frame with `rowid` and `estimate`.
#' @exportS3Method marginaleffects::get_predict
get_predict.flexybayes <- function(
  model,
  newdata = NULL,
  type = "response",
  ...
) {
  .fb_get_predict(model, newdata = newdata, type = type, ...)
}

#' @rdname get_predict.flexybayes
#' @param model A `flexybayes_inla` fit.
#' @exportS3Method marginaleffects::get_predict
get_predict.flexybayes_inla <- function(
  model,
  newdata = NULL,
  type = "response",
  ...
) {
  .fb_get_predict(model, newdata = newdata, type = type, ...)
}


# Called from .onLoad -- marginaleffects >= 0.18 auto-discovers via the
# registered S3 methods, so no explicit registration is needed.
.register_marginaleffects <- function() {
  requireNamespace("marginaleffects", quietly = TRUE)
}
