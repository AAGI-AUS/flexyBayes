# fb_as_draws_simple -- per-fit posterior extraction as a plain named list.
#
# The generic returns a named list whose values are numeric vectors of
# posterior draws, keyed by the engine's own parameter names, with methods
# for `flexybayes_brms` (the brms path) and `flexybayes_inla` (the INLA
# path). It needs no posterior-package dependency. `.fb_canonical_draws()`
# (R/fb_posterior.R) renames these draws to the canonical tokens the
# posterior generics return, and `.fb_genomic_draws()`
# (R/fb_genomic_extras.R) reads them for genomic_summary().

#' Extract per-parameter posterior draws from a model fit
#'
#' S3 generic that extracts a named list of numeric posterior-draw vectors
#' from a fit, keyed by the engine's own parameter names. The methods that
#' reach an active engine are `flexybayes_brms` and `flexybayes_inla`. The
#' `brmsfit` method reads a brms fit inside a flexyBayes fit, a two-part
#' non-linear fit from [fb_refit()] under the names of the linear model.
#' User-defined methods can extend the generic.
#'
#' On an INLA fit the draws come from the fitted approximation: each draw
#' takes one of INLA's hyperparameter configurations, chosen with its
#' posterior weight, and draws the latent field from that configuration's
#' Gaussian approximation, as `INLA::inla.posterior.sample()` does. Every
#' random number comes from R's generator, so `set.seed()` before the call
#' makes the draws repeatable.
#'
#' @param fit A model fit object carrying a posterior, dispatched on by
#'   the methods listed above.
#' @param ... Method-specific arguments: `n_samples` for the INLA method.
#'   The brms method returns every stored draw and refuses `n_samples` with
#'   `method_argument_not_supported`.
#' @returns A named list of numeric vectors, one element per parameter,
#'   each element holding that parameter's posterior draws.
#' @examples
#' # A named list of draw vectors, with no posterior dependency.
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#'   fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
#'                     verbose = FALSE)
#'   dr <- fb_as_draws_simple(fit)
#'   names(dr)
#' }
#' }
#' @export
fb_as_draws_simple <- function(fit, ...) UseMethod("fb_as_draws_simple")

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.flexybayes_inla <- function(fit, n_samples = 1000L, ...) {
  .check_installed(
    "INLA",
    "Package 'INLA' is required to extract draws from a ",
    "flexybayes_inla object."
  )
  if (is.null(fit$inla)) {
    stop("Cannot extract draws: fit$inla is missing.", call. = FALSE)
  }

  draws <- .fb_inla_config_sample(fit, as.integer(n_samples))
  out <- list()
  for (pn in rownames(draws$latent)) {
    out[[pn]] <- draws$latent[pn, ]
  }
  for (pn in colnames(draws$hyperpar)) {
    out[[pn]] <- draws$hyperpar[, pn]
  }
  out
}

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.flexybayes_brms <- function(fit, ...) {
  if ("n_samples" %in% names(list(...))) {
    .fb_refuse_method_argument(
      "fb_as_draws_simple() on a brms fit",
      "n_samples",
      paste0(
        "it returns every draw the sampler stored; the fit's `n_samples` ",
        "and `chains` set how many there are."
      )
    )
  }
  .check_installed(
    "posterior",
    "Package 'posterior' is required to extract draws from a ",
    "flexybayes_brms object."
  )
  if (is.null(fit$brms)) {
    stop("Cannot extract draws: fit$brms is missing.", call. = FALSE)
  }
  m <- tryCatch(
    .fb_brms_draws_matrix(fit$brms),
    error = function(e) {
      stop(
        "posterior::as_draws_matrix() failed on the ",
        "brms fit: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  cols <- colnames(m)
  if (is.null(cols)) {
    cols <- paste0("V", seq_len(ncol(m)))
  }
  setNames(lapply(seq_len(ncol(m)), function(j) as.numeric(m[, j])), cols)
}

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.brmsfit <- function(fit, ...) {
  .check_installed(
    "posterior",
    "Package 'posterior' is required to extract draws from a brmsfit."
  )
  m <- .fb_brms_draws_matrix(fit)
  setNames(lapply(seq_len(ncol(m)), function(j) as.numeric(m[, j])), colnames(m))
}

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.default <- function(fit, ...) {
  stop(
    "fb_as_draws_simple() does not know how to ",
    "extract draws from an object of class ",
    paste(class(fit), collapse = "/"),
    ". Define an `fb_as_draws_simple.<class>` method.",
    call. = FALSE
  )
}
