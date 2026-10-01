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
#' reach an active engine are `flexybayes_brms` and `flexybayes_inla`.
#' User-defined methods can extend the generic.
#'
#' @param fit A model fit object carrying a posterior, dispatched on by
#'   the methods listed above.
#' @param ... Method-specific arguments: `n_samples` for the INLA method.
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

  n_samples <- as.integer(n_samples)
  # Force single-threaded sampling, matching the summary path in
  # methods_inla.R. inla.posterior.sample() otherwise spawns one process
  # per core, which fails under the two-core limit `R CMD check --as-cran`
  # enforces (the worker count is reported as "N simultaneous processes
  # spawned").
  samples <- tryCatch(
    .fb_inla_posterior_sample(fit, n_samples, num.threads = "1:1"),
    error = function(e) {
      stop(
        "INLA::inla.posterior.sample() failed: ",
        conditionMessage(e),
        ". Re-fit with control.compute = list(config = TRUE).",
        call. = FALSE
      )
    }
  )

  if (length(samples) == 0L) {
    stop("INLA::inla.posterior.sample() returned an empty list.", call. = FALSE)
  }

  latent_names <- rownames(samples[[1]]$latent)
  hyperpar_names <- names(samples[[1]]$hyperpar)

  out <- list()
  for (pn in latent_names) {
    out[[pn]] <- vapply(
      samples,
      function(s) as.numeric(s$latent[pn, 1]),
      numeric(1)
    )
  }
  for (pn in hyperpar_names) {
    out[[pn]] <- vapply(
      samples,
      function(s) as.numeric(s$hyperpar[pn]),
      numeric(1)
    )
  }
  out
}

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.flexybayes_brms <- function(fit, ...) {
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
fb_as_draws_simple.default <- function(fit, ...) {
  stop(
    "fb_as_draws_simple() does not know how to ",
    "extract draws from an object of class ",
    paste(class(fit), collapse = "/"),
    ". Define an `fb_as_draws_simple.<class>` method.",
    call. = FALSE
  )
}
