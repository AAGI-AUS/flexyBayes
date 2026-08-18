# as_draws_support.R -- the posterior-package draws view over any fit.
#
# Registers `as_draws()`, `as_draws_df()` and `as_draws_matrix()` for
# flexyBayes fits on the posterior package's own generics via
# @exportS3Method (delayed S3 registration), so posterior stays in
# Suggests. All three read one seam: `.fb_canonical_draws()`, which
# composes the backend draws with the per-backend canonical-name map and
# value transforms. A fit on either engine therefore yields the same
# parameter tokens -- (Intercept), the fixed-effect terms, sigma,
# sd_<group> -- with the variance components on the standard-deviation
# scale, which is the same view triangulate() compares.
#
# as_draws_df() is the constructor; the other two derive from its result,
# so there is one place where a draws object is built.

# --- the shared construction ------------------------------------- #

# .fb_check_posterior() --- the Suggests guard, refusal-first.
#
# posterior is a suggested package: a fit runs, prints and predicts
# without it. The guard names the package and the install command rather
# than letting a bare requireNamespace() failure surface from inside the
# draws extraction.
#
# @noRd
# @keywords internal
.fb_check_posterior <- function(what) {
  if (.fb_posterior_available()) {
    return(invisible(TRUE))
  }
  stop(
    what,
    " builds its result with the posterior package, which is not ",
    "installed. Install it with install.packages(\"posterior\"), or ",
    "read the draws as a named list of numeric vectors with ",
    "fb_as_draws_simple(fit), which needs no extra package.",
    call. = FALSE
  )
}

# .fb_check_n_draws() --- one positive count.
#
# @noRd
# @keywords internal
.fb_check_n_draws <- function(n_draws) {
  ok <- (is.numeric(n_draws) &&
    length(n_draws) == 1L &&
    is.finite(n_draws) &&
    n_draws >= 1L)
  if (!ok) {
    stop(
      "`n_draws` must be a single positive number of draws; got a ",
      class(n_draws)[[1L]],
      " of length ",
      length(n_draws),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# .fb_door_draws_df() --- the canonical draws as a posterior draws_df.
#
# The one construction the three methods share. `n_draws` reaches the
# INLA sampler, which draws that many samples from the fitted
# approximation; the brms path already carries every draw the sampler
# kept and ignores it.
#
# @noRd
# @keywords internal
.fb_door_draws_df <- function(x, n_draws) {
  .fb_check_posterior("as_draws()")
  .fb_check_n_draws(n_draws)

  draws <- .fb_canonical_draws(x, n_samples = as.integer(n_draws))

  if (length(draws) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "fit_lacks_posterior_draws",
      message = paste0(
        "as_draws() cannot build a draws object for this fit: no ",
        "canonical parameter could be read from it. The fit's class is <",
        paste(class(x), collapse = ", "),
        ">. A fit on an active engine carries its posterior under `$inla` ",
        "or `$brms`; check the fit completed (return_code = FALSE)."
      )
    ))
  }

  lens <- vapply(draws, length, integer(1L))
  if (length(unique(lens)) != 1L) {
    stop(
      "as_draws(): the backend returned draw vectors of unequal length (",
      paste(paste0(names(lens), " = ", lens), collapse = ", "),
      "). A draws object needs one value per parameter per draw.",
      call. = FALSE
    )
  }

  posterior::as_draws_df(draws)
}


# --- the three methods -------------------------------------------- #

#' Posterior draws from a flexyBayes fit
#'
#' Methods for the \pkg{posterior} package's `as_draws()`,
#' `as_draws_df()` and `as_draws_matrix()` generics, so a fitted model
#' hands its posterior to the Bayesian workflow ecosystem --
#' \pkg{posterior} summaries, \pkg{bayesplot} displays, anything reading
#' a `draws` object -- in the shape those tools expect.
#'
#' The parameter names are canonical and the same on every engine:
#' `(Intercept)`, one entry per fixed-effect term, `sigma` for the
#' residual standard deviation, and `sd_<group>` for the
#' standard deviation of each random-effect group. Variance components
#' are on the **standard-deviation scale** whatever the engine stored:
#' INLA's precision hyperparameters are transformed on the way out, so
#' `sd_g` is a draw of a standard deviation and not of a precision. This
#' is the same view [triangulate()] compares, and it is built by the same
#' internal seam, so a name in one is a name in the other.
#'
#' @section Where the draws come from:
#'
#' On a **brms** fit the draws are the sampler's own: the kept
#' post-warmup iterations, renamed. `n_draws` is ignored, because the
#' number of draws was fixed when the model was sampled.
#'
#' On an **INLA** fit there is no sampler. The fit is a nested Laplace
#' approximation, and the draws are sampled *from that fitted
#' approximation* by `INLA::inla.posterior.sample()`, reading the
#' configuration store every flexyBayes INLA fit is built with
#' (`control.compute = list(config = TRUE)`). Two consequences follow,
#' and they are the ones Tutorials 01 and 10 teach:
#'
#' \itemize{
#'   \item The **fitting** is deterministic. The Laplace approximation
#'     draws no random numbers, so a `seed` passed to [flexybayes()]
#'     changes nothing about the fit itself.
#'   \item The **sampling** step is not, and it is not reproducible from
#'     `set.seed()` alone. The hyperparameter draws follow R's random
#'     stream, but the latent field -- the intercept and the
#'     fixed-effect terms -- is drawn by INLA's own generator, which
#'     `INLA::inla.posterior.sample()` seeds at random unless its own
#'     `seed` argument is given, and flexyBayes leaves that argument at
#'     its default. `?INLA::inla.posterior.sample` states that
#'     reproducing a sample needs that seed *and* R's RNG state fixed.
#'     Treat the result as one sample of the posterior rather than a
#'     fixed object: Monte-Carlo error in anything computed from it
#'     shrinks with `n_draws`, so raise `n_draws` before reading a small
#'     difference between two sets of draws as a difference between two
#'     posteriors.
#' }
#'
#' @param x A fitted `flexybayes` object from either active engine.
#' @param n_draws Number of posterior draws to return, matching the
#'   default of [fb_as_draws_simple()]. Used on the INLA path, where the
#'   draws are sampled from the fitted approximation; ignored on the brms
#'   path, which returns the draws the sampler kept.
#' @param ... Passed to the underlying draws extraction.
#' @returns A \pkg{posterior} draws object holding one column per
#'   canonical parameter: a `draws_df` from `as_draws_df()` and from
#'   `as_draws()`, a `draws_matrix` from `as_draws_matrix()`.
#'
#' @seealso [fb_as_draws_simple()] for the same draws as a plain named
#'   list, with no \pkg{posterior} dependency; [canonical_names()] for
#'   the name map itself; [prior_summary()] for the priors those draws
#'   were taken under; [summary.flexybayes()], whose `$converge` slot
#'   carries the engine's own convergence diagnostics; [triangulate()]
#'   for a two-engine comparison over the same canonical view.
#'
#' @exportS3Method posterior::as_draws_df
as_draws_df.flexybayes <- function(x, n_draws = 1000L, ...) {
  .fb_door_draws_df(x, n_draws = n_draws)
}

#' @rdname as_draws_df.flexybayes
#' @exportS3Method posterior::as_draws
as_draws.flexybayes <- function(x, n_draws = 1000L, ...) {
  posterior::as_draws(.fb_door_draws_df(x, n_draws = n_draws))
}

#' @rdname as_draws_df.flexybayes
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.flexybayes <- function(x, n_draws = 1000L, ...) {
  posterior::as_draws_matrix(.fb_door_draws_df(x, n_draws = n_draws))
}
