# tidiers.R -- broom-style tidy / glance / augment methods for the
# flexyBayes fit classes.
#
# A flexyBayes fit is a backend-specific S3 list: the greta route returns
# `flexybayes`, the brms route `c("flexybayes_brms", "flexybayes", ...)`,
# and the INLA route `c("flexybayes_inla", "list")`. Downstream code that
# wanted a flat one-row-per-term summary had to reach into backend-specific
# slots, which is fragile and -- worse -- different across backends, so a
# triangulation table that compared greta against INLA was hand-built every
# time. These methods register against the `tidy()` generic (re-exported by
# `broom`), so `broom::tidy(fit)` and `generics::tidy(fit)` both return a
# stable, documented `data.frame` with the canonical `broom` column names
# (`term`, `estimate`, `std.error`, `conf.low`, `conf.high`) regardless of
# which backend produced the fit. Cross-engine tables stop being hand-built.

# ---- Re-export the tidy generic ----------------------------------------
# Following the kernR convention: re-export `generics::tidy` so a user who
# has only flexyBayes loaded can call `tidy(fit)` without attaching broom.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment

# ---- tidy() ------------------------------------------------------------

#' Tidy a flexyBayes fit into a one-row-per-term data frame
#'
#' Turns a flexyBayes fit into a flat, `broom`-style `data.frame`: one row
#' per model term, with stable, documented columns. The method is registered
#' against the `tidy()` generic (re-exported by `broom`), so
#' `broom::tidy(fit)` and `generics::tidy(fit)` both dispatch here.
#'
#' This is the supported accessor for cross-engine summaries. The hub returns
#' backend-specific objects -- `flexybayes` (greta), `flexybayes_brms`
#' (brms), `flexybayes_inla` (INLA) -- whose internal layouts differ. Tidying
#' through this generic yields the same columns across all three, so a
#' greta-versus-INLA triangulation table can be assembled by `rbind`-ing two
#' `tidy()` outputs rather than reaching into each backend's slots by hand.
#'
#' The credible intervals are posterior quantile-based intervals, not
#' frequentist confidence intervals; they are reported in the
#' `broom`-canonical `conf.low` / `conf.high` columns. The `std.error` column
#' carries the posterior standard deviation of each term, again under the
#' `broom`-canonical (dotless to the user, dotted in the column name) label.
#'
#' @param x A flexyBayes fit: `flexybayes` (greta backend) or
#'   `flexybayes_brms` (brms backend, which inherits this method).
#' @param conf.int Logical. Whether to attach credible intervals. Defaults to
#'   `TRUE`.
#' @param conf.level Numeric in `(0, 1)`. The credible level for the
#'   intervals. Defaults to `0.95`.
#' @param effects Character. Which effects to return: `"fixed"` for the
#'   population-level (fixed) coefficients or `"random"` for the
#'   variance-component summary. Defaults to `"fixed"`.
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns A `data.frame` with one row per term and columns:
#'   \describe{
#'     \item{term}{Character. The coefficient or variance-component name.}
#'     \item{estimate}{Numeric. The posterior mean.}
#'     \item{std.error}{Numeric. The posterior standard deviation.}
#'     \item{conf.low}{Numeric. The lower credible bound (present when
#'       `conf.int = TRUE`).}
#'     \item{conf.high}{Numeric. The upper credible bound (present when
#'       `conf.int = TRUE`).}
#'   }
#'   An empty `data.frame` is returned when the requested effects are absent
#'   (for example `effects = "random"` on a fixed-effects-only fit).
#'
#' @seealso [glance.flexybayes()], [augment.flexybayes()],
#'   [tidy.flexybayes_inla()]
#' @examplesIf requireNamespace("generics", quietly = TRUE)
#' \dontrun{
#' fit <- flexybayes(yield ~ env, data = dat, backend = "greta")
#' tidy(fit)
#' tidy(fit, effects = "random")
#' }
#' @export
tidy.flexybayes <- function(
  x,
  conf.int = TRUE,
  conf.level = 0.95,
  effects = c("fixed", "random"),
  ...
) {
  effects <- match.arg(effects)

  if (effects == "fixed") {
    return(.tidy_fixed_flexybayes(x, conf.int, conf.level))
  }
  .tidy_random_flexybayes(x, conf.int)
}

# Fixed-effect tidier shared by the greta / brms fit classes. Reads the
# posterior mean from coef(), the posterior SD from the diagonal of vcov(),
# and -- when asked -- the credible interval from confint().
.tidy_fixed_flexybayes <- function(x, conf.int, conf.level) {
  beta <- stats::coef(x)
  if (length(beta) == 0L) {
    return(.empty_tidy(conf.int))
  }

  v <- stats::vcov(x)
  se <- if (!is.null(v)) sqrt(diag(v)) else rep(NA_real_, length(beta))

  out <- data.frame(
    term = names(beta),
    estimate = unname(beta),
    std.error = unname(se),
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    ci <- stats::confint(x, level = conf.level)
    out$conf.low <- unname(ci[, 1L])
    out$conf.high <- unname(ci[, 2L])
  }

  out
}

# Variance-component tidier shared by the greta / brms fit classes. The
# components live on `fit$extras$variance_comps`, already summarised to
# component / estimate / sd / q2.5 / q97.5 by the backend post-fit code.
.tidy_random_flexybayes <- function(x, conf.int) {
  vc <- x$extras$variance_comps
  if (is.null(vc) || nrow(vc) == 0L) {
    return(.empty_tidy(conf.int))
  }

  out <- data.frame(
    term = vc$component,
    estimate = vc$estimate,
    std.error = vc$sd,
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    out$conf.low <- vc$q2.5
    out$conf.high <- vc$q97.5
  }

  out
}

#' Tidy a per-row INLA fit into a one-row-per-term data frame
#'
#' The INLA backend returns a `flexybayes_inla` object that does not inherit
#' from `flexybayes`, so it needs its own `tidy()` method. The fixed-effect
#' summary is read directly off INLA's `summary.fixed` table, whose `mean`,
#' `sd`, and `0.025quant` / `0.975quant` columns map cleanly onto the
#' `broom`-canonical `estimate`, `std.error`, `conf.low`, and `conf.high`.
#'
#' Because the INLA fixed-effect intervals come from the marginal posteriors
#' INLA has already integrated, the `conf.level` argument is accepted for
#' generic compatibility but only the 95% bounds INLA reports are returned;
#' a one-off message notes this when a different level is requested rather
#' than silently ignoring it.
#'
#' @param x A `flexybayes_inla` fit.
#' @param conf.int Logical. Whether to attach the credible-interval columns.
#'   Defaults to `TRUE`.
#' @param conf.level Numeric in `(0, 1)`. Accepted for generic compatibility;
#'   INLA reports the 95% marginal bounds, so a non-0.95 request is noted.
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns A `data.frame` with one row per fixed-effect term and the columns
#'   `term`, `estimate`, `std.error`, and (when `conf.int = TRUE`)
#'   `conf.low` / `conf.high`.
#'
#' @seealso [tidy.flexybayes()]
#' @examplesIf requireNamespace("INLA", quietly = TRUE)
#' \dontrun{
#' fit <- flexybayes(yield ~ env, data = dat, backend = "inla")
#' tidy(fit)
#' }
#' @export
tidy.flexybayes_inla <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  sf <- x$inla$summary.fixed
  if (is.null(sf) || nrow(sf) == 0L) {
    return(.empty_tidy(conf.int))
  }

  if (isTRUE(conf.int) && !isTRUE(all.equal(conf.level, 0.95))) {
    message(
      "tidy.flexybayes_inla(): INLA reports 95% marginal credible bounds; ",
      "conf.level = ", conf.level, " is ignored. Resample via vcov() for a ",
      "different level."
    )
  }

  out <- data.frame(
    term = rownames(sf),
    estimate = sf[["mean"]],
    std.error = sf[["sd"]],
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    out$conf.low <- sf[["0.025quant"]]
    out$conf.high <- sf[["0.975quant"]]
  }

  rownames(out) <- NULL
  out
}

# Empty-but-well-typed tidy frame, so a fixed-effects-only fit asked for
# random effects (or a degenerate fit) returns the documented columns with
# zero rows rather than an untyped empty data.frame.
.empty_tidy <- function(conf.int) {
  out <- data.frame(
    term = character(0),
    estimate = numeric(0),
    std.error = numeric(0),
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    out$conf.low <- numeric(0)
    out$conf.high <- numeric(0)
  }
  out
}

# ---- glance() ----------------------------------------------------------

#' Glance at a flexyBayes fit
#'
#' Returns a one-row `data.frame` of model-level statistics: the response
#' family and link, the observation and parameter counts, the
#' log-likelihood, the worst convergence diagnostics across monitored
#' parameters, and the wall-clock run time.
#'
#' @param x A flexyBayes fit (`flexybayes` or `flexybayes_brms`).
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns A one-row `data.frame` with `nobs`, `npar`, `logLik`, `family`,
#'   `link`, `chains`, `samples`, `max_rhat`, `min_ess`, and `run_time`.
#'
#' @seealso [tidy.flexybayes()]
#' @examplesIf requireNamespace("generics", quietly = TRUE)
#' \dontrun{
#' glance(fit)
#' }
#' @export
glance.flexybayes <- function(x, ...) {
  mi <- x$extras$model_info
  ll <- stats::logLik(x)

  conv <- x$extras$convergence
  max_rhat <- if (!is.null(conv$gelman)) {
    max(conv$gelman$psrf[, "Point est."], na.rm = TRUE)
  } else {
    NA_real_
  }

  min_ess <- if (!is.null(conv$n_eff)) {
    min(conv$n_eff, na.rm = TRUE)
  } else {
    NA_real_
  }

  data.frame(
    nobs = mi$n_obs,
    npar = mi$n_params,
    logLik = as.numeric(ll),
    family = mi$family,
    link = mi$link,
    chains = x$extras$call_info$chains,
    samples = x$extras$call_info$n_samples,
    max_rhat = round(max_rhat, 3L),
    min_ess = round(min_ess, 0L),
    run_time = round(x$extras$run_time, 1L),
    stringsAsFactors = FALSE
  )
}

# ---- augment() ---------------------------------------------------------

#' Augment a flexyBayes fit with fitted values and residuals
#'
#' Returns the model frame with two observation-level columns added: the
#' posterior-mean fitted value and the response residual.
#'
#' @param x A flexyBayes fit (`flexybayes` or `flexybayes_brms`).
#' @param data Optional `data.frame` to augment. Defaults to the data the
#'   model was fitted to.
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns The supplied (or original) `data.frame` with `.fitted` and
#'   `.resid` columns appended.
#'
#' @seealso [tidy.flexybayes()]
#' @examplesIf requireNamespace("generics", quietly = TRUE)
#' \dontrun{
#' augment(fit)
#' }
#' @export
augment.flexybayes <- function(x, data = NULL, ...) {
  if (is.null(data)) {
    data <- x$glm$data
  }

  out <- data
  out$.fitted <- stats::fitted(x)
  out$.resid <- stats::residuals(x)
  out
}

# ---- INLA: glance() / augment() not available (tidy() only) ------------
# INLA fits do not inherit "flexybayes", so without these the generic
# raises a bare "no applicable method" error. These give an on-brand,
# actionable refusal instead and make methods("glance"/"augment") list INLA.

#' @rdname glance.flexybayes
#' @export
glance.flexybayes_inla <- function(x, ...) {
  stop(
    "glance() is not available for INLA fits (`flexybayes_inla`). ",
    "Use tidy() for a coefficient-level summary, or summary() / ",
    "fb_structured_cov() for variance components.",
    call. = FALSE
  )
}

#' @rdname augment.flexybayes
#' @export
augment.flexybayes_inla <- function(x, data = NULL, ...) {
  stop(
    "augment() is not available for INLA fits (`flexybayes_inla`). ",
    "Use tidy() for a coefficient-level summary, or predict() for ",
    "fitted values.",
    call. = FALSE
  )
}
