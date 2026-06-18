# fb_extreme_value.R -- generalised extreme value (GEV) family + fitter.
#
# The generalised extreme value distribution is the limiting law of block
# maxima (annual maximum rainfall, daily maximum temperature, peak yields).
# It does not fit the single-column GLM-link mixed-model emit path the other
# flexyBayes families use: there is no canonical mean-link parameterisation,
# and the quantity of interest is the (location, scale, shape) triple plus
# the return levels they imply, not a regression coefficient. So GEV is
# provided as a dedicated, self-contained fitter (`fb_gev()`) -- the same
# pattern the genomic fitters `fb_gwas()` / `fb_gblup_cv()` follow -- with a
# matching family-object constructor (`fb_family_gev()`) registered against
# the family system so `family = "gen_extreme_value"` is recognised and the
# user is routed to `fb_gev()` rather than refused generically.
#
# The fitter uses dependency-free maximum likelihood (base `optim`). A
# scalable Bayesian GEV belongs on INLA's native `gev` / `bgev` family
# rather than the greta backend, which has no GEV distribution and whose
# free-form-likelihood idioms do not recover the shape reliably; that route
# is left for a future release rather than wired half-heartedly here.
# References: Coles (2001), *An Introduction to Statistical Modeling of
# Extreme Values*; Jenkinson (1955).

# ---- Family-object constructor -----------------------------------------

#' Generalised extreme value (GEV) family object
#'
#' Constructs the family descriptor for the generalised extreme value
#' distribution, the limiting law of block maxima. The object mirrors the
#' shape of a base `stats::family()` descriptor -- a named list carrying the
#' family name, the parameter names, and the natural support -- so it sits
#' alongside the other flexyBayes families, while signalling that GEV is
#' fitted through the dedicated `fb_gev()` entry point rather than the
#' GLM-link emit path.
#'
#' The GEV is parameterised by a location \eqn{\mu}, a positive scale
#' \eqn{\sigma}, and a shape \eqn{\xi} (the extreme value index). The shape
#' governs the tail: \eqn{\xi > 0} gives the heavy-tailed Frechet type,
#' \eqn{\xi < 0} the bounded Weibull type, and \eqn{\xi \to 0} the Gumbel
#' limit.
#'
#' @returns An object of class `c("fb_family_gev", "fb_family")`: a list with
#'   `family` (the canonical string `"gen_extreme_value"`), `parameters`
#'   (the character vector `c("location", "scale", "shape")`), `n_par` (the
#'   integer `3`), `link` (`"identity"` on the location), and `fitter`
#'   (`"fb_gev"`).
#'
#' @seealso [fb_gev()]
#' @examples
#' fam <- fb_family_gev()
#' fam$parameters
#' @export
fb_family_gev <- function() {
  structure(
    list(
      family = "gen_extreme_value",
      parameters = c("location", "scale", "shape"),
      n_par = 3L,
      link = "identity",
      fitter = "fb_gev"
    ),
    class = c("fb_family_gev", "fb_family")
  )
}

# ---- Density / simulation primitives ------------------------------------

# GEV log-density in the Jenkinson / von Mises parameterisation. `xi` is the
# shape; the Gumbel limit (`abs(xi) < tol`) is handled separately because the
# `(1 + xi z)` term degenerates there. Out-of-support points (where
# `1 + xi z <= 0`) return -Inf rather than NaN so the likelihood is well
# defined across the whole real line.
.dgev_log <- function(y, location, scale, shape, tol = 1e-08) {
  if (scale <= 0) {
    return(rep(-Inf, length(y)))
  }
  z <- (y - location) / scale
  if (abs(shape) < tol) {
    return(-log(scale) - z - exp(-z))
  }
  t <- 1 + shape * z
  out <- rep(-Inf, length(y))
  ok <- t > 0
  t_ok <- t[ok]
  out[ok] <- -log(scale) -
    (1 + 1 / shape) * log(t_ok) -
    t_ok^(-1 / shape)
  out
}

#' Simulate from a generalised extreme value distribution
#'
#' Draws block-maxima observations from a GEV with the given location, scale,
#' and shape, by inverting the GEV cumulative distribution function. The
#' Gumbel limit (\eqn{\xi \to 0}) is handled exactly.
#'
#' @param n Integer. The number of observations to draw.
#' @param location Numeric. The location parameter \eqn{\mu}.
#' @param scale Numeric. The positive scale parameter \eqn{\sigma}.
#' @param shape Numeric. The shape parameter \eqn{\xi}. Defaults to `0`
#'   (Gumbel).
#'
#' @returns A numeric vector of length `n`.
#'
#' @seealso [fb_gev()]
#' @examples
#' set.seed(1)
#' y <- rgev(100L, location = 10, scale = 2, shape = 0.15)
#' summary(y)
#' @export
rgev <- function(n, location, scale, shape = 0) {
  .check_positive_scalar(scale, "scale")
  u <- stats::runif(n)
  if (abs(shape) < 1e-08) {
    return(location - scale * log(-log(u)))
  }
  location + scale * ((-log(u))^(-shape) - 1) / shape
}

# ---- Fitter -------------------------------------------------------------

#' Fit a generalised extreme value distribution to block maxima
#'
#' Estimates the location, scale, and shape of a GEV distribution from a
#' vector of block-maxima observations by maximum likelihood (base `optim()`,
#' dependency-free and deterministic given the data).
#'
#' Maximum likelihood maximises the GEV log-likelihood directly, with the
#' scale optimised on the log scale to keep it positive and the standard
#' errors recovered from the observed-information (Hessian) matrix by the
#' delta method on the scale.
#'
#' Return levels (the level exceeded once per `m` blocks on average) are
#' computed from the fitted parameters and reported on `fit$return_levels`
#' for the requested return periods.
#'
#' A scalable Bayesian GEV belongs on INLA's native `gev` / `bgev` family
#' and is planned for a future release; the greta backend ships no GEV
#' distribution.
#'
#' @param y Numeric vector of block-maxima observations. Must contain at
#'   least four finite values.
#' @param return_periods Numeric vector of return periods (in blocks) at which
#'   to report return levels. Defaults to `c(10, 50, 100)`.
#' @param conf_level Numeric in `(0, 1)`. The interval level for the parameter
#'   summary. Defaults to `0.95`.
#'
#' @returns An object of class `c("fb_gev_fit", "fb_family_fit")`: a list with
#'   `estimates` (a `data.frame` of `term` / `estimate` / `std.error` /
#'   `conf.low` / `conf.high`), `return_levels` (a `data.frame` of
#'   `return_period` / `return_level`), `method`, `n_obs`, and `logLik`.
#'
#' @seealso [fb_family_gev()], [rgev()], [tidy.fb_gev_fit()]
#' @examples
#' set.seed(1)
#' y <- rgev(200L, location = 10, scale = 2, shape = 0.1)
#' fit <- fb_gev(y)
#' fit$estimates
#' @export
fb_gev <- function(
  y,
  return_periods = c(10, 50, 100),
  conf_level = 0.95
) {
  y <- .check_block_maxima(y)

  fit <- .fb_gev_ml(y, conf_level)

  est <- fit$estimates
  pars <- stats::setNames(est$estimate, est$term)
  fit$return_levels <- .gev_return_levels(
    location = pars[["location"]],
    scale = pars[["scale"]],
    shape = pars[["shape"]],
    return_periods = return_periods
  )
  fit$method <- "ml"
  fit$n_obs <- length(y)
  structure(fit, class = c("fb_gev_fit", "fb_family_fit"))
}

# Maximum-likelihood GEV fit. Optimises (location, log-scale, shape) and
# recovers SEs from the Hessian, with the scale SE delta-corrected back to
# the natural scale. The line search occasionally probes out-of-support
# parameters; the penalised negative log-likelihood returns a large finite
# value there, and the NaN warnings that base `log()` raises on those probes
# are suppressed so the fit is quiet.
.fb_gev_ml <- function(y, conf_level) {
  nll <- function(par) {
    loc <- par[1L]
    scale <- exp(par[2L])
    shape <- par[3L]
    ll <- .dgev_log(y, loc, scale, shape)
    if (any(!is.finite(ll))) {
      return(1e10)
    }
    -sum(ll)
  }

  start <- c(stats::median(y), log(stats::IQR(y) / 2 + 0.1), 0.1)
  opt <- suppressWarnings(stats::optim(
    start,
    nll,
    method = "BFGS",
    hessian = TRUE
  ))
  if (!isTRUE(opt$convergence == 0L)) {
    warning(
      "fb_gev(): the maximum-likelihood optimiser did not converge ",
      "(convergence code ", opt$convergence, "). Treat the estimates as ",
      "provisional and inspect the data.",
      call. = FALSE
    )
  }

  loc <- opt$par[1L]
  log_scale <- opt$par[2L]
  scale <- exp(log_scale)
  shape <- opt$par[3L]

  se <- .safe_hessian_se(opt$hessian)
  # Delta method: SE(scale) = SE(log_scale) * d exp(log_scale) = SE * scale.
  se_natural <- c(se[1L], se[2L] * scale, se[3L])

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  estimate <- c(loc, scale, shape)
  estimates <- data.frame(
    term = c("location", "scale", "shape"),
    estimate = estimate,
    std.error = se_natural,
    conf.low = estimate - z * se_natural,
    conf.high = estimate + z * se_natural,
    stringsAsFactors = FALSE
  )

  list(estimates = estimates, logLik = -opt$value)
}

# Return level for return period m: the (1 - 1/m) quantile of the GEV.
.gev_return_levels <- function(location, scale, shape, return_periods) {
  p <- 1 - 1 / return_periods
  yp <- -log(p)
  level <- if (abs(shape) < 1e-08) {
    location - scale * log(yp)
  } else {
    location + scale * (yp^(-shape) - 1) / shape
  }
  data.frame(
    return_period = return_periods,
    return_level = level,
    stringsAsFactors = FALSE
  )
}

# ---- tidy / print -------------------------------------------------------

#' Tidy a GEV fit
#'
#' Returns the GEV parameter summary as a `broom`-style `data.frame`, one row
#' per parameter (`location`, `scale`, `shape`), with the canonical columns.
#'
#' @param x An `fb_gev_fit` object from [fb_gev()].
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns A `data.frame` with `term`, `estimate`, `std.error`, `conf.low`,
#'   and `conf.high`.
#'
#' @seealso [fb_gev()]
#' @examples
#' set.seed(1)
#' fit <- fb_gev(rgev(200L, 10, 2, 0.1))
#' tidy(fit)
#' @export
tidy.fb_gev_fit <- function(x, ...) {
  x$estimates
}

#' Print a GEV fit
#'
#' @param x An `fb_gev_fit` object.
#' @param ... Currently unused; present for generic compatibility.
#' @returns `x`, invisibly.
#' @export
print.fb_gev_fit <- function(x, ...) {
  cat("Generalised extreme value fit  [flexyBayes]\n")
  cat(strrep("-", 50L), "\n")
  cat("  Method :", x$method, "  N =", x$n_obs, "\n")
  est <- x$estimates
  for (i in seq_len(nrow(est))) {
    cat(sprintf(
      "  %-9s %8.3f  (%.3f, %.3f)\n",
      est$term[i],
      est$estimate[i],
      est$conf.low[i],
      est$conf.high[i]
    ))
  }
  if (!is.null(x$return_levels)) {
    cat("  Return levels:\n")
    rl <- x$return_levels
    for (i in seq_len(nrow(rl))) {
      cat(sprintf(
        "    %g-block: %8.3f\n",
        rl$return_period[i],
        rl$return_level[i]
      ))
    }
  }
  cat(strrep("-", 50L), "\n")
  invisible(x)
}

# ---- Shared validation helpers -----------------------------------------

# Validate and clean a block-maxima vector: numeric, finite, length >= 4.
.check_block_maxima <- function(y) {
  if (!is.numeric(y)) {
    stop("`y` must be a numeric vector of block maxima.", call. = FALSE)
  }
  y <- y[is.finite(y)]
  if (length(y) < 4L) {
    stop(
      "`y` must contain at least four finite values to identify the three ",
      "GEV parameters; got ", length(y), ".",
      call. = FALSE
    )
  }
  y
}

# Validate a positive scalar argument.
.check_positive_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0) {
    stop("`", name, "` must be a single positive number.", call. = FALSE)
  }
  invisible(value)
}

# Standard errors from a Hessian, guarding a non-invertible matrix: when the
# Hessian cannot be inverted (flat or singular likelihood) the SEs are NA
# rather than an error, and the caller's intervals widen to NA accordingly.
.safe_hessian_se <- function(hessian) {
  vc <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(vc)) {
    return(rep(NA_real_, ncol(hessian)))
  }
  d <- diag(vc)
  d[d < 0] <- NA_real_
  sqrt(d)
}
