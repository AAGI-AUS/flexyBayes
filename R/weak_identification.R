# weak_identification.R: a warning for a variance component the data
# barely identify.
#
# Since 0.13.0 the default prior on a standard deviation is a half-normal
# with no upper bound, so a component the data cannot identify is reported
# at whatever value the prior and the likelihood settle on. The boundary
# check (R/boundary_collapse.R) catches a component pinned at zero. This
# file catches one that sits high because only the prior holds it there.
#
# The reading is the share of the prior's probability that the posterior's
# 95% interval still spans: F(upper) - F(lower), with F the prior's
# distribution function on the standard deviation. A posterior the data
# have narrowed spans little of the prior. A posterior that is the prior
# spans 95% of it. The share is unchanged by any monotone transform of the
# standard deviation and by any rescaling of the response or a covariate
# that rescales the prior with it, which the package defaults do
# (R/inla_default_priors.R). A null component also spans little: the data
# cut away the prior's upper tail.
#
# Measured on the default half-normal prior, INLA, one thread, 10 seeds per
# design, `y ~ 1 + (1 | g)` with group and residual SD 1 unless stated:
#
#   design                                 share of prior probability
#   3 levels of 10 rows                    0.55 to 0.80
#   random copy of a 6-level fixed factor  0.88
#   4 levels of 10 rows                    0.24 to 0.49
#   5 levels of 10 rows                    0.19 to 0.36
#   8 levels of 10 rows                    0.15 to 0.20
#   30 levels of 8 rows                    0.07
#   5 levels of 6 rows, group SD 0         0.16 to 0.20
#   20 levels of 8 rows, group SD 0        0.05 to 0.07
#
# An uncorrelated slope over 3 levels gave 0.49 to 0.60 and over 20 levels
# 0.05 to 0.06 (3 seeds), the same to 0.001 with the covariate multiplied
# by 1000, and multiplying the response by 1000 left the 3-level share
# unchanged. The threshold 0.4 sits between the highest share of any design
# with five or more levels (0.36) and the lowest at three levels (0.49).
#
# The share needs the prior the engine received and a posterior on the same
# coordinates, so it is computed on INLA fits. A prior this file cannot read
# leaves the component unflagged.
#
# The second reading is the width of the 95% interval over the posterior
# median, both on the standard deviation. Over the same ten seeds
# (inst/extdata/weak_identification_calibration.csv):
#
#   design                                 interval / median
#   3 levels of 10 rows                    4.30 to 5.15
#   random copy of a 6-level fixed factor  3.21 to 3.22
#   4 levels of 10 rows                    2.61 to 3.58
#   5 levels of 10 rows                    1.97 to 2.66
#   8 levels of 10 rows                    1.30 to 1.50
#   30 levels of 8 rows                    0.58 to 0.64
#   5 levels of 6 rows, group SD 0         2.55 to 4.24
#   20 levels of 8 rows, group SD 0        1.56 to 3.19
#
# A null component reads high because its median is near zero. The random
# copy of a fixed factor reads below the threshold, so either reading flags
# a component. Computed on INLA fits only.

# A component whose 95% interval spans at least this share of its prior's
# probability is weakly identified.
.FB_WEAK_IDENTIFICATION_SHARE <- 0.4

# The 99th percentile of the interval-to-median ratio over the ten-seed
# designs with five or more levels
# (data-raw/weak_identification_calibration.R).
.FB_WEAK_IDENTIFICATION_RATIO <- 4.24

# The distribution function, on the standard deviation, of a prior as INLA
# reports it received it, or NULL for a prior this cannot read.
#
# @noRd
# @keywords internal
.fb_inla_prior_sd_cdf <- function(prior, param) {
  if (is.na(prior)) {
    return(NULL)
  }
  num <- "([-+0-9.eE]+)"
  grab <- function(pattern) {
    m <- regmatches(prior, regexec(pattern, prior))[[1L]]
    if (length(m)) as.numeric(m[-1L]) else NULL
  }
  if (startsWith(prior, "expression:")) {
    v <- grab(paste0("expression: *L=", num, "; *U=", num))
    if (!is.null(v)) {
      return(function(x) pmin(pmax((x - v[[1L]]) / (v[[2L]] - v[[1L]]), 0), 1))
    }
    v <- grab(paste0("expression: *U=", num))
    if (!is.null(v)) {
      return(function(x) pmin(x / v[[1L]], 1))
    }
    # The lognormal is on log sd, with location mu.
    v <- grab(paste0("expression: *s=", num, "; *mu=", num))
    if (!is.null(v)) {
      return(function(x) stats::pnorm((log(x) - v[[2L]]) / v[[1L]]))
    }
    v <- grab(paste0("expression: *s=", num))
    if (is.null(v)) {
      return(NULL)
    }
    # The lognormal is quadratic in theta. The half-normal is not.
    if (grepl("theta*theta", prior, fixed = TRUE)) {
      return(function(x) stats::pnorm(log(x) / v[[1L]]))
    }
    return(function(x) 2 * stats::pnorm(x / v[[1L]]) - 1)
  }
  p <- .fb_manifest_param_values(param)
  if (identical(prior, "pc.prec") && length(p) >= 2L) {
    # P(sd > u) = alpha, an exponential on the standard deviation.
    rate <- -log(p[[2L]]) / p[[1L]]
    return(function(x) stats::pexp(x, rate))
  }
  if (identical(prior, "loggamma") && length(p) >= 2L) {
    # A gamma prior on the precision 1 / sd^2.
    return(function(x) {
      stats::pgamma(x^-2, shape = p[[1L]], rate = p[[2L]], lower.tail = FALSE)
    })
  }
  NULL
}

# The share of the prior's probability inside each standard deviation's 95%
# interval, named by canonical component.
#
# The fit's marginals are on the response's units and the prior INLA
# received is on the coordinates it fitted, so each interval is divided by
# the factor that maps those coordinates back.
#
# @noRd
# @keywords internal
.fb_inla_prior_share <- function(object, mapper) {
  fit <- object$inla
  readback <- .fb_manifest_inla_readback(fit)
  rows <- readback$rows[readback$rows$section != "fixed", , drop = FALSE]
  theta_names <- rownames(fit$internal.summary.hyperpar) %||% character(0)
  out <- numeric(0)
  if (!isTRUE(readback$complete) || nrow(rows) != length(theta_names)) {
    return(out)
  }
  marginals <- fit$internal.marginals.hyperpar %||% list()
  for (i in seq_along(theta_names)) {
    nm <- theta_names[[i]]
    if (!startsWith(nm, "Log precision for ")) {
      next
    }
    # INLA lists its hyperparameters in the order of its prior record.
    paired <- switch(
      rows$section[[i]],
      random = endsWith(nm, paste0(" for ", rows$tag[[i]])),
      family = grepl(" observations(\\[[0-9]+\\])?$", nm),
      FALSE
    )
    canonical <- unname(mapper[sub("^Log precision", "Precision", nm)])
    cdf <- .fb_inla_prior_sd_cdf(rows$prior[[i]], rows$param[[i]])
    m <- marginals[[nm]]
    if (
      !paired ||
        length(canonical) != 1L ||
        is.na(canonical) ||
        is.null(cdf) ||
        is.null(m)
    ) {
      next
    }
    factor <- .fb_support_sd_factor(
      nm,
      object$extras$response_standardisation,
      object$extras$covariate_standardisation
    )
    limits <- exp(-INLA::inla.qmarginal(c(0.975, 0.025), m) / 2) / factor
    share <- cdf(limits[[2L]]) - cdf(limits[[1L]])
    if (is.finite(share)) {
      out[[canonical]] <- share
    }
  }
  out
}

# The width of each standard deviation's 95% interval over its posterior
# median, named by component.
#
# @noRd
# @keywords internal
.fb_interval_ratio <- function(vc) {
  medians <- attr(vc, "posterior_median") %||% numeric(0)
  cmp <- as.character(vc$component)
  keep <- cmp %in% names(medians)
  ratio <- (vc$q97.5[keep] - vc$q2.5[keep]) / unname(medians[cmp[keep]])
  ratio <- stats::setNames(ratio, cmp[keep])
  ratio[is.finite(ratio)]
}

# Whether each component is weakly identified: its interval spans at least
# the threshold share of its prior, or is at least the threshold multiple of
# its median. `sd_spatial` is left to the spatial field's own test.
#
# @noRd
# @keywords internal
.fb_weakly_identified <- function(components, prior_share, interval_ratio) {
  read <- function(x) {
    out <- rep(NA_real_, length(components))
    has <- components %in% names(x)
    out[has] <- x[components[has]]
    out
  }
  share <- read(prior_share)
  ratio <- read(interval_ratio)
  startsWith(components, "sd_") &
    components != "sd_spatial" &
    ((!is.na(share) & share >= .FB_WEAK_IDENTIFICATION_SHARE) |
      (!is.na(ratio) & ratio >= .FB_WEAK_IDENTIFICATION_RATIO))
}

# @noRd
# @keywords internal
.fb_weak_identification_reasons <- function(object) {
  vc <- object$extras$variance_comps
  if (!is.data.frame(vc) || nrow(vc) == 0L) {
    return(character(0L))
  }
  share <- attr(vc, "prior_share") %||% numeric(0)
  ratio <- attr(vc, "interval_ratio") %||% numeric(0)
  cmp <- as.character(vc$component)
  hit <- which(.fb_weakly_identified(cmp, share, ratio))
  if (length(hit) == 0L) {
    return(character(0L))
  }
  vapply(
    hit,
    function(i) {
      read <- c(
        if (cmp[[i]] %in% names(share)) {
          sprintf("%.0f%% of its prior's probability", 100 * share[[cmp[[i]]]])
        },
        if (cmp[[i]] %in% names(ratio)) {
          sprintf("%.2f times its median", ratio[[cmp[[i]]]])
        }
      )
      sprintf(
        "%s (95%% interval %.3g to %.3g, %s)",
        cmp[[i]],
        vc$q2.5[[i]],
        vc$q97.5[[i]],
        paste(read, collapse = ", ")
      )
    },
    character(1)
  )
}

# @noRd
# @keywords internal
.fb_warn_weak_identification <- function(object) {
  silent <- getOption(
    "flexyBayes.silence_weak_identification_warning",
    FALSE
  )
  if (isTRUE(silent)) {
    return(invisible(NULL))
  }
  reasons <- .fb_weak_identification_reasons(object)
  if (length(reasons) == 0L) {
    return(invisible(NULL))
  }
  warning(
    "flexyBayes: a variance component in this fit is weakly identified: ",
    paste(reasons, collapse = "; "),
    ". Its 95% interval spans 40% or more of its prior's probability, or ",
    "is 4.24 or more times its posterior median, so the data say little ",
    "about the reported value. ",
    "Common ",
    "causes are a grouping factor with few levels and a random term that ",
    "splits the rows as a fixed term already does. Report the component ",
    "as weakly identified, not as an estimate, state a prior you can ",
    "defend with fb_prior(), or remove the term the data cannot separate. ",
    "summary(fit)$varcomp marks the row `weakly identified`. Silence via ",
    "options(flexyBayes.silence_weak_identification_warning = TRUE).",
    call. = FALSE
  )
  invisible(NULL)
}
