# inla_standardise.R -- Gaussian INLA fits on a standardised response.
#
# A Gaussian fit on INLA runs on y* = (y - m) / s, where m is the mean of the
# response when the model has an intercept and 0 otherwise, and s is its
# standard deviation. Every prior the package or the user states is restated
# on those coordinates. A quantity that carries no stated prior takes INLA's
# own default on the standardised coordinates, which makes that default
# scale-free: the fit of a response in kilograms and the fit of the same
# response in tonnes are one model. The stored result is mapped back where the
# fit is built, and draws are mapped back where they are taken, because INLA's
# samplers read the engine's internal state rather than the stored tables.

.FB_INLA_STANDARDISED_RESPONSE <- "fb_standardised_response"

# --- the standardisation record ----------------------------------------- #

#' The centre and scale a Gaussian INLA fit is run on
#'
#' @param fb An `fb_terms` object.
#' @param data The data the fit is built from.
#' @returns A list with `centre`, `scale` and `n_obs`, the number of finite
#'   responses.
#' @noRd
#' @keywords internal
.fb_inla_standardisation <- function(fb, data) {
  y <- suppressWarnings(as.numeric(.fb_response_values(fb, data)))
  y <- y[is.finite(y)]
  scale <- if (length(y) >= 2L) stats::sd(y) else NA_real_
  centre <- if (isTRUE(fb$intercept) && length(y) > 0L) mean(y) else 0
  .fb_inla_standardisation_record(centre, scale, n_obs = length(y))
}

#' The standardisation of an aggregated Gaussian fit, from its statistics
#'
#' The aggregated route never holds the rows, so the centre and scale are
#' read off the cell sufficient statistics: the total count, the sum of the
#' response and the sum of its squares.
#'
#' @param fb An `fb_terms` object.
#' @param n,s1,s2 The total count, sum and sum of squares of the response.
#' @param n_cells The number of cells the aggregated likelihood is written
#'   over, which is the number of observations INLA sees.
#' @returns The record described for `.fb_inla_standardisation()`.
#' @noRd
#' @keywords internal
.fb_inla_standardisation_from_stats <- function(fb, n, s1, s2, n_cells) {
  resp_var <- if (n > 1) (s2 - s1^2 / n) / (n - 1) else NA_real_
  scale <- if (is.finite(resp_var) && resp_var > 0) sqrt(resp_var) else NA_real_
  centre <- if (isTRUE(fb$intercept) && n > 0) s1 / n else 0
  .fb_inla_standardisation_record(centre, scale, n_obs = n_cells)
}

.fb_inla_standardisation_record <- function(centre, scale, n_obs) {
  # A response with no spread has nothing to divide by; the fit then runs on
  # the response shifted by its mean, which is still an exact transform.
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  list(centre = centre, scale = scale, n_obs = n_obs)
}

# The factor a hyperparameter is mapped back by. Every latent field and the
# residual carry the response's units, so a precision divides by s^2; a
# correlation or a copy's scaling is unchanged.
.fb_inla_hyper_factor <- function(name, spec) {
  if (startsWith(name, "Precision for ")) {
    return(1 / spec$scale^2)
  }
  1
}


# --- stated priors on the standardised coordinates ----------------------- #

#' One variance-component prior restated for a field divided by `k`
#'
#' A field whose standard deviation is sigma on the original scale has
#' standard deviation sigma / k on the standardised one. Each translatable
#' family is closed under that map: a uniform or half-normal scale divides by
#' `k`, a lognormal location moves by `-log(k)`, an exponential rate and a
#' PC upper bound follow sigma, and a log-gamma prior on the precision keeps
#' its shape with the rate divided by `k^2`.
#'
#' @param entry One `priors_to_inla()` entry.
#' @param k The divisor, the response scale.
#' @returns The mapped entry.
#' @noRd
#' @keywords internal
.fb_inla_map_hyper_entry <- function(entry, k) {
  meta <- entry$meta
  family <- meta$family %||% entry$prior %||% ""

  if (identical(family, "uniform")) {
    lower <- meta$lower / k
    upper <- meta$upper / k
    return(list(
      prior = .inla_uniform_sd_expr(upper = upper, lower = lower),
      meta = utils::modifyList(meta, list(lower = lower, upper = upper))
    ))
  }
  if (identical(family, "half_normal")) {
    scale <- meta$scale / k
    return(list(
      prior = .inla_halfnormal_sd_expr(scale),
      meta = utils::modifyList(meta, list(scale = scale))
    ))
  }
  if (identical(family, "lognormal")) {
    location <- (meta$location %||% 0) - log(k)
    return(list(
      prior = .inla_lognormal_sd_expr(meta$scale, location = location),
      meta = utils::modifyList(meta, list(location = location))
    ))
  }
  if (identical(family, "exponential")) {
    rate <- meta$rate * k
    return(list(
      prior = "pc.prec",
      param = c(1 / rate, exp(-1)),
      meta = utils::modifyList(meta, list(rate = rate))
    ))
  }
  if (identical(family, "pc.prec")) {
    return(list(
      prior = "pc.prec",
      param = c(entry$param[[1L]] / k, entry$param[[2L]]),
      meta = meta
    ))
  }
  if (identical(family, "loggamma")) {
    return(list(
      prior = "loggamma",
      param = c(entry$param[[1L]], entry$param[[2L]] / k^2),
      meta = meta
    ))
  }
  stop(
    "flexyBayes internal error: no standardisation rule for the INLA prior ",
    "\"", family, "\". A stated prior that cannot be restated on the ",
    "standardised scale would change the model the fit runs.",
    call. = FALSE
  )
}

#' Every stated INLA precision prior restated on the standardised coordinates
#'
#' @param hyper_ctrl The `priors_to_inla()` list, already checked for keys
#'   no term reads.
#' @param spec The standardisation record.
#' @returns `hyper_ctrl` with each entry divided by the response scale.
#' @noRd
#' @keywords internal
.fb_inla_standardise_hyper <- function(hyper_ctrl, spec) {
  lapply(hyper_ctrl, .fb_inla_map_hyper_entry, k = spec$scale)
}

#' The stated fixed-effect prior restated on the standardised coordinates
#'
#' INLA states the fixed-effect prior as a mean and a precision. On
#' y* = (y - m) / s a slope is beta / s, so its mean divides by s and its
#' precision multiplies by s^2; the intercept is (beta0 - m) / s. Only what
#' the caller stated is mapped; a coefficient with no stated prior keeps
#' INLA's default on the standardised coordinates.
#'
#' @param control_fixed The `control.fixed` list the emit built, possibly
#'   empty.
#' @param spec The standardisation record.
#' @returns The mapped `control.fixed` list.
#' @noRd
#' @keywords internal
.fb_inla_standardise_control_fixed <- function(control_fixed, spec) {
  s <- spec$scale
  m <- spec$centre
  over <- function(v, f) {
    if (is.list(v)) lapply(v, f) else f(v)
  }
  out <- control_fixed
  if (!is.null(control_fixed$mean)) {
    out$mean <- over(control_fixed$mean, function(v) v / s)
  }
  if (!is.null(control_fixed$prec)) {
    out$prec <- over(control_fixed$prec, function(v) v * s^2)
  }
  if (!is.null(control_fixed$mean.intercept)) {
    out$mean.intercept <- (control_fixed$mean.intercept - m) / s
  }
  if (!is.null(control_fixed$prec.intercept)) {
    out$prec.intercept <- control_fixed$prec.intercept * s^2
  }
  out
}


# --- the result on the original scale ------------------------------------ #

# Columns of an INLA summary table that move with the quantity, and the one
# that scales with its magnitude. `kld` and `ID` are left alone.
.FB_INLA_LOCATION_COLS <- c(
  "mean", "0.025quant", "0.5quant", "0.975quant", "mode"
)

.fb_inla_map_table <- function(tab, a, b) {
  if (is.null(tab) || !NROW(tab)) {
    return(tab)
  }
  cols <- intersect(.FB_INLA_LOCATION_COLS, colnames(tab))
  for (col in cols) {
    tab[, col] <- a + b * tab[, col]
  }
  if ("sd" %in% colnames(tab)) {
    tab[, "sd"] <- abs(b) * tab[, "sd"]
  }
  tab
}

.fb_inla_map_marginal <- function(marg, a, b) {
  if (is.null(marg) || !NROW(marg)) {
    return(marg)
  }
  marg[, 1L] <- a + b * marg[, 1L]
  marg[, 2L] <- marg[, 2L] / abs(b)
  marg
}

# One row of a table, mapped by row name.
.fb_inla_map_rows <- function(tab, a_of, b_of) {
  if (is.null(tab) || !NROW(tab)) {
    return(tab)
  }
  for (nm in rownames(tab)) {
    tab[nm, ] <- .fb_inla_map_table(tab[nm, , drop = FALSE], a_of(nm), b_of(nm))
  }
  tab
}

#' An INLA result fitted on a standardised response, on the original scale
#'
#' Every stored reading of the fit is mapped back: fixed effects, random
#' effects, hyperparameters, fitted values and the linear predictor, with
#' their marginals; the marginal likelihood, which loses `n log(s)`; DIC and
#' WAIC, which gain `2 n log(s)` while their effective parameter counts are
#' unchanged; and CPO, which divides by `s`. `$.args`, `$mode` and `$misc`
#' keep the engine's own coordinates, since they are what INLA ran and what
#' its samplers read.
#'
#' @param fit An `inla` result.
#' @param spec The standardisation record the fit was run under.
#' @returns The result on the original scale.
#' @noRd
#' @keywords internal
.fb_inla_back_transform <- function(fit, spec) {
  s <- spec$scale
  m <- spec$centre
  log_s <- log(s)

  # Fixed effects, latent fields and hyperparameters ---------------------

  fixed_a <- function(nm) if (identical(nm, "(Intercept)")) m else 0
  fit$summary.fixed <- .fb_inla_map_rows(
    fit$summary.fixed, fixed_a, function(nm) s
  )
  for (nm in names(fit$marginals.fixed)) {
    fit$marginals.fixed[[nm]] <- .fb_inla_map_marginal(
      fit$marginals.fixed[[nm]], fixed_a(nm), s
    )
  }

  for (nm in names(fit$summary.random)) {
    fit$summary.random[[nm]] <- .fb_inla_map_table(
      fit$summary.random[[nm]], 0, s
    )
    if (!is.null(fit$marginals.random[[nm]])) {
      fit$marginals.random[[nm]] <- lapply(
        fit$marginals.random[[nm]],
        .fb_inla_map_marginal,
        a = 0,
        b = s
      )
    }
  }

  hyper_b <- function(nm) .fb_inla_hyper_factor(nm, spec)
  fit$summary.hyperpar <- .fb_inla_map_rows(
    fit$summary.hyperpar, function(nm) 0, hyper_b
  )
  for (nm in names(fit$marginals.hyperpar)) {
    fit$marginals.hyperpar[[nm]] <- .fb_inla_map_marginal(
      fit$marginals.hyperpar[[nm]], 0, hyper_b(nm)
    )
  }
  # The internal hyperparameters are log precisions, which shift.
  internal_a <- function(nm) {
    if (startsWith(nm, "Log precision for ")) -2 * log_s else 0
  }
  fit$internal.summary.hyperpar <- .fb_inla_map_rows(
    fit$internal.summary.hyperpar, internal_a, function(nm) 1
  )
  for (nm in names(fit$internal.marginals.hyperpar)) {
    fit$internal.marginals.hyperpar[[nm]] <- .fb_inla_map_marginal(
      fit$internal.marginals.hyperpar[[nm]], internal_a(nm), 1
    )
  }

  # Fitted values and the linear predictor (identity link) --------------

  for (slot in c("fitted.values", "linear.predictor")) {
    summary_slot <- paste0("summary.", slot)
    marginal_slot <- paste0("marginals.", slot)
    fit[[summary_slot]] <- .fb_inla_map_table(fit[[summary_slot]], m, s)
    if (!is.null(fit[[marginal_slot]])) {
      fit[[marginal_slot]] <- lapply(
        fit[[marginal_slot]],
        .fb_inla_map_marginal,
        a = m,
        b = s
      )
    }
  }

  # Likelihood-based criteria --------------------------------------------

  if (!is.null(fit$mlik)) {
    fit$mlik[, 1L] <- fit$mlik[, 1L] - spec$n_obs * log_s
  }
  if (is.list(fit$dic)) {
    n_dic <- sum(is.finite(fit$dic$local.dic))
    for (nm in intersect(c("dic", "mean.deviance", "deviance.mean"),
                         names(fit$dic))) {
      fit$dic[[nm]] <- fit$dic[[nm]] + 2 * n_dic * log_s
    }
    if (!is.null(fit$dic$family.dic) && !is.null(fit$dic$family)) {
      n_by_family <- tapply(
        is.finite(fit$dic$local.dic), fit$dic$family, sum
      )
      fit$dic$family.dic <- fit$dic$family.dic +
        2 * as.numeric(n_by_family[seq_along(fit$dic$family.dic)]) * log_s
    }
    if (!is.null(fit$dic$local.dic)) {
      fit$dic$local.dic <- fit$dic$local.dic + 2 * log_s
    }
  }
  if (is.list(fit$waic)) {
    n_waic <- sum(is.finite(fit$waic$local.waic))
    if (!is.null(fit$waic$waic)) {
      fit$waic$waic <- fit$waic$waic + 2 * n_waic * log_s
    }
    if (!is.null(fit$waic$local.waic)) {
      fit$waic$local.waic <- fit$waic$local.waic + 2 * log_s
    }
  }
  if (is.list(fit$cpo) && !is.null(fit$cpo$cpo)) {
    fit$cpo$cpo <- fit$cpo$cpo / s
  }
  if (is.list(fit$po) && !is.null(fit$po$po)) {
    fit$po$po <- fit$po$po / s
  }

  fit
}


# --- draws on the original scale ----------------------------------------- #

#' Posterior draws from a standardised INLA fit, on the original scale
#'
#' The one entry every caller takes draws through.
#' `INLA::inla.posterior.sample()` reads the engine's internal configurations,
#' which stay on the standardised scale whatever the stored tables say, so
#' each latent row and hyperparameter is mapped back here by the same rule the
#' stored result was.
#'
#' @param object A `flexybayes_inla` fit.
#' @param n Number of draws.
#' @param ... Passed to `INLA::inla.posterior.sample()`.
#' @returns The sampler's list of draws, on the original scale.
#' @noRd
#' @keywords internal
.fb_inla_posterior_sample <- function(object, n, ...) {
  draws <- INLA::inla.posterior.sample(n, object$inla, ...)
  spec <- object$extras$response_standardisation
  if (is.null(spec) || !length(draws)) {
    return(draws)
  }

  rows <- rownames(draws[[1L]]$latent)
  prefixes <- sub(":[0-9]+$", "", rows)
  fixed <- rownames(object$inla$summary.fixed) %||% character(0)
  random <- names(object$inla$summary.random) %||% character(0)
  a <- numeric(length(rows))
  for (p in unique(prefixes)) {
    if (p %in% c("Predictor", "(Intercept)")) {
      a[prefixes == p] <- spec$centre
    } else if (!p %in% c(fixed, random)) {
      stop(
        "flexyBayes internal error: the INLA draw row \"", p, "\" has no ",
        "rule for mapping back from the standardised scale.",
        call. = FALSE
      )
    }
  }

  lapply(draws, function(z) {
    z$latent[, 1L] <- a + spec$scale * z$latent[, 1L]
    if (length(z$hyperpar)) {
      z$hyperpar <- z$hyperpar * vapply(
        names(z$hyperpar), .fb_inla_hyper_factor, numeric(1L), spec = spec
      )
    }
    z
  })
}

#' Joint hyperparameter draws from a standardised INLA fit, on the original
#' scale
#'
#' @param object A `flexybayes_inla` fit.
#' @param n Number of draws.
#' @param ... Passed to `INLA::inla.hyperpar.sample()`.
#' @returns A matrix of draws, one column per hyperparameter.
#' @noRd
#' @keywords internal
.fb_inla_hyperpar_sample <- function(object, n, ...) {
  smp <- INLA::inla.hyperpar.sample(n, object$inla, ...)
  spec <- object$extras$response_standardisation
  if (is.null(spec) || !is.matrix(smp)) {
    return(smp)
  }
  for (nm in colnames(smp)) {
    smp[, nm] <- smp[, nm] * .fb_inla_hyper_factor(nm, spec)
  }
  smp
}
