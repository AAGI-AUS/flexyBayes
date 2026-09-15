# fb_posterior -- canonical posterior-draws view over any backend fit.
#
# The reshape makes a canonical posterior interface the single seam every
# output dialect and the backend-conformance battery read from. This file
# seeds that interface with `.fb_canonical_draws()`: it composes the raw
# backend draws (`fb_as_draws_simple()`) with the per-backend canonical-name
# map + value transforms (`canonical_names()`), so a fit on ANY backend
# yields the same canonical parameter tokens ((Intercept), the fixed-effect
# terms, sigma, sd_<group>) with variance components on the SD scale.
#
# INLA's raw `inla.posterior.sample()` output -- fixed effects indexed as
# "(Intercept):1", hyperparameters named "Precision for ...", plus one
# per-observation "Predictor:k" latent node -- is thereby canonicalised the
# same way `triangulate()` already does internally, but as a reusable
# accessor. Testing INLA as a conformance CANDIDATE (rather than only a
# name-mapped reference) needs this.
#
# Internal seed; the public `fb_posterior` interface (native marginals +
# optional joint draws + provenance + diagnostic-validity flags) builds on
# this in a later phase.

# --- canonical draws --------------------------------------------------- #

# .fb_canonical_draws() -- a named list of posterior draw vectors keyed by
# CANONICAL parameter name, for any backend carrying a canonical_names()
# mapper. `drop = TRUE` discards backend-native nuisance nodes (e.g. INLA's
# per-observation "Predictor:k") that map to no canonical parameter. When the
# backend supplies no mapper (an empty map), the raw draws are returned
# unchanged -- they are assumed already canonical.
#
# Value transforms (e.g. INLA precision -> SD) key on the native name and
# apply BEFORE any rename, matching triangulate()'s ordering, so a caller
# never double-applies them.
.fb_canonical_draws <- function(fit, n_samples = 1000L, drop = TRUE) {
  draws <- fb_as_draws_simple(fit, n_samples = n_samples)
  if (!is.list(draws) || is.null(names(draws))) {
    stop(
      ".fb_canonical_draws(): fb_as_draws_simple() must return a named ",
      "list of numeric vectors.",
      call. = FALSE
    )
  }

  cn <- tryCatch(
    canonical_names(fit, drop = drop),
    error = function(e) list(map = character(0), transform = list())
  )
  value_transform <- cn$transform %||% list()
  name_map <- cn$map %||% character(0)

  for (nm in names(value_transform)) {
    if (!is.null(draws[[nm]]) && is.function(value_transform[[nm]])) {
      draws[[nm]] <- value_transform[[nm]](draws[[nm]])
    }
  }

  # No mapper -> the draws are already canonical; return unchanged.
  if (length(name_map) == 0L) {
    return(draws)
  }

  out <- list()
  for (nm in names(draws)) {
    if (nm %in% names(name_map)) {
      out[[as.character(name_map[[nm]])]] <- draws[[nm]]
    } else if (!isTRUE(drop)) {
      out[[nm]] <- draws[[nm]]
    }
  }
  out
}

# --- MCMC convergence diagnostics -------------------------------------- #

# .fb_draws_array_or_null() -- a chain-structured `posterior::draws_array`
# for an MCMC-inference fit, or NULL when the backend is not MCMC. Only an
# MCMC backend carries a native multi-chain draws source (brms's `$brms`
# slot; a future MCMC backend may pre-store a multi-chain `$draws`); a
# deterministic approximation (INLA/Laplace) has none, so R-hat / bulk-ESS
# are undefined for it (v5-audit: INLA samples are not Markov chains).
.fb_draws_array_or_null <- function(fit) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(NULL)
  }
  if (!is.null(fit$brms)) {
    da <- tryCatch(
      posterior::as_draws_array(fit$brms),
      error = function(e) NULL
    )
    if (!is.null(da)) {
      return(da)
    }
  }
  if (inherits(fit$draws, "draws_array") &&
        isTRUE(posterior::nchains(fit$draws) >= 2L)) {
    return(fit$draws)
  }
  NULL
}

# .fb_mcmc_diagnostics() -- the worst-case bulk-ESS + R-hat over a fit's
# parameters, with `applicable = FALSE` for a non-MCMC backend. The
# conformance battery's ESS/R-hat gate reads this: it enforces a floor only
# where `applicable`, and treats a Laplace/INLA fit as exempt rather than
# fabricating an MCMC diagnostic on samples that are not Markov chains.
.fb_mcmc_diagnostics <- function(fit) {
  da <- .fb_draws_array_or_null(fit)
  if (is.null(da)) {
    return(list(
      applicable = FALSE, min_ess_bulk = NA_real_, max_rhat = NA_real_
    ))
  }
  s <- posterior::summarise_draws(
    da,
    ess_bulk = posterior::ess_bulk,
    rhat = posterior::rhat
  )
  list(
    applicable = TRUE,
    min_ess_bulk = suppressWarnings(min(s$ess_bulk, na.rm = TRUE)),
    max_rhat = suppressWarnings(max(s$rhat, na.rm = TRUE))
  )
}
