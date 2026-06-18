# Truth-display helpers --- v0.3.8
#
# Shared helpers for print.flexybayes() / summary.flexybayes() /
# print.flexybayes_aggregated() / summary.flexybayes_aggregated() that
# render the two adjacent "Representation:" and "Engine:" lines. The
# split disambiguates representation regime (exact vs aggregated_exact)
# from inference engine (greta MCMC vs INLA Laplace vs brms / Stan
# HMC) --- closes the 2026-05-25 audit Critical Fix #4.
#
# These helpers also back fb_plan()'s engine_label / representation_label
# slots (see R/fb_plan.R for the plan-side equivalents .engine_label_for()
# and .representation_label_for() --- the fit-side helpers below take a
# materialised fit object plus its $extras$backend_decision trace).

# .repr_label_for_fit() --- the Representation: line content.
# Reads $exactness directly (the canonical slot) and
# attaches the compression ratio for aggregated_exact when N/K >= 2.
# v0.3.10: blocks-format vm/ped fits surface
# the "(block-diagonal, K blocks)" annotation by reading the routing
# trace's representation_plan for block_diagonal entries.
.repr_label_for_fit <- function(fit, bd) {
  exact <- fit$exactness %||% "exact"
  # Approximate fits
  # carry a bracketed [APPROX: <scheme>] badge --- the load-bearing,
  # monochrome-friendly visual distinction marking that the fit is not
  # exact. validate_approximation(fit) reports the realised bias bound.
  if (startsWith(exact, "approximate_")) {
    scheme <- sub("^approximate_", "", exact)
    return(paste0("approximate [APPROX: ", scheme, "]"))
  }
  if (identical(exact, "aggregated_exact")) {
    am <- fit$extras$aggregation_meta
    if (
      !is.null(am) &&
        !is.null(am$N) &&
        !is.null(am$K) &&
        am$K > 0L &&
        am$N / am$K >= 2
    ) {
      return(sprintf("aggregated_exact (compression %.0f:1)", am$N / am$K))
    }
    return("aggregated_exact")
  }
  block_entries <- if (!is.null(bd$representation_plan)) {
    Filter(
      function(rp) {
        identical(
          rp$representation_class %||% NA_character_,
          .representation_class("block_diagonal")
        )
      },
      bd$representation_plan
    )
  } else {
    NULL
  }
  if (length(block_entries)) {
    k_total <- sum(vapply(
      block_entries,
      function(rp) {
        as.integer(rp$block_count %||% NA_integer_)
      },
      integer(1L)
    ))
    if (!is.na(k_total) && k_total > 0L) {
      return(sprintf("exact (block-diagonal, %d blocks)", k_total))
    }
    return("exact (block-diagonal)")
  }
  exact
}

# .engine_label_for_fit() --- the Engine: line content. Maps the
# backend_decision trace's $backend + $path tuple to a human-readable
# inference-engine + approximation-regime phrase.
.engine_label_for_fit <- function(fit, bd) {
  backend <- bd$backend %||% NA_character_
  path <- bd$path %||% NA_character_
  if (is.na(backend)) {
    return("(unknown engine)")
  }
  if (identical(backend, "greta")) {
    return("greta MCMC")
  }
  if (identical(backend, "brms")) {
    return("brms / Stan HMC")
  }
  if (identical(backend, "inla")) {
    if (identical(path, "aggregated_inla")) {
      return("INLA Laplace (aggregated)")
    }
    return("INLA Laplace")
  }
  if (identical(backend, "gretaR")) {
    return("gretaR (R-native; dormant)")
  }
  backend
}
