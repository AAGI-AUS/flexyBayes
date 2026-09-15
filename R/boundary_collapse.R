# boundary_collapse.R -- a variance component pinned at zero says so.
#
# `.fb_spatial_collapse_reasons()` in emit_inla.R catches exactly one
# term type: `sd_spatial` running to its floor against the nugget. Every
# other variance component can do the same thing, and did so silently.
#
# On a collapsed run `summary(fit)$varcomp` does mark the row
# `note = "collapsed"`. What was missing is a warning: a note in a table
# column is easy to read past, and a zero genotype variance read as a
# result rather than as a failed mode is a scientific conclusion.
#
# The threshold is calibrated against a measured floor, not a guess. An
# earlier version of this file put it at 0.05 of the residual SD on the
# strength of one simulation at n = 240, and claimed a genuinely null
# component returns an upper bound near 0.43 of the residual SD. The
# package's own test suite falsified that: `mk_inla_data()` in
# test-emit-inla.R is `y = rnorm(30)` with `g` assigned cyclically -- no
# group effect exists -- and the fit returns an upper bound at 0.025 of
# the residual SD, so the warning fired on seven fixtures whose
# components are null by construction.
#
# Measured 2026-08-28 by a sweep of 112 INLA fits with the group SD set
# to exactly zero, crossing n in {30, 60, 120, 240, 480} with 5, 10 and
# 20 groups, 8 seeds per cell (scratch script, `y ~ 1 + f(g, "iid")`,
# package default prior). The ratio q97.5(sd_g) / sigma is flat in both
# n and the group count -- min 0.0228, 5th percentile 0.0239, median
# 0.0263 -- because it is a floor set by the prior, not by the data. A
# truly null component does not go below about 0.023.
#
# The degenerate mode sits two orders of magnitude lower: `agridat`
# besag.met county C1 returned sd_gen upper bounds of 0.00396 against a
# residual SD near 15.6 (a ratio of 2.5e-4) on one run of three, the
# other two returning 0.0176 and 11.75 from the identical call. The fit
# is bistable and which mode it finds is not reproducible. This is the
# INLA degenerate-mode non-reproducibility the getting-started and
# dispatch vignettes already teach, reaching a plain iid random effect.
#
# The threshold is `.FB_COLLAPSE_FRACTION` (0.01), reused rather than
# redefined. That is the same constant, against the same residual
# reference, that `.fb_varcomp_notes()` uses to mark a row
# `note = "collapsed"` in `summary(fit)$varcomp`, and its comment there
# says the two surfaces agreeing is the point of having both. A separate
# number here would have opened a band where the table says collapsed
# and the warning stays silent. 0.01 sits about 2.3 times below the
# lowest null ratio measured and about 40 times above the degenerate
# one. The floor is a property of the default prior; a user who sets a
# very different prior on the component moves it.
#
# On a collapsed run `summary(fit)$varcomp` does mark the row
# `note = "collapsed"`. What was missing is a warning: a note in a table
# column is easy to read past, and a zero genotype variance read as a
# result rather than as a failed mode is a scientific conclusion.
#
# Two limits on reach, both structural rather than oversight. The
# detector needs a residual SD to compare against, so it reads Gaussian-
# scale fits and stays silent on families that carry no `sigma` row. And
# it needs credible bounds, so it covers the per-row emit paths only:
# the aggregated (streaming) emitters store `variance_comps` as posterior
# means with no quantiles (R/emit_gaussian_aggregated.R,
# R/emit_count_aggregated.R), and there is nothing there to read an upper
# bound from until those emitters compute one.
#
# The detector is deliberately the same shape as the spatial one and
# reads the same canonical table (`extras$variance_comps`, columns
# `component`, `estimate`, `sd`, `q2.5`, `q97.5`, plus a
# `posterior_median` attribute), which both active engines populate. It
# says only what the table shows: an upper credible bound that is a
# negligible fraction of the residual scale. It does not claim to know
# whether the component is truly zero.

# Residual scale for the fit, on the SD scale. Prefer the posterior
# median attribute the canonical table carries; fall back to the table's
# own `sigma` row.
#
# @noRd
# @keywords internal
.fb_residual_scale <- function(vc) {
  medians <- attr(vc, "posterior_median") %||% numeric(0)
  if ("sigma" %in% names(medians)) {
    val <- medians[["sigma"]]
    if (is.numeric(val) && length(val) == 1L && is.finite(val) && val > 0) {
      return(val)
    }
  }
  i <- match("sigma", as.character(vc$component))
  if (!is.na(i) && is.finite(vc$estimate[[i]]) && vc$estimate[[i]] > 0) {
    return(vc$estimate[[i]])
  }
  NA_real_
}

# Components whose upper credible bound is a negligible fraction of the
# residual scale. `sd_spatial` is excluded: the spatial detector already
# reports it, with grid-specific advice this one cannot give, and two
# warnings for one fact is worse than one.
#
# @noRd
# @keywords internal
.fb_boundary_collapse_reasons <- function(object) {
  vc <- object$extras$variance_comps
  if (!is.data.frame(vc) || nrow(vc) == 0L) {
    return(character(0L))
  }
  if (!all(c("component", "q97.5") %in% names(vc))) {
    return(character(0L))
  }
  sigma <- .fb_residual_scale(vc)
  if (!is.finite(sigma)) {
    return(character(0L))
  }

  cmp <- as.character(vc$component)
  candidates <- which(
    startsWith(cmp, "sd_") &
      cmp != "sd_spatial" &
      is.finite(vc$q97.5)
  )
  hit <- candidates[
    vc$q97.5[candidates] < .FB_COLLAPSE_FRACTION * sigma
  ]
  if (length(hit) == 0L) {
    return(character(0L))
  }

  vapply(
    hit,
    function(i) {
      sprintf(
        "%s (upper credible bound %.4g against a residual SD of %.4g)",
        cmp[[i]],
        vc$q97.5[[i]],
        sigma
      )
    },
    character(1)
  )
}

# @noRd
# @keywords internal
.fb_warn_boundary_collapse <- function(object) {
  silent <- getOption(
    "flexyBayes.silence_boundary_collapse_warning", FALSE
  )
  if (isTRUE(silent)) {
    return(invisible(NULL))
  }
  reasons <- .fb_boundary_collapse_reasons(object)
  if (length(reasons) == 0L) {
    return(invisible(NULL))
  }
  warning(
    "flexyBayes: a variance component in this fit is at the boundary -- ",
    paste(reasons, collapse = "; "),
    ". The term carries no signal the data could separate from residual ",
    "noise, and its variance has been absorbed into the residual. The ",
    "convergence block reports a converged mode regardless, because the ",
    "optimiser did converge, to a solution with that term at zero. This ",
    "is not always a property of the data: the same model on a random ",
    "subset of the same trial can return a well-identified component. ",
    "Three routes: give the term an informative prior with fb_prior() so ",
    "it is not left to run to a boundary; fit the same model on the other ",
    "engine for a second reading; or report the component as ",
    "unidentified rather than as zero. Silence via ",
    "options(flexyBayes.silence_boundary_collapse_warning = TRUE).",
    call. = FALSE
  )
  invisible(NULL)
}
