# boundary_collapse.R -- a variance component pinned at zero says so.
#
# `.fb_spatial_collapse_reasons()` in emit_inla.R catches exactly one
# term type: `sd_spatial` running to its floor against the nugget. Every
# other variance component can do the same thing, and did so silently.
#
# Measured 2026-08-28 on `agridat::besag.met` county C1,
# `yield ~ rep, random = ~ gen` on INLA. Three runs of the identical
# call returned an `sd_gen` upper credible bound of 0.00396, 0.0176 and
# 11.75 against a residual SD near 15: the fit is bistable, landing
# either on a degenerate mode with the component pinned at zero or on a
# well-identified one, and which it finds is not reproducible. This is
# the INLA degenerate-mode non-reproducibility the getting-started and
# dispatch vignettes already teach, reaching a plain iid random effect.
#
# On a collapsed run `summary(fit)$varcomp` does mark the row
# `note = "collapsed"`. What was missing is a warning: a note in a table
# column is easy to read past, and a zero genotype variance read as a
# result rather than as a failed mode is a scientific conclusion.
#
# The threshold is calibrated for the degenerate mode, not for a small
# variance. A component that is genuinely zero still returns a posterior
# well away from the boundary -- 240 rows simulated with no group effect
# at all give an upper bound near 0.43 of the residual SD -- so this
# does not fire on an honestly-small component.
#
# The detector is deliberately the same shape as the spatial one and
# reads the same canonical table (`extras$variance_comps`, columns
# `component`, `estimate`, `sd`, `q2.5`, `q97.5`, plus a
# `posterior_median` attribute), which both active engines populate. It
# says only what the table shows: an upper credible bound that is a
# negligible fraction of the residual scale. It does not claim to know
# whether the component is truly zero.

# The same fraction the spatial detector uses, for one reason: a reader
# comparing the two warnings should not have to learn two thresholds.
.FB_BOUNDARY_COLLAPSE_FRACTION <- 0.05

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
    vc$q97.5[candidates] < .FB_BOUNDARY_COLLAPSE_FRACTION * sigma
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
