# inla_restart_probe.R -- is a converged mode the mode, or just a mode?
#
# `.lgm_check_numerical()` reports a fit as numerically confirmed when INLA
# returns `mode.status == 0` and a finite marginal likelihood. Neither
# distinguishes the solution the optimiser found from the solution it should
# have found, and this package has now measured that gap twice.
#
# On `agridat::besag.met` with `random = ~ diag(env):gen`, twenty-four fits
# across four prior settings returned no genetic variance in any environment
# on nineteen of them, and substantial variance on five. ASReml, on the same
# rows, returns 5.84 and 5.45 in two of the six environments, stably. Every
# one of the twenty-four reported `mode.status = 0` and
# `Numerical confirm: PASS`. Tightening the prior by a factor of twenty
# changed nothing, so this is mode selection rather than a prior this package
# chose -- which is what separates it from the factor-analytic specific
# variance and the `iid2d` Wishart, both of which were our own defaults.
#
# The project's own conclusion from those cases is the design here: only the
# marginal likelihood separates two converged modes. So when a component comes
# back pinned at zero, the fit is restarted from a mode where that component
# has real variance, and the two marginal likelihoods are compared. A restart
# that finds a materially better solution means the first one was not the
# answer, and the numerical-confirm verdict says so rather than passing.
#
# The probe costs one extra fit and runs only when something collapsed, so a
# fit with nothing pinned at zero pays nothing.

# How much better a restart's marginal log-likelihood must be before the
# first solution is called inferior. Two log units is the conventional
# "positive evidence" step on the Kass-Raftery scale, and it is far above
# the run-to-run movement of INLA's own integration on a stable fit.
.FB_RESTART_MLIK_MARGIN <- 2

# Ratio of residual variance below which a component counts as pinned at
# zero for this purpose. 1e-4 on the variance scale is the square of the
# package's 1% rule on the SD scale (.FB_COLLAPSE_FRACTION), so the two
# surfaces agree by construction rather than by coincidence.
.FB_RESTART_COLLAPSE_VARIANCE_RATIO <- 1e-4


# Names of the hyperparameters that came back pinned at zero variance.
#
# Read from the raw INLA hyperparameter table rather than from the
# canonicalised variance-component frame, because the probe runs inside the
# emit, before a flexybayes object exists.
#
# @noRd
.fb_inla_collapsed_hyper <- function(fit) {
  hp <- fit$summary.hyperpar
  if (is.null(hp) || !nrow(hp)) {
    return(character(0))
  }
  nm <- rownames(hp)
  resid_i <- which(nm == "Precision for the Gaussian observations")
  if (!length(resid_i)) {
    return(character(0))
  }
  resid_prec <- hp[resid_i[[1L]], "mean"]
  if (!is.finite(resid_prec) || resid_prec <= 0) {
    return(character(0))
  }
  prec_rows <- setdiff(grep("^Precision for ", nm), resid_i)
  if (!length(prec_rows)) {
    return(character(0))
  }
  # The LOW end of the precision interval: a component is pinned only when
  # even its most generous variance is negligible.
  low <- hp[prec_rows, "0.025quant"]
  keep <- is.finite(low) &
    low > resid_prec / .FB_RESTART_COLLAPSE_VARIANCE_RATIO
  nm[prec_rows][keep]
}


# Restart the fit from a mode where the pinned components carry real
# variance, and report whether that solution is better.
#
# Returns NULL when the probe does not apply or could not run -- a probe that
# fails is not evidence of stability, and the caller records the difference.
#
# @noRd
.fb_inla_restart_probe <- function(fit, inla_args, collapsed) {
  if (!length(collapsed)) {
    return(NULL)
  }
  theta <- fit$mode$theta
  if (is.null(theta) || !length(theta) || !all(is.finite(theta))) {
    return(NULL)
  }
  nm <- names(theta)
  if (is.null(nm)) {
    return(NULL)
  }
  # Guarded here rather than left to the caller's tryCatch: a fit with no
  # marginal likelihood has nothing to compare against, and that is a reason
  # to decline the probe, not an error to swallow.
  if (is.null(fit$mlik) || !nrow(fit$mlik)) {
    return(NULL)
  }
  resid_i <- grep("Gaussian observations", nm)
  resid_theta <- if (length(resid_i)) theta[[resid_i[[1L]]]] else stats::median(theta)

  # A pinned component's log precision has run away upward. Start it two log
  # units above the residual's, which is a variance about an eighth of the
  # residual's -- squarely inside the range a real component occupies, and
  # nowhere near the boundary the first fit ran to.
  start <- theta
  collapsed_theta <- sub("^Precision for ", "Log precision for ", collapsed)
  hit <- nm %in% collapsed_theta
  if (!any(hit)) {
    # INLA does not always spell a hyperparameter's theta name as the
    # summary table spells it, so fall back to matching on the term label.
    # This is deliberately loose, and it is safe to be loose: the fallback
    # only decides where the restart STARTS. A start on the wrong entry
    # makes the probe less likely to escape the basin, and cannot put a
    # wrong answer in the result, because which fit is returned is decided
    # by the marginal likelihood afterwards and not by this match.
    hit <- vapply(
      nm,
      function(x) {
        any(vapply(
          collapsed,
          function(c1) grepl(sub("^Precision for ", "", c1), x, fixed = TRUE),
          logical(1)
        ))
      },
      logical(1)
    )
  }
  if (!any(hit)) {
    return(NULL)
  }
  start[hit] <- resid_theta + 2

  args2 <- inla_args
  args2$control.mode <- list(theta = unname(start), restart = TRUE)
  # The restart runs in INLA's external program, and setTimeLimit() cannot
  # bound a fit that runs outside R (cairn
  # 2026-08-14-external-binary-timeouts-need-a-process-watchdog). A restart
  # that hangs therefore hangs the fit. Accepted at this stage because the probe
  # only fires on a fit that already converged once on the same data, so the
  # surface is known to be tractable, and a watchdog belongs with the one that
  # cairn calls for rather than as a second private implementation.
  fit2 <- tryCatch(.inla_call(args2), error = function(e) NULL)
  if (is.null(fit2) || is.null(fit2$mlik) || !nrow(fit2$mlik)) {
    return(NULL)
  }
  m1 <- fit$mlik[1, 1]
  m2 <- fit2$mlik[1, 1]
  if (!is.finite(m1) || !is.finite(m2)) {
    return(NULL)
  }
  list(
    restarted = TRUE,
    mlik_original = m1,
    mlik_restart = m2,
    delta = m2 - m1,
    better = (m2 - m1) > .FB_RESTART_MLIK_MARGIN,
    collapsed = collapsed,
    # The restarted fit itself, so the caller can return it rather than the
    # solution it has just shown to be worse. Carried here and dropped by
    # the caller once adopted, so it never reaches the returned object twice.
    fit = fit2
  )
}


# The reason string a better restart contributes to the numerical-confirm
# gate. Kept separate so the wording is testable without running a fit.
#
# @noRd
.fb_restart_reason <- function(probe) {
  paste0(
    "a restart from a different starting mode found a better solution ",
    "(marginal log-likelihood ",
    format(round(probe$mlik_restart, 3), nsmall = 3),
    " against ",
    format(round(probe$mlik_original, 3), nsmall = 3),
    ", a difference of ",
    format(round(probe$delta, 3), nsmall = 3),
    "), so ",
    if (length(probe$collapsed) == 1L) {
      "the component "
    } else {
      paste0("the ", length(probe$collapsed), " components ")
    },
    "reported at zero variance are an artefact of where the optimiser ",
    "started, not a property of the data"
  )
}


# The note recorded when the restart is adopted.
#
# Adoption is a change to which answer the user gets, so it is stated rather
# than left to be inferred from a number that moved.
#
# @noRd
.fb_restart_adopted_note <- function(probe) {
  paste0(
    "the first optimisation converged to a mode with ",
    format(round(probe$delta, 1), nsmall = 1),
    " less marginal log-likelihood, pinning ",
    length(probe$collapsed),
    if (length(probe$collapsed) == 1L) " component" else " components",
    " at zero variance. A restart from a mode where those components carry ",
    "real variance found the better solution, and that is the fit returned. ",
    "The likelihood surface here has more than one converged mode, so treat ",
    "other quantities from this fit with the same care"
  )
}
