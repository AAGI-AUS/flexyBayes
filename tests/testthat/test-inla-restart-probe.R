# =============================================================================
# The restart probe.
#
# `mode.status == 0` and a finite marginal likelihood do not make a fit
# trustworthy. On `agridat::besag.met` with `random = ~ diag(env):gen`, 24
# fits across four prior settings returned no genetic variance in any
# environment on 19 of them and substantial variance on 5, and all 24
# passed that check. ASReml on
# the same rows returns 5.84 and 5.45 in two of six environments, stably.
#
# Restarting from a mode where the pinned components carry real variance
# finds a solution 65 to 154 marginal-log-likelihood units better.
#
# These tests pin the detector, its threshold, and that a healthy fit is
# neither flagged nor charged for the extra fit.
# =============================================================================

fake_hyper <- function(resid_prec, comp_precs, comp_low = NULL) {
  nm <- c(
    "Precision for the Gaussian observations",
    paste0("Precision for f", seq_along(comp_precs))
  )
  low <- c(
    resid_prec * 0.9,
    if (is.null(comp_low)) comp_precs * 0.9 else comp_low
  )
  m <- data.frame(
    mean = c(resid_prec, comp_precs),
    `0.025quant` = low,
    check.names = FALSE
  )
  rownames(m) <- nm
  m
}

test_that("a component pinned at zero variance is detected from the raw table", {
  # A component whose variance is below 1e-4 of the residual's (the 1% rule
  # on the SD scale, squared) counts as pinned.
  fit <- list(
    summary.hyperpar = fake_hyper(
      resid_prec = 0.0025, # residual variance 400
      comp_precs = c(1e10, 0.01) # first is pinned, second is healthy
    )
  )
  hit <- .fb_inla_collapsed_hyper(fit)
  expect_length(hit, 1L)
  expect_equal(hit, "Precision for f1")
})

test_that("a healthy fit yields no pinned components", {
  fit <- list(
    summary.hyperpar = fake_hyper(
      resid_prec = 0.0025,
      comp_precs = c(0.01, 0.02, 0.005)
    )
  )
  expect_equal(.fb_inla_collapsed_hyper(fit), character(0))
})

test_that("the detector reads the low end of the interval, not the mean", {
  # A component whose mean precision is enormous but whose interval reaches
  # a plausible variance is not pinned: the fit is uncertain, not collapsed.
  fit <- list(
    summary.hyperpar = fake_hyper(
      resid_prec = 0.0025,
      comp_precs = c(1e10),
      comp_low = c(0.01)
    )
  )
  expect_equal(.fb_inla_collapsed_hyper(fit), character(0))
})

test_that("the probe does not run when nothing is pinned", {
  # The cost contract: a fit with no collapsed component pays for no extra
  # fit at all.
  expect_null(.fb_inla_restart_probe(list(), list(), character(0)))
})

test_that("a better restart is judged better only beyond the declared margin", {
  probe_at <- function(delta) {
    list(
      mlik_original = -100,
      mlik_restart = -100 + delta,
      delta = delta,
      collapsed = "Precision for f1"
    )
  }
  # The margin is two log units, the conventional positive-evidence step.
  expect_gt(.FB_RESTART_MLIK_MARGIN, 0)
  expect_equal(.FB_RESTART_MLIK_MARGIN, 2)
  txt <- .fb_restart_reason(probe_at(154.02))
  expect_match(txt, "restart from a different starting mode", fixed = TRUE)
  expect_match(txt, "154.02", fixed = TRUE)
  expect_match(txt, "not a property of the data", fixed = TRUE)
})

test_that("the reason names how many components were affected", {
  one <- .fb_restart_reason(list(
    mlik_original = -1,
    mlik_restart = -1,
    delta = 5,
    collapsed = "Precision for f1"
  ))
  many <- .fb_restart_reason(list(
    mlik_original = -1,
    mlik_restart = -1,
    delta = 5,
    collapsed = paste0("Precision for f", 1:6)
  ))
  expect_match(one, "the component ", fixed = TRUE)
  expect_match(many, "the 6 components ", fixed = TRUE)
})

test_that("a live fit that needed a restart returns the restarted fit", {
  skip_on_cran()
  # No skip_on_ci(): where INLA is installed, this live path runs.
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  e <- new.env()
  utils::data("besag.met", package = "agridat", envir = e)
  d <- e$besag.met
  d$env <- factor(d$county)
  d$gen <- factor(d$gen)
  d <- d[!is.na(d$yield), ]
  d <- d[order(d$env, d$gen), ]
  rownames(d) <- NULL
  fit <- suppressMessages(suppressWarnings(flexybayes(
    yield ~ env,
    random = ~ diag(env):gen,
    data = d,
    backend = "inla",
    verbose = FALSE
  )))

  probe <- fit$num_check$restart_probe
  v <- summary(fit)$varcomp
  sigma <- v$estimate[v$component == "sigma"]
  gen <- v$estimate[startsWith(v$component, "sd_gen__")]
  pinned <- sum(gen < 0.01 * sigma)

  # The structure is bistable, so this asserts the invariant, not a
  # particular outcome. Whichever mode the first optimisation reached, the
  # fit that comes back must be the better one and must say whether it was
  # swapped.
  if (isTRUE(fit$num_check$restart_adopted)) {
    expect_false(is.null(probe))
    expect_true(isTRUE(probe$better))
    expect_gt(probe$delta, .FB_RESTART_MLIK_MARGIN)
    expect_true(nzchar(fit$num_check$restart_adopted_note))
    # The point of adopting: the components the first mode pinned at zero
    # are not pinned in the fit the user receives.
    expect_lt(pinned, length(gen))
  } else {
    # No adoption means the first mode was not beaten, so nothing was
    # silently swapped and the probe either did not run or found no better
    # solution.
    expect_true(is.null(probe) || !isTRUE(probe$better))
  }
})

test_that("the adopted fit is the restarted fit, not the original", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  e <- new.env()
  utils::data("besag.met", package = "agridat", envir = e)
  d <- e$besag.met
  d$env <- factor(d$county)
  d$gen <- factor(d$gen)
  d <- d[!is.na(d$yield), ]
  d <- d[order(d$env, d$gen), ]
  rownames(d) <- NULL
  fit <- suppressMessages(suppressWarnings(flexybayes(
    yield ~ env,
    random = ~ diag(env):gen,
    data = d,
    backend = "inla",
    verbose = FALSE
  )))
  probe <- fit$num_check$restart_probe
  skip_if(
    !isTRUE(fit$num_check$restart_adopted),
    "this run's first optimisation was not beaten"
  )
  # The returned object must carry the restart's marginal likelihood, not
  # the one it was shown to be worse than. Without this the adoption could
  # record itself while still handing back the original fit.
  expect_equal(fit$inla$mlik[1, 1], probe$mlik_restart, tolerance = 1e-8)
  expect_gt(fit$inla$mlik[1, 1], probe$mlik_original)
  # And the probe record must not smuggle a second copy of the fit into the
  # returned object.
  expect_null(probe$fit)
})
