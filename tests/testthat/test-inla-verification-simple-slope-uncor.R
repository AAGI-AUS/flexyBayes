# ADR 0020 §Decision 5 -- three-arbitrator INLA verification gate
# for uncorrelated random slopes (x || g).
#
# Per the design spec "INLA mapping verification policy" + ADR 0020
# §Decision 5: the INLA mapping is registered if and only if the
# three-arbitrator verification test passes. The arbitrators were
# INLA, a native engine since withdrawn entirely (see NEWS.md, 0.9.3),
# and lme4 (this is a gaussian-identity term, so rule 1 of §3.4
# applies).
#
# Verification criterion (as designed, no longer runnable -- see
# below): simulate J = 20 groups with a known intercept variance and
# a known slope variance; fit flexybayes(..., backend = "inla"), the
# withdrawn engine's equivalent, and lme4::lmer(... REML = FALSE).
# Assert that INLA's posterior for both sd_<g> and sd_<x>_<g> matches
# both peers within W_1 <= 0.20 * tau_true on a 4000-draw chain.
#
# Outcome handling. If verification passes, we write
# inst/extdata/inla-verification/simple_slope_uncor.rds with
# pass = TRUE; emit_inla() then admits the (x || g) term. If
# verification fails, we leave the artefact absent (or write
# pass = FALSE) and emit_inla() refuses with the
# flexybayes_inla_simple_slope_uncor_deferred condition.
#
# v0.2.6 ship state. The verification artefact is NOT generated at
# ship time; the artefact was intended to be generated only by an
# explicit local rehearsal. For the v0.2.6 release the INLA mapping
# refuses; (x || g) fits via brms (backend = "auto" already falls
# back there automatically).
#
# 0.9.3 withdrawal (see NEWS.md): the native engine that was one of
# the three arbitrators is gone entirely, not quarantined -- unlike
# the earlier quarantine framing, re-entry (should it ever be
# proposed) is a fresh implementation, not a repair of retained code,
# so there is no rehearsal design left to keep "on file for when the
# arbitrator set is revisited". The three-arbitrator rehearsal test
# that used to live here (fitting the withdrawn engine's arm) is
# deleted rather than kept unreachable; the criterion would need a
# fresh two-arbitrator (or re-admitted-engine) design before any
# INLA mapping for this term class can be registered. Tests 1 and 3
# below are unaffected -- they pin the structural refusal/deferral
# behaviour, which does not depend on which engines the criterion
# names.
#
# This test file is skipped via skip_if_not_installed("INLA");
# acceptance is structural (the file exists, and the artefact state
# is what emit_inla() consults).

mk_inla_verification_fixture <- function(
  seed = 20260523L,
  J = 20L,
  n_per = 12L,
  beta = 0.5,
  sd_int = 25,
  sd_slope = 4,
  sigma_e = 20
) {
  set.seed(seed)
  g <- factor(rep(seq_len(J), each = n_per))
  x <- rep(seq_len(n_per) - 1L, times = J)
  u_int <- rnorm(J, sd = sd_int)
  u_slope <- rnorm(J, sd = sd_slope)
  y <- 250 +
    beta * x +
    u_int[as.integer(g)] +
    u_slope[as.integer(g)] * x +
    rnorm(length(x), sd = sigma_e)
  list(
    data = data.frame(y = y, x = x, g = g),
    tau_int = sd_int,
    tau_slope = sd_slope,
    sigma_e = sigma_e
  )
}


# ---------------------------------------------------------------- #
# (1) Artefact existence determines emit_inla() admission           #
# ---------------------------------------------------------------- #

test_that("emit_inla() consults the verification artefact for (x || g) admission", {
  # This subtest is run regardless of INLA installation; it asserts
  # the policy contract directly. The .check_inla_verification_*
  # helper is the gate -- it consults system.file() for the
  # artefact path. If the artefact is absent (the v0.2.6 ship
  # state), the gate refuses with the structured condition.
  artefact_path <- system.file(
    "extdata",
    "inla-verification",
    "simple_slope_uncor.rds",
    package = "flexyBayes"
  )
  if (!nzchar(artefact_path) || !file.exists(artefact_path)) {
    err <- tryCatch(
      flexyBayes:::.check_inla_verification_simple_slope_uncor(),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_inla_simple_slope_uncor_deferred")
    expect_identical(err$deferral_target, "a future release")
    expect_identical(err$workaround, "backend = \"brms\"")
  } else {
    art <- readRDS(artefact_path)
    if (isTRUE(art$pass)) {
      expect_invisible(
        flexyBayes:::.check_inla_verification_simple_slope_uncor()
      )
    } else {
      err <- tryCatch(
        flexyBayes:::.check_inla_verification_simple_slope_uncor(),
        error = function(e) e
      )
      expect_s3_class(err, "flexybayes_inla_simple_slope_uncor_deferred")
    }
  }
})


# ---------------------------------------------------------------- #
# (2) Refusal message names the brms workaround                     #
# ---------------------------------------------------------------- #

test_that("INLA verification refusal points users to brms as the workaround", {
  err <- tryCatch(
    flexyBayes:::.check_inla_verification_simple_slope_uncor(),
    error = function(e) e
  )
  # Skip cleanly if the artefact happens to be present + passing
  # (e.g., a future local rehearsal generated it). The ship-state
  # behaviour is the failing branch.
  if (!inherits(err, "flexybayes_inla_simple_slope_uncor_deferred")) {
    testthat::skip(
      "INLA verification artefact present + pass -- nothing to refuse"
    )
  }
  msg <- conditionMessage(err)
  expect_true(grepl("backend = \"brms\"", msg, fixed = TRUE))
  expect_true(grepl("three-arbitrator verification test", msg, fixed = TRUE))
  expect_true(grepl("future release", msg))
})
