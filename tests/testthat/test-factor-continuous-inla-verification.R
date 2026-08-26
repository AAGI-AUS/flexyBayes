# INLA mapping verification test for the factor_numeric_interaction
# term class (ADR 0019 §Decision 5 / the design spec three-arbitrator
# decision rule).
#
# Policy: the INLA mapper for the new term class registers IF AND
# ONLY IF this verification test demonstrates faithful posterior
# agreement with both Arbitrator B and lme4::lmer on a gaussian-
# identity fixture, per the v2 benchmark agreement contract. Failure
# to verify (INLA not installed, INLA verification mismatch, INLA fit
# error) means the on-disk verification artefact stays absent or
# carries `pass = FALSE`, and the lgm_gate() check
# .lgm_check_factor_numeric_interaction_inla_verified() refuses
# INLA dispatch on the factor:continuous indexed-slope term class
# with a deferral message.
#
# 0.9.3: Arbitrator B was the native engine withdrawn entirely this
# release (see NEWS.md). The verification below cannot currently be
# re-run as designed -- rebuilding the three-arbitrator criterion
# around the two active engines (e.g. substituting brms for
# Arbitrator B) is a methodology decision for whoever takes up the
# rebuild, not a mechanical substitution this withdrawal makes on its
# own authority. The main test therefore records the historical
# methodology (kept verbatim below, with the withdrawn engine's name
# generalised to "Arbitrator B") and skips rather than attempting a
# fit through a backend that raises unknown_backend. The
# gate's refusal (`.lgm_check_factor_numeric_interaction_inla_verified()`
# returning `pass = FALSE`) is exercised directly in
# test-factor-continuous-emit.R's Subtest 10 and in the two companion
# tests below, both of which are unaffected by the withdrawal.
#
# Tier: Tier 2 (full devtools::test()). Skipped on CRAN and while
# Arbitrator B has no active-engine implementation.

# Where the host-local verification artefact may be written. An installed
# library and the R CMD check directory are both out of bounds -- CRAN
# policy allows a package's tests the session temporary directory and
# nothing else -- so the source tree is used only when the suite is
# actually running from it, which is what a real DESCRIPTION and inst/
# two levels up mean. Every other host writes under tempdir(), so the
# developer rehearsal hook still works for the length of the session and
# leaves nothing behind after it.
.fnia_artefact_dir <- function() {
  if (
    file.exists(file.path("..", "..", "DESCRIPTION")) &&
      dir.exists(file.path("..", "..", "inst"))
  ) {
    return(file.path("..", "..", "inst", "extdata", "inla-verification"))
  }
  file.path(tempdir(), "flexybayes-inla-verification")
}

# Where to read it back from. Reading an installed copy is fine, so that
# is tried first; otherwise the artefact is wherever this session wrote
# it.
.fnia_artefact_path <- function() {
  installed <- system.file(
    "extdata",
    "inla-verification",
    "factor_numeric_interaction.rds",
    package = "flexyBayes",
    mustWork = FALSE
  )
  if (nzchar(installed)) {
    return(installed)
  }
  file.path(.fnia_artefact_dir(), "factor_numeric_interaction.rds")
}

test_that("INLA mapping for factor_numeric_interaction passes 3-arbitrator gate", {
  skip(paste0(
    "Arbitrator B (the native engine withdrawn entirely in 0.9.3, see ",
    "NEWS.md) has no active-engine implementation; the three-arbitrator ",
    "criterion cannot be re-run until it is rebuilt around the two ",
    "active engines. See the design-record comment below for the ",
    "fixture and methodology a rebuild would start from."
  ))

  # -- Design record (not executed): the pre-0.9.3 fixture + methodology --
  #
  # Fixture: gaussian-identity, 3-level factor x continuous interaction,
  # n = 240, seed 20260523L, beta_true = ((Intercept)=0, fb=0.3, fc=-0.2,
  # x=0.8, fb:x=0.5, fc:x=-0.4). Tolerance per coefficient:
  # 0.10 * max(0.05, |beta_true|) + 0.05 (design spec §3.4).
  #
  # Arbitrator A -- lme4::lm (gaussian-identity reduces to lm on a
  # fixed-effects-only model); REML SEs serve as the tolerance unit.
  #
  # Arbitrator B -- the withdrawn engine's indexed emit (the v0.2.6
  # "all-levels + global intercept" parameterisation eta = mu_atg +
  # tau_f[f_id] + beta_x * x + slope_dev_f_x[f_id] * x), converted to
  # treatment-coded canonical slots via (Intercept) = mu_atg +
  # tau_f[1,1]; fb = tau_f[2,1] - tau_f[1,1]; fc = tau_f[3,1] -
  # tau_f[1,1]; x = beta_x; fb:x = slope_dev_f_x_raw[1,1]; fc:x =
  # slope_dev_f_x_raw[2,1] -- taking the mean of the linear combination
  # on the draws, not the mean of each block separately, since the two
  # MCMC parameter blocks are not individually identified. A rebuilt
  # Arbitrator B on brms would read canonical coefficients directly via
  # coef(fit) (see test-factor-continuous-emit.R's brms rewrite for the
  # pattern) rather than reconstructing them from raw parameter names.
  #
  # Arbitrator C -- INLA via the native `f:x` formula syntax, INLA
  # default priors, posterior mean of summary.fixed.
  #
  # Verdict: INLA's posterior mean must agree with BOTH Arbitrator B and
  # lme4 within tolerance for every coefficient. The verification
  # artefact (list(timestamp, R_version, INLA_version, flexyBayes_v,
  # fixture, agreement, pass), written to .fnia_artefact_dir()) records
  # the outcome for `.lgm_check_factor_numeric_interaction_inla_verified()`
  # to read when a developer opts in via
  # options(flexyBayes.dev_inla_verification_artefacts = TRUE).
})

# Companion check: the artefact -> gate pickup, in isolation from the
# INLA fit so it is fast.
#
# Since 0.9.0 the pickup is a developer rehearsal hook and nothing more.
# The artefact directory is excluded from the build and the option below
# is off on every shipped surface, so an artefact that reports
# `pass = TRUE` lifts the gate only for a developer who asked for it by
# hand. The default-off half of that contract is asserted in
# test-inla-verification-artefact-policy.R, which is where the P0-2
# reproducibility policy lives.

test_that("the artefact lifts the gate only under the developer option", {
  art_path <- .fnia_artefact_path()
  if (!file.exists(art_path)) {
    testthat::skip("verification artefact not yet produced on this host")
  }

  rec <- readRDS(art_path)
  if (!isTRUE(rec$pass)) {
    testthat::skip(
      "verification artefact records pass = FALSE; gate refuses INLA"
    )
  }

  withr::local_options(
    list(flexyBayes.dev_inla_verification_artefacts = TRUE)
  )
  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    fixed_terms = list(
      list(
        type = "factor",
        var = "f",
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f"
      ),
      list(type = "continuous", var = "x", label = "x"),
      list(
        type = "factor_numeric_interaction",
        factor = "f",
        continuous = "x",
        vars = c("f", "x"),
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f:x"
      )
    ),
    random_terms = list(),
    residual_terms = list(list(type = "units")),
    priors = list(legacy = TRUE),
    source = "brms"
  )
  r10 <- flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(fb)
  expect_true(r10$pass)
})

# Trivial-pass guard: on IRs without a factor_numeric_interaction
# term the §3.4 gate is a no-op (the verification artefact's
# existence is irrelevant when the term class is absent).

test_that("verification gate is no-op without factor_numeric_interaction term", {
  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    fixed_terms = list(list(type = "continuous", var = "x", label = "x")),
    random_terms = list(),
    residual_terms = list(list(type = "units")),
    priors = list(legacy = TRUE),
    source = "brms"
  )
  r <- flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(fb)
  expect_true(r$pass)
})
