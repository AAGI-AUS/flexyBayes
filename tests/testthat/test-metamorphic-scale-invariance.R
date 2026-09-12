# =============================================================================
# Scale invariance: multiplying the response by a constant.
#
# Multiplying y by c is an exact transformation. Every variance scales by c^2,
# so a variance *ratio* and a *correlation* are invariant, and a fit that does
# not reproduce them is answering a different question from the one written.
#
# This file exists because of finding O04. `us(env):gen` on INLA returns a
# genetic correlation of 0.147 where ASReml returns 0.930, and the thing that
# made it undeniable was not the ASReml comparison -- it was multiplying the
# response by ten and watching the correlation move a long way. The single
# value -0.622 that stood here is withdrawn: eight identical rescaled calls
# span -0.894 to +0.686 and none converged, so the rescaled arm has no one
# number to quote. What reproduces is the FAILURE of invariance, and that is
# what the tripwire below asserts. That check
# needs no licence, no second engine and no reference implementation. It took
# under a minute to write and it had never been written, which is the whole
# lesson: a metamorphic check is cheaper than an oracle and sometimes finds
# more.
#
# Two kinds of test here, and the difference matters.
#
#   * A GATE on the structures the ASReml oracle has validated. These assert
#     invariance and must pass. If one starts failing, a fit that was right
#     has become wrong.
#
#   * A TRIPWIRE on the structure that is known to fail it. It asserts the
#     defect, with its measurement and its finding named. A characterisation
#     test is not an endorsement: when `us()` is repaired this test fails, and
#     the failure message says so and says which record to close. The
#     alternative -- an expectation that the defect is absent, skipped until
#     someone fixes it -- is a silent skip, and a silent skip is how a known
#     defect becomes an unknown one.
# =============================================================================

.msi_ratio <- function(fit, numerator, denominator) {
  vc <- summary(fit)$varcomp
  val <- function(nm) {
    row <- vc[vc$component == nm, , drop = FALSE]
    if (!nrow(row)) {
      stop("no variance component named ", nm, call. = FALSE)
    }
    as.numeric(row$estimate[[1L]])
  }
  val(numerator) / val(denominator)
}


test_that("a variance ratio is invariant to rescaling the response (iid intercept)", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  set.seed(20260911L)
  n_g <- 20L
  n_r <- 5L
  d <- data.frame(
    g = factor(rep(paste0("g", seq_len(n_g)), each = n_r)),
    x = stats::rnorm(n_g * n_r)
  )
  u <- stats::rnorm(n_g, sd = 1.5)
  d$y <- 2 + 0.8 * d$x + u[as.integer(d$g)] + stats::rnorm(nrow(d), sd = 1)

  fit_1 <- flexybayes(y ~ x, random = ~ g, data = d,
                      backend = "inla", verbose = FALSE)
  d10 <- d
  d10$y <- 10 * d10$y
  fit_10 <- flexybayes(y ~ x, random = ~ g, data = d10,
                       backend = "inla", verbose = FALSE)

  r1 <- .msi_ratio(fit_1, "sd_g", "sigma")
  r10 <- .msi_ratio(fit_10, "sd_g", "sigma")

  # A tenfold rescaling of a well-identified Gaussian model: the ratio is the
  # same quantity, and the tolerance is the posterior's, not the algebra's.
  expect_equal(r10, r1, tolerance = 0.05)

  # And the slope scales exactly by the factor, which catches a rescaling
  # applied to one side of the model and not the other.
  b1 <- summary(fit_1)$fixed
  b10 <- summary(fit_10)$fixed
  expect_equal(
    b10$estimate[b10$term == "x"] / 10,
    b1$estimate[b1$term == "x"],
    tolerance = 0.05
  )
})


test_that("a variance ratio is invariant to rescaling the response (nested strata)", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  skip_on_cran()

  # The oracle's own substrate for the nested row (scenario 02), so the
  # metamorphic gate runs on exactly the data the ASReml comparison
  # validated -- BLUP agreement r = 1.0000 -- rather than on a simulation
  # chosen to be easy.
  e <- new.env()
  utils::data("john.alpha", package = "agridat", envir = e)
  d <- e$john.alpha
  d$gen <- factor(d$gen)
  d$rep <- factor(d$rep)
  d$block <- factor(d$block)

  fit_1 <- flexybayes(yield ~ gen, random = ~ rep + rep:block, data = d,
                      backend = "inla", verbose = FALSE)
  d10 <- d
  d10$yield <- 10 * d10$yield
  fit_10 <- flexybayes(yield ~ gen, random = ~ rep + rep:block, data = d10,
                       backend = "inla", verbose = FALSE)

  # Read `rep:block`, not `rep`. The choice is the oracle's own identification
  # rule (compare.R: ASReml se/estimate < 0.5), applied to a metamorphic gate.
  # `rep` has three levels here, so its variance component is weakly
  # identified and its ratio wanders at the rescaled mode -- measured over
  # three refits each: sd_rep/sigma 3.167, 3.168, 3.179 unscaled against
  # 3.249, 3.122, 2.737 rescaled, a 19 per cent spread that fails a 10 per
  # cent tolerance intermittently. `rep:block` has eighteen cells and is
  # reproducible to four digits at both scales: 0.8953, 0.8953, 0.8952
  # against 0.9058, 0.9077, 0.9060.
  #
  # The first version of this gate read `sd_rep` and was flaky. The fix is a
  # better-identified quantity, not a wider tolerance: a metamorphic gate is
  # only as reproducible as the component it reads, and SPEC.md is explicit
  # that there is no outcome in which a tolerance moves.
  expect_equal(
    .msi_ratio(fit_10, "sd_rep:block", "sigma"),
    .msi_ratio(fit_1, "sd_rep:block", "sigma"),
    tolerance = 0.10
  )
})


# --- the tripwire's own rule ------------------------------------------ #

# .o04_verdict() --- what the O04 tripwire concludes from one pair of fits.
#
# Extracted so it can be exercised on values a real fit will not produce on
# demand. A disjunction asserted inline is not falsifiable by breaking the
# implementation, because the surviving branch carries it; the failure mode
# that matters here is the one where BOTH branches go quiet at once, and
# that is only reachable by supplying the values.
.o04_verdict <- function(ok_10, c_1, c_10, tol = 0.05) {
  if (is.na(ok_10)) {
    return("unverified")
  }
  if (isFALSE(ok_10)) {
    return("nonconverged")
  }
  if (abs(c_10 - c_1) > tol) {
    return("not_invariant")
  }
  "repaired"
}


test_that("us() on INLA fails scale invariance -- tripwire for finding O04", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  skip_on_cran()
  skip_on_ci()

  e <- new.env()
  utils::data("besag.met", package = "agridat", envir = e)
  d <- e$besag.met
  d$env <- factor(d$county)
  d$gen <- factor(d$gen)
  d <- d[!is.na(d$yield), ]
  d <- d[d$env %in% levels(d$env)[seq_len(2L)], ]
  d$env <- droplevels(d$env)
  d$gen <- droplevels(d$gen)

  fit_of <- function(dat) {
    flexybayes(yield ~ env, random = ~ us(env):gen, data = dat,
               backend = "inla", verbose = FALSE)
  }
  correlation_of <- function(fit) {
    vc <- summary(fit)$varcomp
    row <- vc[grepl("^cor_", vc$component), , drop = FALSE]
    if (!nrow(row)) {
      stop("no correlation component in the summary", call. = FALSE)
    }
    as.numeric(row$estimate[[1L]])
  }

  # What is asserted here, and what deliberately is not.
  #
  # Measured 2026-09-11 across identical calls on these rows:
  #
  #   unscaled,  n = 4   0.1448 to 0.1474          -- tight
  #   rescaled x10, n = 8  -0.894 to +0.686, 0/8 converged
  #
  # So the two halves of this defect have to be asserted differently. The
  # unscaled correlation is reproducible and is nowhere near ASReml's 0.930,
  # and that is the finding a repair would move. The rescaled fit is not a
  # wrong number -- it is not a converged fit at all, and the correlation it
  # returns spans nearly the whole of [-1, 1] across identical calls.
  #
  # Two earlier versions of this test asserted a numeric property of the
  # rescaled fit -- first a gap above 0.3, then a gap above the oracle's
  # declared 0.10 correlation tolerance -- and both were flaky, because no
  # numeric property of a non-convergent fit is assertable. That is the same
  # error corrected elsewhere this session: quoting a digit a refit
  # contradicts. What is stable about the rescaled arm is categorical, so
  # that is what is asserted.
  fit_1 <- fit_of(d)
  c_1 <- correlation_of(fit_1)

  d10 <- d
  d10$yield <- 10 * d10$yield
  fit_10 <- fit_of(d10)

  # The stable half, and the one a repair would move.
  expect_true(
    c_1 < 0.5,
    info = paste0(
      "the unscaled genetic correlation is ", signif(c_1, 3),
      " against ASReml's 0.930 on these rows. If it has risen toward 0.930 ",
      "this is a repair: close finding O04, re-run oracle scenarios ",
      "06_us_dim2 and 07_us_dim3, reclassify the INLA us() row in ",
      "R/translation_registry.R, update the capability-matrix note and ",
      "KNOWN_ISSUES open issue 1, and replace this tripwire with the ",
      "invariance gate the two tests above use."
    )
  )

  # The rescaled half. Rescaling the response is an exact transformation and
  # must not decide whether the model fits at all.
  # This assertion used to read isTRUE(summary(fit_10)$converge$ok) against a
  # field that did not exist. isTRUE(NULL) is FALSE, so it passed for every
  # fit -- clean, failed, and diagnostic-free were indistinguishable to it
  # (audit A25). The field is real now and tri-state.
  #
  # The rule is a disjunction, because the defect is a failure of SCALE
  # INVARIANCE and not a categorical claim that the rescaled fit never
  # converges: during the 11 September audit one rescaled fit reported mode
  # status 0 with a passing confirmation while moving the correlation from
  # about 0.135 to -0.780.
  #
  # A disjunction cannot be falsified by breaking one branch -- the other
  # simply carries it -- so the rule lives in .o04_verdict() and is
  # falsified directly in the test below this one.
  verdict <- .o04_verdict(
    ok_10 = summary(fit_10)$converge$ok,
    c_1 = c_1,
    c_10 = correlation_of(fit_10)
  )

  expect_false(
    identical(verdict, "unverified"),
    info = paste0(
      "the rescaled us() fit reports an UNVERIFIED convergence state, which ",
      "is neither outcome this tripwire distinguishes. Read the diagnostics ",
      "before going further; NA must not pass as though it were a result."
    )
  )
  expect_false(
    identical(verdict, "repaired"),
    info = paste0(
      "the tenfold-rescaled us() fit both converged AND reproduced the ",
      "unscaled correlation. Rescaling the response is exact, so that is ",
      "what a correct fit does: finding O04 may be repaired. Re-run oracle ",
      "scenarios 06_us_dim2 and 07_us_dim3, reclassify the INLA us() row, ",
      "and replace this tripwire with the invariance gate above."
    )
  )
})



test_that("the O04 tripwire rule can reach all four of its verdicts", {
  # Non-convergence, whatever the correlations did.
  expect_identical(.o04_verdict(FALSE, 0.135, 0.135), "nonconverged")
  expect_identical(.o04_verdict(FALSE, 0.135, -0.780), "nonconverged")

  # Converged but not invariant: the audit's own counterexample.
  expect_identical(.o04_verdict(TRUE, 0.135, -0.780), "not_invariant")

  # Converged AND invariant. This is the state that must make the gate red,
  # because it is what a repaired us() looks like and the tripwire's job is
  # to notice a repair rather than to keep asserting a defect that is gone.
  expect_identical(.o04_verdict(TRUE, 0.135, 0.135), "repaired")
  expect_identical(.o04_verdict(TRUE, 0.930, 0.9315), "repaired")

  # Unverified never passes as a result.
  expect_identical(.o04_verdict(NA, 0.135, -0.780), "unverified")
})


test_that("the O04 rule is sensitive at its stated tolerance", {
  # Just inside and just outside 0.05, so the boundary is asserted rather
  # than assumed.
  expect_identical(.o04_verdict(TRUE, 0.100, 0.149), "repaired")
  expect_identical(.o04_verdict(TRUE, 0.100, 0.151), "not_invariant")
})
