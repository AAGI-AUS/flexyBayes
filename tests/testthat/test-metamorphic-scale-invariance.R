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
# response by ten and watching the correlation move to -0.622. That check
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
  expect_false(
    isTRUE(summary(fit_10)$converge$ok),
    info = paste0(
      "the tenfold-rescaled us() fit now reports convergence. Measured ",
      "0 of 8 identical calls converging on 2026-09-11, so this is a ",
      "change worth reading rather than a flake -- see the note above."
    )
  )
})
