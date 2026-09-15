# =============================================================================
# Scale invariance: multiplying the response by a constant.
#
# Multiplying y by c is an exact transformation. Every variance scales by c^2,
# so a variance *ratio* and a *correlation* are invariant, and a fit that does
# not reproduce them is answering a different question from the one written.
#
# This file exists because of finding O04. `us(env):gen` on INLA returned a
# genetic correlation of 0.147 where ASReml returns 0.930, and the thing that
# made it undeniable was not the ASReml comparison -- it was multiplying the
# response by ten and watching the correlation move a long way. That check
# needs no licence, no second engine and no reference implementation, and a
# metamorphic check is cheaper than an oracle and sometimes finds more.
#
# The two structures that failed it -- `us(f):g` and the correlated random
# slope `(x | g)` -- were watched here by tripwires until 2026-09-13, when
# estimated correlation matrices moved to brms and INLA began refusing both.
# What remains is:
#
#   * a GATE on the structures the ASReml oracle has validated. These assert
#     invariance and must pass. If one starts failing, a fit that was right
#     has become wrong.
#
#   * the refusal of the two structures, at the original scale and rescaled,
#     so that neither scale can put them back on INLA.
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

  # Under the default uniform bound INLA's interval for `rep` (three levels)
  # runs past the bound, 0.17 to 3.65 against 3.05, and the support audit
  # refuses the fit at both scales (8 of 8 repeats each). A half-normal
  # prior on each standard deviation, with scale the standard deviation of
  # the modelled response, has no bound to cross and moves with the response
  # exactly as the default does, so the gate still measures invariance.
  .msi_hn_prior <- function(data) {
    s <- stats::sd(data$yield)
    fb_prior(
      sigma ~ half_normal(scale = s),
      sd(group = "rep") ~ half_normal(scale = s),
      sd(group = "rep:block") ~ half_normal(scale = s)
    )
  }
  fit_1 <- flexybayes(yield ~ gen, random = ~ rep + rep:block, data = d,
                      prior = .msi_hn_prior(d),
                      backend = "inla", verbose = FALSE)
  d10 <- d
  d10$yield <- 10 * d10$yield
  fit_10 <- flexybayes(yield ~ gen, random = ~ rep + rep:block, data = d10,
                       prior = .msi_hn_prior(d10),
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


# --- the two structures that failed invariance ------------------------ #

test_that("us() refuses on INLA at the original scale and rescaled (finding O04)", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")

  e <- new.env()
  utils::data("besag.met", package = "agridat", envir = e)
  d <- e$besag.met
  d$env <- factor(d$county)
  d$gen <- factor(d$gen)
  d <- d[!is.na(d$yield), ]
  d <- d[d$env %in% levels(d$env)[seq_len(2L)], ]
  d$env <- droplevels(d$env)
  d$gen <- droplevels(d$gen)

  for (k in c(1, 10)) {
    dk <- d
    dk$yield <- k * dk$yield
    expect_error(
      suppressMessages(flexybayes(
        yield ~ env, random = ~ us(env):gen, data = dk,
        backend = "inla", return_code = TRUE, verbose = FALSE
      )),
      class = "flexybayes_refusal_correlated_covariance_requires_brms",
      info = paste0("response multiplied by ", k)
    )
  }
})


test_that("(x | g) refuses on INLA at the original scale and rescaled (S01 D1)", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")

  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  d <- e$sleepstudy

  for (k in c(1, 10)) {
    dk <- d
    dk$Reaction <- k * dk$Reaction
    expect_error(
      suppressMessages(flexybayes(
        Reaction ~ Days + (Days | Subject), data = dk,
        backend = "inla", return_code = TRUE, verbose = FALSE
      )),
      class = "flexybayes_refusal_correlated_covariance_requires_brms",
      info = paste0("response multiplied by ", k)
    )
  }
})
