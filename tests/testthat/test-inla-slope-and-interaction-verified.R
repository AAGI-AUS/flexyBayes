# The verification that replaces the retired runtime artefact gate.
#
# Until 0.10.0.9000 two INLA mappings: the uncorrelated random slope
# `(x || g)` and the factor-by-numeric fixed interaction `y ~ f * x`, were
# refused at run time unless an on-disk artefact recorded that a
# "three-arbitrator" comparison had passed. One of those arbitrators was an
# engine withdrawn at 0.9.3, so the criterion could not be re-run and the
# refusal became permanent: a working capability blocked by the memory of a
# withdrawn dependency.
#
# The criterion's purpose was "the INLA mapping agrees with an independent
# reference". That is a test, not a runtime check, and this is it. lme4 fits
# the same models by REML and never sees this package's parser or emit, so it
# is a genuine oracle instead of a restatement of our own arithmetic.

.slope_fixture <- function(seed = 7L, n_g = 25L, n_per = 12L) {
  set.seed(seed)
  s1 <- n_g * n_per
  tab <- data.frame(
    g = factor(rep(seq_len(n_g), each = n_per)),
    x = rnorm(s1)
  )
  i1 <- as.integer(tab$g)
  b0 <- rnorm(n_g, 0, 1.0)
  b1 <- rnorm(n_g, 0, 0.5)
  tab$y <- 2 + 0.8 * tab$x + b0[i1] + b1[i1] * tab$x + rnorm(s1, 0, 0.3)
  tab
}

test_that("an uncorrelated random slope fits on INLA and matches lme4 REML", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .slope_fixture()
  fit <- suppressMessages(flexybayes(
    y ~ x + (x || g),
    data = tab,
    family = "gaussian",
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  m <- lme4::lmer(y ~ x + (1 | g) + (0 + x | g), data = tab, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(m))
  ref_int <- vc$sdcor[vc$grp == "g" & vc$var1 == "(Intercept)"][1]
  ref_slp <- vc$sdcor[vc$grp == "g.1" | (vc$grp == "g" & vc$var1 == "x")][1]
  ref_res <- vc$sdcor[vc$grp == "Residual"][1]

  # Posterior means of variance components sit a little above REML point
  # estimates because the posteriors are right-skewed, so the comparison is
  # on the ratio. A mistranslated mapping misses by far more than 25%.
  expect_equal(unname(got[["sigma"]]), unname(ref_res), tolerance = 0.15)
  expect_true(any(grepl("^sd_", names(got))))
  sds <- got[grepl("^sd_", names(got))]
  expect_equal(length(sds), 2L)
  # The generated index name must not reach the label. The first working
  # version reported sd_g_for_slope_x here.
  expect_false(any(grepl("_for_slope_", names(got), fixed = TRUE)))
  expect_true("sd_x_g" %in% names(got))
  expect_equal(
    unname(sort(sds)),
    unname(sort(c(ref_int, ref_slp))),
    tolerance = 0.25
  )
})

test_that("a factor-by-numeric fixed interaction fits on INLA and matches lm", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  set.seed(11L)
  s1 <- 300L
  tab <- data.frame(
    f = factor(sample(c("a", "b", "c"), s1, replace = TRUE)),
    x = rnorm(s1)
  )
  tab$y <- 1 + as.numeric(tab$f) * tab$x + rnorm(s1, 0, 0.3)
  tab$g <- factor(rep(seq_len(30L), length.out = s1))

  fit <- suppressMessages(flexybayes(
    y ~ f * x,
    random = ~g,
    data = tab,
    family = "gaussian",
    backend = "inla",
    verbose = FALSE
  ))
  cf <- fit$inla$summary.fixed
  ref <- coef(lm(y ~ f * x, data = tab))

  # The fixed-effect design is the same matrix in both, so with a weak prior
  # the posterior means should sit essentially on the least-squares
  # estimates. This is the check the retired artefact stood in for.
  for (nm in c("x", "fb:x", "fc:x")) {
    expect_true(nm %in% rownames(cf))
    expect_equal(unname(cf[nm, "mean"]), unname(ref[[nm]]), tolerance = 0.05)
  }
})

test_that("the retired artefact machinery is gone, not merely bypassed", {
  # A gate that no longer fires but still exists invites its own return.
  expect_false(exists(
    ".check_inla_verification_simple_slope_uncor",
    envir = asNamespace("flexyBayes"),
    inherits = FALSE
  ))
  expect_false(exists(
    ".lgm_check_factor_numeric_interaction_inla_verified",
    envir = asNamespace("flexyBayes"),
    inherits = FALSE
  ))
  expect_false(exists(
    ".fb_dev_verification_artefacts_enabled",
    envir = asNamespace("flexyBayes"),
    inherits = FALSE
  ))
  expect_identical(
    system.file("extdata", "inla-verification", package = "flexyBayes"),
    ""
  )
})
