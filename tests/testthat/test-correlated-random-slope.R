# Correlated random slope on INLA and brms: (x | g).
#
# The model is the uncorrelated slope with the correlation freed -- one
# intercept SD, one slope SD, and the correlation between them. It used to
# refuse at ingest, before any engine was chosen, on the ground that no
# active backend represented it. INLA does: the two blocks of deviations
# are ONE iid2d field over an index laid out intercepts-then-slopes, and
# the slope block enters the linear predictor as a fixed-scaling copy of
# that same field weighted by x.
#
# The correlation is the whole point of the class, so the tests do not stop
# at "a correlation came back". They check that a NEGATIVE correlation
# comes back negative, and they check the whole triple against lme4's REML
# fit of the same model -- an oracle this package did not write.

.cs_fixture <- function(seed = 2026L, n_g = 60L, n_i = 12L,
                        sd_int = 1.2, sd_slope = 0.6, rho = 0.6,
                        sd_resid = 0.35, beta_x = 1.5) {
  set.seed(seed)
  cov_is <- rho * sd_int * sd_slope
  mat_s <- matrix(c(sd_int^2, cov_is, cov_is, sd_slope^2), 2L, 2L)
  mat_u <- matrix(rnorm(n_g * 2L), ncol = 2L) %*% chol(mat_s)
  tab <- expand.grid(obs = seq_len(n_i), g = factor(seq_len(n_g)))
  tab$x <- runif(nrow(tab), -2, 2)
  gi <- as.integer(tab$g)
  tab$y <- 5 + beta_x * tab$x + mat_u[gi, 1L] + mat_u[gi, 2L] * tab$x +
    rnorm(nrow(tab), 0, sd_resid)
  tab
}

test_that("the gate accepts (x | g) for INLA", {
  tab <- .cs_fixture(n_g = 12L, n_i = 4L)
  fb <- fb_from_brms(y ~ x + (x | g), data = tab, family = "gaussian")
  expect_true("lgm_compatible" %in% lgm_gate(fb)$capabilities)
})

test_that("the emitted formula is one iid2d field plus a copy", {
  tab <- .cs_fixture(n_g = 12L, n_i = 4L)
  fb <- fb_from_brms(y ~ x + (x | g), data = tab, family = "gaussian")
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                collapse = " "))
  # One field of length 2 * n_g: intercepts then slopes.
  expect_match(txt, "f(g_cor_x_id, model = \"iid2d\", n = 24)", fixed = TRUE)
  # The slope block is the SAME field, weighted by x -- a copy, not a
  # second set of deviations.
  expect_match(txt, "f(g_cor_x_slope_id, x, copy = \"g_cor_x_id\")",
               fixed = TRUE)
})

test_that("the two index names are distinct and derived from the term", {
  # The iid2d field spans 1..n_g for the intercepts and n_g+1..2n_g for the
  # slopes, so the two blocks need two index columns whose names cannot
  # collide with each other or with a user column. That the BLOCKING is
  # right -- group i's intercept paired with group i's slope, not some
  # other group's -- is not visible in a name or a formula; it is what the
  # negative-correlation test below falsifies, since a mispaired field
  # returns a correlation near zero whatever the truth is.
  expect_identical(.fb_inla_cor_slope_index("g", "x"), "g_cor_x_id")
  expect_identical(.fb_inla_cor_slope_slope_index("g", "x"),
                   "g_cor_x_slope_id")
  # A group or slope name R would not accept as a symbol is made safe.
  expect_identical(.fb_inla_cor_slope_index("site name", "x"),
                   "site.name_cor_x_id")
})

test_that("a generated index name colliding with a user column refuses", {
  tab <- .cs_fixture(n_g = 10L, n_i = 3L)
  tab$g_cor_x_id <- 1L
  fb <- fb_from_brms(y ~ x + (x | g), data = tab, family = "gaussian")
  expect_error(
    emit_inla(fb, data = tab, verbose = FALSE),
    "already has one"
  )
})

test_that("INLA recovers both standard deviations and the correlation", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .cs_fixture()
  fit <- suppressMessages(flexybayes(
    y ~ x + (x | g),
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  expect_equal(as.numeric(fit$inla$mode$mode.status), 0)

  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)
  expect_equal(unname(got[["sd_g"]]), 1.2, tolerance = 0.20)
  expect_equal(unname(got[["sd_x_g"]]), 0.6, tolerance = 0.20)
  expect_equal(unname(got[["cor_g__Intercept__x"]]), 0.6, tolerance = 0.20)
  expect_equal(unname(got[["sigma"]]), 0.35, tolerance = 0.10)

  # The interval on the correlation covers the truth and excludes zero --
  # this class exists to estimate that number, so a wide interval centred
  # anywhere would not be a recovery.
  row <- v[v$component == "cor_g__Intercept__x", ]
  expect_lt(row$conf.low, 0.6)
  expect_gt(row$conf.high, 0.6)
  expect_gt(row$conf.low, 0)
})

test_that("a negative correlation comes back negative", {
  # Crossover between the intercept and the slope. A positive-correlation
  # fixture cannot distinguish a working sign from a structure that only
  # ever returns positive numbers.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .cs_fixture(seed = 77L, rho = -0.7)
  fit <- suppressMessages(flexybayes(
    y ~ x + (x | g),
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  row <- v[v$component == "cor_g__Intercept__x", ]
  expect_lt(row$estimate, 0)
  expect_equal(row$estimate, -0.7, tolerance = 0.25)
  expect_lt(row$conf.high, 0)
})

test_that("INLA agrees with lme4's REML fit of the same model", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .cs_fixture()
  fit <- suppressMessages(flexybayes(
    y ~ x + (x | g),
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  ref <- lme4::lmer(y ~ x + (x | g), data = tab)
  vc <- as.data.frame(lme4::VarCorr(ref))
  ref_sd_int <- vc$sdcor[vc$grp == "g" & vc$var1 == "(Intercept)" &
                           is.na(vc$var2)]
  ref_sd_slope <- vc$sdcor[vc$grp == "g" & vc$var1 == "x" & is.na(vc$var2)]
  ref_cor <- vc$sdcor[!is.na(vc$var2)]

  got <- stats::setNames(summary(fit)$varcomp$estimate,
                         summary(fit)$varcomp$component)
  expect_equal(unname(got[["sd_g"]]), ref_sd_int, tolerance = 0.15)
  expect_equal(unname(got[["sd_x_g"]]), ref_sd_slope, tolerance = 0.15)
  expect_equal(unname(got[["cor_g__Intercept__x"]]), ref_cor,
               tolerance = 0.15)
})

test_that("no generated index name reaches a summary label", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .cs_fixture(n_g = 25L, n_i = 6L)
  fit <- suppressMessages(flexybayes(
    y ~ x + (x | g),
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  comps <- summary(fit)$varcomp$component
  expect_true(all(c("sd_g", "sd_x_g", "cor_g__Intercept__x") %in% comps))
  expect_false(any(grepl("_cor_x_id|Rho1:2|component 1", comps)))
})

test_that("an explicit SD prior on (x | g) refuses rather than being dropped", {
  # iid2d's hyperparameters are a single joint Wishart on theta1, with
  # theta2 and theta3 carrying prior = "none". There is no
  # standard-deviation slot, and splicing one in makes INLA drop it with a
  # warning -- the silent prior loss this package fixed in the sectioned
  # residual.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .cs_fixture(n_g = 15L, n_i = 4L)
  expect_error(
    suppressMessages(flexybayes(
      y ~ x + (x | g),
      data = tab,
      prior = fb_prior(sd(group = "g") ~ pc(upper = 3, prob = 0.05)),
      backend = "inla",
      verbose = FALSE
    )),
    "Wishart"
  )
})

test_that("the reported prior is the one that reached the fit", {
  # The package's own SD default is DECLARED for this group, and INLA has
  # nowhere to put it, so the fit runs under INLA's Wishart. Reporting the
  # declared default in the prior column would say the fit used a prior it
  # never saw.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .cs_fixture(n_g = 20L, n_i = 5L)
  fit <- suppressMessages(flexybayes(
    y ~ x + (x | g),
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  cell <- function(cmp) v$prior[v$component == cmp]
  expect_identical(cell("sd_g"), "engine default")
  expect_identical(cell("sd_x_g"), "engine default")
  # The residual DOES take the declared prior, so it must still say so --
  # otherwise the fix is just a blanket relabelling.
  expect_false(identical(cell("sigma"), "engine default"))
  expect_true(all(c("sd_g", "sd_x_g") %in% names(prior_summary(fit)$not_applied)))
})
