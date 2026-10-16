# Correlated random slope: (x | g).
#
# The model is the uncorrelated slope with the correlation freed, one
# intercept SD, one slope SD, and the correlation between them. It is fitted
# on brms, which writes `(x | g)` verbatim.
#
# INLA refuses it by name. The INLA lowering is still in the emitter, the
# two blocks of deviations as ONE iid2d field over an index laid out
# intercepts-then-slopes, the slope block a fixed-scaling copy weighted by x,
# and its structural tests stay below, but the gate refuses the structure
# (correlated_covariance_requires_brms): under the joint Wishart prior that
# field carries, the slope variance on lme4::sleepstudy came back a third
# above REML's and the correlation moved when the response was rescaled.
#
# The correlation is the whole point of the class, so the recovery test does
# not stop at "a correlation came back". It checks that a NEGATIVE
# correlation comes back negative, and it checks the whole triple against
# lme4's REML fit of the same model, an oracle this package did not write.

.cs_fixture <- function(
  seed = 2026L,
  n_g = 60L,
  n_i = 12L,
  sd_int = 1.2,
  sd_slope = 0.6,
  rho = 0.6,
  sd_resid = 0.35,
  beta_x = 1.5
) {
  set.seed(seed)
  cov_is <- rho * sd_int * sd_slope
  mat_s <- matrix(c(sd_int^2, cov_is, cov_is, sd_slope^2), 2L, 2L)
  mat_u <- matrix(rnorm(n_g * 2L), ncol = 2L) %*% chol(mat_s)
  tab <- expand.grid(obs = seq_len(n_i), g = factor(seq_len(n_g)))
  tab$x <- runif(nrow(tab), -2, 2)
  gi <- as.integer(tab$g)
  tab$y <- 5 +
    beta_x * tab$x +
    mat_u[gi, 1L] +
    mat_u[gi, 2L] * tab$x +
    rnorm(nrow(tab), 0, sd_resid)
  tab
}

test_that("the gate refuses (x | g) for INLA, by its own rule", {
  tab <- .cs_fixture(n_g = 12L, n_i = 4L)
  fb <- fb_from_brms(y ~ x + (x | g), data = tab, family = "gaussian")
  res <- lgm_gate(fb)
  expect_s3_class(res, "lgm_refusal")
  expect_identical(res$primary_rule, "correlated_covariance_requires_brms")
  # One reason, not two: the emit-allowlist check does not also call the
  # structure an emit gap.
  rules <- vapply(res$failures, function(f) f$rule_id, character(1L))
  expect_false("random_term_type_inla" %in% rules)
  expect_match(res$failures[[1L]]$reason, "`(x | g)`", fixed = TRUE)
})

test_that("the retained INLA emit is one iid2d field plus a copy", {
  tab <- .cs_fixture(n_g = 12L, n_i = 4L)
  fb <- fb_from_brms(y ~ x + (x | g), data = tab, family = "gaussian")
  txt <- gsub(
    "\\s+",
    " ",
    paste(deparse(.build_inla_formula(fb)), collapse = " ")
  )
  # One field of length 2 * n_g: intercepts then slopes.
  expect_match(txt, "f(g_cor_x_id, model = \"iid2d\", n = 24", fixed = TRUE)
  # The slope block is the SAME field, weighted by x, a copy, not a
  # second set of deviations.
  expect_match(
    txt,
    "f(g_cor_x_slope_id, x, copy = \"g_cor_x_id\")",
    fixed = TRUE
  )
})

test_that("the two index names are distinct and derived from the term", {
  # The iid2d field spans 1..n_g for the intercepts and n_g+1..2n_g for the
  # slopes, so the two blocks need two index columns whose names cannot
  # collide with each other or with a user column.
  expect_identical(.fb_inla_cor_slope_index("g", "x"), "g_cor_x_id")
  expect_identical(.fb_inla_cor_slope_slope_index("g", "x"), "g_cor_x_slope_id")
  # A group or slope name R would not accept as a symbol is made safe.
  expect_identical(
    .fb_inla_cor_slope_index("site name", "x"),
    "site.name_cor_x_id"
  )
})

test_that("a generated index name colliding with a user column refuses", {
  # emit_inla() checks that INLA is installed before it reads the data.
  skip_if_not_installed("INLA")
  tab <- .cs_fixture(n_g = 10L, n_i = 3L)
  tab$g_cor_x_id <- 1L
  fb <- fb_from_brms(y ~ x + (x | g), data = tab, family = "gaussian")
  expect_error(
    emit_inla(fb, data = tab, verbose = FALSE),
    "already has one"
  )
})

test_that("brms recovers a negative correlation and agrees with lme4's REML fit", {
  # Crossover between the intercept and the slope. A positive-correlation
  # fixture cannot distinguish a working sign from a structure that only
  # ever returns positive numbers.
  skip_if_not_installed("brms")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .cs_fixture(seed = 77L, n_g = 40L, n_i = 8L, rho = -0.7)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ x + (x | g),
    data = tab,
    backend = "auto",
    chains = 2L,
    n_samples = 800L,
    warmup = 400L,
    seed = 20260913L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  # `auto` routes the structure to brms instead of refusing it.
  expect_s3_class(fit, "flexybayes_brms")

  # brms's own labels: sd_<g>, sd_<g>__<x> and cor_<g>__Intercept__<x>. The
  # correlation row is the point of the structure. It used to be missing.
  v <- summary(fit)$varcomp
  expect_true(all(c("sd_g", "sd_g__x", "cor_g__Intercept__x") %in% v$component))

  got <- stats::setNames(v$estimate, v$component)
  row <- v[v$component == "cor_g__Intercept__x", ]
  expect_lt(row$estimate, 0)
  expect_lt(row$conf.high, 0)
  expect_equal(unname(got[["sd_g"]]), 1.2, tolerance = 0.25)
  expect_equal(unname(got[["sd_g__x"]]), 0.6, tolerance = 0.25)
  expect_equal(unname(got[["cor_g__Intercept__x"]]), -0.7, tolerance = 0.25)

  ref <- lme4::lmer(y ~ x + (x | g), data = tab)
  vc <- as.data.frame(lme4::VarCorr(ref))
  ref_sd_int <- vc$sdcor[
    vc$grp == "g" & vc$var1 == "(Intercept)" & is.na(vc$var2)
  ]
  ref_sd_slope <- vc$sdcor[vc$grp == "g" & vc$var1 == "x" & is.na(vc$var2)]
  ref_cor <- vc$sdcor[!is.na(vc$var2)]
  expect_equal(unname(got[["sd_g"]]), ref_sd_int, tolerance = 0.15)
  expect_equal(unname(got[["sd_g__x"]]), ref_sd_slope, tolerance = 0.15)
  expect_equal(unname(got[["cor_g__Intercept__x"]]), ref_cor, tolerance = 0.15)
})
