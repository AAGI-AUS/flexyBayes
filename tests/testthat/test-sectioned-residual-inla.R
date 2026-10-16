# Heterogeneous residual by factor level on INLA: dsum(~ units | f).
#
# This refused with the claim that "INLA folds residual variance into the
# likelihood and has no separate residual-structure surface". The
# multiple-likelihood interface is exactly that surface: a matrix response
# with family = rep("gaussian", k) carries one observation precision per
# section. The claim was ours, not INLA's.
#
# Together with the interaction random effect this is the full
# multi-environment-trial model, the one ASReml is actually used for, and
# both halves refused on INLA before 0.10.0.9000.

.sr_fixture <- function(
  seed = 31L,
  n_gen = 30L,
  n_rep = 6L,
  sd_resid = c(0.25, 0.80, 1.60),
  sd_gen = 0.7
) {
  set.seed(seed)
  n_env <- length(sd_resid)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(n_rep)
  )
  g <- rnorm(n_gen, 0, sd_gen)
  tab$y <- 5 +
    g[as.integer(tab$gen)] +
    rnorm(nrow(tab), 0, sd_resid[as.integer(tab$env)])
  tab
}

test_that("the gate accepts a sectioned residual for INLA", {
  tab <- .sr_fixture(n_gen = 8L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~gen,
    residual = ~ dsum(~ units | env),
    data = tab,
    family = "gaussian"
  )
  res <- lgm_gate(fb)
  expect_true("lgm_compatible" %in% res$capabilities)
})

test_that("INLA recovers one residual SD per section", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  sd_resid <- c(0.25, 0.80, 1.60)
  tab <- .sr_fixture(sd_resid = sd_resid)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~gen,
    residual = ~ dsum(~ units | env),
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  for (i in seq_along(sd_resid)) {
    nm <- paste0("sigma_", i)
    expect_true(nm %in% names(got))
    expect_equal(unname(got[[nm]]), sd_resid[[i]], tolerance = 0.20)
  }
  expect_equal(unname(got[["sd_gen"]]), 0.7, tolerance = 0.20)
  # The scalar sigma is gone, there is no single residual to report.
  expect_false("sigma" %in% names(got))
  # And the internal matrix-response column must not surface anywhere.
  expect_false(any(grepl("fb_sectioned_response", names(got), fixed = TRUE)))
})

test_that("the baseline MET model fits: interaction random effect and dsum", {
  # Both halves refused on INLA before 0.10.0.9000, so this combination was
  # unreachable there. It is the model ASReml is bought for.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  set.seed(61L)
  n_gen <- 50L
  n_env <- 3L
  sd_resid <- c(0.25, 0.7, 1.5)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(5L)
  )
  g <- rnorm(n_gen, 0, 0.9)
  ge <- rnorm(n_gen * n_env, 0, 0.5)
  i_ge <- as.integer(tab$gen) + (as.integer(tab$env) - 1L) * n_gen
  tab$y <- 5 +
    g[as.integer(tab$gen)] +
    ge[i_ge] +
    rnorm(nrow(tab), 0, sd_resid[as.integer(tab$env)])

  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env),
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  expect_equal(unname(got[["sd_gen"]]), 0.9, tolerance = 0.20)
  expect_equal(unname(got[["sd_gen:env"]]), 0.5, tolerance = 0.25)
  for (i in seq_along(sd_resid)) {
    expect_equal(
      unname(got[[paste0("sigma_", i)]]),
      sd_resid[[i]],
      tolerance = 0.20
    )
  }
})

test_that("a level-conditioned residual still refuses", {
  # at(f, level):units conditions the residual on selected levels, which is
  # a different model from one variance per level and has no
  # multiple-likelihood lowering.
  tab <- .sr_fixture(n_gen = 8L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~gen,
    residual = ~ at(env, "1"):units,
    data = tab,
    family = "gaussian"
  )
  res <- lgm_gate(fb)
  expect_false("lgm_compatible" %in% res$capabilities)
})

test_that("a residual prior reaches every section instead of being dropped", {
  # The first working version passed `control_family[[1L]]` per section,
  # the `hyper` element on its own. INLA warned that 'prec' was not a valid
  # control.family name and REMOVED it, so the residual prior was silently
  # discarded and the fit ran on INLA's default. The fit still "worked",
  # which is why this needs a test instead of an eye.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .sr_fixture(n_gen = 10L, n_rep = 3L)
  expect_no_warning(
    suppressMessages(flexybayes(
      y ~ env,
      random = ~gen,
      residual = ~ dsum(~ units | env),
      data = tab,
      prior = fb_prior(sigma ~ pc(upper = 2, prob = 0.05)),
      backend = "inla",
      verbose = FALSE
    ))
  )
})
