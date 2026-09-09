# Factor-analytic genotype-by-environment covariance on INLA: fa(f, k):g.
#
# ASReml's headline multi-environment-trial model, and the one flexyBayes
# declared no engine could fit. INLA has no latent model called "factor
# analytic", which is what the old refusal said; it does not need one. Each
# common factor is a single latent field entered once per level of the outer
# factor as a SCALED COPY, and `copy=`'s scaling parameter IS the loading --
# so the loadings arrive with a posterior, which REML cannot give.
#
# Factor j is anchored at level j (loading fixed at 1) and entered only from
# level j onward: the zero upper triangle ASReml uses to identify the loading
# matrix against rotation. An anchored loading is therefore absent from the
# summary rather than reported as an estimate of 1.

.fa_fixture <- function(seed = 51L, n_gen = 80L, n_rep = 6L,
                        loadings = c(1.0, 1.6, 0.5, -0.9),
                        sd_f = 0.9, sd_spec = 0.25, sd_resid = 0.30) {
  set.seed(seed)
  n_env <- length(loadings)
  vec_f <- rnorm(n_gen, 0, sd_f)
  mat_u <- outer(vec_f, loadings) +
    matrix(rnorm(n_gen * n_env, 0, sd_spec), ncol = n_env)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(n_rep)
  )
  tab$y <- 5 + mat_u[cbind(as.integer(tab$gen), as.integer(tab$env))] +
    rnorm(nrow(tab), 0, sd_resid)
  tab
}

test_that("the gate accepts fa(f, k):g for INLA", {
  tab <- .fa_fixture(n_gen = 12L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    family = "gaussian"
  )
  expect_true("lgm_compatible" %in% lgm_gate(fb)$capabilities)
})

test_that("the emitted formula carries scaled copies and specific variances", {
  tab <- .fa_fixture(n_gen = 12L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    family = "gaussian"
  )
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                collapse = " "))
  # The anchor is an ordinary iid field.
  expect_match(txt, "f(gen_fa1_env_1, model = \"iid\"", fixed = TRUE)
  # Every other level is a scaled copy of it -- that scaling is the loading.
  expect_match(txt, "copy = \"gen_fa1_env_1\", fixed = FALSE", fixed = TRUE)
  # And one specific variance per level.
  expect_match(txt, "f(gen_psi_env_4, model = \"iid\"", fixed = TRUE)
})

test_that("fa() recovers the loadings, including a negative one", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  loadings <- c(1.0, 1.6, 0.5, -0.9)
  tab <- .fa_fixture(loadings = loadings)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  # Loadings 2..4; the first is anchored at 1 and not estimated.
  for (e in 2:4) {
    nm <- paste0("fa1_env", e)
    expect_true(nm %in% names(got))
    expect_equal(unname(got[[nm]]), loadings[[e]], tolerance = 0.15)
  }
  # A crossover environment must come back NEGATIVE, not merely small.
  expect_lt(unname(got[["fa1_env4"]]), 0)

  expect_equal(unname(got[["sd_fa1"]]), 0.9, tolerance = 0.20)
  expect_equal(unname(got[["sigma"]]), 0.30, tolerance = 0.15)
  expect_true("sd_gen__env1" %in% names(got))
})

test_that("fa(f, 2) fits two factors with the anchoring pattern", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  set.seed(77L)
  n_gen <- 100L
  n_env <- 5L
  mat_l <- cbind(
    c(1, 1.4, 0.6, -0.8, 1.1),
    c(0, 1, 0.7, 0.4, -0.5)
  )
  mat_u <- cbind(rnorm(n_gen, 0, 0.9), rnorm(n_gen, 0, 0.6)) %*% t(mat_l) +
    matrix(rnorm(n_gen * n_env, 0, 0.25), ncol = n_env)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(5L)
  )
  tab$y <- 5 + mat_u[cbind(as.integer(tab$gen), as.integer(tab$env))] +
    rnorm(nrow(tab), 0, 0.3)

  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ fa(env, 2):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  # mode.status is INLA's own convergence oracle; it is numeric, not integer.
  expect_equal(as.numeric(fit$inla$mode$mode.status), 0)
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  for (e in 2:5) {
    expect_equal(unname(got[[paste0("fa1_env", e)]]), mat_l[e, 1],
                 tolerance = 0.20)
  }
  # Factor 2 is anchored at level 2, so it is estimated from level 3 on.
  for (e in 3:5) {
    expect_equal(unname(got[[paste0("fa2_env", e)]]), mat_l[e, 2],
                 tolerance = 0.25)
  }
  expect_false("fa2_env1" %in% names(got))
})

test_that("the specific variance does not collapse under the package prior", {
  # Left on INLA's default log-gamma precision, the specific variance at the
  # level carrying the largest loading collapses: a true 0.25 came back at
  # 0.007, with the interval excluding the truth. The package's own SD prior
  # reaches these terms via .fb_default_prior_targets() and holds them off
  # the boundary. This is the package's contribution on this structure, so it
  # gets a test rather than a comment.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  tab <- .fa_fixture()
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  spec <- v$estimate[grepl("^sd_gen__env", v$component)]
  expect_length(spec, 4L)
  expect_true(all(spec > 0.02),
              info = paste("specific SDs:", paste(round(spec, 4),
                                                  collapse = ", ")))
})

test_that("generated fa index names do not reach a summary label", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .fa_fixture(n_gen = 30L, n_rep = 3L)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  comps <- summary(fit)$varcomp$component
  expect_false(any(grepl("_psi_|_fa1_env_|Beta for", comps)))
})
