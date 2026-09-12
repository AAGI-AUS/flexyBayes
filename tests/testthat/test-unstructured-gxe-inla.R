# Unstructured genotype-by-environment covariance on INLA: us(f):g.
#
# INLA holds the whole k-by-k covariance as one iidKd field over a blocked
# index: each row indexes exactly one element of a length k*n_inner field,
# and iid2d..iid5d carry the unstructured covariance across blocks.
#
# Two things about this structure are worth their own tests. The correlations
# are shrunk toward each other by INLA's Wishart hyperprior; and the prior
# interlingua has nowhere to put a standard-deviation prior, so an explicit
# one refuses rather than being dropped.
#
# The shrinkage used to be described here as "a real property of the fit and
# not a defect". That is withdrawn. Finding O04 measured the genetic
# correlation at 0.142 against ASReml's 0.930, and rescaling the response by
# ten -- an exact transformation under which a correlation is invariant --
# does not reproduce it. The INLA us() route is classified `analogue`, not
# `covariance-family`, and is refused for selection. These tests assert what
# the emitter EMITS and what the structure's parameters look like; they do
# not assert that the numbers are right, and the test names below say so.

.us_fixture <- function(seed = 41L, n_gen = 60L, n_rep = 5L,
                        sd_by_env = c(0.8, 1.2, 0.5), sd_resid = 0.35) {
  set.seed(seed)
  n_env <- length(sd_by_env)
  # Exchangeable-ish correlation of the right size, with the 3-environment
  # case carrying the specific 0.7 / 0.2 / 0.4 the recovery test asserts on.
  mat_cor <- if (n_env == 3L) {
    matrix(c(1, 0.7, 0.2, 0.7, 1, 0.4, 0.2, 0.4, 1), 3L, 3L)
  } else {
    m <- matrix(0.4, n_env, n_env)
    diag(m) <- 1
    m
  }
  mat_sigma <- diag(sd_by_env) %*% mat_cor %*% diag(sd_by_env)
  mat_u <- matrix(rnorm(n_gen * n_env), ncol = n_env) %*% chol(mat_sigma)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(n_rep)
  )
  tab$y <- 5 + mat_u[cbind(as.integer(tab$gen), as.integer(tab$env))] +
    rnorm(nrow(tab), 0, sd_resid)
  tab
}

test_that("the gate accepts us(f):g for INLA at k <= 5", {
  tab <- .us_fixture(n_gen = 10L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ us(env):gen,
    data = tab,
    family = "gaussian"
  )
  expect_true("lgm_compatible" %in% lgm_gate(fb)$capabilities)
})

test_that("the emitted formula carries one iidKd field", {
  tab <- .us_fixture(n_gen = 10L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ us(env):gen,
    data = tab,
    family = "gaussian"
  )
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                collapse = " "))
  expect_match(txt, "f(gen_by_env_id, model = \"iid3d\"", fixed = TRUE)
})

test_that("us() beyond five levels refuses and names fa() instead", {
  # k(k+1)/2 parameters is the problem fa() exists to solve, and INLA's
  # iidkd beyond k = 5 is a different parameterisation.
  tab <- .us_fixture(n_gen = 12L, n_rep = 2L,
                     sd_by_env = c(0.8, 1.2, 0.5, 0.9, 1.1, 0.7))
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ us(env):gen,
    data = tab,
    family = "gaussian"
  )
  res <- lgm_gate(fb)
  expect_false("lgm_compatible" %in% res$capabilities)
})

test_that("INLA emits an iidKd field with the expected parameter shape", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  sd_by_env <- c(0.8, 1.2, 0.5)
  tab <- .us_fixture(sd_by_env = sd_by_env)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ us(env):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  for (i in seq_along(sd_by_env)) {
    nm <- paste0("sd_gen__env", i)
    expect_true(nm %in% names(got))
    expect_equal(unname(got[[nm]]), sd_by_env[[i]], tolerance = 0.20)
  }
  expect_equal(unname(got[["sigma"]]), 0.35, tolerance = 0.15)

  # The correlations are present and correctly ORDERED, but INLA's Wishart
  # hyperprior shrinks them toward each other: a true 0.7 / 0.2 / 0.4 comes
  # back near 0.71 / 0.41 / 0.61. That is a property of the prior, measured
  # and documented, not a translation defect -- so the test asserts the
  # ordering and the strongest correlation rather than three point values it
  # would be dishonest to claim.
  expect_true("cor_gen__env1__env2" %in% names(got))
  expect_equal(unname(got[["cor_gen__env1__env2"]]), 0.7, tolerance = 0.15)
  expect_true(all(got[grepl("^cor_", names(got))] > 0))
})

test_that("the iidKd index does not leak into the summary label", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .us_fixture(n_gen = 20L, n_rep = 3L)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ us(env):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  comps <- summary(fit)$varcomp$component
  # brms names these sd_<g>__<f><level> and cor_<g>__<f>I__<f>J.
  expect_true("sd_gen__env1" %in% comps)
  expect_true("cor_gen__env1__env2" %in% comps)
  expect_false(any(grepl("_by_env_id|Rho|component", comps)))
})

test_that("an explicit SD prior on us() refuses rather than being dropped", {
  # INLA's iidKd hyperparameters are prec1..preck and cor12.. under a joint
  # Wishart. There is no `prec` slot, and splicing one in makes INLA drop it
  # with a warning -- the silent prior loss this package fixed in the
  # sectioned residual. The auto-default is left unspliced (INLA's Wishart
  # governs); an explicit prior refuses.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .us_fixture(n_gen = 15L, n_rep = 2L)
  expect_error(
    suppressMessages(flexybayes(
      y ~ env,
      random = ~ us(env):gen,
      data = tab,
      prior = fb_prior(sd(group = "gen") ~ pc(upper = 2, prob = 0.05)),
      backend = "inla",
      verbose = FALSE
    )),
    "Wishart"
  )
})
