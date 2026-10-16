# Unstructured genotype-by-environment covariance: us(f):g.
#
# Fitted on brms as `(0 + f | g)`, one standard deviation per level of `f`
# and every correlation between levels.
#
# INLA refuses it by name. The INLA lowering is still in the emitter, the
# whole k-by-k covariance as one iidKd field over a blocked index, each row
# indexing one element of a length k*n_inner field, and its structural
# tests stay below, but the gate refuses the structure
# (correlated_covariance_requires_brms). Finding O04 measured the genetic
# correlation on that route at 0.142 against ASReml's 0.930, rescaling the
# response by ten did not reproduce it, and under an identical Wishart prior
# Stan returns what INLA returns: the joint prior couples the variances and
# the correlations, which is why the structure moved to brms instead of
# being tuned on INLA.

.us_fixture <- function(
  seed = 41L,
  n_gen = 60L,
  n_rep = 5L,
  sd_by_env = c(0.8, 1.2, 0.5),
  sd_resid = 0.35
) {
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
  tab$y <- 5 +
    mat_u[cbind(as.integer(tab$gen), as.integer(tab$env))] +
    rnorm(nrow(tab), 0, sd_resid)
  tab
}

test_that("the gate refuses us(f):g for INLA, by its own rule", {
  tab <- .us_fixture(n_gen = 10L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ us(env):gen,
    data = tab,
    family = "gaussian"
  )
  res <- lgm_gate(fb)
  expect_s3_class(res, "lgm_refusal")
  expect_identical(res$primary_rule, "correlated_covariance_requires_brms")
  rules <- vapply(res$failures, function(f) f$rule_id, character(1L))
  expect_false("random_term_type_inla" %in% rules)
  expect_match(res$failures[[1L]]$reason, "`us(env):gen`", fixed = TRUE)
})

test_that("the retained INLA emit carries one iidKd field", {
  tab <- .us_fixture(n_gen = 10L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ us(env):gen,
    data = tab,
    family = "gaussian"
  )
  txt <- gsub(
    "\\s+",
    " ",
    paste(deparse(.build_inla_formula(fb)), collapse = " ")
  )
  expect_match(txt, "f(gen_by_env_id, model = \"iid3d\"", fixed = TRUE)
})

test_that("us() beyond five levels is refused on INLA and names fa() as well", {
  # k(k+1)/2 parameters is the problem fa() exists to solve, and INLA's
  # iidkd beyond k = 5 is a different parameterisation. The order check still
  # fires beside the correlated-covariance rule, so the refusal also points
  # at the reduced-rank structure.
  tab <- .us_fixture(
    n_gen = 12L,
    n_rep = 2L,
    sd_by_env = c(0.8, 1.2, 0.5, 0.9, 1.1, 0.7)
  )
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ us(env):gen,
    data = tab,
    family = "gaussian"
  )
  res <- lgm_gate(fb)
  expect_false("lgm_compatible" %in% res$capabilities)
  rules <- vapply(res$failures, function(f) f$rule_id, character(1L))
  expect_true(all(
    c("correlated_covariance_requires_brms", "us_gxe_order_inla") %in% rules
  ))
})

test_that("brms fits us(f):g with one SD per level and the correlations", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  sd_by_env <- c(0.8, 1.2, 0.5)
  tab <- .us_fixture(n_gen = 40L, n_rep = 3L, sd_by_env = sd_by_env)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ env,
    random = ~ us(env):gen,
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

  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)
  # brms names these sd_<g>__<f><level> and cor_<g>__<f>I__<f>J.
  for (i in seq_along(sd_by_env)) {
    nm <- paste0("sd_gen__env", i)
    expect_true(nm %in% names(got), info = nm)
    expect_equal(unname(got[[nm]]), sd_by_env[[i]], tolerance = 0.25)
  }
  expect_true("cor_gen__env1__env2" %in% names(got))
  expect_equal(unname(got[["cor_gen__env1__env2"]]), 0.7, tolerance = 0.25)
  expect_equal(unname(got[["sigma"]]), 0.35, tolerance = 0.15)
  expect_false(any(grepl("_by_env_id|Rho|component", v$component)))
})
