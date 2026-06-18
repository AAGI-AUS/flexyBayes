# Boundary test: at a sample size where BOTH the per-row and the
# aggregated paths are feasible, the fitted posterior must agree. The
# aggregated likelihood is algebraically identical to the per-row
# likelihood, so under matched priors INLA's deterministic Laplace
# engine recovers the same fixed-effect coefficients and variance
# components to numerical precision. greta is cross-checked against the
# exact INLA fit within Monte-Carlo tolerance.

# Posterior-summary extractors that read the raw INLA fit on both the
# per-row and the aggregated fit objects (both carry `$inla`).
.eq_inla_beta <- function(fit) {
  s <- fit$inla$summary.fixed
  stats::setNames(s$mean, rownames(s))
}
.eq_inla_hyper <- function(fit) {
  h <- fit$inla$summary.hyperpar
  stats::setNames(h$mean, rownames(h))
}

.eq_make_data <- function(family, n = 6e4L, ng = 30L, seed = 909L) {
  set.seed(seed)
  env_eff <- c(0, 1.2, -0.7, 0.5)
  geno_eff <- stats::rnorm(ng, 0, 0.6)
  env <- sample(seq_along(env_eff), n, replace = TRUE)
  geno <- sample(seq_len(ng), n, replace = TRUE)
  eta <- env_eff[env] + geno_eff[geno]
  y <- switch(
    family,
    gaussian = eta + stats::rnorm(n, 0, 1.5),
    poisson = stats::rpois(n, exp(eta - 0.5)),
    binomial = stats::rbinom(n, 1, stats::plogis(eta - 1))
  )
  data.frame(env = factor(env), geno = factor(geno), y = y)
}

.eq_prior <- function(family) {
  if (identical(family, "gaussian")) {
    fb_prior(
      sigma ~ uniform(lower = 0, upper = 10),
      sd(group = "geno") ~ uniform(lower = 0, upper = 10)
    )
  } else {
    fb_prior(sd(group = "geno") ~ uniform(lower = 0, upper = 10))
  }
}

for (fam in c("gaussian", "binomial", "poisson")) {
  test_that(paste0("INLA: per-row == streamed-aggregated [", fam, "]"), {
    skip_if_not_installed("INLA")
    withr::local_options(flexyBayes.silence_uniform_inla_approx = TRUE)
    df <- .eq_make_data(fam)
    pr <- .eq_prior(fam)

    f_row <- flexybayes(
      y ~ env,
      random = ~geno,
      data = df,
      family = fam,
      backend = "inla",
      aggregate = FALSE,
      prior = pr,
      verbose = FALSE
    )
    f_str <- flexybayes_stream(
      y ~ env,
      random = ~geno,
      source = df,
      family = fam,
      backend = "inla",
      chunk_rows = 1.5e4,
      prior = pr,
      verbose = FALSE
    )

    br <- .eq_inla_beta(f_row)
    bs <- .eq_inla_beta(f_str)
    expect_equal(bs[names(br)], br, tolerance = 1e-4)

    hr <- .eq_inla_hyper(f_row)
    hs <- .eq_inla_hyper(f_str)
    expect_equal(hs[names(hr)], hr, tolerance = 5e-3)

    # The aggregated fit must report a real compression on this design.
    expect_lt(f_str$extras$aggregation_meta$compression, 0.05)
  })
}

for (fam in c("binomial", "poisson")) {
  test_that(paste0("in-memory aggregate=TRUE == per-row [", fam, "]"), {
    skip_if_not_installed("INLA")
    df <- .eq_make_data(fam, n = 3e4L, ng = 20L, seed = 314L)
    pr <- .eq_prior(fam)

    f_row <- flexybayes(
      y ~ env,
      random = ~geno,
      data = df,
      family = fam,
      backend = "inla",
      aggregate = FALSE,
      prior = pr,
      verbose = FALSE
    )
    f_agg <- flexybayes(
      y ~ env,
      random = ~geno,
      data = df,
      family = fam,
      backend = "inla",
      aggregate = TRUE,
      prior = pr,
      verbose = FALSE
    )

    expect_identical(f_agg$exactness, "aggregated_exact")
    expect_identical(backend_decision(f_agg)$path, "aggregated_count")

    br <- .eq_inla_beta(f_row)
    ba <- .eq_inla_beta(f_agg)
    expect_equal(ba[names(br)], br, tolerance = 1e-4)
    hr <- .eq_inla_hyper(f_row)
    ha <- .eq_inla_hyper(f_agg)
    expect_equal(ha[names(hr)], hr, tolerance = 5e-3)
  })
}

test_that("a non-Bernoulli binomial response is not auto-aggregated", {
  skip_if_not_installed("INLA")
  set.seed(8L)
  n <- 2000L
  df <- data.frame(
    env = factor(sample(1:3, n, replace = TRUE)),
    y = stats::rbinom(n, 5L, 0.4)
  ) # counts in 0..5
  # aggregate = TRUE must refuse (cannot recover per-row trials).
  expect_error(
    flexybayes(
      y ~ env,
      data = df,
      family = "poisson",
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    ),
    NA # poisson on these counts is fine and *is* aggregatable
  )
  expect_error(
    suppressWarnings(
      flexybayes(
        y ~ env,
        data = df,
        family = "binomial",
        backend = "inla",
        aggregate = TRUE,
        verbose = FALSE
      )
    ),
    "0/1|Bernoulli|trial"
  )
})

test_that("greta streamed-aggregated matches exact INLA within MC error", {
  skip_on_ci()
  skip_if_not_installed("INLA")
  skip_if_greta_backend_unusable()
  withr::local_options(flexyBayes.silence_uniform_inla_approx = TRUE)

  df <- .eq_make_data("gaussian", n = 3e4L, ng = 20L, seed = 111L)
  pr <- .eq_prior("gaussian")

  f_inla <- flexybayes_stream(
    y ~ env,
    random = ~geno,
    source = df,
    family = "gaussian",
    backend = "inla",
    prior = pr,
    verbose = FALSE
  )
  set.seed(2024L)
  f_greta <- flexybayes_stream(
    y ~ env,
    random = ~geno,
    source = df,
    family = "gaussian",
    backend = "greta",
    n_samples = 600L,
    warmup = 600L,
    chains = 2L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )

  bi <- coef(f_inla)
  bg <- coef(f_greta)
  # Two engines, one exact aggregated likelihood: agree within a loose
  # Monte-Carlo tolerance on the fixed effects.
  expect_equal(bg[names(bi)], bi, tolerance = 0.1)
})
