# =============================================================================
# What the engines themselves do, checked against the engines.
#
# A refusal, a capability note or a help page that says what INLA or brms
# can or cannot do is a claim about someone else's software. Each such claim
# the package keeps is checked here by calling the engine directly, not
# through flexyBayes, so a new engine release that changes the fact fails a
# test instead of leaving a false sentence in a message. Claims about what
# flexyBayes itself fits are tested through the public verbs elsewhere.
#
# The capability table (R/capability_matrix.R) names blocks of this file
# by their description in its `inla_basis` and `brms_basis` fields, and
# test-capability-matrix.R checks that each named block exists and calls
# the engine. Rename a block here and that test fails.
#
# brms is checked through brms::stancode(), which writes the Stan program
# without compiling it. INLA is checked with small fits on one thread.
# =============================================================================

.ec_inla_quiet <- function(expr) {
  INLA::inla.setOption(num.threads = "1:1")
  out <- NULL
  utils::capture.output(
    out <- tryCatch(suppressWarnings(expr), error = function(e) e),
    type = "output"
  )
  out
}

# --- INLA ------------------------------------------------------------ #

test_that("INLA's likelihood list has no hurdle or zero-inflated gamma", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  lik <- names(INLA::inla.models()$likelihood)
  expect_false(any(grepl("hurdle", lik)))
  expect_false(any(grepl("^zeroinflated.*gamma", lik)))
  # The families the package's own messages say INLA does carry.
  expect_true(all(c("gamma", "lognormal", "tweedie", "gev", "bgev") %in% lik))
})

test_that("INLA refuses one variable as both a fixed effect and an f() index", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(1L)
  d <- data.frame(y = stats::rnorm(30L), x = rep(seq_len(10L), 3L))
  res <- .ec_inla_quiet(INLA::inla(y ~ x + f(x, model = "iid"), data = d))
  expect_s3_class(res, "error")
})

test_that("INLA cannot name a data column whose name contains a colon", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(2L)
  d <- data.frame(y = stats::rnorm(30L), z = stats::rnorm(30L))
  names(d)[[2L]] <- "a2:b2"
  res <- .ec_inla_quiet(INLA::inla(y ~ `a2:b2`, data = d))
  expect_s3_class(res, "error")
  expect_match(conditionMessage(res), "a2:b2", fixed = TRUE)
})

test_that("INLA's correlated iid models carry a Wishart prior", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  latent <- INLA::inla.models()$latent
  for (m in c("iid2d", "iid3d", "iid4d", "iid5d", "iidkd")) {
    priors <- vapply(
      latent[[m]]$hyper,
      function(h) as.character(h$prior),
      character(1L)
    )
    expect_true(grepl("^wishart", priors[[1L]]), label = m)
  }
})

test_that("INLA's exchangeable group model fits one variance and one shared correlation", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(3L)
  n_gen <- 40L
  n_env <- 4L
  cor_env <- matrix(0.6, n_env, n_env)
  diag(cor_env) <- 1
  u <- matrix(stats::rnorm(n_gen * n_env), n_gen) %*% chol(cor_env)
  d <- expand.grid(gen = seq_len(n_gen), env = seq_len(n_env), rep = 1:2)
  d$y <- u[cbind(d$gen, d$env)] + stats::rnorm(nrow(d), sd = 0.5)
  fit <- .ec_inla_quiet(INLA::inla(
    y ~ 1 + f(
      gen,
      model = "iid",
      group = env,
      control.group = list(model = "exchangeable")
    ),
    data = d
  ))
  expect_false(inherits(fit, "error"))
  hyper <- rownames(fit$summary.hyperpar)
  expect_true(any(grepl("^Precision for gen", hyper)))
  expect_true(any(grepl("^GroupRho for gen", hyper)))
})

test_that("INLA builds a trait-by-genotype covariance from a known-covariance field and a scaled copy", {
  # Trait 1 carries a field u with known covariance A; trait 2 carries
  # beta * u plus an independent field v with the same A. The trait
  # covariance is then [s1^2, beta s1^2; beta s1^2, beta^2 s1^2 + s2^2]
  # crossed with A, which is the unstructured trait covariance.
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(4L)
  n_gen <- 30L
  a <- crossprod(matrix(stats::rnorm(n_gen * n_gen), n_gen)) / n_gen +
    diag(n_gen)
  q <- solve(a)
  sig <- matrix(c(1, 0.7, 0.7, 1.5), 2L)
  z <- t(chol(a)) %*% matrix(stats::rnorm(n_gen * 2L), n_gen) %*% chol(sig)
  d <- expand.grid(gen = seq_len(n_gen), trait = 1:2, rep = 1:3)
  d$y <- z[cbind(d$gen, d$trait)] + stats::rnorm(nrow(d), sd = 0.4)
  d$t2 <- as.numeric(d$trait == 2L)
  d$ia <- ifelse(d$trait == 1L, d$gen, NA)
  d$ic <- ifelse(d$trait == 2L, d$gen, NA)
  d$ib <- ifelse(d$trait == 2L, d$gen, NA)
  fit <- .ec_inla_quiet(INLA::inla(
    y ~ 1 + t2 +
      f(ia, model = "generic0", Cmatrix = q, values = seq_len(n_gen)) +
      f(ic, copy = "ia", hyper = list(beta = list(fixed = FALSE))) +
      f(ib, model = "generic0", Cmatrix = q, values = seq_len(n_gen)),
    data = d
  ))
  expect_false(inherits(fit, "error"))
  expect_true("Beta for ic" %in% rownames(fit$summary.hyperpar))
})

test_that("INLA's ar model takes a second-order autoregression", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(5L)
  d <- data.frame(t = seq_len(100L))
  d$y <- as.numeric(stats::arima.sim(list(ar = c(0.5, 0.3)), 100L)) +
    stats::rnorm(100L, sd = 0.3)
  fit <- .ec_inla_quiet(INLA::inla(
    y ~ 1 + f(t, model = "ar", order = 2L),
    data = d
  ))
  expect_false(inherits(fit, "error"))
  expect_true("PACF2 for t" %in% rownames(fit$summary.hyperpar))
})

test_that("INLA takes a probit link for the binomial family", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(6L)
  d <- data.frame(x = stats::rnorm(200L))
  d$y <- stats::rbinom(200L, 1L, stats::pnorm(0.3 + 0.8 * d$x))
  fit <- .ec_inla_quiet(INLA::inla(
    y ~ x,
    family = "binomial",
    Ntrials = 1,
    data = d,
    control.family = list(link = "probit")
  ))
  expect_false(inherits(fit, "error"))
  expect_true("x" %in% rownames(fit$summary.fixed))
})

# --- brms ------------------------------------------------------------ #

.ec_stancode <- function(...) {
  code <- tryCatch(brms::stancode(...), error = function(e) e)
  if (inherits(code, "error")) {
    return(code)
  }
  strsplit(as.character(code), "\n", fixed = TRUE)[[1L]]
}

.ec_genotypes <- function() {
  set.seed(7L)
  gen <- factor(rep(paste0("g", seq_len(12L)), each = 3L))
  a <- crossprod(matrix(stats::rnorm(144L), 12L)) / 12 + diag(12L)
  dimnames(a) <- list(levels(gen), levels(gen))
  list(gen = gen, a = a)
}

test_that("brms writes a trait-by-genotype covariance with a known genotype covariance", {
  skip_on_cran()
  skip_if_not_installed("brms")
  g <- .ec_genotypes()
  # Long format: one response column and a trait factor, as flexyBayes
  # data are laid out.
  d_long <- data.frame(
    gen = rep(g$gen, 2L),
    trait = factor(rep(c("t1", "t2"), each = 36L)),
    y = stats::rnorm(72L)
  )
  code <- .ec_stancode(
    y ~ 0 + trait + (0 + trait | gr(gen, cov = A)),
    data = d_long,
    data2 = list(A = g$a)
  )
  expect_false(inherits(code, "error"))
  expect_true(any(grepl("Lcov_1", code, fixed = TRUE)))
  expect_true(any(grepl("cholesky_factor_corr[M_1] L_1", code, fixed = TRUE)))
  # Wide format: one column per trait.
  d_wide <- data.frame(
    gen = g$gen,
    y1 = stats::rnorm(36L),
    y2 = stats::rnorm(36L)
  )
  code_wide <- .ec_stancode(
    brms::bf(brms::mvbind(y1, y2) ~ 1 + (1 | p | gr(gen, cov = A))) +
      brms::set_rescor(FALSE),
    data = d_wide,
    data2 = list(A = g$a)
  )
  expect_false(inherits(code_wide, "error"))
  expect_true(any(grepl("Lcov_1", code_wide, fixed = TRUE)))
})

test_that("brms has no tweedie, compound Poisson or zero-inflated gamma family", {
  skip_on_cran()
  skip_if_not_installed("brms")
  for (fam in c("tweedie", "compound_poisson", "zero_inflated_gamma")) {
    expect_error(brms::brmsfamily(fam), label = fam)
  }
  # The two families the package's messages point to instead.
  expect_s3_class(brms::brmsfamily("hurdle_gamma"), "brmsfamily")
  expect_s3_class(brms::brmsfamily("gen_extreme_value"), "brmsfamily")
})

test_that("brms's mi() term takes a gamma response and refuses a Poisson one", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(8L)
  d <- data.frame(
    x = stats::rnorm(30L),
    count = c(NA, stats::rpois(29L, 3)),
    positive = c(NA, exp(stats::rnorm(29L)))
  )
  expect_s3_class(
    .ec_stancode(count | mi() ~ x, data = d, family = stats::poisson()),
    "error"
  )
  for (fam in list(stats::Gamma(link = "log"), brms::lognormal())) {
    code <- .ec_stancode(positive | mi() ~ x, data = d, family = fam)
    expect_false(inherits(code, "error"), label = fam$family)
  }
})

test_that("brms has an ar() term for a one-dimensional residual autocorrelation", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(9L)
  d <- data.frame(t = seq_len(30L), y = stats::rnorm(30L))
  code <- .ec_stancode(y ~ 1 + ar(time = t, p = 1), data = d)
  expect_false(inherits(code, "error"))
  expect_true(any(grepl("Kar", code, fixed = TRUE)))
})

test_that("brms fits an mgcv smooth s(x) in its own formula", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("mgcv")
  set.seed(10L)
  d <- data.frame(x = stats::runif(40L), y = stats::rnorm(40L))
  code <- .ec_stancode(y ~ s(x), data = d)
  expect_false(inherits(code, "error"))
  expect_true(any(grepl("sds_1", code, fixed = TRUE)))
})

test_that("brms puts an LKJ prior with eta 1 on group-level correlations by default", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(11L)
  d <- data.frame(
    g = factor(rep(seq_len(10L), each = 6L)),
    f = factor(rep(seq_len(3L), 20L)),
    y = stats::rnorm(60L)
  )
  correlated <- .ec_stancode(y ~ 1 + (0 + f | g), data = d)
  expect_true(any(grepl("lkj_corr_cholesky_lpdf(L_1 | 1)", correlated,
    fixed = TRUE)))
  uncorrelated <- .ec_stancode(y ~ 1 + (0 + f || g), data = d)
  expect_false(any(grepl("lkj_corr_cholesky", uncorrelated, fixed = TRUE)))
})

test_that("brms refuses a coefficient prior beside a class-wide horseshoe or R2D2", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(12L)
  d <- data.frame(
    y = stats::rnorm(40L),
    p1 = stats::rnorm(40L),
    x1 = stats::rnorm(40L),
    x2 = stats::rnorm(40L)
  )
  for (special in c("horseshoe(1)", "R2D2()")) {
    res <- .ec_stancode(
      y ~ p1 + x1 + x2,
      data = d,
      prior = c(
        brms::set_prior(special),
        brms::set_prior("normal(0, 1)", coef = "p1")
      )
    )
    expect_s3_class(res, "error")
    expect_match(conditionMessage(res), "special priors", fixed = TRUE)
  }
})

test_that("brms's weights() multiplies each observation's log-likelihood", {
  skip_on_cran()
  skip_if_not_installed("brms")
  set.seed(13L)
  d <- data.frame(y = stats::rnorm(40L), x = stats::rnorm(40L),
    w = stats::runif(40L))
  code <- .ec_stancode(y | weights(w) ~ x, data = d)
  expect_true(any(grepl("target += weights[n] * (normal_lpdf(", code,
    fixed = TRUE)))
})
