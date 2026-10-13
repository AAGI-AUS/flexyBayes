# Class guard DC-01: a value read on the wrong scale.
#
# Every standard deviation an INLA fit reports is checked three ways, on
# fixtures that cover a near-zero component, a random slope on a covariate
# in large units, and a Poisson slope:
#
# 1. units: multiplying the response by 1000 multiplies every reported SD
#    quantile by 1000, and dividing a covariate by 1000 multiplies its
#    slope SD by 1000 and leaves the rest alone (DEF-107, DEF-133, DEF-149,
#    DEF-153). The default priors are stated in the data's units, so both
#    changes leave the model the same;
# 2. surfaces: the medians and 95% limits in the variance table agree with
#    fb_variance_draws() and vpredict() (DEF-001, DEF-018);
# 3. engine: on the Gaussian fixture, sigma and the well-identified SD are
#    within 10% of lme4's REML estimates, which share no code with the
#    package (DEF-001).

.g01_fixtures <- function() {
  set.seed(6)
  d1 <- expand.grid(env = factor(1:6), gen = factor(1:20), blk = factor(1:4))
  d1$eb <- interaction(d1$env, d1$blk)
  # eb has true SD 0.
  d1$y <- 2 + stats::rnorm(20L)[d1$gen] + stats::rnorm(nrow(d1))

  set.seed(20261008L)
  d2 <- data.frame(g = factor(rep(sprintf("g%02d", 1:30), each = 8L)))
  d2$x <- stats::rnorm(nrow(d2))
  b <- stats::rnorm(30L, 0, 0.4)[as.integer(d2$g)]
  u <- stats::rnorm(30L, 0, 0.8)[as.integer(d2$g)]
  d2$y <- 1 + (0.5 + b) * d2$x + u + stats::rnorm(nrow(d2), 0, 0.5)
  d2$n <- stats::rpois(nrow(d2), exp(0.4 + (0.2 + b / 2) * d2$x + u / 2))

  list(
    gaussian_near_zero = list(
      args = list(y ~ env, random = ~ gen + eb, data = d1),
      response = "y", covariate = NA_character_
    ),
    gaussian_slope = list(
      args = list(y ~ x + (x || g), data = d2),
      response = "y", covariate = "x"
    ),
    poisson_slope = list(
      args = list(n ~ x + (x || g), data = d2, family = "poisson"),
      response = NA_character_, covariate = "x"
    )
  )
}

.g01_fit <- function(args) {
  suppressWarnings(suppressMessages(do.call(
    flexybayes, c(args, list(backend = "inla", verbose = FALSE))
  )))
}

# The untransformed fits, fitted once for the three blocks.
.g01_cache <- new.env(parent = emptyenv())
.g01_base_fit <- function(nm) {
  if (is.null(.g01_cache[[nm]])) {
    .g01_cache[[nm]] <- .g01_fit(.g01_fixtures()[[nm]]$args)
  }
  .g01_cache[[nm]]
}

# Median and 95% limits of every SD the fit reports, by component.
.g01_quantiles <- function(fit) {
  vc <- fit$extras$variance_comps
  med <- attr(vc, "posterior_median")
  vc <- vc[match(names(med), vc$component), ]
  out <- rbind(q2.5 = vc$q2.5, q50 = unname(med), q97.5 = vc$q97.5)
  colnames(out) <- names(med)
  out
}

.g01_rel <- function(a, b) max(abs(a / b - 1))

test_that("DC-01 guard: reported SDs follow the units of the response and covariates", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  # Two identical calls differ by up to 3.5% on the near-zero component's
  # lower limit. A value on the wrong scale is off by a factor.
  off <- character(0)
  for (nm in names(.g01_fixtures())) {
    fx <- .g01_fixtures()[[nm]]
    q <- .g01_quantiles(.g01_base_fit(nm))
    if (!is.na(fx$response)) {
      a <- fx$args
      a$data[[fx$response]] <- 1000 * a$data[[fx$response]]
      q_y <- .g01_quantiles(.g01_fit(a))
      err <- .g01_rel(q_y[, colnames(q)], 1000 * q)
      if (err > 0.1) off <- c(off, sprintf("%s, response x 1000: %.2g", nm, err))
    }
    if (!is.na(fx$covariate)) {
      a <- fx$args
      a$data[[fx$covariate]] <- a$data[[fx$covariate]] / 1000
      q_x <- .g01_quantiles(.g01_fit(a))[, colnames(q)]
      slope <- grepl(fx$covariate, colnames(q), fixed = TRUE)
      expect_true(any(slope), label = paste(nm, "reports a slope SD"))
      expected <- q
      expected[, slope] <- q[, slope] * 1000
      err <- .g01_rel(q_x, expected)
      if (err > 0.1) off <- c(off, sprintf("%s, covariate / 1000: %.2g", nm, err))
    }
  }
  expect_identical(off, character(0))
})

test_that("DC-01 guard: the variance table, the draws and vpredict() agree", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  off <- character(0)
  for (nm in names(.g01_fixtures())) {
    fit <- .g01_base_fit(nm)
    q <- .g01_quantiles(fit)
    draws <- fb_variance_draws(fit, 20000L, seed = 1L)
    dq <- vapply(colnames(q), function(k) {
      stats::quantile(draws[[k]], c(0.025, 0.5, 0.975), names = FALSE)
    }, numeric(3L))
    # Improved draws match the medians to about 1e-9 and the limits to
    # under 1e-3 at 20,000 draws (DEF-002).
    if (.g01_rel(dq[2L, ], q["q50", ]) > 1e-4) off <- c(off, paste(nm, "draw medians"))
    if (.g01_rel(dq[-2L, ], q[-2L, ]) > 2e-3) off <- c(off, paste(nm, "draw limits"))
    specs <- lapply(colnames(q), function(k) {
      stats::as.formula(paste0("v", match(k, colnames(q)), " ~ `", k, "`"))
    })
    set.seed(1)
    vp <- do.call(vpredict, c(list(fit), specs, list(n_draws = 20000L)))
    if (.g01_rel(vp$median, q["q50", ]) > 1e-4) off <- c(off, paste(nm, "vpredict medians"))
  }
  expect_identical(off, character(0))
})

test_that("DC-01 guard: well-identified SDs are near lme4's on the Gaussian fixture", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  fx <- .g01_fixtures()$gaussian_near_zero
  q <- .g01_quantiles(.g01_base_fit("gaussian_near_zero"))
  ref <- as.data.frame(lme4::VarCorr(lme4::lmer(
    y ~ env + (1 | gen) + (1 | eb), data = fx$args$data
  )))
  expect_lt(abs(q["q50", "sigma"] / ref$sdcor[ref$grp == "Residual"] - 1), 0.1)
  expect_lt(abs(q["q50", "sd_gen"] / ref$sdcor[ref$grp == "gen"] - 1), 0.1)
})
