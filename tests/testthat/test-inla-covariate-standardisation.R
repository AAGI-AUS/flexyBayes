# =============================================================================
# INLA fits on standardised covariates.
#
# A fit that carries no response standardisation runs on numeric covariates
# divided by their standard deviations, with every stated prior restated on
# those coordinates and every reading mapped back where it is stored, draws
# where they are taken (R/inla_covariate_standardisation.R). Three kinds of
# evidence, none of which reads its expected value off the code under test:
#
#   * the eligibility rule and the restatement algebra against values written
#     out by hand here;
#   * the package's fit against a raw INLA call on the standardised covariate
#     with the same priors written out in raw units, the coefficient divided
#     back by hand;
#   * metamorphic invariance: multiplying a covariate by a constant changes no
#     model, so a coefficient in the covariate's own units must not move.
#
# Tolerances come from the spread of repeated identical fits, measured on INLA
# 26.8.7 and quoted where each is used. At the base of this unit the same
# fits moved 1.8e-2 (Poisson) and 1.7e-2 (Bernoulli) between a covariate and
# the same covariate multiplied by a thousand, against a run-to-run spread of
# 2.0e-7 over five identical fits per cell.
# =============================================================================

.cs_count_data <- function(family, seed = 20260916L) {
  set.seed(seed)
  n_g <- 30L
  n_r <- 8L
  d <- data.frame(g = factor(rep(sprintf("g%02d", seq_len(n_g)), each = n_r)))
  d$x <- stats::rnorm(nrow(d))
  d$z <- stats::rnorm(nrow(d), 3, 2)
  u <- stats::rnorm(n_g, 0, 0.7)
  eta <- 0.4 + 0.5 * d$x + u[as.integer(d$g)]
  d$y <- if (identical(family, "poisson")) {
    stats::rpois(nrow(d), exp(eta))
  } else {
    stats::rbinom(nrow(d), 1L, stats::plogis(eta))
  }
  d
}

.cs_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

.cs_code <- function(...) {
  .cs_quiet(flexybayes(..., backend = "inla", return_code = TRUE,
                       verbose = FALSE))
}

.cs_component <- function(fit, component) {
  v <- summary(fit)$varcomp
  v$estimate[v$component == component][[1L]]
}


# --- what is standardised ------------------------------------------------- #

test_that("a fit is standardised on its covariates only when its response is not", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .cs_count_data("poisson")

  count <- .cs_code(y ~ x, random = ~ g, data = d, family = "poisson")
  expect_equal(count$covariate_standardisation$divisor[["x"]], stats::sd(d$x))
  expect_equal(count$covariate_standardisation$coef_divisor[["x"]],
               stats::sd(d$x))
  expect_null(count$response_standardisation)

  # A Gaussian identity-link fit is already one model in its response's
  # units, and standardising both would restate every prior twice.
  gauss <- .cs_code(z ~ x, random = ~ g, data = d)
  expect_equal(gauss$response_standardisation$scale, stats::sd(d$z))
  expect_null(gauss$covariate_standardisation)
})

test_that("a covariate a formula transforms is left exactly as recorded", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .cs_count_data("poisson")
  d$x <- abs(d$x) + 1

  # log(x / s) is not log(x) divided by anything, so the coefficient could
  # not be mapped back and the column is protected instead.
  expect_null(
    .cs_code(y ~ log(x), random = ~ g, data = d,
             family = "poisson")$covariate_standardisation
  )
  # The same column reached through a transform anywhere in the model
  # protects it, even where it also enters untransformed.
  expect_null(
    .cs_code(y ~ x + log(x), random = ~ g, data = d,
             family = "poisson")$covariate_standardisation
  )
  # An offset is not a covariate.
  d$off <- log(d$x)
  expect_null(
    .cs_code(y ~ offset(off), random = ~ g, data = d,
             family = "poisson")$covariate_standardisation
  )
})

test_that("an interaction's coefficient carries the product of its covariates' divisors", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .cs_count_data("poisson")
  d$f <- factor(rep(c("a", "b"), length.out = nrow(d)))

  numeric_pair <- .cs_code(y ~ x * z, random = ~ g, data = d,
                           family = "poisson")
  cd <- numeric_pair$covariate_standardisation$coef_divisor
  expect_equal(cd[["x"]], stats::sd(d$x))
  expect_equal(cd[["z"]], stats::sd(d$z))
  expect_equal(cd[["x:z"]], stats::sd(d$x) * stats::sd(d$z))

  # A factor-by-numeric interaction: every indexed slope is in the
  # covariate's units, and the factor's own contrast is not.
  mixed <- .cs_code(y ~ f * x, random = ~ g, data = d, family = "poisson")
  cd <- mixed$covariate_standardisation$coef_divisor
  expect_equal(cd[["x"]], stats::sd(d$x))
  expect_equal(cd[["fb:x"]], stats::sd(d$x))
  expect_false("fb" %in% names(cd))
})


# --- the restatement algebra ---------------------------------------------- #

test_that("control.fixed is restated per coefficient, and a scalar expands", {
  cstd <- list(
    coef_divisor = c(x = 4, `f2:x` = 4),
    slope_divisor = numeric(0),
    slope_key_divisor = numeric(0)
  )
  expect_identical(
    .fb_inla_covariate_standardise_control_fixed(list(), cstd),
    list()
  )

  # A coefficient multiplied by 4 has its prior mean multiplied by 4 and its
  # precision divided by 16. `f2:x` carries no row of its own, so it takes
  # the default, restated by its own divisor.
  per_coef <- .fb_inla_covariate_standardise_control_fixed(
    list(mean = list(x = 2, default = 0), prec = list(x = 1 / 16, default = 1)),
    cstd
  )
  expect_equal(per_coef$mean$x, 8)
  expect_equal(per_coef$prec$x, 1 / 256)
  expect_equal(per_coef$mean[["f2:x"]], 0)
  expect_equal(per_coef$prec[["f2:x"]], 1 / 16)

  # A scalar is a statement about every coefficient at once, and the
  # standardised ones no longer share a value, so it expands. The intercept
  # keeps its own entry untouched: the standardisation does not centre.
  scalar <- .fb_inla_covariate_standardise_control_fixed(
    list(mean = 0, prec = 1 / 100, mean.intercept = 0, prec.intercept = 1 / 100),
    cstd
  )
  expect_equal(scalar$prec$x, 1 / 1600)
  expect_equal(scalar$prec$default, 1 / 100)
  expect_equal(scalar$mean$x, 0)
  expect_equal(scalar$prec.intercept, 1 / 100)
  expect_equal(scalar$mean.intercept, 0)
})

test_that("an uncorrelated slope prior is restated on the standardised covariate", {
  cstd <- list(
    coef_divisor = numeric(0),
    slope_divisor = c(g_for_slope_x = 5),
    slope_key_divisor = c(x_g = 5)
  )
  hyper <- list(
    x_g = list(
      prior = .inla_halfnormal_sd_expr(2),
      meta = list(family = "half_normal", scale = 2)
    ),
    g = list(
      prior = .inla_halfnormal_sd_expr(3),
      meta = list(family = "half_normal", scale = 3)
    )
  )
  out <- .fb_inla_covariate_standardise_hyper(hyper, cstd)
  # A standard deviation in response units per covariate unit multiplies by
  # the divisor; a random intercept carries no covariate and does not move.
  expect_equal(out$x_g$meta$scale, 10)
  expect_equal(out$g$meta$scale, 3)

  # And the precision of that field multiplies by the square of the divisor,
  # which is what maps it back.
  expect_equal(
    .fb_inla_covariate_hyper_factor("Precision for g_for_slope_x", cstd), 25
  )
  expect_equal(.fb_inla_covariate_hyper_factor("Precision for g", cstd), 1)
  expect_equal(
    .fb_inla_covariate_hyper_factor(
      "Precision for the Gaussian observations", cstd
    ),
    1
  )
})


# --- the fit against a raw INLA call -------------------------------------- #

test_that("the stored result equals a raw INLA fit on the standardised covariate", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  d <- .cs_count_data("poisson")
  s <- stats::sd(d$x)
  fit <- .cs_quiet(flexybayes(
    y ~ x,
    random = ~ g,
    data = d,
    family = "poisson",
    prior = fb_prior(sd(group = "g") ~ half_normal(scale = 2)),
    prior_fixed_sd = 10,
    backend = "inla",
    verbose = FALSE
  ))

  # The same model, written out in raw units on the standardised covariate:
  # a coefficient multiplied by s has prior SD 10 s, so precision 1 / (10 s)^2,
  # while the intercept keeps 1 / 10^2 because nothing is centred.
  dx <- d
  dx$x <- d$x / s
  raw <- INLA::inla(
    y ~ 1 + x + f(g, model = "iid", hyper = list(
      prec = list(prior = .inla_halfnormal_sd_expr(2))
    )),
    data = dx,
    family = "poisson",
    control.fixed = list(
      mean = list(x = 0, default = 0),
      prec = list(x = 1 / (10 * s)^2, default = 1 / 100),
      mean.intercept = 0,
      prec.intercept = 1 / 100
    ),
    control.compute = list(config = TRUE, return.marginals = TRUE,
                           dic = TRUE, waic = TRUE)
  )
  got <- fit$inla

  # The two calls are the same computation, and INLA has no random stream:
  # over five identical fits at the base the largest relative run-to-run
  # range on these quantities was 2.0e-7. Each tolerance below sits well
  # above that and orders of magnitude below any mapping error, which would
  # move the coefficient by the covariate's standard deviation itself.
  expect_equal(got$summary.fixed["x", "mean"],
               raw$summary.fixed["x", "mean"] / s, tolerance = 1e-6)
  expect_equal(got$summary.fixed["x", "sd"],
               raw$summary.fixed["x", "sd"] / s, tolerance = 1e-6)
  expect_equal(got$summary.fixed["(Intercept)", "mean"],
               raw$summary.fixed["(Intercept)", "mean"], tolerance = 1e-6)
  expect_equal(got$summary.random$g$mean, raw$summary.random$g$mean,
               tolerance = 1e-6)

  # The linear predictor is the same number on both coordinate systems, so
  # the likelihood-based quantities are carried across untouched.
  expect_equal(got$summary.fitted.values$mean,
               raw$summary.fitted.values$mean, tolerance = 1e-6)
  expect_lt(abs(got$mlik[1L, 1L] - raw$mlik[1L, 1L]), 1e-4)
  expect_lt(abs(got$dic$dic - raw$dic$dic), 1e-3)
  expect_lt(abs(got$waic$waic - raw$waic$waic), 1e-3)

  # The fit records what it ran on, and INLA's own argument record keeps the
  # standardised coordinates.
  cstd <- fit$extras$covariate_standardisation
  expect_equal(cstd$divisor[["x"]], s)
  expect_equal(stats::sd(fit$inla$.args$data$x), 1, tolerance = 1e-12)
})


# --- metamorphic invariance ----------------------------------------------- #

test_that("a non-Gaussian coefficient does not move when its covariate is rescaled", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  # Multiplying a covariate by a constant is an exact transformation: the
  # coefficient divides by that constant and the model is the one written.
  # The tolerance is 1e-5, about fifty times the run-to-run spread measured
  # at the base (2.0e-7) and three orders of magnitude below the departure
  # measured there (1.8e-2 on Poisson, 1.7e-2 on Bernoulli at x * 1000).
  for (family in c("poisson", "binomial")) {
    base <- .cs_count_data(family)
    cells <- vapply(c(0.001, 1, 1000), function(cc) {
      d <- base
      d$x <- cc * d$x
      fit <- .cs_quiet(flexybayes(
        y ~ x, random = ~ g, data = d, family = family,
        backend = "inla", verbose = FALSE
      ))
      stats::coef(fit)[["x"]] * cc
    }, numeric(1L))
    expect_lt(diff(range(cells)) / abs(mean(cells)), 1e-5, label = family)
  }
})

test_that("an uncorrelated slope comes back in the covariate's own units", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  base <- .cs_count_data("poisson")
  cells <- lapply(c(1, 1000), function(cc) {
    d <- base
    d$x <- cc * d$x
    fit <- .cs_quiet(flexybayes(
      y ~ x + (x || g), data = d, family = "poisson",
      backend = "inla", verbose = FALSE
    ))
    # The slope field and its standard deviation are in response units per
    # covariate unit, so both carry the factor back out.
    c(sd_slope = .cs_component(fit, "sd_x_g") * cc,
      sd_int = .cs_component(fit, "sd_g"),
      coef_x = stats::coef(fit)[["x"]] * cc,
      ranef_slope = stats::sd(fit$inla$summary.random[[
        paste0("g_for_slope_x")
      ]]$mean) * cc)
  })
  cells <- do.call(rbind, cells)
  for (q in colnames(cells)) {
    expect_lt(diff(range(cells[, q])) / abs(mean(cells[, q])), 1e-4, label = q)
  }
})

test_that("draws are taken on the original scale", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  d <- .cs_count_data("poisson")
  d$x <- 1000 * d$x
  fit <- .cs_quiet(flexybayes(
    y ~ x, random = ~ g, data = d, family = "poisson",
    backend = "inla", verbose = FALSE
  ))
  sf <- fit$inla$summary.fixed

  # INLA's sampler reads the standardised internal state; unmapped, these
  # draws would sit a thousandfold away from the summary they belong to.
  draws <- fb_as_draws_simple(fit, n_samples = 2000L)
  x_draws <- draws[["x:1"]]
  expect_lt(abs(mean(x_draws) - sf["x", "mean"]),
            5 * sf["x", "sd"] / sqrt(length(x_draws)))
  expect_equal(stats::sd(x_draws), sf["x", "sd"], tolerance = 0.1)

  v <- vcov(fit, n_samples = 2000L)
  expect_equal(sqrt(diag(v)), sf$sd, tolerance = 0.1, ignore_attr = TRUE)
})

test_that("vpredict comes back in the covariate's own units", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  # vpredict() is the only caller of the joint hyperparameter sampler, so it
  # is the only route on which that mapping shows. A slope standard deviation
  # is in response units per covariate unit; unmapped, it would sit a
  # thousandfold from the table it belongs to at x * 1000.
  #
  # The tolerance is the sampler's own. Over five repeated draw sets on one
  # fit, 2,000 draws each, the estimate ranged over 3.3 per cent of its mean
  # at x and 2.9 per cent at x * 1000, and stayed within 2.0 per cent of the
  # table value. 0.1 is about three times that spread, and a units error on
  # this path is a factor of a thousand.
  cells <- vapply(c(1, 1000), function(cc) {
    d <- .cs_count_data("poisson")
    d$x <- cc * d$x
    fit <- .cs_quiet(flexybayes(
      y ~ x + (x || g), data = d, family = "poisson",
      backend = "inla", verbose = FALSE
    ))
    vp <- .cs_quiet(vpredict(fit, r ~ sd_x_g, n_draws = 2000L))
    expect_equal(vp$estimate, .cs_component(fit, "sd_x_g"),
                 tolerance = 0.1, label = paste0("x * ", cc))
    vp$estimate * cc
  }, numeric(1L))
  expect_lt(diff(range(cells)) / abs(mean(cells)), 0.1)
})
