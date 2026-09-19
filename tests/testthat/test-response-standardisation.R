# =============================================================================
# Gaussian INLA fits on a standardised response.
#
# Every Gaussian INLA fit runs on y* = (y - m) / s, with stated priors
# restated on those coordinates and the result mapped back where it is
# stored, draws where they are taken. Three kinds of evidence, none of which
# reads its expected value off the code under test:
#
#   * the prior map against base R's own distribution functions;
#   * the package's fit against a raw INLA call at a response scale where the
#     raw call is sound, with the same priors written out in raw units here;
#   * metamorphic invariance: multiplying the response or a slope covariate by
#     a constant changes no model, so a scale-free summary must not move.
# =============================================================================

.rs_sleepstudy <- function(multiplier = 1, days_multiplier = 1) {
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  d <- e$sleepstudy
  d$Reaction <- multiplier * d$Reaction
  d$Days <- days_multiplier * d$Days
  d
}

.rs_varcomp <- function(fit) {
  v <- summary(fit)$varcomp
  stats::setNames(v$estimate, v$component)
}

# An INLA expression prior, evaluated in R. The expressions use a small
# C-like syntax; the ones tested here need only the statement separators and
# `return( ... )` rewritten.
.rs_eval_expression <- function(expr, theta) {
  body <- sub("^expression:", "", expr)
  body <- gsub("return\\((.*)\\);\\s*$", "\\1", body)
  env <- new.env()
  env$theta <- theta
  eval(parse(text = gsub(";", "\n", body)), envir = env)
}


# --- the prior map ------------------------------------------------------- #

test_that("each prior family maps to the distribution of sigma / k", {
  k <- 7.3
  p <- c(0.05, 0.25, 0.5, 0.9, 0.99)

  # Uniform on (L, U).
  u <- .fb_inla_map_hyper_entry(
    list(prior = .inla_uniform_sd_expr(12, 2),
         meta = list(family = "uniform", lower = 2, upper = 12)),
    k
  )
  expect_equal(stats::qunif(p, u$meta$lower, u$meta$upper),
               stats::qunif(p, 2, 12) / k)

  # Half-normal with scale 4.
  hn <- .fb_inla_map_hyper_entry(
    list(prior = .inla_halfnormal_sd_expr(4),
         meta = list(family = "half_normal", scale = 4)),
    k
  )
  expect_equal(stats::qnorm((1 + p) / 2, sd = hn$meta$scale),
               stats::qnorm((1 + p) / 2, sd = 4) / k)

  # Lognormal on sigma with log-scale sd 1.5.
  ln <- .fb_inla_map_hyper_entry(
    list(prior = .inla_lognormal_sd_expr(1.5),
         meta = list(family = "lognormal", scale = 1.5)),
    k
  )
  expect_equal(stats::qlnorm(p, ln$meta$location, ln$meta$scale),
               stats::qlnorm(p, 0, 1.5) / k)

  # Exponential with rate 0.2, carried as a PC prior.
  ex <- .fb_inla_map_hyper_entry(
    list(prior = "pc.prec", param = c(5, exp(-1)),
         meta = list(family = "exponential", rate = 0.2)),
    k
  )
  pc_rate <- -log(ex$param[[2L]]) / ex$param[[1L]]
  expect_equal(stats::qexp(p, pc_rate), stats::qexp(p, 0.2) / k)

  # PC prior P(sigma > 3) = 0.05, which is exponential on sigma.
  pc <- .fb_inla_map_hyper_entry(list(prior = "pc.prec", param = c(3, 0.05)), k)
  expect_equal(stats::qexp(p, -log(pc$param[[2L]]) / pc$param[[1L]]),
               stats::qexp(p, -log(0.05) / 3) / k)

  # Log-gamma on the precision: tau ~ Gamma(1, 5e-5), and sigma / k has
  # precision k^2 tau.
  lg <- .fb_inla_map_hyper_entry(list(prior = "loggamma", param = c(1, 5e-5)), k)
  expect_equal(stats::qgamma(p, lg$param[[1L]], rate = lg$param[[2L]]),
               k^2 * stats::qgamma(p, 1, rate = 5e-5))
})

test_that("the mapped lognormal expression is the density of theta under it", {
  # theta = -2 log(sigma), so p(theta) = dlnorm(sigma) * sigma / 2.
  location <- -log(7.3)
  scale <- 1.5
  expr <- .inla_lognormal_sd_expr(scale, location = location)
  theta <- c(-6, -2, 0, 1.7, 5)
  sigma <- exp(-theta / 2)
  expect_equal(
    .rs_eval_expression(expr, theta),
    log(stats::dlnorm(sigma, location, scale) * sigma / 2)
  )
  # And at location 0 the expression is the one the package always emitted.
  expect_identical(.inla_lognormal_sd_expr(scale, location = 0),
                   .inla_lognormal_sd_expr(scale))
})

test_that("an unmapped prior family refuses rather than passing through", {
  expect_error(
    .fb_inla_map_hyper_entry(list(prior = "gaussian", param = c(0, 1)), 2),
    "no standardisation rule"
  )
})

test_that("control.fixed is restated only where a prior was stated", {
  spec <- list(centre = 300, scale = 50)
  expect_identical(.fb_inla_standardise_control_fixed(list(), spec), list())

  scalar <- .fb_inla_standardise_control_fixed(
    list(mean = 0, prec = 1 / 49, mean.intercept = 0, prec.intercept = 1 / 49),
    spec
  )
  expect_equal(scalar$mean, 0)
  expect_equal(scalar$prec, 50^2 / 49)
  expect_equal(scalar$mean.intercept, -300 / 50)
  expect_equal(scalar$prec.intercept, 50^2 / 49)

  per_coef <- .fb_inla_standardise_control_fixed(
    list(mean = list(x = 2, default = 0), prec = list(x = 1 / 16, default = 1)),
    spec
  )
  expect_equal(per_coef$mean, list(x = 2 / 50, default = 0))
  expect_equal(per_coef$prec, list(x = 50^2 / 16, default = 50^2))
  expect_null(per_coef$mean.intercept)
})


# --- the fit against a raw INLA call at a safe scale --------------------- #

test_that("the stored result equals a raw INLA fit under the same priors", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()

  # A response scale where the raw call is sound: sleepstudy divided by ten,
  # standard deviation 5.6. The priors are written out in raw units here.
  d <- .rs_sleepstudy(0.1)
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- suppressWarnings(flexybayes(
    Reaction ~ Days + (1 | Subject),
    data = d,
    prior = fb_prior(
      sd(group = "Subject") ~ half_normal(scale = 5),
      sigma ~ half_normal(scale = 5)
    ),
    prior_fixed_sd = 100,
    backend = "inla",
    verbose = FALSE
  ))
  raw <- INLA::inla(
    Reaction ~ 1 + Days + f(Subject, model = "iid", hyper = list(
      prec = list(prior = .inla_halfnormal_sd_expr(5))
    )),
    data = d,
    family = "gaussian",
    control.family = list(hyper = list(
      prec = list(prior = .inla_halfnormal_sd_expr(5))
    )),
    control.fixed = list(mean = 0, prec = 1e-4, mean.intercept = 0,
                         prec.intercept = 1e-4),
    control.compute = list(config = TRUE, return.marginals = TRUE,
                           dic = TRUE, waic = TRUE)
  )
  got <- fit$inla

  # Neither INLA call is bit-reproducible, so the tolerances come from the
  # spread of the two against each other over twenty identical repeats on
  # this data (INLA 26.8.7). Largest relative difference: fixed-effect mean
  # 1.9e-6, fixed-effect sd 9.0e-3, hyperparameter mean 5.2e-3, random effect
  # 1.1e-3, fitted value 4.1e-4, effective parameters 1.8e-3. Largest absolute
  # difference: log marginal likelihood 0.0037, DIC 0.074, WAIC 0.036. Each
  # tolerance sits about three times above its measured largest value; an
  # identity error in the mapping moves these quantities by orders of
  # magnitude more.
  expect_equal(got$summary.fixed$mean, raw$summary.fixed$mean,
               tolerance = 1e-3)
  expect_equal(got$summary.fixed$sd, raw$summary.fixed$sd, tolerance = 0.03)
  expect_equal(got$summary.hyperpar$mean, raw$summary.hyperpar$mean,
               tolerance = 0.02)
  expect_equal(got$summary.random$Subject$mean, raw$summary.random$Subject$mean,
               tolerance = 1e-2)
  expect_equal(got$summary.fitted.values$mean, raw$summary.fitted.values$mean,
               tolerance = 2e-3)
  expect_lt(abs(got$mlik[1L, 1L] - raw$mlik[1L, 1L]), 0.05)
  expect_lt(abs(got$dic$dic - raw$dic$dic), 0.25)
  expect_lt(abs(got$waic$waic - raw$waic$waic), 0.15)
  expect_equal(got$dic$p.eff, raw$dic$p.eff, tolerance = 1e-2)

  # The fit records what it ran on, and INLA's own argument record keeps the
  # standardised coordinates.
  st <- fit$extras$response_standardisation
  expect_equal(st$centre, mean(d$Reaction))
  expect_equal(st$scale, stats::sd(d$Reaction))
  expect_equal(
    mean(fit$inla$.args$data[[.FB_INLA_STANDARDISED_RESPONSE]]), 0,
    tolerance = 1e-12
  )
})


# --- metamorphic invariance ---------------------------------------------- #

test_that("a default-prior fit does not move when the response is rescaled", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  cells <- lapply(c(0.01, 1, 1000), function(k) {
    fit <- suppressWarnings(flexybayes(
      Reaction ~ Days + (1 | Subject), data = .rs_sleepstudy(k),
      backend = "inla", verbose = FALSE
    ))
    g <- .rs_varcomp(fit)
    c(ratio = g[["sd_Subject"]] / g[["sigma"]],
      days = coef(fit)[["Days"]] / k,
      status = fit$inla$mode$mode.status)
  })
  cells <- do.call(rbind, cells)
  expect_true(all(cells[, "status"] == 0))
  expect_lt(diff(range(cells[, "ratio"])) / mean(cells[, "ratio"]), 1e-3)
  expect_lt(diff(range(cells[, "days"])) / mean(cells[, "days"]), 1e-3)
})

test_that("an uncorrelated slope does not move under response or covariate scale", {
  # The slope covariate is not rescaled: on a standardised response the slope
  # field held its ratio to the residual at Days x 1, 1000 and 100000, and a
  # version that divided Days by its standard deviation killed the INLA engine
  # at x 100000. Days shrunk a thousandfold moved the fit with and without that
  # division, so no cell below that scale is asserted here.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  cells <- lapply(list(c(1, 1), c(100, 1), c(1, 1000), c(100, 1000)), function(k) {
    fit <- suppressWarnings(flexybayes(
      Reaction ~ Days + (Days || Subject),
      data = .rs_sleepstudy(k[[1L]], k[[2L]]),
      backend = "inla", verbose = FALSE
    ))
    g <- .rs_varcomp(fit)
    c(slope = g[["sd_Days_Subject"]] * k[[2L]] / g[["sigma"]],
      intercept = g[["sd_Subject"]] / g[["sigma"]])
  })
  cells <- do.call(rbind, cells)
  expect_lt(diff(range(cells[, "slope"])) / mean(cells[, "slope"]), 1e-3)
  expect_lt(diff(range(cells[, "intercept"])) / mean(cells[, "intercept"]), 1e-3)
})


# --- draws ----------------------------------------------------------------- #

test_that("draws are taken on the original scale", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  fit <- suppressWarnings(flexybayes(
    Reaction ~ Days + (1 | Subject), data = .rs_sleepstudy(1000),
    backend = "inla", verbose = FALSE
  ))
  sf <- fit$inla$summary.fixed

  # INLA's sampler reads the standardised internal state; unmapped, these
  # draws would sit near zero with a spread a thousandth of the summary's.
  draws <- fb_as_draws_simple(fit, n_samples = 2000L)
  days <- draws[["Days:1"]]
  expect_lt(abs(mean(days) - sf["Days", "mean"]),
            5 * sf["Days", "sd"] / sqrt(length(days)))
  expect_equal(stats::sd(days), sf["Days", "sd"], tolerance = 0.1)

  sigma <- 1 / sqrt(draws[["Precision for the Gaussian observations"]])
  expect_equal(mean(sigma), .rs_varcomp(fit)[["sigma"]], tolerance = 0.05)

  v <- vcov(fit, n_samples = 2000L)
  expect_equal(sqrt(diag(v)), sf$sd, tolerance = 0.1, ignore_attr = TRUE)

  vp <- vpredict(fit, r ~ sd_Subject / sigma)
  g <- .rs_varcomp(fit)
  expect_equal(vp$estimate, g[["sd_Subject"]] / g[["sigma"]], tolerance = 0.1)
})


# --- the aggregated route and the record ---------------------------------- #

# [removed] "the aggregated route matches the per-row fixed effects at large
# scale" -- an aggregated Gaussian fit now refuses by name, because the
# aggregated Gaussian route cannot carry the residual prior. The comparison
# of its standardised fit with the per-row one returns with that composition.

test_that("only a Gaussian identity-link fit is standardised", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(4L)
  d <- data.frame(g = factor(rep(letters[1:8], each = 10L)))
  d$y <- stats::rpois(nrow(d), 20)
  d$z <- stats::rnorm(nrow(d), 5, 2)

  pois <- suppressMessages(flexybayes(y ~ 1, random = ~ g, data = d,
                                      family = "poisson", backend = "inla",
                                      return_code = TRUE, verbose = FALSE))
  expect_null(pois$response_standardisation)
  expect_match(deparse(pois$formula)[[1L]], "^y ~")

  gauss <- suppressMessages(flexybayes(z ~ 1, random = ~ g, data = d,
                                       backend = "inla", return_code = TRUE,
                                       verbose = FALSE))
  expect_equal(gauss$response_standardisation$scale, stats::sd(d$z))
  expect_match(deparse(gauss$formula)[[1L]],
               paste0("^", .FB_INLA_STANDARDISED_RESPONSE, " ~"))
})
