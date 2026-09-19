# =============================================================================
# A missing response on the brms backend.
#
# brms drops rows whose response is NA. Silently: a 48-row dataset with six
# missing responses reached Stan as N = 42, with a warning easy to miss and no
# trace in the returned object. That made `na_action = "augment"` on brms mean
# complete-case deletion -- the argument promising the design is preserved,
# delivering the opposite.
#
# The tests below assert on the NUMBER OF ROWS THE ENGINE ACTUALLY USED,
# because that is where the defect lived. Every fit returned a valid object of
# the right class with sensible estimates throughout; nothing about the wrapper
# would have revealed it.
#
# brms's `mi()` addition term samples the missing values instead. brms accepts
# it for more families than are admitted here, though not for poisson or
# bernoulli ("Argument 'mi' is not supported for family 'poisson(log)'"). The
# gaussian and lognormal families get augmentation, each checked by a live fit,
# and every other family gets a refusal rather than silent deletion.
# =============================================================================

.brms_miss_data <- function(family = "gaussian", n_missing = 6L, seed = 1L) {
  set.seed(seed)
  d <- expand.grid(col_num = 1:6, row_num = 1:8)
  d$blk <- factor(rep(1:6, length.out = nrow(d)))
  b <- stats::rnorm(6L, 0, 0.5)
  eta <- 0.5 + b[as.integer(d$blk)]
  d$y <- switch(
    family,
    gaussian = eta + stats::rnorm(nrow(d)),
    poisson = stats::rpois(nrow(d), exp(eta)),
    binomial = stats::rbinom(nrow(d), 1L, stats::plogis(eta))
  )
  d$y[sample.int(nrow(d), n_missing)] <- NA
  d
}

# The lognormal fixture. Eight blocks of twelve rather than six of eight, so
# the block standard deviation is estimable, with ten of the 96 responses
# missing completely at random.
.brms_miss_lognormal_data <- function(seed = 11L) {
  set.seed(seed)
  d <- expand.grid(col_num = 1:8, row_num = 1:12)
  d$blk <- factor(rep(1:8, length.out = nrow(d)))
  b <- stats::rnorm(8L, 0, 0.4)
  d$y <- exp(1 + b[as.integer(d$blk)] + stats::rnorm(nrow(d), 0, 0.5))
  d$y[sample.int(nrow(d), 10L)] <- NA
  d
}

test_that("brms samples a missing Gaussian response rather than dropping it", {
  skip_if_not_installed("brms")
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_data("gaussian")
  fit <- suppressMessages(suppressWarnings(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    backend = "brms", na_action = "augment",
    n_samples = 200, warmup = 200, chains = 2,
    verbose = FALSE, mcmc_verbose = FALSE
  )))
  # The whole point: every row reaches the engine, including the six whose
  # response is missing. Before the fix this was 42.
  expect_equal(stats::nobs(fit$brms), nrow(d))
})

test_that("the sampled missing responses reach summary(fit)$missing", {
  # The device is only half-delivered if the fit carries the posterior and
  # no accessor reports it. One row per design cell that lost its response,
  # which is exactly the gap between the two observation counts.
  skip_if_not_installed("brms")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_data("gaussian")
  fit <- suppressMessages(suppressWarnings(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    backend = "brms", na_action = "augment",
    n_samples = 200, warmup = 200, chains = 2,
    verbose = FALSE, mcmc_verbose = FALSE
  )))

  captured <- utils::capture.output(print(s <- summary(fit)))
  mv <- s$missing
  expect_s3_class(mv, "data.frame")
  expect_identical(
    names(mv),
    c("row", "estimate", "std.error", "conf.low", "conf.high")
  )
  expect_identical(
    nrow(mv),
    stats::nobs(fit) - stats::nobs(fit, type = "observed")
  )
  # The rows named are the rows whose response was missing, and the
  # accessor and the summary slot are one construction.
  expect_identical(sort(mv$row), sort(which(is.na(d$y))))
  expect_identical(stats::coef(fit, what = "missing"), mv)
  expect_true(all(is.finite(mv$estimate)))
  expect_true(all(mv$std.error > 0))
  expect_true(all(mv$conf.low <= mv$conf.high))
})

test_that("a complete brms fit reports a typed zero-row missing table", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_data("gaussian", n_missing = 0L)
  fit <- suppressMessages(suppressWarnings(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    backend = "brms", na_action = "augment",
    n_samples = 200, warmup = 200, chains = 2,
    verbose = FALSE, mcmc_verbose = FALSE
  )))
  mv <- stats::coef(fit, what = "missing")
  expect_s3_class(mv, "data.frame")
  expect_identical(nrow(mv), 0L)
  expect_true(is.integer(mv$row))
})

test_that("the emitted brms formula marks the response as partially missing", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_data("gaussian")
  code <- suppressMessages(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    backend = "brms", na_action = "augment", return_code = TRUE
  ))
  txt <- paste(as.character(code), collapse = "\n")
  # brms generates a Yl parameter vector for the missing responses.
  expect_match(txt, "Ymi|Jmi|Yl")
})

test_that("a complete response does not get the mi() treatment", {
  # The guard must not alter models that were never affected.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_data("gaussian", n_missing = 0L)
  code <- suppressMessages(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    backend = "brms", na_action = "augment", return_code = TRUE
  ))
  txt <- paste(as.character(code), collapse = "\n")
  expect_false(grepl("Ymi|Jmi", txt))
})

test_that("a missing non-Gaussian response is refused, not deleted", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  for (fam in c("poisson", "binomial")) {
    expect_error(
      suppressMessages(flexybayes(
        fixed = y ~ 1, random = ~ blk, residual = ~ units,
        data = .brms_miss_data(fam), family = fam,
        backend = "brms", na_action = "augment", return_code = TRUE
      )),
      class = "flexybayes_brms_cannot_augment_nongaussian",
      label = fam
    )
  }
})

test_that("a missing lognormal response is sampled and matches deletion", {
  # Under ignorability the parameter posterior is the same whether the
  # missing cells are augmented or omitted, so the complete-case fit is the
  # oracle for the augmented one. The identity cannot see a route that
  # deletes under the augmentation name, since deletion targets the same
  # posterior, so the carried rows are asserted first.
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_lognormal_data()
  fit_one <- function(na_action) {
    suppressMessages(flexybayes(
      y ~ 1 + (1 | blk), data = d, family = brms::lognormal(),
      # brms centres its default intercept prior on the response median, and
      # at zero once mi() is in the formula. Stating the prior keeps the two
      # fits apart in their rows only.
      prior = fb_prior(b("(Intercept)") ~ normal(1, 2.5)),
      backend = "brms", na_action = na_action,
      n_samples = 1000, warmup = 1000, chains = 2, seed = 1001L,
      # Stan's default of 0.8 left occasional divergent transitions on this
      # design, in both representations.
      control = list(adapt_delta = 0.95),
      verbose = FALSE, mcmc_verbose = FALSE
    ))
  }
  fit_aug <- expect_no_warning(fit_one("augment"))
  fit_omit <- expect_no_warning(fit_one("omit"))

  expect_equal(stats::nobs(fit_aug$brms), nrow(d))
  mv <- summary(fit_aug)$missing
  expect_identical(nrow(mv), sum(is.na(d$y)))
  expect_identical(sort(mv$row), sort(which(is.na(d$y))))
  # The imputed values are draws of the response, so they lie on its support.
  expect_true(all(mv$conf.low > 0))

  par_vec <- c("b_Intercept", "sd_blk__Intercept", "sigma")
  summ <- function(fit) {
    posterior::summarise_draws(
      posterior::as_draws_array(fit$brms, variable = par_vec),
      "mean",
      "sd"
    )
  }
  s_aug <- summ(fit_aug)
  s_omit <- summ(fit_omit)
  expect_identical(s_aug$variable, s_omit$variable)
  # Four Monte Carlo units each: the spread of that summary over eight seeds
  # of each representation at these settings, combined across the two fits.
  tol_mean <- c(b_Intercept = 0.0575, sd_blk__Intercept = 0.0455,
    sigma = 0.00583)
  tol_sd <- c(b_Intercept = 0.0435, sd_blk__Intercept = 0.0560,
    sigma = 0.00488)
  for (i in seq_along(s_aug$variable)) {
    p <- s_aug$variable[[i]]
    expect_lte(abs(s_aug$mean[[i]] - s_omit$mean[[i]]), tol_mean[[p]],
      label = paste("posterior mean difference,", p))
    expect_lte(abs(s_aug$sd[[i]] - s_omit$sd[[i]]), tol_sd[[p]],
      label = paste("posterior sd difference,", p))
  }
})

test_that("a missing Student-t response is refused by name on brms", {
  # brms accepts mi() for student(), but emit_brms() has no translation for
  # it: admitting the family here ended in an untyped emitter error. The
  # brms::student() family object now stops earlier, at the family gate, as
  # the string does. `t`, the INLA spelling, reaches no brms route and is
  # refused by the augmentation gate itself.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_lognormal_data()
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1 + (1 | blk), data = d, family = brms::student(),
      backend = "brms", na_action = "augment", return_code = TRUE
    )),
    class = "flexybayes_refusal_unsupported_family"
  )
  expect_error(
    .fb_brms_missing_response_mode(list(response = "y", family = "t"), d),
    class = "flexybayes_refusal_brms_cannot_augment_nongaussian"
  )
})

test_that("na_action = 'omit' still deletes on brms, for any family", {
  # Deletion remains available and correct; what is refused is deletion
  # wearing the name of augmentation.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .brms_miss_data("poisson")
  expect_no_error(suppressMessages(flexybayes(
    fixed = y ~ 1, random = ~ blk, residual = ~ units, data = d,
    family = "poisson", backend = "brms", na_action = "omit",
    return_code = TRUE
  )))
})

test_that("INLA carries a missing response for every family", {
  # The route the refusal points at has to work, or the message is unhelpful.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  for (fam in c("gaussian", "poisson", "binomial")) {
    expect_no_error(
      suppressMessages(suppressWarnings(flexybayes(
        fixed = y ~ 1, random = ~ blk, residual = ~ units,
        data = .brms_miss_data(fam), family = fam,
        backend = "inla", na_action = "augment", verbose = FALSE
      ))),
      message = fam
    )
  }
})
