# =============================================================================
# What update() carries, and what it refuses to guess.
#
# update() re-issues every recorded argument to flexybayes(). An argument
# the fit never recorded would come back as its default, a different model
# under the same name. Every route records the full set, including the
# missing-response policy. A fit whose record is short is refused by name,
# not re-fitted with defaults.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


.uc_data <- function(seed = 5150L, n = 60L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b"), length.out = n)),
    g = factor(rep(letters[1:5], length.out = n)),
    # A second grouping factor the model does not use, so a re-fit has a
    # term to ADD that is not already the fixed effect. `~ g + f` would
    # put `f` on both sides of `y ~ f`, which is refused at plan time.
    h = factor(rep(LETTERS[1:4], length.out = n))
  )
  b_g <- stats::rnorm(5L, sd = 0.8)[as.integer(d$g)]
  d$y <- 1 + 0.6 * (d$f == "b") + b_g + stats::rnorm(n, sd = 0.5)
  d
}

.uc_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

# The nineteen fields a re-fit needs, in the order the per-row emits write
# them. The last three are tested in test-update-engine-carry.R.
.UC_CALL_FIELDS <- c(
  "fixed", "random", "residual", "data_name", "family", "link",
  "known_matrices", "weights", "n_samples", "warmup", "chains", "seed",
  "control", "prior_fixed_sd", "prior_vc_sd", "na_action", "backend",
  "aggregate", "verbose"
)


# ---------------------------------------------------------------- #
# 1. A default Gaussian call (F-1)                                  #
# ---------------------------------------------------------------- #
#
# The aggregated Gaussian route cannot carry the residual prior, so a
# default Gaussian call is fitted per row. The test that compared the
# aggregated Gaussian record with the per-row one would now compare a
# per-row record with itself. The aggregated record is compared on the
# count route in section 4.

test_that("a default Gaussian call is fitted per row, with the full record", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_data(), verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_false(inherits(fit, "flexybayes_aggregated"))
  expect_identical(names(fit$extras$call_info), .UC_CALL_FIELDS)
})

test_that("update() re-fits a default Gaussian fit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  d <- .uc_data()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = d, verbose = FALSE, mcmc_verbose = FALSE
  ))

  # y carries no h effect, so the added sd_h is weakly identified.
  refit <- withr::with_options(
    list(flexyBayes.silence_weak_identification_warning = TRUE),
    suppressMessages(stats::update(fit, random = ~ g + h))
  )
  expect_s3_class(refit, "flexybayes")
  expect_identical(refit$extras$call_info$random, ~ g + h)
  # The arguments the record grew are carried, not re-defaulted.
  expect_identical(refit$extras$call_info$na_action, "augment")
  expect_identical(
    refit$extras$call_info$known_matrices,
    fit$extras$call_info$known_matrices
  )
  expect_identical(names(refit$extras$call_info), .UC_CALL_FIELDS)
})


# ---------------------------------------------------------------- #
# 2. The missing-response policy (F-2)                              #
# ---------------------------------------------------------------- #

test_that("an omit fit re-fits as omit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  d <- .uc_data()
  d$y[c(4L, 22L)] <- NA
  # A grouping factor that is not a fixed effect, so the override adds a
  # well-posed random term. A factor on both sides of the model is
  # confounded, a case tested in test-summary-asreml-shape.R.
  d$h <- factor(rep(paste0("H", 1:4), length.out = nrow(d)))

  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = d, backend = "inla", aggregate = FALSE,
    na_action = "omit", verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(fit$extras$call_info$na_action, "omit")
  expect_identical(fit$extras$na_action$na_action, "omit")
  expect_identical(stats::nobs(fit), 58L)

  # y carries no h effect, so the added sd_h is weakly identified.
  refit <- withr::with_options(
    list(flexyBayes.silence_weak_identification_warning = TRUE),
    suppressMessages(stats::update(fit, random = ~ g + h))
  )
  expect_identical(refit$extras$call_info$na_action, "omit")
  expect_identical(refit$extras$na_action$na_action, "omit")
  # The re-fit is computed on the same rows, not on the augmented design.
  expect_identical(stats::nobs(refit), 58L)
})

test_that("an augment fit re-fits as augment", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  d <- .uc_data()
  d$y[c(4L, 22L)] <- NA

  # A well-posed added term, for the reason given in the omit test above.
  d$h <- factor(rep(paste0("H", 1:4), length.out = nrow(d)))

  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = d, backend = "inla", aggregate = FALSE,
    na_action = "augment", verbose = FALSE, mcmc_verbose = FALSE
  ))
  # y carries no h effect, so the added sd_h is weakly identified.
  refit <- withr::with_options(
    list(flexyBayes.silence_weak_identification_warning = TRUE),
    suppressMessages(stats::update(fit, random = ~ g + h))
  )
  expect_identical(refit$extras$call_info$na_action, "augment")
  expect_identical(stats::nobs(refit), 60L)
  expect_identical(as.integer(stats::nobs(refit, type = "observed")), 58L)
})


# ---------------------------------------------------------------- #
# 3. The refusal stays reachable                                    #
# ---------------------------------------------------------------- #

test_that("a record with no recorded policy refuses by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_data(), backend = "inla",
    aggregate = FALSE, verbose = FALSE, mcmc_verbose = FALSE
  ))
  fit$extras$call_info$na_action <- NULL

  err <- tryCatch(
    stats::update(fit, random = ~ g + h),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_update_call_not_reconstructable")
  expect_match(conditionMessage(err), "na_action", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# 4. The aggregated count route                                     #
# ---------------------------------------------------------------- #
#
# An ordinary Poisson or binomial call auto-routes to the aggregation layer,
# which records the same fields as the per-row route.

.uc_count_data <- function(seed = 4242L, n = 90L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    g = factor(rep(letters[1:6], length.out = n)),
    # Crossed with `f`, so the re-fit can add it without aliasing `f`. With
    # `f`'s own period, `h` relabelled `f`: its variance was unidentified,
    # its mode sat on the default bound, and INLA put the group precision
    # outside that bound in 8 of 24 re-fits.
    h = factor(rep(LETTERS[1:3], each = 3L, length.out = n))
  )
  b_g <- stats::rnorm(6L, sd = 0.3)[as.integer(d$g)]
  eta <- 1.2 + 0.4 * (d$f == "b") - 0.3 * (d$f == "c") + b_g
  d$y <- stats::rpois(n, exp(eta))
  d
}

test_that("the count fixture crosses f and h", {
  d <- .uc_count_data()

  expect_true(all(table(d$f, d$h) > 0L))
})

test_that("a default Poisson call auto-routes to the aggregated emit", {
  # The premise of the next two tests.
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_count_data(), family = "poisson",
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
})

test_that("the aggregated count record carries the same fields", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_count_data(), family = "poisson",
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(names(fit$extras$call_info), .UC_CALL_FIELDS)
})

test_that("update() re-fits a default aggregated Poisson fit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_count_data(), family = "poisson",
    verbose = FALSE, mcmc_verbose = FALSE
  ))

  # y carries no h effect, so the added sd_h is weakly identified.
  refit <- withr::with_options(
    list(flexyBayes.silence_weak_identification_warning = TRUE),
    suppressMessages(stats::update(fit, random = ~ g + h))
  )
  expect_s3_class(refit, "flexybayes")
  expect_identical(refit$extras$call_info$random, ~ g + h)
  expect_identical(refit$extras$call_info$family, "poisson")
  expect_identical(refit$extras$call_info$na_action, "augment")
  expect_identical(names(refit$extras$call_info), .UC_CALL_FIELDS)
})

test_that("a short count record still refuses by name", {
  # The refusal is the correct answer for a record that genuinely lacks
  # a field, and completing the emit must not special-case an older
  # object into a silent default.
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_count_data(), family = "poisson",
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  fit$extras$call_info$known_matrices <- NULL

  err <- tryCatch(
    stats::update(fit, random = ~ g + h),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_update_call_not_reconstructable")
  expect_match(conditionMessage(err), "known_matrices", fixed = TRUE)
})

# =============================================================================
# The unnamed stats::update() idiom, `update(fit, . ~ . + z)`, is refused
# by name. Overrides are matched to the recorded arguments by name, so an
# unnamed one would otherwise be dropped without a word.
# =============================================================================

test_that("update() refuses the unnamed stats::update() replacement-formula idiom", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, data = .uc_data(), backend = "inla", verbose = FALSE
  ))
  orig_fixed <- fit$extras$call_info$fixed

  err <- tryCatch(
    stats::update(fit, . ~ . + g),
    condition = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_update_unnamed_argument_not_supported"
  )
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "update(fit, . ~ . + z)", fixed = TRUE)
  expect_match(conditionMessage(err), "fixed = ", fixed = TRUE)

  # An identity re-fit keeps the recorded formula.
  identity_refit <- suppressMessages(stats::update(fit))
  expect_identical(identity_refit$extras$call_info$fixed, orig_fixed)
})

test_that("update() with a mix of one unnamed and one named argument still refuses", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, data = .uc_data(), backend = "inla", verbose = FALSE
  ))
  err <- tryCatch(
    stats::update(fit, . ~ . + g, verbose = FALSE),
    condition = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_update_unnamed_argument_not_supported"
  )
})

test_that("update() with only named arguments is unaffected by the unnamed-argument guard", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, data = .uc_data(), backend = "inla", verbose = FALSE
  ))
  # In this fixture `hB + hD` equals `fb`, so `f + h` needs fixed_alias.
  refit <- suppressMessages(
    stats::update(fit, fixed = y ~ f + h, fixed_alias = "drop")
  )
  expect_s3_class(refit, "flexybayes")
  expect_identical(refit$extras$fixed_alias$dropped$column, "hB")
  expect_identical(refit$extras$call_info$fixed, y ~ f + h)
  expect_true(any(grepl("^h", summary(refit)$fixed$term)))
})
