# =============================================================================
# What update() carries, and what it refuses to guess.
#
# update() re-issues every recorded argument to flexybayes(). An argument
# the fit never recorded comes back as its default, so a short record means
# a re-fit that is a different model under the same name. Two gaps closed
# here:
#
#   * the aggregated Gaussian route recorded eleven of the sixteen fields,
#     so an ordinary `flexybayes(y ~ f, random = ~ g, data = d)` -- which
#     auto-routes there -- refused update() outright;
#   * no route recorded the missing-response policy, so an `omit` fit would
#     have re-fitted as the default `augment`, on a different set of rows
#     from the one it was named after.
#
# The refusal has to stay reachable. A fit whose record is genuinely short
# is refused by name rather than re-fitted with defaults filled in.
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
# them. Read off emit_inla() rather than restated, so a change to one
# record is a change to the contract all three are held to. The last
# three joined at 0.9.1: without them a re-fit took the formal defaults
# for the engine, the representation and the reporting, so an identity
# update() could come back on another engine (see
# test-update-engine-carry.R).
.UC_CALL_FIELDS <- c(
  "fixed", "random", "residual", "data_name", "family", "link",
  "known_matrices", "weights", "n_samples", "warmup", "chains", "seed",
  "control", "prior_fixed_sd", "prior_vc_sd", "na_action", "backend",
  "aggregate", "verbose"
)


# ---------------------------------------------------------------- #
# 1. The aggregated route (F-1)                                     #
# ---------------------------------------------------------------- #

test_that("a default Gaussian call auto-routes to the aggregated emit", {
  # The premise of the next test. If this stops holding, the update() test
  # below stops exercising the aggregated record and silently passes on
  # the per-row one.
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_data(), verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
})

test_that("the aggregated record carries the same fields as the others", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_data(), verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_identical(names(fit$extras$call_info), .UC_CALL_FIELDS)

  per_row <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .uc_data(), backend = "inla",
    aggregate = FALSE, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(
    names(fit$extras$call_info),
    names(per_row$extras$call_info)
  )
})

test_that("update() re-fits a default aggregated Gaussian fit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .uc_silence()
  d <- .uc_data()
  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = d, verbose = FALSE, mcmc_verbose = FALSE
  ))

  refit <- suppressMessages(stats::update(fit, random = ~ g + h))
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
  # A grouping factor that is not already a fixed effect, so the override
  # below adds a well-posed random term. This test's subject is the
  # missing-response policy travelling with the re-fit, and the term the
  # override adds is incidental to that -- it only has to be a term. It
  # used to be `f`, which is also the fixed effect, and a factor on both
  # sides of the model is confounded: once the prior policy started
  # pricing the added term's SD, INLA's marginals for it came back
  # degenerate and the varcomp table could not be built. That behaviour
  # is correct and is pinned on purpose in
  # test-summary-asreml-shape.R, not carried here as a side effect.
  d$h <- factor(rep(paste0("H", 1:4), length.out = nrow(d)))

  fit <- suppressMessages(flexybayes(
    y ~ f, random = ~g, data = d, backend = "inla", aggregate = FALSE,
    na_action = "omit", verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(fit$extras$call_info$na_action, "omit")
  expect_identical(fit$extras$na_action$na_action, "omit")
  expect_identical(stats::nobs(fit), 58L)

  refit <- suppressMessages(stats::update(fit, random = ~ g + h))
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
  refit <- suppressMessages(stats::update(fit, random = ~ g + h))
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
# 4. The aggregated COUNT route (F-3)                               #
# ---------------------------------------------------------------- #
#
# The Gaussian aggregated record was completed in Phase 2; the count one
# was not, so an ordinary Poisson or binomial call -- which auto-routes
# to the same aggregation layer -- still refused update(). The repair is
# the one applied there: five formals, five recorded fields, and the same
# five passed through dispatch.

.uc_count_data <- function(seed = 4242L, n = 90L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    g = factor(rep(letters[1:6], length.out = n)),
    # See .uc_data(): a term the re-fit can add without aliasing `f`.
    h = factor(rep(LETTERS[1:3], length.out = n))
  )
  b_g <- stats::rnorm(6L, sd = 0.3)[as.integer(d$g)]
  eta <- 1.2 + 0.4 * (d$f == "b") - 0.3 * (d$f == "c") + b_g
  d$y <- stats::rpois(n, exp(eta))
  d
}

test_that("a default Poisson call auto-routes to the aggregated emit", {
  # The premise of the next two tests, asserted rather than assumed.
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

  refit <- suppressMessages(stats::update(fit, random = ~ g + h))
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
