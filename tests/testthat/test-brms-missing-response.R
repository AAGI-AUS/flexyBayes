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
# brms's `mi()` addition term samples the missing values instead, but it is
# Gaussian-only -- `bf(y | mi() ~ ...)` with family = poisson raises
# "Argument 'mi' is not supported for family 'poisson(log)'". So Gaussian gets
# augmentation and everything else gets a refusal, rather than silent deletion.
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

test_that("brms samples a missing Gaussian response rather than dropping it", {
  skip_if_not_installed("brms")
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
