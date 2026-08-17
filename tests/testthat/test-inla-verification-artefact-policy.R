# test-inla-verification-artefact-policy.R -- an unshipped file may not
# decide what the package fits.
#
# Until 0.9.0 the INLA gate for the factor:numeric indexed interaction
# consulted an .rds under inst/extdata/inla-verification/. That directory
# is untracked and was not build-excluded, so a locally built tarball
# carried a file a clean clone of the same commit did not, and the two
# fitted different models. The 2026-08-16 adversarial review found it as
# P0-2.
#
# The policy since: the artefacts are a developer rehearsal hook. They
# are excluded from the build, and they are read only when
# `flexyBayes.dev_inla_verification_artefacts` is switched on by hand.
# The tests below prove that a PASSING artefact on disk changes nothing
# with the option at its default.

.fni_fixture <- function(n = 120L, seed = 20260816L) {
  set.seed(seed)
  d <- data.frame(
    x = stats::rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n))
  )
  d$y <- 1 + 0.5 * d$x + as.numeric(d$f) * 0.3 * d$x + stats::rnorm(n, 0, 0.3)
  d
}

# An IR carrying exactly one factor_numeric_interaction term, so the
# term-class-specific gate check is the thing under test.
.fni_ir <- function() {
  flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    fixed_terms = list(
      list(
        type = "factor",
        var = "f",
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f"
      ),
      list(type = "continuous", var = "x", label = "x"),
      list(
        type = "factor_numeric_interaction",
        factor = "f",
        continuous = "x",
        vars = c("f", "x"),
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f:x"
      )
    ),
    random_terms = list(),
    residual_terms = list(list(type = "units")),
    priors = list(legacy = TRUE),
    source = "brms"
  )
}

# Run the gate check against a temporary artefact that reports
# `pass = TRUE`, at a given setting of the developer option. The path
# helper is stubbed rather than the package tree written to, so the probe
# is identical whether the suite runs from source or from an install.
# Returns a list of the two things under test: what the artefact reader
# saw, and what the gate decided.
.probe_with_passing_artefact <- function(dev_enabled) {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "factor_numeric_interaction.rds")
  saveRDS(list(pass = TRUE, fixture = "test double"), path)
  testthat::local_mocked_bindings(
    .factor_numeric_interaction_verification_path = function() path,
    .package = "flexyBayes"
  )
  withr::local_options(
    list(flexyBayes.dev_inla_verification_artefacts = dev_enabled)
  )
  list(
    artefact_reads_pass =
      flexyBayes:::.factor_numeric_interaction_inla_verified(),
    gate =
      flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(
        .fni_ir()
      )
  )
}


# ---------------------------------------------------------------- #
# (1) The gate refuses regardless of what is on disk                 #
# ---------------------------------------------------------------- #

test_that("the gate refuses factor:numeric on INLA with no artefact", {
  fb <- .fni_ir()
  withr::local_options(
    list(flexyBayes.dev_inla_verification_artefacts = NULL)
  )
  res <- flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(fb)
  expect_false(isTRUE(res$pass))
})

test_that("a passing artefact does not lift the refusal at the default", {
  probe <- .probe_with_passing_artefact(dev_enabled = NULL)
  # The artefact itself reads as passing ...
  expect_true(probe$artefact_reads_pass)
  # ... and the gate refuses anyway, because shipped behaviour never
  # consults it.
  expect_false(isTRUE(probe$gate$pass))
})

test_that("the developer option is the only thing that admits the mapping", {
  probe <- .probe_with_passing_artefact(dev_enabled = TRUE)
  expect_true(probe$artefact_reads_pass)
  expect_true(isTRUE(probe$gate$pass))
})

test_that("the gate is a trivial pass without a factor:numeric term", {
  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    fixed_terms = list(list(type = "continuous", var = "x", label = "x")),
    random_terms = list(),
    residual_terms = list(list(type = "units")),
    priors = list(legacy = TRUE),
    source = "brms"
  )
  res <- flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(fb)
  expect_true(isTRUE(res$pass))
})


# ---------------------------------------------------------------- #
# (2) The refusal message names brms, never greta                    #
# ---------------------------------------------------------------- #

test_that("the factor:numeric refusal recommends brms and not greta", {
  skip_if_not_installed("INLA")
  d <- .fni_fixture()
  withr::local_options(
    list(flexyBayes.dev_inla_verification_artefacts = NULL)
  )
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ f * x,
      data = d,
      backend = "inla",
      return_code = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal")
  expect_s3_class(
    err,
    "flexybayes_lgm_factor_numeric_interaction_inla_verified"
  )
  msg <- conditionMessage(err)
  expect_match(msg, "backend = \"brms\"", fixed = TRUE)
  expect_false(grepl("backend = \"greta\"", msg, fixed = TRUE))
})

test_that("the same model fits on brms", {
  skip_if_not_installed("brms")
  d <- .fni_fixture()
  emitted <- tryCatch(
    suppressMessages(flexybayes(
      y ~ f * x,
      data = d,
      backend = "brms",
      return_code = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_false(inherits(emitted, "error"))
})


# ---------------------------------------------------------------- #
# (3) The (x || g) emit follows the same policy                      #
# ---------------------------------------------------------------- #

test_that("the (x || g) INLA deferral is host-independent", {
  withr::local_options(
    list(flexyBayes.dev_inla_verification_artefacts = NULL)
  )
  err <- tryCatch(
    flexyBayes:::.check_inla_verification_simple_slope_uncor(),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_inla_simple_slope_uncor_deferred")
  expect_match(
    conditionMessage(err),
    "host-independent",
    fixed = TRUE
  )
})


# ---------------------------------------------------------------- #
# (4) The artefact directory is excluded from the build              #
# ---------------------------------------------------------------- #

test_that(".Rbuildignore excludes the verification artefact directory", {
  path <- testthat::test_path("..", "..", ".Rbuildignore")
  skip_if(
    !file.exists(path),
    ".Rbuildignore is not in the installed tree"
  )
  patterns <- readLines(path, warn = FALSE)
  expect_true("^inst/extdata/inla-verification/" %in% patterns)
})
