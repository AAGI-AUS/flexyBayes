# Tests for the family response-domain gate (U-FAM-SUPPORT, 2026-09-16).
#
# Before this fix, family = "lognormal" on a response with values at or
# below zero passed .resolve_family() (b979c5d admitted lognormal to the
# roster) and reached brms::brm(), which stopped with the untyped
# condition "Family 'lognormal' requires response greater than 0.",
# a plain simpleError, not a flexyBayes refusal. gamma, hurdle_gamma and
# beta carry the same restricted-domain likelihoods and reached the same
# unrefused failure, on brms as an untyped stop() and on INLA as an
# unrelated inla_program_failed diagnosis naming neither the family nor
# the offending row.
#
# .fb_refuse_family_response_domain() (R/na_action.R) now checks the
# observed response against .fb_family_response_domain()
# (R/family_traits.R) inside .fb_apply_na_action(), before any engine
# runs and before the plan / review_code / return_code short-circuits,
# so every route through flexybayes() (and its sugar: fb(), fb_brms(),
# fb_inla()) gets the same typed refusal (family_response_domain_violation).

# ---------------------------------------------------------------- #
# 1. Unit tests of the guard function                                #
# ---------------------------------------------------------------- #

test_that(".fb_refuse_family_response_domain() is a no-op for a family with no domain restriction", {
  fb <- list(family = "gaussian", response = "y")
  d <- data.frame(y = c(-5, 0, 5, NA))
  expect_null(flexyBayes:::.fb_refuse_family_response_domain(fb, d))
})

test_that(".fb_refuse_family_response_domain() is a no-op for a compliant response, per family", {
  cases <- list(
    list(family = "lognormal", y = c(0.1, 1, 5, 100)),
    list(family = "gamma", y = c(0.01, 1, 5, 100)),
    list(family = "hurdle_gamma", y = c(0, 0.01, 1, 5)),
    list(family = "beta", y = c(0.001, 0.5, 0.999))
  )
  for (cs in cases) {
    fb <- list(family = cs$family, response = "y")
    d <- data.frame(y = cs$y)
    expect_null(
      flexyBayes:::.fb_refuse_family_response_domain(fb, d),
      label = cs$family
    )
  }
})

test_that(".fb_refuse_family_response_domain() treats a missing response as na_action's business, not a domain violation", {
  fb <- list(family = "lognormal", response = "y")
  d <- data.frame(y = c(1, 2, NA, NA))
  expect_null(flexyBayes:::.fb_refuse_family_response_domain(fb, d))
})

test_that(".fb_refuse_family_response_domain() refuses lognormal on a non-positive response, naming the column and value", {
  fb <- list(family = "lognormal", response = "y")
  d <- data.frame(y = c(1, 2, -0.5, 3))
  err <- tryCatch(
    flexyBayes:::.fb_refuse_family_response_domain(fb, d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
  expect_s3_class(err, "flexybayes_family_response_domain_refusal")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "family = \"lognormal\"", fixed = TRUE)
  expect_match(conditionMessage(err), "`y`", fixed = TRUE)
  expect_match(conditionMessage(err), "row 3: -0.5", fixed = TRUE)
  expect_match(conditionMessage(err), "greater than 0", fixed = TRUE)
  expect_match(conditionMessage(err), "gaussian", fixed = TRUE)
  expect_identical(err$column, "y")
  expect_identical(err$offending_rows, 3L)
})

test_that(".fb_refuse_family_response_domain() refuses gamma at exactly zero, naming hurdle_gamma", {
  fb <- list(family = "gamma", response = "y")
  d <- data.frame(y = c(1, 2, 0, 3))
  err <- tryCatch(
    flexyBayes:::.fb_refuse_family_response_domain(fb, d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
  expect_match(conditionMessage(err), "hurdle_gamma", fixed = TRUE)
})

test_that(".fb_refuse_family_response_domain() refuses gamma on a negative response, naming gaussian (not hurdle_gamma)", {
  fb <- list(family = "gamma", response = "y")
  d <- data.frame(y = c(1, 2, -3, 3))
  err <- tryCatch(
    flexyBayes:::.fb_refuse_family_response_domain(fb, d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
  expect_match(conditionMessage(err), "gaussian", fixed = TRUE)
  expect_false(grepl("hurdle_gamma", conditionMessage(err), fixed = TRUE))
})

test_that(".fb_refuse_family_response_domain() refuses hurdle_gamma only on a negative response, not on zero", {
  fb <- list(family = "hurdle_gamma", response = "y")
  ok <- data.frame(y = c(0, 0, 1, 2))
  expect_null(flexyBayes:::.fb_refuse_family_response_domain(fb, ok))

  bad <- data.frame(y = c(0, -1, 1, 2))
  err <- tryCatch(
    flexyBayes:::.fb_refuse_family_response_domain(fb, bad),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
  expect_identical(err$offending_rows, 2L)
  expect_match(conditionMessage(err), "0 or greater", fixed = TRUE)
})

test_that(".fb_refuse_family_response_domain() refuses beta at either boundary of the open unit interval", {
  fb <- list(family = "beta", response = "y")
  err_zero <- tryCatch(
    flexyBayes:::.fb_refuse_family_response_domain(
      fb,
      data.frame(y = c(0.2, 0.5, 0))
    ),
    error = function(e) e
  )
  err_one <- tryCatch(
    flexyBayes:::.fb_refuse_family_response_domain(
      fb,
      data.frame(y = c(0.2, 0.5, 1))
    ),
    error = function(e) e
  )
  expect_s3_class(
    err_zero,
    "flexybayes_refusal_family_response_domain_violation"
  )
  expect_s3_class(
    err_one,
    "flexybayes_refusal_family_response_domain_violation"
  )
  expect_match(
    conditionMessage(err_zero),
    "strictly between 0 and 1",
    fixed = TRUE
  )
  expect_match(
    conditionMessage(err_one),
    "strictly between 0 and 1",
    fixed = TRUE
  )
})

test_that(".fb_refuse_family_response_domain() finds the violation, ignoring NA rows in between", {
  fb <- list(family = "lognormal", response = "y")
  d <- data.frame(y = c(1, NA, -2, NA, 3))
  err <- tryCatch(
    flexyBayes:::.fb_refuse_family_response_domain(fb, d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
  expect_identical(err$offending_rows, 3L)
})

# ---------------------------------------------------------------- #
# 2. Live entry points, refuses before .dispatch_backend() runs    #
# ---------------------------------------------------------------- #

.domain_case_data <- function(family, seed = 20260916L, n = 20L) {
  set.seed(seed)
  y <- switch(
    family,
    lognormal = c(stats::rlnorm(n - 1L), -0.5),
    gamma = c(stats::rgamma(n - 1L, 2, 1), -1),
    hurdle_gamma = c(stats::rgamma(n - 2L, 2, 1), 0, -1),
    beta = c(stats::rbeta(n - 1L, 2, 3), 1)
  )
  data.frame(y = y, x = stats::rnorm(length(y)))
}

test_that("flexybayes() refuses a domain-violating response before the engine runs, on every listed family", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .dispatch_backend = function(...) stop("the engine was reached")
  )
  for (fam in c("lognormal", "gamma", "hurdle_gamma", "beta")) {
    for (engine in c("brms", "inla")) {
      d <- .domain_case_data(fam)
      err <- tryCatch(
        suppressMessages(flexybayes(
          y ~ x,
          data = d,
          family = fam,
          backend = engine,
          verbose = FALSE
        )),
        error = function(e) e
      )
      expect_s3_class(
        err,
        "flexybayes_refusal_family_response_domain_violation"
      )
    }
  }
})

test_that("fb() refuses a domain-violating lognormal response before the engine runs", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .dispatch_backend = function(...) stop("the engine was reached")
  )
  d <- .domain_case_data("lognormal")
  err <- tryCatch(
    suppressWarnings(suppressMessages(fb(
      y ~ x,
      data = d,
      family = "lognormal",
      backend = "brms",
      verbose = FALSE
    ))),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
})

test_that("fb_brms() refuses a domain-violating lognormal response before the engine runs", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .dispatch_backend = function(...) stop("the engine was reached")
  )
  d <- .domain_case_data("lognormal")
  err <- tryCatch(
    suppressWarnings(suppressMessages(fb_brms(
      y ~ x,
      data = d,
      family = "lognormal",
      verbose = FALSE
    ))),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
})

test_that("fb_inla() refuses a domain-violating lognormal response before the engine runs", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .dispatch_backend = function(...) stop("the engine was reached")
  )
  d <- .domain_case_data("lognormal")
  err <- tryCatch(
    suppressWarnings(suppressMessages(fb_inla(
      y ~ x,
      data = d,
      family = "lognormal",
      verbose = FALSE
    ))),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
})

test_that("flexybayes() with return_code = TRUE also refuses before any code is generated", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .domain_case_data("lognormal")
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ x,
      data = d,
      family = "lognormal",
      backend = "brms",
      return_code = TRUE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_response_domain_violation")
})

test_that("a missing response under na_action = \"augment\" reaches the engine unrefused", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .dispatch_backend = function(...) stop("the engine was reached")
  )
  set.seed(20260916L)
  d <- data.frame(y = stats::rlnorm(20L), x = stats::rnorm(20L))
  d$y[c(3L, 7L)] <- NA
  expect_error(
    suppressMessages(flexybayes(
      y ~ x,
      data = d,
      family = "lognormal",
      backend = "brms",
      na_action = "augment",
      verbose = FALSE
    )),
    "the engine was reached",
    fixed = TRUE
  )
})

# ---------------------------------------------------------------- #
# 3. Out of scope, by construction, and proven so                    #
# ---------------------------------------------------------------- #

test_that("flexybayes_stream() cannot reach the domain check: its family roster excludes every restricted-domain family", {
  set.seed(20260916L)
  d <- data.frame(y = c(stats::rlnorm(19L), -0.5), x = stats::rnorm(20L))
  err <- tryCatch(
    suppressMessages(flexybayes_stream(
      y ~ x,
      source = d,
      family = "lognormal",
      backend = "inla",
      verbose = FALSE
    )),
    error = function(e) e
  )
  # Refused by the existing family-roster gate, never reaching
  # .fb_refuse_family_response_domain(), flexybayes_stream() only
  # aggregates gaussian / binomial / poisson, none of which carries a
  # restricted response domain.
  expect_s3_class(err, "flexybayes_refusal_aggregation_route_unavailable")
  expect_false(inherits(err, "flexybayes_family_response_domain_refusal"))
})

test_that("fb_from_brms() and fb_from_asreml() build the IR without checking the response domain", {
  # Ingest adapters build the intermediate representation only. They
  # never reach an engine on their own (no exported verb fits an fb_terms
  # IR directly), so the response-domain check, like the analogous
  # .refuse_non_binary_binomial_response(), lives downstream, at
  # .fb_apply_na_action() inside flexybayes(), not in the ingest layer.
  d <- data.frame(y = c(stats::rlnorm(19L), -0.5), x = stats::rnorm(20L))

  ir_brms <- fb_from_brms(formula = y ~ x, data = d, family = "lognormal")
  expect_s3_class(ir_brms, "fb_terms")
  expect_identical(ir_brms$family, "lognormal")

  ir_asreml <- fb_from_asreml(fixed = y ~ x, data = d, family = "lognormal")
  expect_s3_class(ir_asreml, "fb_terms")
  expect_identical(ir_asreml$family, "lognormal")
})
