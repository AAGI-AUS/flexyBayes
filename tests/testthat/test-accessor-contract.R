# Tests for the accessor contract (D7 / WP-D, ORCHESTRA fitness audit
# 2026-08-25 section 3.1-3.4). Downstream tools -- the AAGI ORCHESTRA
# workspace among them -- read a fitted object purely through exported
# accessors, never through engine-specific internals, so the shape of
# those accessors is itself a contract worth pinning: column names,
# column types, and (separately) the invariant that every registered
# refusal condition's class vector carries a marker ending in `_refusal`
# (the shared gate predicate `is_orchestra_decline()` keys on exactly
# this pattern, `ORCHESTRA_dev/integration/refusal_contract.R:99`).
#
# One INLA fit and one brms fit, each built once and read from several
# angles, per the item's own scope. Live-engine tests skip cleanly when
# the engine isn't installed and never run on CRAN.

# ---------------------------------------------------------------- #
# On an INLA fit                                                   #
# ---------------------------------------------------------------- #

test_that("accessor contract: pinned column names/types on an INLA fit", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  set.seed(20260826L)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- suppressMessages(flexybayes(
    y ~ x + (1 | g),
    data = d,
    backend = "inla",
    seed = 99L,
    verbose = FALSE
  ))

  s <- summary(fit)

  # summary(fit)$varcomp
  expect_named(
    s$varcomp,
    c("component", "estimate", "std.error", "conf.low", "conf.high",
      "prior", "note"),
    ignore.order = FALSE
  )
  expect_type(s$varcomp$component, "character")
  expect_type(s$varcomp$estimate, "double")
  expect_type(s$varcomp$std.error, "double")
  expect_type(s$varcomp$conf.low, "double")
  expect_type(s$varcomp$conf.high, "double")
  expect_type(s$varcomp$prior, "character")
  expect_true(nrow(s$varcomp) >= 1L)

  # summary(fit)$fixed
  expect_named(
    s$fixed,
    c("term", "estimate", "std.error", "conf.low", "conf.high"),
    ignore.order = FALSE
  )
  expect_type(s$fixed$term, "character")
  expect_type(s$fixed$estimate, "double")
  expect_type(s$fixed$std.error, "double")
  expect_identical(nrow(s$fixed), 2L) # (Intercept), x

  # predict(fit, classify = )
  pr <- predict(fit, classify = "x")
  expect_s3_class(pr, "fb_predict_classify")
  expect_named(
    pr,
    c("x", "estimate", "std.error", "conf.low", "conf.high"),
    ignore.order = FALSE
  )
  expect_true(all(vapply(pr, is.numeric, logical(1))))

  # the backend-decision accessor
  bd <- backend_decision(fit)
  expect_named(
    bd,
    c("backend", "path", "gate_checks", "reason", "preflight_summary",
      "representation_plan", "rejected_routes", "routing_policy_version"),
    ignore.order = FALSE
  )
  expect_identical(bd$backend, "inla")
  expect_type(bd$path, "character")
  expect_type(bd$reason, "character")
  expect_type(bd$routing_policy_version, "character")

  # the seed
  expect_identical(fit$extras$call_info$seed, 99L)
})

# ---------------------------------------------------------------- #
# On a brms fit                                                    #
# ---------------------------------------------------------------- #

test_that("accessor contract: pinned column names/types on a brms fit", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(20260826L)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ x + (1 | g),
    data = d,
    backend = "brms",
    n_samples = 300,
    warmup = 300,
    chains = 1,
    seed = 123L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )))

  s <- summary(fit)

  # summary(fit)$varcomp -- identical shape to the INLA fit above; the
  # accessor contract is that a caller need not know which engine ran.
  expect_named(
    s$varcomp,
    c("component", "estimate", "std.error", "conf.low", "conf.high",
      "prior", "note"),
    ignore.order = FALSE
  )
  expect_type(s$varcomp$component, "character")
  expect_type(s$varcomp$estimate, "double")
  expect_type(s$varcomp$conf.low, "double")
  expect_type(s$varcomp$conf.high, "double")

  # summary(fit)$fixed
  expect_named(
    s$fixed,
    c("term", "estimate", "std.error", "conf.low", "conf.high"),
    ignore.order = FALSE
  )
  expect_identical(nrow(s$fixed), 2L)

  # predict(fit, classify = )
  pr <- predict(fit, classify = "x")
  expect_s3_class(pr, "fb_predict_classify")
  expect_named(
    pr,
    c("x", "estimate", "std.error", "conf.low", "conf.high"),
    ignore.order = FALSE
  )
  expect_true(all(vapply(pr, is.numeric, logical(1))))

  # the backend-decision accessor -- same field set as the INLA fit
  bd <- backend_decision(fit)
  expect_named(
    bd,
    c("backend", "path", "gate_checks", "reason", "preflight_summary",
      "representation_plan", "rejected_routes", "routing_policy_version"),
    ignore.order = FALSE
  )
  expect_identical(bd$backend, "brms")
  expect_type(bd$path, "character")

  # the seed
  expect_identical(fit$extras$call_info$seed, 123L)
})

# ---------------------------------------------------------------- #
# Every registered refusal condition's class ends in `_refusal`    #
# ---------------------------------------------------------------- #

test_that("every registered refusal condition carries a class ending in '_refusal'", {
  # No engine needed: exercises the shared constructor
  # .fb_refusal_condition() across the whole registry, which is what
  # every typed refusal in the package is built through. This is the
  # gate-facing invariant ORCHESTRA's is_orchestra_decline() keys on
  # (a class matching `_(refusal|abstention)$`).
  codes <- ls(flexyBayes:::.refusal_registry)
  expect_true(length(codes) > 0L)

  ends_in_refusal <- function(classes) any(grepl("_refusal$", classes))

  offenders <- character(0)
  for (code in codes) {
    cond <- flexyBayes:::.fb_refusal_condition(
      reason_code = code,
      message = "accessor-contract probe"
    )
    if (!ends_in_refusal(class(cond))) {
      offenders <- c(offenders, code)
    }
  }
  expect_identical(offenders, character(0))
})
