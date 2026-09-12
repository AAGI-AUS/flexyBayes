# =============================================================================
# fb_decision(): the decision rule, validated at construction.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 3.3 and 6.
#
# The whole reason this object exists rather than four arguments on
# fb_select() is that its parts are only jointly meaningful, so the tests are
# mostly about the combinations that are NOT meaningful. Each refusal is
# fired on the input that should produce it, and each is checked to be silent
# on the input that should not -- a validator that refused everything would
# pass a file that only ever fed it bad input.
# =============================================================================


test_that("the default is a validated projection on elpd at one standard error", {
  d <- fb_decision()
  expect_s3_class(d, "fb_decision")
  expect_identical(d$method, "projection")
  expect_identical(d$metric, "elpd")
  expect_identical(d$rule, "one_se")
  expect_identical(d$search, "forward")
  expect_true(d$validate_search)
  expect_null(d$tolerance)
  expect_null(d$effect_size)
  expect_null(d$nterms_max)
})


test_that("the element names are the ones the projection reads", {
  # .fb_sel_decision_spec() copies by name, so a renamed element would be
  # silently dropped and the default used in its place.
  d <- fb_decision(rule = "tolerance", tolerance = 0.25, search = "L1",
                   nterms_max = 3L, validate_search = FALSE)
  spec <- .fb_sel_decision_spec(d)
  expect_identical(spec$rule, "tolerance")
  expect_identical(spec$tolerance, 0.25)
  expect_identical(spec$search, "L1")
  expect_identical(spec$nterms_max, 3L)
  expect_false(spec$validate_search)
})


# --- the rule that has no number ---------------------------------------- #

test_that("rule = \"tolerance\" without a tolerance refuses at construction", {
  expect_error(
    fb_decision(rule = "tolerance"),
    class = "flexybayes_refusal_selection_decision_invalid"
  )
  expect_error(fb_decision(rule = "tolerance"), "nats per resampling unit")
  # And the negative control: with one, it builds.
  expect_s3_class(fb_decision(rule = "tolerance", tolerance = 0), "fb_decision")
})


test_that("a tolerance handed to a rule that ignores it refuses", {
  # The failure this closes is a silent one: a caller who writes
  # `tolerance = 0.1` and leaves `rule` at its default gets a
  # one-standard-error answer and no indication that the number they chose
  # was never read.
  expect_error(fb_decision(tolerance = 0.1), "does not read it")
  expect_error(fb_decision(rule = "min", tolerance = 0.1), "does not read it")
})


test_that("a non-finite or negative tolerance refuses", {
  expect_error(fb_decision(rule = "tolerance", tolerance = -1), "non-negative")
  expect_error(fb_decision(rule = "tolerance", tolerance = NA_real_),
               "non-negative")
  expect_error(fb_decision(rule = "tolerance", tolerance = c(1, 2)),
               "non-negative")
})


# --- the threshold decision --------------------------------------------- #

test_that("posterior_threshold without an effect size refuses, and says why", {
  expect_error(
    fb_decision(method = "posterior_threshold"),
    class = "flexybayes_refusal_selection_decision_invalid"
  )
  # The refusal states the reason a zero default would be wrong rather than
  # only that the argument is missing.
  expect_error(
    fb_decision(method = "posterior_threshold"),
    "no mass at zero"
  )
  expect_s3_class(
    fb_decision(method = "posterior_threshold", effect_size = 0.2),
    "fb_decision"
  )
})


test_that("an effect size handed to a method that ignores it refuses", {
  expect_error(fb_decision(effect_size = 0.2), "does not read it")
  expect_error(fb_decision(method = "full", effect_size = 0.2),
               "does not read it")
})


test_that("the threshold is a probability and is held to it", {
  expect_error(fb_decision(threshold = 0), "strictly inside")
  expect_error(fb_decision(threshold = 1), "strictly inside")
  expect_error(fb_decision(threshold = 1.5), "strictly inside")
  expect_s3_class(fb_decision(threshold = 0.9), "fb_decision")
})


test_that("nterms_max must be a whole number of terms", {
  expect_error(fb_decision(nterms_max = 0), "at least 1")
  expect_error(fb_decision(nterms_max = 2.5), "whole number")
  expect_identical(fb_decision(nterms_max = 4)$nterms_max, 4L)
})


test_that("validate_search is a flag, and the fast mode says what it is", {
  expect_error(fb_decision(validate_search = NA), "TRUE or FALSE")
  out <- paste(
    utils::capture.output(print(fb_decision(validate_search = FALSE))),
    collapse = " "
  )
  expect_match(out, "optimistic", fixed = TRUE)
  expect_match(out, "NOT cross-validated", fixed = TRUE)
})


test_that("print states the rule rather than only the method", {
  out <- paste(
    utils::capture.output(print(fb_decision(method = "posterior_threshold",
                                            effect_size = 0.3,
                                            threshold = 0.8))),
    collapse = " "
  )
  expect_match(out, "P(|beta| > 0.3) > 0.8", fixed = TRUE)
  expect_match(out, "not an inclusion probability", fixed = TRUE)

  full <- paste(
    utils::capture.output(print(fb_decision(method = "full"))),
    collapse = " "
  )
  expect_match(full, "no selection", fixed = TRUE)
})


test_that("the vocabulary of every axis is closed", {
  expect_error(fb_decision(method = "stepwise"))
  expect_error(fb_decision(metric = "aic"))
  expect_error(fb_decision(rule = "bic"))
  expect_error(fb_decision(search = "backward"))
})
