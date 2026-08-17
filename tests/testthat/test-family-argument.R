# test-family-argument.R -- the `family` argument accepts what an R user
# reaches for.
#
# `flexybayes()` documents `family` as a string, but `stats::binomial()`
# is the first thing an R user tries, and until 0.9.0 it produced the raw
# base-R error "the condition has length > 1" from `.resolve_family()`'s
# scalar `if`. A family object now resolves; anything that is neither a
# name nor a usable family object refuses typed.

test_that("a stats family object resolves to its family and link", {
  expect_identical(
    flexyBayes:::.resolve_family(stats::binomial(), NULL),
    list(family = "binomial", link = "logit")
  )
  expect_identical(
    flexyBayes:::.resolve_family(stats::gaussian(), NULL),
    list(family = "gaussian", link = "identity")
  )
  expect_identical(
    flexyBayes:::.resolve_family(stats::poisson(), NULL),
    list(family = "poisson", link = "log")
  )
})

test_that("a family object's own link is honoured, not overwritten", {
  # binomial(link = "probit") must not silently become the logit the
  # string path defaults to.
  expect_identical(
    flexyBayes:::.resolve_family(stats::binomial(link = "probit"), NULL)$link,
    "probit"
  )
  # Gamma() names the inverse link. `family = "gamma"` defaults to log.
  # The two spellings therefore mean different models, and the object
  # spelling is read as written.
  expect_identical(
    flexyBayes:::.resolve_family(stats::Gamma(), NULL)$link,
    "inverse"
  )
  expect_identical(
    flexyBayes:::.resolve_family("gamma", NULL)$link,
    "log"
  )
})

test_that("a family generator passed unevaluated resolves", {
  expect_identical(
    flexyBayes:::.resolve_family(stats::binomial, NULL),
    list(family = "binomial", link = "logit")
  )
})

test_that("the existing string path is unchanged", {
  expect_identical(
    flexyBayes:::.resolve_family("gaussian", NULL),
    list(family = "gaussian", link = "identity")
  )
  expect_identical(
    flexyBayes:::.resolve_family("BINOMIAL", "probit"),
    list(family = "binomial", link = "probit")
  )
})

test_that("a malformed family object refuses typed", {
  bogus <- structure(list(nonsense = TRUE), class = "family")
  err <- tryCatch(
    flexyBayes:::.resolve_family(bogus, NULL),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_argument_not_recognised")
  expect_s3_class(err, "flexybayes_family_argument_refusal")
  expect_match(conditionMessage(err), "binomial()", fixed = TRUE)
})

test_that("a non-family, non-string family argument refuses typed", {
  for (bad in list(1L, c("gaussian", "binomial"), NULL, list(a = 1))) {
    err <- tryCatch(
      flexyBayes:::.resolve_family(bad, NULL),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_family_argument_not_recognised")
  }
})

test_that("a family object contradicting `link` refuses rather than picking", {
  err <- tryCatch(
    flexyBayes:::.resolve_family(stats::binomial(), "probit"),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_family_argument_not_recognised")
  expect_match(conditionMessage(err), "name different links", fixed = TRUE)
})

test_that("a family object agreeing with `link` is accepted", {
  expect_identical(
    flexyBayes:::.resolve_family(stats::binomial(), "logit"),
    list(family = "binomial", link = "logit")
  )
})


# ---------------------------------------------------------------- #
# End to end through the entry point                                 #
# ---------------------------------------------------------------- #

test_that("family = binomial() plans and emits like the string spelling", {
  skip_if_not_installed("INLA")
  set.seed(20260816L)
  n <- 200L
  d <- data.frame(
    a = factor(sample(c("a1", "a2"), n, replace = TRUE)),
    b = factor(sample(c("b1", "b2"), n, replace = TRUE))
  )
  d$y <- stats::rbinom(n, 1L, 0.4)

  p_obj <- suppressMessages(flexybayes(
    y ~ a + b,
    data = d,
    family = stats::binomial(),
    backend = "inla",
    plan = TRUE
  ))
  p_str <- suppressMessages(flexybayes(
    y ~ a + b,
    data = d,
    family = "binomial",
    backend = "inla",
    plan = TRUE
  ))
  expect_s3_class(p_obj, "fb_plan")
  expect_identical(p_obj$backend_chosen, p_str$backend_chosen)
  expect_identical(p_obj$gate_outcome, p_str$gate_outcome)
})
