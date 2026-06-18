# Tests for fb_terms — flexyBayes intermediate representation (IR)
#
# Phase 0.A: skeleton class, constructor, validator, predicate,
# accessors, print and format S3 methods. No ingest logic exercised
# here; ingest tests land in phase 0.B (test-fb-from-asreml.R).

# ---------------------------------------------------------------- #
# Constructor                                                      #
# ---------------------------------------------------------------- #

test_that("new_fb_terms() constructs a valid object with defaults", {
  obj <- new_fb_terms(
    response = "yield",
    family = gaussian(),
    source = "asreml"
  )
  expect_s3_class(obj, "fb_terms")
  expect_s3_class(obj, "list")
  expect_identical(obj$response, "yield")
  expect_true(obj$intercept)
  expect_length(obj$fixed_terms, 0L)
  expect_length(obj$random_terms, 0L)
  expect_length(obj$rcov_terms, 0L)
  expect_length(obj$addition_terms, 0L)
  expect_identical(obj$source, "asreml")
  expect_null(obj$priors)
  expect_null(obj$link)
  expect_identical(obj$capabilities, character())
})

test_that("new_fb_terms() accepts brms source and character family", {
  obj <- new_fb_terms(
    response = "y",
    family = "gaussian",
    source = "brms"
  )
  expect_identical(obj$source, "brms")
  expect_identical(obj$family, "gaussian")
})

test_that("new_fb_terms() preserves complex term descriptors", {
  fixed <- list(
    list(type = "factor", var = "treatment", n_levels = 3L),
    list(type = "continuous", var = "age")
  )
  random <- list(
    list(type = "fa_gxe", outer = "env", inner = "geno", k = 2L)
  )
  rcov <- list(
    list(type = "at_units", var = "env")
  )
  obj <- new_fb_terms(
    response = "yield",
    family = gaussian(),
    fixed_terms = fixed,
    random_terms = random,
    rcov_terms = rcov,
    source = "asreml"
  )
  expect_identical(obj$fixed_terms, fixed)
  expect_identical(obj$random_terms, random)
  expect_identical(obj$rcov_terms, rcov)
})

# ---------------------------------------------------------------- #
# Validator                                                        #
# ---------------------------------------------------------------- #

test_that("validator rejects malformed `response`", {
  expect_error(
    new_fb_terms(
      response = c("a", "b"),
      family = gaussian(),
      source = "asreml"
    ),
    "`response` must be a non-empty length-1 character"
  )
  expect_error(
    new_fb_terms(response = "", family = gaussian(), source = "asreml"),
    "`response` must be a non-empty length-1 character"
  )
  expect_error(
    new_fb_terms(
      response = NA_character_,
      family = gaussian(),
      source = "asreml"
    ),
    "`response` must be a non-empty length-1 character"
  )
})

test_that("validator rejects malformed `family`", {
  expect_error(
    new_fb_terms(
      response = "y",
      family = list("not a family"),
      source = "asreml"
    ),
    "`family` must be a `family` object"
  )
  expect_error(
    new_fb_terms(
      response = "y",
      family = c("gaussian", "binomial"),
      source = "asreml"
    ),
    "`family` must be a `family` object"
  )
})

test_that("validator rejects malformed `intercept`", {
  expect_error(
    new_fb_terms(
      response = "y",
      family = gaussian(),
      intercept = "yes",
      source = "asreml"
    ),
    "`intercept` must be TRUE or FALSE"
  )
  expect_error(
    new_fb_terms(
      response = "y",
      family = gaussian(),
      intercept = c(TRUE, FALSE),
      source = "asreml"
    ),
    "`intercept` must be TRUE or FALSE"
  )
  expect_error(
    new_fb_terms(
      response = "y",
      family = gaussian(),
      intercept = NA,
      source = "asreml"
    ),
    "`intercept` must be TRUE or FALSE"
  )
})

test_that("validator rejects malformed term descriptors", {
  expect_error(
    new_fb_terms(
      response = "y",
      family = gaussian(),
      fixed_terms = list("not a list element"),
      source = "asreml"
    ),
    "fixed_terms\\[\\[1\\]\\]"
  )
  expect_error(
    new_fb_terms(
      response = "y",
      family = gaussian(),
      random_terms = list(list(no_type_field = TRUE)),
      source = "asreml"
    ),
    "random_terms\\[\\[1\\]\\]"
  )
  expect_error(
    new_fb_terms(
      response = "y",
      family = gaussian(),
      fixed_terms = list(list(type = "")),
      source = "asreml"
    ),
    "non-empty character `type`"
  )
})

test_that("validator rejects malformed `source`", {
  expect_error(
    new_fb_terms(response = "y", family = gaussian(), source = "invalid"),
    "should be one of"
  )
})

test_that("validator catches mutation that violates invariants", {
  obj <- new_fb_terms(response = "y", family = gaussian(), source = "asreml")
  obj$fixed_terms <- list(list(no_type_field = TRUE))
  expect_error(validate_fb_terms(obj), "fixed_terms\\[\\[1\\]\\]")
})

# ---------------------------------------------------------------- #
# Predicate                                                        #
# ---------------------------------------------------------------- #

test_that("is_fb_terms() recognises the class", {
  obj <- new_fb_terms(response = "y", family = gaussian(), source = "brms")
  expect_true(is_fb_terms(obj))
  expect_false(is_fb_terms(list()))
  expect_false(is_fb_terms(data.frame()))
  expect_false(is_fb_terms(NULL))
})

# ---------------------------------------------------------------- #
# Accessors                                                        #
# ---------------------------------------------------------------- #

test_that("accessors return their corresponding fields", {
  fixed <- list(
    list(type = "factor", var = "g"),
    list(type = "continuous", var = "x")
  )
  random <- list(list(type = "simple", var = "site"))
  obj <- new_fb_terms(
    response = "y",
    family = gaussian(),
    link = "log",
    intercept = FALSE,
    fixed_terms = fixed,
    random_terms = random,
    capabilities = c("lgm_compatible", "no_distributional"),
    source = "brms"
  )
  expect_identical(fb_response(obj), "y")
  expect_identical(fb_link(obj), "log")
  expect_false(fb_intercept(obj))
  expect_identical(fb_fixed_terms(obj), fixed)
  expect_identical(fb_random_terms(obj), random)
  expect_identical(
    fb_capabilities(obj),
    c("lgm_compatible", "no_distributional")
  )
  expect_identical(fb_source(obj), "brms")
})

# ---------------------------------------------------------------- #
# Print and format S3 methods                                      #
# ---------------------------------------------------------------- #

test_that("print.fb_terms() produces a readable multi-line summary", {
  obj <- new_fb_terms(
    response = "yield",
    family = gaussian(),
    fixed_terms = list(list(type = "factor", var = "treatment")),
    random_terms = list(list(type = "simple", var = "site")),
    source = "asreml"
  )
  out <- capture.output(print(obj))
  expect_true(any(grepl("fb_terms", out)))
  expect_true(any(grepl("source:.*asreml", out)))
  expect_true(any(grepl("response:.*yield", out)))
  expect_true(any(grepl("fixed:.*1", out)))
  expect_true(any(grepl("random:.*1", out)))
  expect_true(any(grepl("not yet evaluated by lgm_gate", out)))
})

test_that("print.fb_terms() displays evaluated capabilities when set", {
  obj <- new_fb_terms(
    response = "y",
    family = gaussian(),
    capabilities = c("lgm_compatible", "no_distributional"),
    source = "brms"
  )
  out <- capture.output(print(obj))
  expect_true(any(grepl("lgm_compatible.*no_distributional", out)))
  expect_false(any(grepl("not yet evaluated", out)))
})

test_that("format.fb_terms() returns a one-line summary", {
  obj <- new_fb_terms(
    response = "yield",
    family = gaussian(),
    fixed_terms = list(
      list(type = "factor", var = "treatment"),
      list(type = "continuous", var = "x")
    ),
    random_terms = list(list(type = "simple", var = "site")),
    source = "asreml"
  )
  s <- format(obj)
  expect_length(s, 1L)
  expect_match(s, "fb_terms")
  expect_match(s, "asreml")
  expect_match(s, "yield ~")
  expect_match(s, "2 fixed")
  expect_match(s, "1 random")
})
