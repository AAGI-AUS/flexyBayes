# Tests for the Stage 3A model-level aggregation plan (audit P1.5,
# 2026-05-23). Per-term `aggregated_likelihood_candidate` flags
# remain the dispatcher contract through v0.3.1; this plan is the
# v0.3.2 supersedure path. At v0.3.0.9000 the plan is attached to
# `.fb_preflight()` results for diagnostic visibility -- no dispatch
# consumer yet.

# Hand-build minimal <fb_terms> IR (avoids the brms walker's pre-
# ingest refusals for some scope tests).
.test_make_plan_ir <- function(
  family = "gaussian",
  link = "identity",
  fixed_terms = list(),
  random_terms = list(),
  n = 1000L
) {
  structure(
    list(
      response = "y",
      family = family,
      link = link,
      intercept = TRUE,
      fixed_terms = fixed_terms,
      random_terms = random_terms,
      rcov_terms = list(),
      addition_terms = list(),
      priors = list(),
      data_summary = list(n = n),
      capabilities = character(),
      source = "test"
    ),
    class = c("fb_terms", "list")
  )
}


# ---------------------------------------------------------------- #
# Eligibility cases                                                 #
# ---------------------------------------------------------------- #

test_that("eligible: y ~ f + (1 | g) with K/N below threshold", {
  set.seed(1L)
  N <- 1000L
  df <- data.frame(
    y = rnorm(N),
    f = factor(sample(letters[1:5], N, replace = TRUE)),
    g = factor(sample.int(20L, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = df)
  ds <- .fb_dataset(df)
  plan <- .fb_aggregation_plan(fb, ds)

  expect_s3_class(plan, "fb_aggregation_plan")
  expect_true(plan$eligible)
  expect_identical(plan$reason_codes, character(0L))
  # K_est = L_f * L_g = 5 * 20 = 100; N = 1000; ratio = 0.1.
  expect_identical(plan$K_est, 100L)
  expect_identical(plan$N, N)
  expect_equal(plan$compression_est, 0.1)
  expect_false(plan$requires_materialisation)
  expect_length(plan$cell_key_terms, 2L)
})

test_that("eligible at the productivity boundary (exactly K/N = threshold)", {
  # Construct so K_est = N * threshold = 500 (5 * 100 levels).
  set.seed(2L)
  N <- 1000L
  df <- data.frame(
    y = rnorm(N),
    f = factor(sample(letters[1:5], N, replace = TRUE)),
    g = factor(sample.int(100L, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = df)
  plan <- .fb_aggregation_plan(fb, .fb_dataset(df))
  expect_true(plan$eligible) # boundary is inclusive (<=)
  expect_equal(plan$compression_est, 0.5)
})


# ---------------------------------------------------------------- #
# Out-of-scope cases                                                #
# ---------------------------------------------------------------- #

test_that("out of scope: non-aggregatable family (gamma)", {
  set.seed(3L)
  df <- data.frame(
    y = rgamma(50L, 2),
    g = factor(sample(1:5, 50L, replace = TRUE))
  )
  fb <- fb_from_asreml(y ~ g, data = df, family = "gamma")
  plan <- .fb_aggregation_plan(fb, .fb_dataset(df))

  expect_false(plan$eligible)
  expect_true("non_aggregatable_family" %in% plan$reason_codes)
})

test_that("in scope: binomial and poisson with factor cell keys", {
  set.seed(7L)
  for (fam in c("binomial", "poisson")) {
    y <- if (identical(fam, "binomial")) {
      rbinom(200L, 1L, 0.4)
    } else {
      rpois(200L, 2)
    }
    df <- data.frame(y = y, g = factor(sample(1:5, 200L, replace = TRUE)))
    fb <- fb_from_asreml(y ~ 1, random = ~g, data = df, family = fam)
    plan <- .fb_aggregation_plan(fb, .fb_dataset(df))
    expect_true(plan$eligible, info = fam)
  }
})

test_that("out of scope: poisson with a continuous fixed effect", {
  set.seed(3L)
  df <- data.frame(y = rpois(50L, 2), x = rnorm(50L))
  fb <- fb_from_brms(y ~ x, data = df, family = "poisson")
  plan <- .fb_aggregation_plan(fb, .fb_dataset(df))

  expect_false(plan$eligible)
  expect_true("continuous_cell_key_data_dependent" %in% plan$reason_codes)
})

test_that("out of scope: smooth (s(x)) fixed term", {
  fb <- .test_make_plan_ir(
    fixed_terms = list(list(type = "smooth_mgcv", var = "x", k = 10L))
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  plan <- .fb_aggregation_plan(fb, ds)
  expect_false(plan$eligible)
  expect_true("smooth_term_not_aggregatable" %in% plan$reason_codes)
})

test_that("out of scope: uncorrelated random slope (x || g)", {
  set.seed(4L)
  df <- data.frame(
    y = rnorm(100L),
    x = rnorm(100L),
    g = factor(rep(letters[1:5], 20L))
  )
  fb <- fb_from_brms(y ~ (x || g), data = df)
  plan <- .fb_aggregation_plan(fb, .fb_dataset(df))
  expect_false(plan$eligible)
  expect_true("random_slope_not_aggregatable" %in% plan$reason_codes)
})

test_that("out of scope: structured-covariance random term (audit P2.8)", {
  fb <- .test_make_plan_ir(
    n = 100L,
    random_terms = list(list(type = "vm", var = "g", var_n = 10L))
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 100L,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list(g = letters[1:10])
  )
  plan <- .fb_aggregation_plan(fb, ds)
  expect_false(plan$eligible)
  expect_true(
    "structured_random_not_aggregatable" %in%
      plan$reason_codes
  )
})


# ---------------------------------------------------------------- #
# Compression edge cases                                            #
# ---------------------------------------------------------------- #

test_that("data-dependent: numeric fixed term forces requires_materialisation", {
  set.seed(5L)
  df <- data.frame(
    y = rnorm(100L),
    x = rnorm(100L),
    g = factor(rep(letters[1:5], 20L))
  )
  fb <- fb_from_brms(y ~ x + (1 | g), data = df)
  plan <- .fb_aggregation_plan(fb, .fb_dataset(df))

  expect_false(plan$eligible)
  expect_identical(plan$reason_codes, "continuous_cell_key_data_dependent")
  expect_true(plan$requires_materialisation)
  expect_true(is.na(plan$K_est))
})

test_that("unproductive: K/N above threshold (factor levels close to N)", {
  # Hand-built IR with explicit var_n = 80 -- fb_from_brms re-factors
  # group columns via factor() which drops unused levels, so a real-
  # data path cannot guarantee 80 levels on 100 rows.
  fb <- .test_make_plan_ir(
    n = 100L,
    random_terms = list(list(type = "simple", var = "g", var_n = 80L))
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 100L,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(80L)))
  )
  plan <- .fb_aggregation_plan(fb, ds)

  # 80 levels on 100 rows -> K/N = 0.8 > 0.5 threshold.
  expect_false(plan$eligible)
  expect_identical(plan$reason_codes, "compression_unproductive")
  expect_identical(plan$K_est, 80L)
  expect_equal(plan$compression_est, 0.8)
})

test_that("level count unresolvable: metadata-only dataset without dictionaries", {
  fb <- .test_make_plan_ir(
    fixed_terms = list(list(type = "factor", var = "f"))
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", f = "factor"),
    dictionaries = list() # no dictionary for `f`
  )
  plan <- .fb_aggregation_plan(fb, ds)

  expect_false(plan$eligible)
  expect_identical(plan$reason_codes, "compression_level_count_unresolvable")
  expect_true(plan$requires_materialisation)
  expect_true(is.na(plan$K_est))
})


# ---------------------------------------------------------------- #
# Preflight integration                                             #
# ---------------------------------------------------------------- #

test_that(".fb_preflight() attaches the aggregation_plan slot", {
  set.seed(7L)
  N <- 200L
  df <- data.frame(
    y = rnorm(N),
    g = factor(sample.int(10L, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ (1 | g), data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  expect_s3_class(pf$aggregation_plan, "fb_aggregation_plan")
  expect_true(pf$aggregation_plan$eligible)
  expect_identical(pf$aggregation_plan$K_est, 10L)
})


# ---------------------------------------------------------------- #
# Print method                                                      #
# ---------------------------------------------------------------- #

test_that("print.fb_aggregation_plan: eligible + reason_codes rendering", {
  set.seed(8L)
  N <- 1000L
  df <- data.frame(
    y = rnorm(N),
    f = factor(sample(letters[1:5], N, replace = TRUE)),
    g = factor(sample.int(20L, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = df)
  plan <- .fb_aggregation_plan(fb, .fb_dataset(df))

  out <- capture.output(print(plan))
  expect_true(any(grepl("<fb_aggregation_plan>", out, fixed = TRUE)))
  expect_true(any(grepl("eligible = TRUE", out, fixed = TRUE)))
  expect_true(any(grepl("cell_key_terms:", out, fixed = TRUE)))
  expect_true(any(grepl("K_est = 100", out, fixed = TRUE)))
  expect_true(any(grepl("compression_est = 0.100", out, fixed = TRUE)))
})

test_that("print.fb_aggregation_plan: NA fields render cleanly", {
  fb <- .test_make_plan_ir(
    fixed_terms = list(list(type = "smooth_mgcv", var = "x", k = 10L))
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  plan <- .fb_aggregation_plan(fb, ds)

  out <- capture.output(print(plan))
  expect_true(any(grepl("eligible = FALSE", out, fixed = TRUE)))
  expect_true(any(grepl("smooth_term_not_aggregatable", out, fixed = TRUE)))
  expect_true(any(grepl("K_est = NA", out, fixed = TRUE)))
  expect_true(any(grepl("compression_est = NA", out, fixed = TRUE)))
})


# ---------------------------------------------------------------- #
# Internal-only contract                                            #
# ---------------------------------------------------------------- #

test_that(".fb_aggregation_plan is internal -- no exported binding", {
  ns <- asNamespace("flexyBayes")
  exp <- getNamespaceExports(ns)
  expect_false(".fb_aggregation_plan" %in% exp)
  expect_false("fb_aggregation_plan" %in% exp)
  expect_true(exists(".fb_aggregation_plan", envir = ns, inherits = FALSE))
})
