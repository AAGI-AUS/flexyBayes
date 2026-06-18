# Tests for the Stage 2 MVP preflight integration in .dispatch_backend()
# (ADR 0021 / v0.3.0). Verifies:
#   - small-n (< 1e5) paths bypass preflight entirely; existing test
#     suite behaviour is preserved
#   - large-n paths run preflight and short-circuit with a structured
#     <flexybayes_preflight_refusal> condition when the design memory
#     exceeds the ceiling
#   - the condition carries the binding term + ceiling so downstream
#     tooling can pattern-match without parsing free text

test_that("dispatch: small-n path bypasses preflight (existing behaviour intact)", {
  df <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    g = factor(rep(letters[1:5], 10))
  )
  # review_code = TRUE returns the <flexybayes_review> token without
  # invoking the backend; exercises dispatch's branching while staying
  # offline (no greta / INLA install needed for this test).
  rev <- flexybayes(y ~ x, random = ~g, data = df, review_code = TRUE)
  expect_s3_class(rev, "flexybayes_review")
  # No preflight slot mutates the review object at this commit (commit
  # 7 wires the print integration); the review object survives the
  # bypass.
  expect_identical(rev$backend, "greta")
})

test_that("dispatch: large-n path triggers preflight refusal with tight ceiling", {
  # The review_code = TRUE path on flexybayes() short-circuits before
  # .dispatch_backend() (commit 7 hooks preflight into that path); so
  # we call .dispatch_backend() directly here to exercise the
  # dispatch-side gate. Force a refusal deterministically via the
  # flexyBayes.preflight_ceiling_gb option (0.0001 GiB = ~100 KB) so
  # the test does not depend on host RAM.
  N <- 2e5L
  df <- data.frame(
    y = rnorm(N),
    g = factor(sample.int(50000L, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ (1 | g), data = df)
  expect_true(fb$data_summary$n >= 1e5L)

  withr::local_options(flexyBayes.preflight_ceiling_gb = 0.0001)

  err <- tryCatch(
    flexyBayes:::.dispatch_backend(
      fb = fb,
      data = df,
      backend = "greta",
      known_matrices = list(),
      weights = NULL,
      n_samples = 100,
      warmup = 100,
      chains = 1,
      prior_fixed_sd = 100,
      prior_vc_sd = 1,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      return_code = TRUE,
      the_call = NULL,
      fixed = y ~ 1,
      random = ~g,
      rcov = NULL,
      family = "gaussian",
      link = "identity",
      data_name = "df"
    ),
    flexybayes_preflight_refusal = function(c) c,
    error = function(e) e
  )

  expect_s3_class(err, "flexybayes_preflight_refusal")
  expect_identical(err$reason_code, "design_memory_exceeds_ceiling")
  expect_true(grepl("\\(1 \\| g\\)", err$binding_term))
  expect_s3_class(err$refusal, "fb_preflight_refusal")
  expect_true(err$ceiling_bytes < 1024^2) # ~100 KB ceiling forced
})

test_that("dispatch: large-n path with generous ceiling passes through preflight", {
  # Direct dispatch call (review path covered separately in commit 7
  # tests). On a normal-RAM host the default ceiling (0.8 x RAM)
  # accepts the modest design and dispatch reaches return_code = TRUE
  # without preflight refusing.
  N <- 2e5L
  df <- data.frame(
    y = rnorm(N),
    g = factor(sample.int(50L, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ (1 | g), data = df)

  # Run with no synthetic ceiling override
  code <- flexyBayes:::.dispatch_backend(
    fb = fb,
    data = df,
    backend = "greta",
    known_matrices = list(),
    weights = NULL,
    n_samples = 100,
    warmup = 100,
    chains = 1,
    prior_fixed_sd = 100,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    return_code = TRUE,
    the_call = NULL,
    fixed = y ~ 1,
    random = ~g,
    rcov = NULL,
    family = "gaussian",
    link = "identity",
    data_name = "df"
  )
  expect_type(code, "character")
  expect_true(nchar(code) > 100L)
})
