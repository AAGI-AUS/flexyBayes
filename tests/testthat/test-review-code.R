# Tests for the `review_code = TRUE` option on flexybayes().
#
# Contract: review_code = TRUE returns a `<flexybayes_review>`
# deferred-execution token instead of fitting. cat_code() round-trips
# the generated code string; proceed() advances into a fit (class
# `flexybayes`); a second proceed() returns the cached fit. The R-RNG
# snapshot captured at construction is restored on the first proceed()
# so that `set.seed(s); fit_direct <- flexybayes(...)` and
# `set.seed(s); rev <- flexybayes(..., review_code = TRUE);
# fit_via <- proceed(rev)` agree on the posterior point estimates.
# Mutual-exclusion guard: `return_code = TRUE` + `review_code = TRUE`
# raises a structured refusal.
#
# Review mode only ever produces Stan code via brms (see NEWS.md, 0.9.3 --
# the native engine this mode also once supported was withdrawn entirely;
# `backend = "auto"` under review_code now resolves to brms unconditionally,
# since it is the only remaining code-producing engine). brms is required
# for the proceed-into-fit assertions; the field / class / mutual-exclusion
# checks run without it.

.check_brms_review_prereqs <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
}

mk_review_data <- function() {
  set.seed(42)
  n <- 30L
  data.frame(
    yield = rnorm(n, mean = 100, sd = 10),
    env = factor(rep(1:3, length.out = n)),
    geno = factor(rep(1:5, length.out = n))
  )
}


# ---------------------------------------------------------------- #
# (a) review_code = TRUE returns a <flexybayes_review> with the    #
#     expected fields and class.                                   #
# ---------------------------------------------------------------- #

test_that("review_code = TRUE returns a flexybayes_review with the expected fields", {
  skip_if_not_installed("brms")
  d <- mk_review_data()
  rev <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE
  )

  expect_s3_class(rev, "flexybayes_review")
  expect_true(is.environment(rev))

  # Slots
  expect_type(rev$code, "character")
  expect_length(rev$code, 1L)
  expect_true(nzchar(rev$code))
  expect_identical(rev$backend, "stan_via_brms")
  expect_s3_class(rev$ir, "fb_terms")
  expect_identical(rev$data_name, "d")
  expect_true(is.call(rev$call))
  expect_false(isTRUE(rev$proceeded))
  expect_null(rev$fit)
  expect_true(is.list(rev$proceed_args))
})


# ---------------------------------------------------------------- #
# (b) print() emits the two-line summary.                          #
# ---------------------------------------------------------------- #

test_that("print(<flexybayes_review>) emits the two-line summary", {
  skip_if_not_installed("brms")
  d <- mk_review_data()
  rev <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE
  )

  out <- utils::capture.output(print(rev))
  expect_length(out, 2L)
  expect_match(out[[1]], "^<flexybayes_review>")
  expect_match(out[[1]], "backend=stan_via_brms")
  expect_match(out[[2]], "cat_code", fixed = TRUE)
  expect_match(out[[2]], "proceed", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# (c) cat_code() round-trips the code string.                      #
# ---------------------------------------------------------------- #

test_that("cat_code(rev) emits and returns the stored code", {
  skip_if_not_installed("brms")
  d <- mk_review_data()
  rev <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE
  )

  emitted <- utils::capture.output(returned <- cat_code(rev))
  expect_identical(returned, rev$code)
  expect_true(any(grepl("generated with brms", emitted, fixed = TRUE)))
})


# ---------------------------------------------------------------- #
# (d) proceed() returns a flexybayes-class fit; cache + (e) second #
#     call returns the same fit identically.                       #
# ---------------------------------------------------------------- #

test_that("proceed(rev) advances into a fit and caches the result", {
  .check_brms_review_prereqs()
  d <- mk_review_data()

  set.seed(20260521L)
  # Structural test (does proceed() advance into a fit and cache it?),
  # not numeric agreement. The 200-draw / chains=1 budget produced a
  # genuine R-hat warning plus a residual low-ESS warning that carries
  # flexyBayes's own convergence message rather than rstan's -- neither
  # is caught by .muffle_ess_warnings() (by design: it does not widen
  # to cover a real non-convergence signal), so the budget is raised
  # (with a higher adapt_delta) until the fit genuinely converges. The
  # RNG-snapshot-at-construction contract this test is about (see file
  # header) is independent of n_samples/warmup/adapt_delta.
  rev <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    n_samples = 1000,
    warmup = 1000,
    chains = 1,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE
  ))
  # Seed slot is populated from the user's set.seed() above.
  expect_false(is.null(rev$seed))
  expect_type(rev$seed, "integer")

  fit1 <- suppressMessages(proceed(rev))

  expect_s3_class(fit1, "flexybayes")
  expect_true(isTRUE(rev$proceeded))
  expect_identical(rev$fit, fit1)

  # Second proceed() returns the cached fit. Silence the one-time
  # cached-fit note via the documented option.
  prev_opt <- options(flexyBayes.silence_review_cached_note = TRUE)
  on.exit(options(prev_opt), add = TRUE)
  fit2 <- proceed(rev)
  expect_identical(fit2, fit1)
})


# ---------------------------------------------------------------- #
# (f) Mutual-exclusion: return_code + review_code raises a         #
#     structured refusal.                                          #
# ---------------------------------------------------------------- #

test_that("return_code + review_code raises a structured refusal", {
  # The mutual-exclusion guard fires before any backend work, so this
  # does not need brms installed.
  d <- mk_review_data()
  err <- tryCatch(
    flexybayes(
      yield ~ env,
      random = ~geno,
      data = d,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      return_code = TRUE,
      review_code = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_code_flags_mutually_exclusive")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "mutually exclusive")
})


# ---------------------------------------------------------------- #
# (g) Seed snapshot lets two review objects -- built at the same   #
#     outer seed -- reproduce the same code string. The R-RNG      #
#     state captured at construction is what matters for the       #
#     proceed-cache contract; Stan-side non-determinism on the      #
#     brms chain is not part of the snapshot guarantee.            #
# ---------------------------------------------------------------- #

test_that("review objects built at the same outer seed carry identical code + seed", {
  skip_if_not_installed("brms")
  d <- mk_review_data()

  set.seed(20260521L)
  rev_a <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE
  )

  set.seed(20260521L)
  rev_b <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE
  )

  expect_identical(rev_a$code, rev_b$code)
  expect_identical(rev_a$seed, rev_b$seed)
})


# ---------------------------------------------------------------- #
# (h) Session-level default option opts call sites into review     #
#     mode without changing call shape.                            #
# ---------------------------------------------------------------- #

test_that("flexyBayes.review_code_default opts a call into review mode", {
  skip_if_not_installed("brms")
  d <- mk_review_data()

  prev_opt <- options(flexyBayes.review_code_default = TRUE)
  on.exit(options(prev_opt), add = TRUE)

  rev <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  expect_s3_class(rev, "flexybayes_review")
})


# ---------------------------------------------------------------- #
# (i) Audit P1.1 (2026-05-24): session-level review_code default   #
#     must not bypass the unsupported-backend refusal. Eight       #
#     focused cases pin the order: option resolution happens BEFORE#
#     the backend guard for both flexybayes() and fb_brms(), and   #
#     the refusal fires for both the explicit argument path and    #
#     the option-default path.                                     #
# ---------------------------------------------------------------- #

test_that("flexybayes(): explicit review_code = TRUE, backend = 'inla' refuses", {
  d <- mk_review_data()
  expect_error(
    flexybayes(
      yield ~ env,
      random = ~geno,
      data = d,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      review_code = TRUE,
      backend = "inla"
    ),
    'supported with backend = "brms"',
    fixed = TRUE
  )
})

# ADR 0031 Q1 (0.9.3: the native engine this mode once also supported was
# withdrawn entirely -- see NEWS.md): review_code with backend = "auto"
# resolves to brms (the only remaining code-producing engine) rather than
# refusing. Explicit backend = "inla" still refuses (tests above/below).
test_that("flexybayes(): review_code = TRUE, backend = 'auto' resolves to brms", {
  skip_if_not_installed("brms")
  d <- mk_review_data()
  rev <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE,
    backend = "auto"
  )
  expect_s3_class(rev, "flexybayes_review")
})

test_that("flexybayes(): option-driven review default, backend = 'inla' refuses", {
  d <- mk_review_data()
  prev_opt <- options(flexyBayes.review_code_default = TRUE)
  on.exit(options(prev_opt), add = TRUE)
  expect_error(
    flexybayes(
      yield ~ env,
      random = ~geno,
      data = d,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      backend = "inla"
    ),
    'supported with backend = "brms"',
    fixed = TRUE
  )
})

test_that("flexybayes(): option-driven review default, backend = 'auto' resolves to brms", {
  d <- mk_review_data()
  prev_opt <- options(flexyBayes.review_code_default = TRUE)
  on.exit(options(prev_opt), add = TRUE)
  rev <- flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    backend = "auto"
  )
  expect_s3_class(rev, "flexybayes_review")
})

test_that("fb_brms(): explicit review_code = TRUE, backend = 'inla' refuses", {
  d <- mk_review_data()
  expect_error(
    fb(
      yield ~ env + (1 | geno),
      data = d,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      review_code = TRUE,
      backend = "inla"
    ),
    'supported with',
    fixed = TRUE
  )
})

test_that("fb_brms(): review_code = TRUE, backend = 'auto' resolves to brms (ADR 0031 Q1)", {
  d <- mk_review_data()
  rev <- fb(
    yield ~ env + (1 | geno),
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    review_code = TRUE,
    backend = "auto"
  )
  expect_s3_class(rev, "flexybayes_review")
})

test_that("fb_brms(): option-driven review default, backend = 'inla' refuses", {
  d <- mk_review_data()
  prev_opt <- options(flexyBayes.review_code_default = TRUE)
  on.exit(options(prev_opt), add = TRUE)
  expect_error(
    fb(
      yield ~ env + (1 | geno),
      data = d,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      backend = "inla"
    ),
    'supported with',
    fixed = TRUE
  )
})

test_that("fb_brms(): option-driven review default, backend = 'auto' resolves to brms", {
  d <- mk_review_data()
  prev_opt <- options(flexyBayes.review_code_default = TRUE)
  on.exit(options(prev_opt), add = TRUE)
  rev <- fb(
    yield ~ env + (1 | geno),
    data = d,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    backend = "auto"
  )
  expect_s3_class(rev, "flexybayes_review")
})
