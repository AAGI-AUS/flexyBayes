# =============================================================================
# Two spellings of one model reach the same verdict.
#
# `y ~ x + (1 | g:h)` and `y ~ x, random = ~ g:h` are the same model in two
# grammars, and the README says the grammar "is detected from the call
# shape", a promise of surface-independent behaviour. On INLA they
# behaved differently: the ASReml spelling was stopped by `lgm_gate()` with
# a typed refusal naming the engine that does fit it, and the bar spelling
# was read as a simple random intercept on a group named "g:h", passed the
# gate, reached the INLA binary, and failed there with a message about the
# emitted formula instead of about the user's model.
#
# A refusal contract that holds on one surface and not the other is worse
# than no contract, because the user learns to trust it. The fix is at the
# ingest layer instead of at the gate: both surfaces are classified by the
# same walker, so the two descriptors are identical and every downstream
# gate sees one model.
# =============================================================================

.gsp_data <- function() {
  set.seed(20260819L)
  g <- factor(paste0("g", rep(seq_len(10L), times = 6L)))
  f <- factor(rep(c("e1", "e2", "e3"), each = 20L))
  h <- factor(rep(c("b1", "b2", "b3", "b4"), times = 15L))
  x <- stats::rnorm(60L)
  eta <- 0.6 +
    0.5 * x +
    stats::rnorm(10L, 0, 0.6)[as.integer(g)] +
    stats::rnorm(3L, 0, 0.4)[as.integer(f)]
  data.frame(
    y = eta + stats::rnorm(60L, 0, 0.8),
    x = x,
    g = g,
    f = f,
    h = h
  )
}

test_that("the two grammars build the same descriptor for a nested term", {
  d <- .gsp_data()
  asreml_ir <- fb_from_asreml(y ~ x, random = ~ g:h, data = d)
  bar_ir <- suppressMessages(fb_from_brms(y ~ x + (1 | g:h), data = d))
  expect_identical(
    bar_ir$random_terms[[1L]],
    asreml_ir$random_terms[[1L]]
  )
  expect_identical(bar_ir$random_terms[[1L]]$type, "nested")
})

test_that("a three-way bar interaction lands on the same combo descriptor", {
  d <- .gsp_data()
  asreml_ir <- fb_from_asreml(y ~ x, random = ~ g:h:f, data = d)
  bar_ir <- suppressMessages(fb_from_brms(y ~ x + (1 | g:h:f), data = d))
  expect_identical(
    bar_ir$random_terms[[1L]],
    asreml_ir$random_terms[[1L]]
  )
  expect_identical(bar_ir$random_terms[[1L]]$type, "combo")
})

test_that("both grammars fit the nested term on INLA, and agree", {
  # Until 0.10.0.9000 this asserted that both grammars REFUSED, and the
  # parity it checked was parity of a refusal. INLA fits a nested random
  # effect as an iid field over the interaction index, so the refusal was
  # ours and it is gone. The parity the test exists for is unchanged and is
  # now checked on the fit itself, which is a stronger statement: the
  # brms-style and ASReml-style spellings of one model must produce the same
  # components with the same values.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .gsp_data()
  bar <- suppressMessages(flexybayes(
    y ~ x + (1 | g:h),
    data = d,
    family = "gaussian",
    backend = "inla",
    verbose = FALSE
  ))
  asreml <- suppressMessages(flexybayes(
    y ~ x,
    random = ~ g:h,
    data = d,
    family = "gaussian",
    backend = "inla",
    verbose = FALSE
  ))
  expect_identical(class(bar), class(asreml))

  v_bar <- summary(bar)$varcomp
  v_asr <- summary(asreml)$varcomp
  expect_identical(v_bar$component, v_asr$component)
  expect_true("sd_g:h" %in% v_bar$component)
  # The two spellings build the same IR and the same INLA formula, and the
  # estimates agree to about 1e-6 relative instead of to the bit: INLA's
  # optimiser does not land on identical numbers from the two term orderings
  # the grammars produce. 1e-4 is far tighter than any modelling difference
  # and still fails if the grammars ever diverge in substance.
  expect_equal(v_bar$estimate, v_asr$estimate, tolerance = 1e-4)
})

test_that("both grammars still emit the nested term on brms", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .gsp_data()
  bar <- suppressMessages(flexybayes(
    y ~ x + (1 | g:h),
    data = d,
    family = "gaussian",
    backend = "brms",
    return_code = TRUE
  ))
  asreml <- suppressMessages(flexybayes(
    y ~ x,
    random = ~ g:h,
    data = d,
    family = "gaussian",
    backend = "brms",
    return_code = TRUE
  ))
  expect_type(as.character(bar), "character")
  expect_identical(as.character(bar), as.character(asreml))
})

test_that("the numeric-in-a-random-interaction guard fires on both surfaces", {
  d <- .gsp_data()
  bar <- tryCatch(
    suppressMessages(fb_from_brms(y ~ x + (1 | f:x), data = d)),
    condition = function(c) c
  )
  asreml <- tryCatch(
    fb_from_asreml(y ~ x, random = ~ f:x, data = d),
    condition = function(c) c
  )
  expect_s3_class(
    bar,
    "flexybayes_refusal_numeric_variable_in_random_interaction"
  )
  expect_identical(class(bar), class(asreml))
})
