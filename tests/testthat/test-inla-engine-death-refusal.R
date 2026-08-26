# test-inla-engine-death-refusal.R -- C2 / S2 / FS-25.
#
# R/emit_inla.R's do.call(INLA::inla, inla_args) used to catch any error
# and re-wrap it as `stop("INLA fit failed: ", conditionMessage(e))` --
# untyped, no reason code, no remedy. FS-25 measured a real instance: the
# fit ran for 41.6 minutes, then failed with a raw pass-through of "The
# inla-program exited with an error. Unless you interrupted it yourself,
# ...". This file exercises the fix: the call now routes through
# .inla_call(), and the tryCatch around it classifies the failure --
# a recognised engine-death message pattern becomes the typed refusal
# `inla_program_failed` carrying the design size, the largest verified
# size, the binding term, and the remedies; anything else (flexyBayes's
# own code) propagates unwrapped.
#
# The mocking approach is the one the spec names directly:
# testthat::local_mocked_bindings(.inla_call = ...) substitutes a
# controlled failure without needing a real 41-minute INLA crash.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

mk_inla_data <- function() {
  set.seed(2026L)
  n <- 30L
  data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = factor(rep(1:5, length.out = n))
  )
}

mk_fb_for_inla <- function(
  data,
  random_terms = list(list(type = "simple", var = "g"))
) {
  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    link = "identity",
    fixed_terms = list(list(type = "continuous", var = "x")),
    random_terms = random_terms,
    residual_terms = list(list(type = "units")),
    source = "asreml"
  )
  flexyBayes:::lgm_gate(fb)
}


# ---------------------------------------------------------------- #
# .inla_is_program_death_message() -- the pattern gate               #
# ---------------------------------------------------------------- #

test_that(".inla_is_program_death_message() matches the observed FS-25 text and its siblings, and nothing else", {
  expect_true(flexyBayes:::.inla_is_program_death_message(
    paste0(
      "The inla-program exited with an error. Unless you interrupted ",
      "it yourself, ..."
    )
  ))
  expect_true(flexyBayes:::.inla_is_program_death_message(
    "inla.core.safe: something went wrong internally"
  ))
  expect_true(flexyBayes:::.inla_is_program_death_message(
    "the fitting process exited unexpectedly"
  ))
  # Own-code / unrelated errors must NOT match.
  expect_false(flexyBayes:::.inla_is_program_death_message(
    "object 'bogus_var' not found"
  ))
  expect_false(flexyBayes:::.inla_is_program_death_message(
    "argument \"formula\" is missing, with no default"
  ))
  expect_false(flexyBayes:::.inla_is_program_death_message(NA_character_))
  expect_false(flexyBayes:::.inla_is_program_death_message(character(0)))
})


# ---------------------------------------------------------------- #
# .inla_largest_verified() -- reads the S4 ceilings artefact          #
# ---------------------------------------------------------------- #

test_that(".inla_largest_verified() reports the largest COMPLETED per-row fit on record", {
  # The 2026-08-22 ceilings study has no completed flexyBayes fit (its
  # 911,808-row rung was preflight refused; the record was corrected on
  # 2026-08-26), so the largest verified size is the boundary study's
  # 1,000,000-row per-row INLA fit, whose latent-field size was not recorded.
  res <- flexyBayes:::.inla_largest_verified()
  expect_false(is.null(res))
  expect_identical(res$n, 1e6)
  expect_true(is.na(res$random_effects))
  expect_identical(res$run_date, "2026-08-19")
})


# ---------------------------------------------------------------- #
# Engine death -> typed refusal, every field                        #
# ---------------------------------------------------------------- #

test_that("mocked engine death becomes inla_program_failed with every field", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)

  testthat::local_mocked_bindings(
    .inla_call = function(...) {
      stop(paste0(
        "The inla-program exited with an error. Unless you interrupted ",
        "it yourself, this is a bug, please report it to <...>."
      ))
    },
    .package = "flexyBayes"
  )

  expect_error(
    flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE),
    class = "flexybayes_refusal_inla_program_failed"
  )

  err <- tryCatch(
    flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE),
    flexybayes_refusal_inla_program_failed = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_inla_program_failed")
  expect_s3_class(err, "flexybayes_inla_engine_refusal")
  expect_s3_class(err, "flexybayes_refusal")
  expect_identical(err$reason_code, "inla_program_failed")

  # Design size.
  expect_identical(err$n_rows, 30)
  expect_false(is.na(err$n_latent_effects_est))
  expect_identical(err$n_latent_effects_est, 5) # g has 5 levels

  # Binding term (largest single random-effect contributor).
  expect_identical(err$binding_term, "g")
  expect_identical(err$binding_n_levels, 5)

  # Largest verified size, read from the C2 CSV row.
  expect_identical(err$largest_verified_n, 1e6)
  expect_true(is.na(err$largest_verified_random_effects))

  # The engine's own message is carried, and the rendered message
  # states the remedies.
  expect_match(err$engine_message, "inla-program exited", fixed = TRUE)
  expect_match(err$message, "inla-program exited", fixed = TRUE)
  expect_match(err$message, "Remedies", fixed = TRUE)
  expect_match(err$message, "1,000,000", fixed = TRUE)
  expect_match(err$message, "g (5 levels)", fixed = TRUE)
})

test_that("a NULL/empty INLA result with no R error also becomes inla_program_failed", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)

  testthat::local_mocked_bindings(
    .inla_call = function(...) NULL,
    .package = "flexyBayes"
  )

  err <- tryCatch(
    flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE),
    flexybayes_refusal_inla_program_failed = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_inla_program_failed")
  expect_identical(err$reason_code, "inla_program_failed")
  expect_match(err$engine_message, "NULL", fixed = TRUE)
})

test_that("a non-death engine error is NOT wrapped -- flexyBayes's own code propagates unchanged", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)

  testthat::local_mocked_bindings(
    .inla_call = function(...) stop("object 'totally_bogus_var' not found"),
    .package = "flexyBayes"
  )

  err <- tryCatch(
    flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE),
    error = function(e) e
  )
  expect_false(inherits(err, "flexybayes_refusal_inla_program_failed"))
  expect_false(inherits(err, "flexybayes_refusal"))
  expect_match(conditionMessage(err), "totally_bogus_var", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# Sanity: the real fixture still fits cleanly (unmocked)             #
# ---------------------------------------------------------------- #

test_that("the shared fixture still fits cleanly with the real INLA engine (regression guard)", {
  skip_if_no_inla()
  d <- mk_inla_data()
  fb <- mk_fb_for_inla(d)
  fit <- flexyBayes:::emit_inla(fb = fb, data = d, verbose = FALSE)
  expect_s3_class(fit, "flexybayes_inla")
})
