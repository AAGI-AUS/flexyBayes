# Telling a broken compiler apart from a broken model.
#
# The brms path compiles a Stan program. When the C++ toolchain is missing
# or mismatched, the compiler writes its own message to the console and
# what reaches R is a downstream symptom -- `invalid connection`. Reported
# verbatim that sends the reader looking for a fault in their model, and
# the R Journal review of 0.10.0 arrived carrying exactly that reading.
#
# The diagnosis is not a pattern match against error strings from compilers
# this machine cannot produce. It is a probe: compile a three-line program
# and see. So both branches are tested, and the branch that matters is the
# NEGATIVE one -- a healthy toolchain must not be blamed for a fit that
# failed for some other reason.

test_that("the probe compiles a three-line program where a toolchain exists", {
  skip_on_cran()
  skip_if_not_installed("pkgbuild")
  skip_if_not(isTRUE(pkgbuild::has_build_tools(debug = FALSE)),
              "no C++ build tools on this machine")
  expect_true(.fb_cxx_toolchain_works())
})

test_that("a healthy toolchain is not blamed for an unrelated failure", {
  # The falsification. A model that fails for its own reasons must come
  # back with the engine's own words and NO toolchain claim, or the
  # diagnostic would send every brms user to reinstall a compiler that
  # works.
  testthat::local_mocked_bindings(.fb_cxx_toolchain_works = function() TRUE)
  err <- tryCatch(
    .fb_stop_brms_failure("the model has more parameters than observations"),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_false(inherits(err, "flexybayes_refusal"))
  expect_match(conditionMessage(err), "more parameters than observations")
  expect_false(grepl("toolchain", conditionMessage(err)))
  expect_false(grepl("xcode-select", conditionMessage(err), fixed = TRUE))
})

test_that("a broken toolchain is named, with the repair and the way round", {
  testthat::local_mocked_bindings(.fb_cxx_toolchain_works = function() FALSE)
  err <- tryCatch(
    .fb_stop_brms_failure("invalid connection"),
    error = function(e) e
  )
  msg <- conditionMessage(err)
  # It says where the fault is NOT, because that is what the reader has
  # been guessing at.
  expect_match(msg, "not in the model", fixed = TRUE)
  # It says how to check, and how to repair on the platform this shows up
  # on most.
  expect_match(msg, "check_build_tools", fixed = TRUE)
  expect_match(msg, "xcode-select --install", fixed = TRUE)
  # It names the engine that needs no compiler at all.
  expect_match(msg, "backend = \"inla\"", fixed = TRUE)
  # And it still carries what the engine actually said.
  expect_match(msg, "invalid connection", fixed = TRUE)
})

test_that("the toolchain failure is a registered, typed refusal", {
  testthat::local_mocked_bindings(.fb_cxx_toolchain_works = function() FALSE)
  err <- tryCatch(
    .fb_stop_brms_failure("invalid connection"),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal")
  expect_identical(err$reason_code, "brms_stan_toolchain_compile_failure")
  expect_true(
    "brms_stan_toolchain_compile_failure" %in% fb_refusals()$reason_code
  )
})

test_that("a probe that cannot run at all does not accuse the toolchain", {
  # NA means the probe itself could not be carried out (no writable temp
  # directory, say). An answer that was never obtained is not evidence of
  # a fault, so the engine's own message stands.
  testthat::local_mocked_bindings(.fb_cxx_toolchain_works = function() NA)
  err <- tryCatch(
    .fb_stop_brms_failure("invalid connection"),
    error = function(e) e
  )
  expect_false(inherits(err, "flexybayes_refusal"))
  expect_match(conditionMessage(err), "brms::brm() failed: invalid connection",
               fixed = TRUE)
})
