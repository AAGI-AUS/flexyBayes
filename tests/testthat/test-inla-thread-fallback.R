# A host may report no core count. INLA turns that into `NA:1`, emits a parser
# failure and silently chooses another plan. The package supplies `1:1` only
# for that missing-host-information case; valid defaults and explicit caller
# settings remain outside this fallback.

test_that("an unknown host core count has a valid conservative INLA fallback", {
  expect_identical(.fb_inla_threads_if_cores_unknown(NA_integer_), "1:1")
  expect_identical(.fb_inla_threads_if_cores_unknown(integer()), "1:1")
  expect_identical(.fb_inla_threads_if_cores_unknown(0L), "1:1")
  expect_identical(.fb_inla_threads_if_cores_unknown(-1L), "1:1")
  expect_null(.fb_inla_threads_if_cores_unknown(1L))
  expect_null(.fb_inla_threads_if_cores_unknown(8L))
})

test_that("this host never hands INLA the malformed NA thread spelling", {
  value <- .fb_inla_threads_if_cores_unknown()
  if (is.null(value)) {
    expect_true(parallel::detectCores(logical = TRUE) >= 1L)
  } else {
    expect_identical(value, "1:1")
  }
  expect_false(identical(value, "NA:1"))
})
