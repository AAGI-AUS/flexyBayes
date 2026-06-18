# fb_backend_status() -- the backend-readiness diagnostic. It is read-only and
# runs without any backend installed, so these checks are environment-
# independent (no skip gates needed).

test_that("fb_backend_status() returns the documented shape", {
  st <- fb_backend_status()
  expect_s3_class(st, "fb_backend_status")
  expect_s3_class(st, "data.frame")
  expect_setequal(st$backend, c("greta", "INLA", "brms"))
  expect_identical(names(st), c("backend", "installed", "usable", "note"))
  expect_type(st$installed, "logical")
  expect_type(st$usable, "logical")
  expect_type(st$note, "character")
  expect_false(anyNA(st$installed))
  expect_false(anyNA(st$usable))
})

test_that("fb_backend_status() invariant: usable implies installed", {
  st <- fb_backend_status()
  # A backend cannot be usable unless its package is installed.
  expect_true(all(st$usable[!st$installed] == FALSE))
  # Every note carries an actionable, non-empty message.
  expect_true(all(nzchar(st$note)))
})

test_that("print.fb_backend_status() renders a table and returns invisibly", {
  st <- fb_backend_status()
  expect_output(print(st), "flexyBayes backend readiness")
  expect_identical(withVisible(print(st))$visible, FALSE)
})
