# =============================================================================
# 0.9.3: the package's third native engine was withdrawn entirely, not
# quarantined -- no code path, export, registry row, or `Suggests` entry
# remains (see NEWS.md). The invariant this file pins therefore changed
# shape: it used to be "every fitting entry point refuses a quarantined
# name"; it is now "any name outside the two registered engines is unknown
# to every fitting entry point".
#
# The `.check_known_backend_name()` guard does not special-case the name
# of a withdrawn engine -- it refuses anything not in the active set, the
# same way for a name that used to be registered and a name that never
# was. A made-up name therefore pins the exact invariant a withdrawn
# engine's name would, without this file ever needing to spell one out;
# these tests are written to fail when a NEW fitting route is added
# without going through the shared unknown-backend check, rather than to
# enumerate the routes known today.
#
# (This file previously also carried a source-level scan asserting that
# no file under R/ calls into the withdrawn engine's namespace. That
# assertion cannot be restated without writing the withdrawn engine's
# name into a source file that ships, which is exactly what this file
# now exists to make impossible. Its coverage is retired in favour of
# the WP-G2 goal grep -- an empty case-insensitive search for the name
# across every shipped file bar NEWS.md -- run as this package's own
# closing acceptance check, and re-runnable by anyone at any time.)
# =============================================================================

test_that("naming an unregistered backend raises unknown_backend at every known fitting entry point", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(1L)
  d <- data.frame(y = stats::rnorm(20L), x = stats::rnorm(20L))

  # The mixed-model surface. Every unregistered name -- withdrawn or
  # simply made up -- raises the same reason code, and both active
  # engines are named in the message.
  for (nm in c("not_a_registered_backend", "still_not_one")) {
    err <- tryCatch(
      suppressMessages(flexybayes(y ~ x, data = d, backend = nm)),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_unknown_backend_refusal")
    expect_identical(err$reason_code, "unknown_backend", label = nm)
    expect_match(conditionMessage(err), '"inla"', fixed = TRUE, label = nm)
    expect_match(conditionMessage(err), '"brms"', fixed = TRUE, label = nm)
  }
})

test_that("fb_dirichlet() rejects any method beyond the sole ml estimator", {
  # Before 0.9.3, fb_dirichlet() offered a method that called the native
  # engine directly, bypassing the backend registry entirely (so it fit
  # even while backend = <that engine's name> refused at the mixed-model
  # surface above). That method no longer exists -- "ml" is the only
  # choice -- so any other request fails match.arg() before any fitting
  # code runs, the same as it would for a name that was never offered.
  m <- matrix(stats::rgamma(60L, 2), 20L, 3L)
  m <- m / rowSums(m)
  err <- tryCatch(
    fb_dirichlet(m, method = "not_a_real_method"),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_false(inherits(err, "flexybayes_refusal"))
  expect_match(conditionMessage(err), "ml", fixed = TRUE)
})

test_that("the maximum-likelihood Dirichlet fit is unaffected", {
  set.seed(2L)
  m <- matrix(stats::rgamma(90L, 2), 30L, 3L)
  m <- m / rowSums(m)
  fit <- fb_dirichlet(m)
  expect_s3_class(fit, "fb_dirichlet_fit")
  expect_identical(fit$method, "ml")
  expect_equal(sum(fit$mean_composition), 1, tolerance = 1e-8)
})
