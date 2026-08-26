# Tests for the C4 producer: fb_log_posterior().
#
# The native engine that produced a real C4 log-density was withdrawn
# entirely in 0.9.3 (see NEWS.md), and neither active engine (INLA,
# brms) can produce a faithful one either -- INLA's posterior is a
# deterministic Laplace/grid approximation with no evaluable natural-
# scale log-density, and brms's Stan log-density lives on a version-
# fragile unconstrained scale. fb_log_posterior() therefore abstains
# unconditionally with a classed `fb_c4_unavailable` condition
# regardless of what `fit` carries; these tests pin that abstention
# for every reachable dispatch path. Run everywhere; no backend
# required.

# ---- abstain paths (no backend required) ------------------------------------

test_that("the brms backend abstains with a classed, informative condition", {
  brms_fit <- structure(
    list(),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  cnd <- tryCatch(fb_log_posterior(brms_fit), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
  expect_match(conditionMessage(cnd), "brms backend")
})

test_that("the INLA backend abstains with a classed, informative condition", {
  inla_fit <- structure(list(), class = c("flexybayes_inla", "list"))
  cnd <- tryCatch(fb_log_posterior(inla_fit), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
  expect_match(conditionMessage(cnd), "INLA")
})

test_that("a non-flexyBayes object abstains via the default method", {
  cnd <- tryCatch(fb_log_posterior(list(a = 1)), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
})

test_that("a bare flexybayes-classed fit with no engine-specific method abstains via the default", {
  # fb_log_posterior() has no bare `flexybayes` method -- only `.default`,
  # `.flexybayes_brms`, and `.flexybayes_inla` -- so an object carrying
  # only the parent class (e.g. an engine class this version does not
  # recognise) falls through to `.default` rather than crashing.
  unrecognised <- structure(list(), class = c("flexybayes", "list"))
  cnd <- tryCatch(fb_log_posterior(unrecognised), error = function(e) e)
  expect_s3_class(cnd, "fb_c4_unavailable")
})
