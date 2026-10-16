# =============================================================================
# Parts of a brms::bf() object that the brms ingest does not carry.
#
# The ingest reads the formula inside bf() and nothing else. Before 0.14.0 a
# formula for a distributional parameter (`sigma ~ x`), a fixed value for one
# (`sigma = 2`) and a family named inside bf() were dropped, and both engines
# fitted a different model with no message. A non-linear formula
# (`nl = TRUE`) stopped with the untyped "object 'a' not found". Each now
# refuses by name, on every public route into a fit: fb() and flexybayes()
# on each backend, fb_plan(), update(), the deprecated wrappers, and
# fb_select() on a base whose recorded model carries such a part.
# =============================================================================

fbbf_data <- function() {
  set.seed(21L)
  d <- expand.grid(rep = factor(1:2), gen = factor(1:6), env = factor(1:3))
  d$x <- stats::rnorm(nrow(d))
  d$z <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d), 10, 2)
  d$count <- stats::rpois(nrow(d), 5)
  d
}

fbbf_class <- "flexybayes_refusal_brms_ingest_feature_unsupported"

# Every route from a formula to a fit, for one brms::bf() object. The brms
# fit uses a short chain: it only runs if the refusal does not come first.
fbbf_routes <- function(f, d, family = "gaussian") {
  list(
    fb_auto = function() fb(f, data = d, family = family, verbose = FALSE),
    fb_inla = function() {
      fb(f, data = d, family = family, backend = "inla", verbose = FALSE)
    },
    fb_brms = function() {
      fb(
        f,
        data = d,
        family = family,
        backend = "brms",
        chains = 1L,
        n_samples = 200L,
        warmup = 100L,
        seed = 1L,
        verbose = FALSE
      )
    },
    fb_brms_code = function() {
      fb(
        f,
        data = d,
        family = family,
        backend = "brms",
        return_code = TRUE,
        verbose = FALSE
      )
    },
    flexybayes_inla = function() {
      flexybayes(
        f,
        data = d,
        family = family,
        backend = "inla",
        verbose = FALSE
      )
    },
    fb_plan = function() fb_plan(f, data = d, family = family),
    fb_from_brms = function() {
      suppressWarnings(fb_from_brms(f, data = d, family = family))
    },
    fb_inla_wrapper = function() {
      suppressWarnings(fb_inla(f, data = d, family = family, verbose = FALSE))
    },
    fb_brms_wrapper = function() {
      suppressWarnings(fb_brms(
        f,
        data = d,
        family = family,
        return_code = TRUE,
        verbose = FALSE
      ))
    }
  )
}

fbbf_expect_refused <- function(f, d, pattern, family = "gaussian") {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  routes <- fbbf_routes(f, d, family)
  for (nm in names(routes)) {
    err <- tryCatch(routes[[nm]](), error = function(e) e)
    expect_s3_class(err, fbbf_class)
    if (inherits(err, "error")) {
      expect_match(conditionMessage(err), pattern, fixed = TRUE, info = nm)
    } else {
      fail(paste0("route ", nm, " returned a value, not a refusal"))
    }
  }
}


test_that("a formula for a distributional parameter refuses on every route", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  d <- fbbf_data()
  fbbf_expect_refused(brms::bf(y ~ x + (1 | gen), sigma ~ x), d, "sigma ~ x")
})


test_that("a fixed value for a distributional parameter refuses on every route", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  d <- fbbf_data()
  fbbf_expect_refused(brms::bf(y ~ x + (1 | gen), sigma = 2), d, "sigma = 2")
})


test_that("a brms non-linear formula refuses by name on every route", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  d <- fbbf_data()
  fbbf_expect_refused(
    brms::bf(y ~ a * exp(x), a ~ 1 + (1 | gen), nl = TRUE),
    d,
    "non-linear"
  )
})


test_that("a family inside brms::bf() that differs from `family` refuses", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  d <- fbbf_data()
  fbbf_expect_refused(
    brms::bf(count ~ x + (1 | gen), family = stats::poisson()),
    d,
    "poisson"
  )
})


test_that("a family inside brms::bf() that agrees with `family` is admitted", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- fbbf_data()
  expect_s3_class(
    fb_plan(
      brms::bf(count ~ x + (1 | gen), family = stats::poisson()),
      data = d,
      family = "poisson"
    ),
    "fb_plan"
  )
  expect_s3_class(
    fb_plan(brms::bf(y ~ x + (1 | gen), family = stats::gaussian()), data = d),
    "fb_plan"
  )
})


test_that("update() into a formula with an unfitted part refuses", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- fbbf_data()
  base <- fb(
    brms::bf(y ~ x + (1 | gen)),
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  expect_s3_class(base, "flexybayes")
  expect_error(
    update(base, fixed = brms::bf(y ~ x + (1 | gen), sigma ~ x)),
    class = fbbf_class
  )
})
