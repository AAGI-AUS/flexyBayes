# ---------------------------------------------------------------- #
# Hub tidy() coverage -- the broom-style tidy generic must return    #
# the same canonical columns across the backend fit classes the hub  #
# returns (the bare `flexybayes` method, reached via inheritance by  #
# any engine with no more specific override -- currently brms; and   #
# INLA's own `flexybayes_inla`), so a cross-engine table can be       #
# assembled by rbind() rather than hand-built. The bare-method       #
# coverage is exercised in test-methods via the mock fit; this file  #
# adds the generics-dispatch check and the INLA-class coverage on a  #
# real (skip-guarded) INLA fit.                                      #
# ---------------------------------------------------------------- #

.canonical_tidy_cols <- c(
  "term",
  "estimate",
  "std.error",
  "conf.low",
  "conf.high"
)

# A minimal bare-`flexybayes`-class fit stub, self-contained in this file
# so the test does not depend on the make_mock_flexybayes() helper defined
# inside test-methods.R. It carries just enough of the contract for
# coef() / vcov() / confint() / variance_comps to drive
# tidy.flexybayes() -- the bare method, reached via inheritance by any
# engine with no more specific override.
.tidy_hub_mock_bare <- function() {
  glm_obj <- list(coefficients = c("(Intercept)" = 50.1))
  attr(glm_obj, "posterior_vcov") <- matrix(
    0.5, 1L, 1L,
    dimnames = list("(Intercept)", "(Intercept)")
  )
  class(glm_obj) <- c("flexybayes_glm", "glm", "lm")

  vc <- data.frame(
    component = c("sigma_g", "sigma_e"),
    estimate = c(1.0, 2.0),
    sd = c(0.3, 0.2),
    q2.5 = c(0.5, 1.5),
    q97.5 = c(1.5, 2.5),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      glm = glm_obj,
      extras = list(
        variance_comps = vc,
        parse_info = list(
          fixed = list(response = "y", intercept = TRUE, terms = list())
        )
      )
    ),
    class = "flexybayes"
  )
}

test_that("the tidy generic is re-exported and dispatches by class", {
  # generics::tidy is re-exported, so a flexyBayes-only session can call
  # tidy() without attaching broom.
  expect_true(is.function(tidy))
  expect_true("tidy.flexybayes" %in% as.character(utils::methods("tidy")))
  expect_true("tidy.flexybayes_inla" %in% as.character(utils::methods("tidy")))
})

test_that("tidy() on a bare-flexybayes fit dispatches via the generic", {
  # conf.int = FALSE: confint.flexybayes() (the bare parent method, which
  # a mock of this class dispatches to) refuses unconditionally with
  # fit_lacks_posterior_draws -- it no longer reads a draws slot off any
  # engine, so no hand-built mock can satisfy it. The conf.int = TRUE
  # default path (which a real brms fit reaches via its own
  # confint.flexybayes_brms() override, not this one) is covered on a
  # live fit below.
  fit <- .tidy_hub_mock_bare()
  td <- generics::tidy(fit, conf.int = FALSE)
  expect_true(is.data.frame(td))
  expect_true(all(c("term", "estimate", "std.error") %in% names(td)))
})

test_that("tidy() on a live brms fit returns the full canonical column set", {
  skip_on_cran()
  skip_if_not_installed("brms")

  set.seed(12L)
  dat <- data.frame(y = rnorm(30L, 2 + 0.5 * (1:30), 1), x = 1:30)
  fit <- flexybayes(
    y ~ x,
    data = dat,
    backend = "brms",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    seed = 20260427L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )

  # Default conf.int = TRUE dispatches through confint.flexybayes_brms(),
  # the live override reading actual posterior draws -- the case the
  # bare-mock test above cannot exercise.
  td <- tidy(fit)
  expect_true(is.data.frame(td))
  expect_true(all(.canonical_tidy_cols %in% names(td)))
  expect_true(all(is.finite(td$conf.low)))
  expect_true(all(is.finite(td$conf.high)))
  expect_true(all(td$conf.low <= td$conf.high))
})

test_that("tidy(effects = 'random') returns the variance components", {
  fit <- .tidy_hub_mock_bare()
  td <- tidy(fit, effects = "random")
  expect_true(all(
    c("term", "estimate", "std.error") %in% names(td)
  ))
  expect_equal(nrow(td), 2L)
})

test_that("tidy.flexybayes_inla returns the canonical columns", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  set.seed(11L)
  dat <- data.frame(
    y = rnorm(60L, 2 + 0.5 * (1:60), 1),
    x = 1:60
  )
  fit <- tryCatch(
    flexybayes(y ~ x, data = dat, backend = "inla", verbose = FALSE),
    error = function(e) {
      testthat::skip(paste("INLA fit unavailable:", conditionMessage(e)))
    }
  )
  skip_if_not(inherits(fit, "flexybayes_inla"))

  td <- tidy(fit)
  expect_true(is.data.frame(td))
  expect_true(all(.canonical_tidy_cols %in% names(td)))
  # An intercept-and-slope model tidies to two terms.
  expect_gte(nrow(td), 2L)
  # The slope estimate is positive and its interval is finite.
  slope <- td[td$term == "x", ]
  expect_equal(nrow(slope), 1L)
  expect_true(is.finite(slope$estimate))
  expect_lt(slope$conf.low, slope$conf.high)
})

test_that("tidy.flexybayes_inla on an empty fixed-summary returns 0 rows", {
  # A degenerate INLA fit shell with no fixed effects tidies to an empty
  # but well-typed frame rather than erroring.
  shell <- structure(
    list(inla = list(summary.fixed = NULL)),
    class = c("flexybayes_inla", "list")
  )
  td <- tidy(shell)
  expect_true(is.data.frame(td))
  expect_equal(nrow(td), 0L)
  expect_true(all(.canonical_tidy_cols %in% names(td)))
})

test_that("glance() / augment() on an INLA fit refuse with an informative error", {
  # INLA fits (`flexybayes_inla`) support tidy() but not glance()/augment().
  # The methods raise an actionable error (pointing to tidy() / summary())
  # rather than dispatching to the bare method's implementation, and they
  # appear in methods() so dispatch is explicit rather than a bare "no
  # applicable method".
  expect_true(
    "glance.flexybayes_inla" %in% as.character(utils::methods("glance"))
  )
  expect_true(
    "augment.flexybayes_inla" %in% as.character(utils::methods("augment"))
  )

  shell <- structure(list(), class = c("flexybayes_inla", "list"))
  expect_error(generics::glance(shell), "not available for INLA")
  expect_error(generics::augment(shell), "not available for INLA")
})

# ---------------------------------------------------------------- #
# effects = "random" answers on both engines                        #
# ---------------------------------------------------------------- #
#
# `tidy.flexybayes_inla()` had no `effects` argument, so the request was
# absorbed by `...` and the FIXED-effect table came back under a call that
# asked for variance components. Not an error -- a wrong answer, on the
# engine that is the package's default route.

test_that("tidy(effects = 'random') returns variance components on INLA", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  set.seed(31L)
  dat <- data.frame(
    g = factor(rep(seq_len(8L), each = 10L)),
    x = stats::rnorm(80L)
  )
  dat$y <- 2 + 0.5 * dat$x +
    rep(stats::rnorm(8L, 0, 0.8), each = 10L) +
    stats::rnorm(80L, 0, 0.5)
  fit <- suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = dat,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  ))

  fixed <- tidy(fit)
  random <- tidy(fit, effects = "random")

  # Same columns on both requests, per the broom contract.
  expect_true(all(.canonical_tidy_cols %in% names(random)))
  # And genuinely different rows -- the defect was that they were equal.
  expect_false(identical(fixed$term, random$term))
  expect_true("sigma" %in% random$term)
  expect_true("sd_g" %in% random$term)
  expect_false(any(random$term %in% c("(Intercept)", "x")))
  expect_true(all(is.finite(random$estimate)))

  # The two engines answer the same question with the same shape.
  expect_identical(names(tidy(fit, effects = "fixed")), names(fixed))
})

test_that("an unrecognised effects value stops rather than falling through", {
  shell <- structure(
    list(inla = list(summary.fixed = NULL)),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  expect_error(tidy(shell, effects = "variance"), "arg")
  expect_error(
    generics::tidy(
      structure(list(), class = c("flexybayes", "list")),
      effects = "variance"
    ),
    "arg"
  )
})

test_that("tidy(effects = 'random') works on an aggregated fit", {
  # The aggregated emits write a list of posterior means under
  # `variance_comps`, so `nrow()` of it is NULL and the raw read errored.
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  set.seed(32L)
  dat <- data.frame(
    g = factor(rep(seq_len(10L), each = 50L)),
    f1 = factor(rep(c("a", "b"), length.out = 500L))
  )
  dat$y <- stats::rnorm(500L)
  fit <- suppressMessages(fb(
    y ~ f1 + (1 | g),
    data = dat,
    backend = "inla",
    aggregate = "auto",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  skip_if_not(inherits(fit, "flexybayes_aggregated"))

  random <- tidy(fit, effects = "random")
  expect_true(all(.canonical_tidy_cols %in% names(random)))
  expect_gt(nrow(random), 0L)
})
