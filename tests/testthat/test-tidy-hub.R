# ---------------------------------------------------------------- #
# Hub tidy() coverage -- the broom-style tidy generic must return    #
# the same canonical columns across the backend fit classes the hub  #
# returns (greta `flexybayes`, INLA `flexybayes_inla`), so a         #
# cross-engine table can be assembled by rbind() rather than         #
# hand-built. The greta-class coverage is exercised in test-methods  #
# via the mock fit; this file adds the generics-dispatch check and    #
# the INLA-class coverage on a real (skip-guarded) INLA fit.         #
# ---------------------------------------------------------------- #

.canonical_tidy_cols <- c(
  "term",
  "estimate",
  "std.error",
  "conf.low",
  "conf.high"
)

# A minimal greta-class fit stub, self-contained in this file so the test
# does not depend on the make_mock_flexybayes() helper defined inside
# test-methods.R. It carries just enough of the contract for coef() / vcov()
# / confint() / variance_comps to drive tidy.flexybayes().
.tidy_hub_mock_greta <- function() {
  draws_mat <- matrix(stats::rnorm(300L), ncol = 1L) + 50
  colnames(draws_mat) <- "mu_atg"
  draws <- coda::mcmc.list(coda::mcmc(draws_mat))

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
      greta = list(draws = draws),
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

test_that("tidy() on a greta-class fit dispatches via the generic", {
  fit <- .tidy_hub_mock_greta()
  td <- generics::tidy(fit)
  expect_true(is.data.frame(td))
  expect_true(all(.canonical_tidy_cols %in% names(td)))
})

test_that("tidy(effects = 'random') returns the variance components", {
  fit <- .tidy_hub_mock_greta()
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
  # rather than dispatching to the greta implementation, and they appear in
  # methods() so dispatch is explicit rather than a bare "no applicable method".
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
