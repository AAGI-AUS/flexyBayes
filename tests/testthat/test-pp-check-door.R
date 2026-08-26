# =============================================================================
# pp_check() on a flexyBayes fit, and the plot() type of the same name.
#
# A posterior predictive check simulates datasets from the fitted model
# and shows them against the observed response. Until 0.9.1 the package
# claimed exactly that in the plot() documentation and drew something
# else -- observed values against fitted values, plus a density of the
# fitted values -- on every backend that carried both, while an INLA fit
# got a message with no reason code that no caller could catch.
#
# Both halves are fixed here, and the tests hold both:
#
#   * on a brms fit the check is the real one, delegated to
#     brms::pp_check(), whose default overlays the densities of
#     replicated datasets on the density of the response;
#   * on a fit with no predictive draws, BOTH entry points raise the same
#     catchable refusal, and the refusal names the residual displays the
#     fit can draw;
#   * no path draws anything under the pp_check name that is not a
#     posterior predictive check, so the documentation sentence is true
#     wherever it applies.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

.ppc_data <- function(seed = 20260817L, n = 60L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b"), length.out = n)),
    g = factor(rep(letters[1:5], length.out = n)),
    x = stats::rnorm(n)
  )
  b_g <- stats::rnorm(5L, sd = 0.8)[as.integer(d$g)]
  d$y <- 1 + 0.6 * (d$f == "b") + 0.4 * d$x + b_g + stats::rnorm(n, sd = 0.5)
  d
}

.ppc_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

.ppc_cache <- new.env(parent = emptyenv())

.ppc_inla_fit <- function() {
  if (is.null(.ppc_cache$inla)) {
    .ppc_cache$inla <- suppressMessages(fb(
      y ~ f + x,
      random = ~g,
      data = .ppc_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .ppc_cache$inla
}

.ppc_brms_fit <- function() {
  if (is.null(.ppc_cache$brms)) {
    .ppc_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x + (1 | g),
      data = .ppc_data(),
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .ppc_cache$brms
}


# ---------------------------------------------------------------- #
# 1. The real check, on a fit that can answer it                    #
# ---------------------------------------------------------------- #

test_that("pp_check() on a brms fit draws the real check", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bayesplot")
  skip_on_cran()
  skip_on_ci()
  .ppc_silence()
  fit <- .ppc_brms_fit()

  p <- suppressMessages(suppressWarnings(bayesplot::pp_check(fit)))

  expect_s3_class(p, "ggplot")
  # Replicated datasets, not point estimates: brms's default check is
  # the density overlay, one line per replicate plus the observed one.
  expect_gt(length(p$layers), 0L)
})

test_that("pp_check() forwards its arguments to brms", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bayesplot")
  skip_on_cran()
  skip_on_ci()
  .ppc_silence()
  fit <- .ppc_brms_fit()

  p <- suppressMessages(suppressWarnings(
    bayesplot::pp_check(fit, type = "hist", ndraws = 4L)
  ))
  expect_s3_class(p, "ggplot")

  # An unknown type is brms's own refusal, reached because the arguments
  # travel: the method does not silently drop what it was given.
  expect_error(
    suppressMessages(bayesplot::pp_check(fit, type = "not_a_ppc_type")),
    "valid ppc type"
  )
})

test_that("plot(type = 'pp_check') forwards to the same check on brms", {
  skip_if_not_installed("brms")
  skip_if_not_installed("bayesplot")
  skip_on_cran()
  skip_on_ci()
  .ppc_silence()
  fit <- .ppc_brms_fit()

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- suppressMessages(suppressWarnings(plot(fit, type = "pp_check")))

  expect_s3_class(p, "ggplot")
})


# ---------------------------------------------------------------- #
# 2. The refusal, on a fit that cannot                              #
# ---------------------------------------------------------------- #

test_that("pp_check() on an INLA fit refuses by name", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("bayesplot")
  skip_on_cran()
  skip_on_ci()
  .ppc_silence()
  fit <- .ppc_inla_fit()

  err <- tryCatch(bayesplot::pp_check(fit), condition = function(e) e)

  expect_s3_class(err, "flexybayes_refusal_pp_check_requires_predictive_draws")
  msg <- conditionMessage(err)
  expect_match(msg, "plot(fit, type = \"residuals\")", fixed = TRUE)
  expect_match(msg, "plot(fit, type = \"variogram\")", fixed = TRUE)
  expect_match(msg, "replicated datasets", fixed = TRUE)
})

test_that("plot(type = 'pp_check') raises the same catchable condition", {
  # The 0.9.0 behaviour was message() plus invisible(NULL): no reason
  # code, nothing to catch, and a caller checking the return value could
  # not tell a declined display from a drawn one.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .ppc_silence()
  fit <- .ppc_inla_fit()

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot(fit, type = "pp_check"),
    class = "flexybayes_refusal_pp_check_requires_predictive_draws"
  )

  # Same condition from both doors, message included.
  a <- tryCatch(plot(fit, type = "pp_check"), condition = function(e) e)
  b <- tryCatch(
    flexyBayes:::pp_check.flexybayes(fit),
    condition = function(e) e
  )
  expect_identical(conditionMessage(a), conditionMessage(b))
})

test_that("the legacy observed-versus-fitted panel is gone, not renamed", {
  # A fit carrying a response and fitted values but no predictive draws
  # is exactly the case that used to draw the old panel. It refuses now,
  # so nothing is drawn under the pp_check name that is not a check.
  shape <- structure(
    list(
      glm = structure(
        list(y = rnorm(10L), fitted.values = rnorm(10L)),
        class = c("flexybayes_glm", "list")
      )
    ),
    class = c("flexybayes", "list")
  )

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot(shape, type = "pp_check"),
    class = "flexybayes_refusal_pp_check_requires_predictive_draws"
  )
})


# ---------------------------------------------------------------- #
# 3. The Suggests guards                                            #
# ---------------------------------------------------------------- #

test_that("brms absent refuses with a written message", {
  local_mocked_bindings(
    .fb_brms_available = function() FALSE,
    .package = "flexyBayes"
  )
  fake <- structure(
    list(brms = structure(list(), class = "brmsfit")),
    class = c("flexybayes_brms", "flexybayes")
  )
  err <- tryCatch(
    flexyBayes:::pp_check.flexybayes(fake),
    condition = function(e) e
  )
  msg <- conditionMessage(err)
  expect_match(msg, "pp_check()", fixed = TRUE)
  expect_match(msg, "install.packages(\"brms\")", fixed = TRUE)
  expect_false(grepl("there is no package called", msg, fixed = TRUE))
})

test_that("bayesplot absent refuses with a written message", {
  # Defensive rather than reachable in an ordinary install -- brms
  # imports bayesplot, so a fit that can be delegated normally arrives
  # with both. The guard exists because plot(fit, type = "pp_check")
  # reaches the delegation without loading the generic's home package.
  local_mocked_bindings(
    .fb_bayesplot_available = function() FALSE,
    .package = "flexyBayes"
  )
  fake <- structure(
    list(brms = structure(list(), class = "brmsfit")),
    class = c("flexybayes_brms", "flexybayes")
  )
  err <- tryCatch(
    flexyBayes:::pp_check.flexybayes(fake),
    condition = function(e) e
  )
  msg <- conditionMessage(err)
  expect_match(msg, "install.packages(\"bayesplot\")", fixed = TRUE)
  expect_match(msg, "plot(fit, type = \"residuals\")", fixed = TRUE)
  expect_false(grepl("there is no package called", msg, fixed = TRUE))
})


# ---------------------------------------------------------------- #
# 4. The taxonomy                                                   #
# ---------------------------------------------------------------- #

test_that("the refusal code is registered", {
  entries <- ls(flexyBayes:::.refusal_registry)
  expect_true("pp_check_requires_predictive_draws" %in% entries)

  entry <- flexyBayes:::.lookup_refusal("pp_check_requires_predictive_draws")
  expect_identical(entry$since_version, "0.9.1")
  expect_match(entry$message_template, "replicated datasets", fixed = TRUE)

  tab <- fb_refusals(reason_code = "pp_check_requires_predictive_draws")
  expect_identical(nrow(tab), 1L)
})
