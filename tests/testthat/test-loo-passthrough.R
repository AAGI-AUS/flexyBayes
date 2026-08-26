# =============================================================================
# loo() on a flexyBayes fit.
#
# Approximate leave-one-out cross-validation needs the log-likelihood of
# each observation at each posterior draw. A Stan sampler stores it; a
# nested Laplace approximation never computed it. So the method does one
# of three things, and which one is decided by what the fit holds rather
# than by its class name:
#
#   * brms engine -- pass straight through to brms::loo() and return
#     loo's own object, diagnostics included;
#   * a fit with a posterior but no pointwise log-likelihood -- refuse as
#     loo_requires_sampler_draws, naming the information criteria the fit
#     DOES carry so the alternative is in the message;
#   * a fit with no posterior at all -- refuse as the sibling code
#     fit_lacks_posterior_draws, so the two states are told apart by the
#     condition class.
#
# The slot paths in the message are asserted against a live fit, not
# recited: the per-row INLA emit asks for WAIC and DIC at fit time and
# the aggregated emits do not, and a message that names a value the fit
# does not carry is worse than one that names none.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

.loo_data <- function(seed = 20260817L, n = 60L) {
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

.loo_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

.loo_cache <- new.env(parent = emptyenv())

.loo_inla_fit <- function() {
  if (is.null(.loo_cache$inla)) {
    .loo_cache$inla <- suppressMessages(fb(
      y ~ f + x,
      random = ~g,
      data = .loo_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .loo_cache$inla
}

.loo_brms_fit <- function() {
  if (is.null(.loo_cache$brms)) {
    .loo_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x + (1 | g),
      data = .loo_data(),
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .loo_cache$brms
}

# Two constructed shapes for the refusal routing. Neither needs an
# engine: what the method reads is which slots are present.
.loo_shape_no_posterior <- function() {
  structure(list(extras = list()), class = c("flexybayes_brms", "flexybayes"))
}

.loo_shape_draws_only <- function() {
  # `.fb_refuse_loo()`'s has_posterior check reads $inla or a top-level
  # $draws -- a third branch reading the withdrawn native engine's own
  # draws slot was removed with the engine at 0.9.3 (see NEWS.md);
  # nest the fixture's draws directly under $draws so it still
  # represents "has posterior, no log-likelihood" rather than "no
  # posterior at all".
  structure(
    list(draws = 1),
    class = c("flexybayes", "list")
  )
}


# ---------------------------------------------------------------- #
# 1. The brms passthrough                                           #
# ---------------------------------------------------------------- #

test_that("loo() on a brms fit returns loo's own object", {
  skip_if_not_installed("brms")
  skip_if_not_installed("loo")
  skip_on_cran()
  skip_on_ci()
  .loo_silence()
  fit <- .loo_brms_fit()

  out <- suppressWarnings(loo::loo(fit))

  expect_s3_class(out, "loo")
  expect_true(all(c("estimates", "pointwise", "diagnostics") %in% names(out)))
  expect_true(is.finite(out$estimates["elpd_loo", "Estimate"]))

  # The passthrough is a passthrough: the same object brms would have
  # returned from the inner fit, not a re-derivation.
  direct <- suppressWarnings(brms::loo(fit$brms))
  expect_equal(
    out$estimates["elpd_loo", "Estimate"],
    direct$estimates["elpd_loo", "Estimate"]
  )
})


# ---------------------------------------------------------------- #
# 2. The INLA refusal names what the fit does carry                 #
# ---------------------------------------------------------------- #

test_that("loo() on an INLA fit refuses by name and names the WAIC slot", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("loo")
  skip_on_cran()
  skip_on_ci()
  .loo_silence()
  fit <- .loo_inla_fit()

  err <- tryCatch(loo::loo(fit), condition = function(e) e)

  expect_s3_class(err, "flexybayes_refusal_loo_requires_sampler_draws")
  msg <- conditionMessage(err)
  expect_match(msg, "fit$inla$waic$waic", fixed = TRUE)
  expect_match(msg, "fit$inla$dic$dic", fixed = TRUE)
  expect_match(msg, "Laplace", fixed = TRUE)
  expect_match(msg, "log-likelihood", fixed = TRUE)

  # The named slots hold what the message says they hold: one number
  # each, and the message quotes the fit's own values.
  expect_true(is.numeric(fit$inla$waic$waic))
  expect_length(fit$inla$waic$waic, 1L)
  expect_true(is.numeric(fit$inla$dic$dic))
  expect_length(fit$inla$dic$dic, 1L)
  expect_match(msg, format(fit$inla$waic$waic, digits = 6L), fixed = TRUE)
  expect_match(msg, format(fit$inla$dic$dic, digits = 6L), fixed = TRUE)
})

test_that("an aggregated fit is not told it carries criteria it does not", {
  # The aggregated emits fit on cell-level sufficient statistics and ask
  # INLA for neither WAIC nor DIC (control.compute), so the same message
  # written from memory rather than from the fit would name two slots
  # that are empty.
  skip_if_not_installed("loo")
  agg <- structure(
    list(inla = list()),
    class = c("flexybayes_aggregated", "flexybayes_inla", "flexybayes")
  )

  err <- tryCatch(loo::loo(agg), condition = function(e) e)

  expect_s3_class(err, "flexybayes_refusal_loo_requires_sampler_draws")
  msg <- conditionMessage(err)
  expect_match(msg, "carries no information criterion either", fixed = TRUE)
  expect_match(msg, "aggregate = FALSE", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# 3. The two refusals are told apart by state, not by class name    #
# ---------------------------------------------------------------- #

test_that("a fit with no posterior reaches the sibling refusal", {
  skip_if_not_installed("loo")
  err <- tryCatch(
    loo::loo(.loo_shape_no_posterior()),
    condition = function(e) e
  )

  expect_s3_class(err, "flexybayes_refusal_fit_lacks_posterior_draws")
  expect_false(inherits(err, "flexybayes_refusal_loo_requires_sampler_draws"))
  expect_match(
    conditionMessage(err),
    "loo_requires_sampler_draws",
    fixed = TRUE
  )
})

test_that("a fit carrying draws but no log-likelihood reaches the other", {
  skip_if_not_installed("loo")
  err <- tryCatch(
    loo::loo(.loo_shape_draws_only()),
    condition = function(e) e
  )

  expect_s3_class(err, "flexybayes_refusal_loo_requires_sampler_draws")
  expect_false(inherits(err, "flexybayes_refusal_fit_lacks_posterior_draws"))
  expect_match(
    conditionMessage(err),
    "fit_lacks_posterior_draws",
    fixed = TRUE
  )
})


# ---------------------------------------------------------------- #
# 4. The Suggests guard on the passthrough                          #
# ---------------------------------------------------------------- #

test_that("brms absent refuses with a written message, not a namespace error", {
  # Exercised without uninstalling brms, by mocking the named
  # availability predicate the guard consults.
  local_mocked_bindings(
    .fb_brms_available = function() FALSE,
    .package = "flexyBayes"
  )
  fake <- structure(
    list(brms = structure(list(), class = "brmsfit")),
    class = c("flexybayes_brms", "flexybayes")
  )
  err <- tryCatch(
    flexyBayes:::loo.flexybayes(fake),
    condition = function(e) e
  )
  msg <- conditionMessage(err)
  expect_match(msg, "loo()", fixed = TRUE)
  expect_match(msg, "install.packages(\"brms\")", fixed = TRUE)
  expect_false(grepl("there is no package called", msg, fixed = TRUE))
})


# ---------------------------------------------------------------- #
# 5. The taxonomy                                                   #
# ---------------------------------------------------------------- #

test_that("the refusal code is registered and cross-references its sibling", {
  entries <- ls(flexyBayes:::.refusal_registry)
  expect_true("loo_requires_sampler_draws" %in% entries)

  new_entry <- flexyBayes:::.lookup_refusal("loo_requires_sampler_draws")
  old_entry <- flexyBayes:::.lookup_refusal("fit_lacks_posterior_draws")

  # AD-4: the two are siblings, and each says so, so the registry does
  # not read as one reason registered twice.
  expect_match(new_entry$description, "fit_lacks_posterior_draws", fixed = TRUE)
  expect_match(
    old_entry$description,
    "loo_requires_sampler_draws",
    fixed = TRUE
  )
  expect_identical(new_entry$since_version, "0.9.1")
  expect_match(new_entry$message_template, "fit$inla$waic$waic", fixed = TRUE)

  tab <- fb_refusals(reason_code = "loo_requires_sampler_draws")
  expect_identical(nrow(tab), 1L)
})
