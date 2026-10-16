# The reading functions a package that extends flexyBayes uses in place of
# the package's internals (?fb_extension), and the two-part brms read side.

fbext_quiet <- function(env = parent.frame()) {
  withr::local_options(
    list(
      flexyBayes.inla_threads = "1:1",
      flexyBayes.silence_default_prior_note = TRUE
    ),
    .local_envir = env
  )
}

fbext_data <- function() {
  withr::with_seed(20261007L, {
    d <- data.frame(
      site = factor(rep(sprintf("S%02d", 1:8), each = 6)),
      x = stats::rnorm(48),
      w = rep(c(1, 0.5), 24),
      v = stats::runif(48, 0.05, 0.2)
    )
    d$y <- 1 + 0.5 * d$x + stats::rnorm(8)[d$site] + stats::rnorm(48)
    d
  })
}

test_that("the extension API is version 1", {
  expect_identical(fb_extension_api(), 1L)
})

test_that("fb_response_values() reads a column or evaluates an expression", {
  d <- data.frame(y = c(1, 4, 9))
  expect_identical(fb_response_values("y", d), c(1, 4, 9))
  expect_identical(fb_response_values("sqrt(y)", d), c(1, 2, 3))
  expect_null(fb_response_values("missing_column", d))
  expect_null(fb_response_values(NULL, d))
})

test_that("fb_offset_values() sums the offsets and refuses one that is not finite", {
  d <- data.frame(a = c(1, 2, 4), b = c(0, 1, 0))
  expect_null(fb_offset_values(character(0), d))
  expect_equal(
    fb_offset_values(c("offset(log(a))", "offset(b)"), d),
    log(d$a) + d$b
  )
  d$a[[2L]] <- 0
  expect_error(
    fb_offset_values("offset(log(a))", d),
    class = "flexybayes_refusal_fixed_offset_invalid"
  )
  expect_error(
    fb_offset_values("offset(log(area))", d, newdata = TRUE),
    class = "flexybayes_refusal_newdata_offset_missing"
  )
})

test_that("fb_engine_settings() returns recognised settings and refuses others", {
  expect_identical(fb_engine_settings(NULL), list())
  expect_identical(
    fb_engine_settings(list(chains = 2L, adapt_delta = 0.95)),
    list(chains = 2L, adapt_delta = 0.95)
  )
  expect_error(
    fb_engine_settings(list(chains = 2L, iter = 10L)),
    "unrecognised setting 'iter'",
    fixed = TRUE
  )
  expect_error(fb_engine_settings(list(2L)), "named list")
  engine <- suppressWarnings(fb_engine("brms", chains = 2L))
  expect_identical(fb_engine_settings(engine), list(chains = 2L))
})

test_that("fb_fixed_design() gives the columns, their terms and ASReml's order", {
  d <- data.frame(
    a = factor(c("a1", "a2", "a3", "a1", "a2", "a3")),
    b = factor(c("b1", "b1", "b1", "b2", "b2", "b2")),
    x = c(1, 3, 2, 5, 4, 6)
  )
  design <- fb_fixed_design(~ x + a:b, d)
  expect_identical(design$X, stats::model.matrix(~ x + a:b, d))
  expect_identical(
    design$column_terms,
    c("(Intercept)", "x", rep("a:b", 6L))
  )
  # Within a:b, descending levels with the first factor slowest.
  ordered <- colnames(design$X)[design$asreml_order]
  expect_identical(ordered[1:2], c("(Intercept)", "x"))
  expect_identical(
    ordered[3:8],
    c("aa3:bb2", "aa3:bb1", "aa2:bb2", "aa2:bb1", "aa1:bb2", "aa1:bb1")
  )

  # A transform learned on the training rows alone.
  scaled <- fb_fixed_design(~ scale(x), d, train_rows = 1:3)
  expect_equal(
    unname(scaled$X[, 2L]),
    (d$x - mean(d$x[1:3])) / stats::sd(d$x[1:3])
  )
})

test_that("draw names of a two-part fit lose the part prefix brms writes", {
  nm <- c(
    "b_a_Intercept",
    "b_a_regionR2",
    "b_b_x01",
    "b_b_a_soil",
    "sdb_b_x01",
    "hs_global_b",
    "sd_cluster__a_Intercept",
    "r_cluster__a[C01,Intercept]",
    "sigma"
  )
  expect_identical(
    .fb_nl_canonical_names(nm, c(.FB_NL_PROTECTED_PART, .FB_NL_BLOCK_PART)),
    c(
      "b_Intercept",
      "b_regionR2",
      "b_x01",
      "b_a_soil",
      "sdb_x01",
      "hs_global",
      "sd_cluster__Intercept",
      "r_cluster[C01,Intercept]",
      "sigma"
    )
  )
})

test_that("fb_as_draws_simple() reads a brmsfit through the two-part names", {
  testthat::local_mocked_bindings(
    .fb_brms_draws_matrix = function(brmsfit) {
      matrix(1:4, ncol = 2L, dimnames = list(NULL, c("b_x", "sigma")))
    },
    .check_installed = function(...) invisible(TRUE)
  )
  out <- fb_as_draws_simple(structure(list(), class = "brmsfit"))
  expect_identical(out, list(b_x = c(1, 2), sigma = c(3, 4)))
})

test_that("the reading functions return what the fit was fitted with", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbext_quiet()
  d <- fbext_data()
  fit <- flexybayes(
    y ~ x,
    random = ~site,
    data = d,
    weights = d$w,
    backend = "inla",
    verbose = FALSE
  )
  terms <- fb_model_terms(fit)
  expect_identical(terms, fit$extras$fb_terms %||% fit$fb)
  expect_identical(terms$response, "y")
  expect_equal(fb_model_weights(fit), d$w)
  expect_equal(fb_model_weights(terms), d$w)
  expect_null(fb_known_variance(fit))

  verdict <- fb_convergence(fit)
  expect_identical(verdict$engine, "inla")
  expect_identical(verdict$ok, summary(fit)$converge$ok)

  draws <- fb_inla_posterior_sample(fit, 3L, seed = 1L, num.threads = "1:1")
  expect_length(draws, 3L)
  # On the scale of the data: the draws of the linear predictor sit around
  # the response, not around zero.
  predictor <- vapply(
    draws,
    function(z) {
      mean(z$latent[grepl("^Predictor:", rownames(z$latent)), 1L])
    },
    numeric(1L)
  )
  expect_lt(max(abs(predictor - mean(d$y))), 1)
  expect_error(fb_inla_posterior_sample(list(), 1L), "fit on INLA")

  plain <- flexybayes(
    y ~ x,
    random = ~site,
    data = d,
    known_variance = d$v,
    backend = "inla",
    verbose = FALSE
  )
  expect_null(fb_model_weights(plain))
  known <- fb_known_variance(plain)
  expect_equal(known$values, d$v)
  expect_null(known$matrix)
  expect_equal(plain$data[[known$column]], d$v)
  expect_true(known$inla_index %in% names(plain$inla$summary.random))
})

test_that("fb_check_refit() refuses a fit recorded from a formula it did not fit", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  fbext_quiet()
  d <- fbext_data()
  fit <- flexybayes(
    y ~ x,
    random = ~site,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  expect_identical(fb_check_refit(fit), fit)
  fit$extras$call_info$fixed <- brms::bf(y ~ x + (1 | site), sigma ~ x)
  err <- expect_error(
    fb_check_refit(fit, caller = "myverb(): "),
    class = "flexybayes_brms_ingest_refusal"
  )
  expect_match(conditionMessage(err), "^myverb\\(\\): ")
})
