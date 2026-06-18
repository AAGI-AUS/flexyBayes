# Tests for fb_greta() -- direct-greta-model entry to flexyBayes.
#
# Contract (ADR 0012 in design-notes-priors-brms-lgm.md):
# - fb_greta(model) returns a c("flexybayes",
#   "flexybayes_direct_greta", "list") object;
# - the IR is built post-hoc from the greta model graph with
#   source = "greta" and intercept = NA;
# - coef() returns canonical-named posterior means of the target
#   parameters;
# - print() emits the direct-greta header line;
# - predict() requires a user-supplied predictor function;
# - return_code = TRUE and review_code = TRUE both raise structured
#   refusals under the v0.2 contract (deferred to v0.3);
# - non-greta_model inputs raise structural refusals;
# - backend_decision(fit) returns the trivial direct-entry trace;
# - triangulate() works against another flexybayes fit (via the
#   parent-class fb_as_draws_simple method).
#
# Greta is required for the fitting subtests; class / refusal /
# IR-shape subtests run without it via a stub model object.


# Build a tiny greta_model for the integration subtests. Returns a
# list with `model` and `data` for downstream assertions.
mk_greta_toy <- function() {
  set.seed(20260522L)
  n <- 30L
  dat <- data.frame(
    x = stats::rnorm(n),
    y = stats::rnorm(n)
  )
  greta::as_data
  y <- greta::as_data(dat$y)
  x <- greta::as_data(dat$x)
  b0 <- greta::normal(0, 10)
  b1 <- greta::normal(0, 10)
  sigma <- greta::uniform(0, 5)
  mu <- b0 + b1 * x
  greta::distribution(y) <- greta::normal(mu, sigma)
  list(
    model = greta::model(b0, b1, sigma),
    data = dat
  )
}


# ---------------------------------------------------------------- #
# (a) Non-greta_model input raises a structural refusal             #
# ---------------------------------------------------------------- #

test_that("fb_greta() refuses a non-model input", {
  # v0.5.0: fb_greta() is the greta engine pin -- it accepts a formula,
  # a native greta_model, or a greta-source IR (fb_from_greta()). A bare
  # string is none of those, so it errors as an invalid model
  # specification (the `model = ` argument is remapped to the universal
  # model slot for call-compatibility).
  err <- tryCatch(
    fb_greta(model = "not a greta_model"),
    error = function(e) conditionMessage(e)
  )
  expect_true(nzchar(err))
  expect_true(grepl("formula|model", err, ignore.case = TRUE))
})


# ---------------------------------------------------------------- #
# (b) review_code / return_code structured refusals (v0.2)         #
# ---------------------------------------------------------------- #

test_that("fb_greta() refuses review_code on a native graph", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  # v0.5.0: a native greta graph carries no flexyBayes-generated code,
  # so the inspect-then-fit token does not apply -- the native-greta
  # dispatch refuses review_code with a clear message.
  err <- tryCatch(
    fb_greta(toy$model, data = toy$data, review_code = TRUE),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("review_code", err, fixed = TRUE))
  expect_true(grepl("native greta", err, fixed = TRUE))
})

test_that("fb_greta() defers return_code with a documented refusal", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  err <- tryCatch(
    fb_greta(toy$model, data = toy$data, return_code = TRUE),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("return_code", err, fixed = TRUE))
})

test_that("fb_greta() refuses return_code + review_code combination", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  err <- tryCatch(
    fb_greta(
      toy$model,
      data = toy$data,
      return_code = TRUE,
      review_code = TRUE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("mutually exclusive", err, fixed = TRUE))
})


# ---------------------------------------------------------------- #
# (c) IR ingest path: fb_from_greta() shape                          #
# ---------------------------------------------------------------- #

test_that("fb_from_greta() builds an IR with source = 'greta'", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  fb <- suppressMessages(
    fb_from_greta(toy$model, data = toy$data)
  )
  expect_s3_class(fb, "fb_terms")
  expect_identical(fb$source, "greta")
  expect_true(is.na(fb$intercept))
  expect_length(fb$fixed_terms, 0L)
  expect_length(fb$random_terms, 0L)
  expect_true(!is.null(fb$greta_meta))
  expect_setequal(names(fb$greta_meta$arrays), c("b0", "b1", "sigma"))
  expect_setequal(names(fb$greta_meta$canonical_map), c("b0", "b1", "sigma"))
  expect_equal(fb$greta_meta$n_data, 30L)
  expect_identical(fb$greta_meta$likelihood, "normal")
  expect_identical(fb$family, "gaussian")
})


# ---------------------------------------------------------------- #
# (d) canonical_names: validation + auto-fill                       #
# ---------------------------------------------------------------- #

test_that("fb_from_greta() validates canonical_names membership", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  err <- tryCatch(
    fb_from_greta(
      toy$model,
      data = toy$data,
      canonical_names = c(nonexistent = "X")
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("nonexistent", err, fixed = TRUE))
  expect_true(grepl("target", err))
})

test_that("fb_from_greta() applies canonical_names where supplied", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  fb <- suppressMessages(
    fb_from_greta(
      toy$model,
      data = toy$data,
      canonical_names = c(b0 = "(Intercept)", b1 = "Days")
    )
  )
  expect_identical(fb$greta_meta$canonical_map[["b0"]], "(Intercept)")
  expect_identical(fb$greta_meta$canonical_map[["b1"]], "Days")
  # Unmapped target falls back to verbatim.
  expect_identical(fb$greta_meta$canonical_map[["sigma"]], "sigma")
})


# ---------------------------------------------------------------- #
# (e) Full fit: class + coef + print + backend_decision             #
# ---------------------------------------------------------------- #

test_that("fb_greta() returns a flexybayes_direct_greta object", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  # v0.5.0: a canonical-name map is attached at IR-build time via
  # fb_from_greta(); the greta pin then fits the IR.
  ir <- suppressMessages(fb_from_greta(
    toy$model,
    data = toy$data,
    canonical_names = c(b0 = "(Intercept)", b1 = "Days", sigma = "sigma_e")
  ))
  fit <- suppressMessages(fb_greta(
    ir,
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  expect_s3_class(fit, "flexybayes")
  expect_s3_class(fit, "flexybayes_direct_greta")
  expect_true(inherits(fit$glm, "flexybayes_glm"))
  expect_true(!is.null(fit$greta$model))
  expect_true(!is.null(fit$greta$draws))

  beta <- coef(fit)
  expect_true(is.numeric(beta))
  expect_setequal(names(beta), c("(Intercept)", "Days", "sigma_e"))

  # print() smoke-test: no error; output contains the entry line.
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("direct-greta entry", out)))
  expect_true(any(grepl("\\(Intercept\\) -> .+|Canonical map", out)))

  # backend_decision() trivial trace.
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "greta")
  expect_identical(bd$path, "direct_entry")
  expect_true(grepl("bypasses lgm_gate", bd$reason))
})


# ---------------------------------------------------------------- #
# (f) predict() requires a user-supplied predictor                  #
# ---------------------------------------------------------------- #

test_that("predict.flexybayes_direct_greta() requires predictor", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  fit <- suppressMessages(fb_greta(
    toy$model,
    data = toy$data,
    n_samples = 50L,
    warmup = 50L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  err <- tryCatch(
    predict(fit, newdata = toy$data),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("predictor", err, fixed = TRUE))

  # With a valid predictor, predict() returns nrow(newdata) values.
  f <- function(theta, newdata) {
    theta[["b0"]] + theta[["b1"]] * newdata$x
  }
  preds <- predict(fit, newdata = toy$data, predictor = f, n_draws = 25L)
  expect_length(preds, nrow(toy$data))
  expect_true(all(is.finite(preds)))
})


# ---------------------------------------------------------------- #
# (g) triangulate() between two fb_greta() fits compares draws      #
# ---------------------------------------------------------------- #
# Sanity check that fb_greta() composes with triangulate(). Uses
# two fits of the same toy model with different seeds; the
# triangulation should report close agreement on common parameters.

test_that("triangulate() runs across two fb_greta() fits", {
  skip_if_no_greta()
  toy <- mk_greta_toy()
  fit_a <- suppressMessages(fb_greta(
    toy$model,
    data = toy$data,
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_b <- suppressMessages(fb_greta(
    toy$model,
    data = toy$data,
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  tri <- triangulate(fit_a, fit_b)
  expect_true(!is.null(tri))
  # Triangulate output is a list; some implementations carry a class.
  expect_true(is.list(tri))
})
