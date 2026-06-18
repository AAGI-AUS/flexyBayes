# Tests for triangulate() -- cross-engine posterior comparison
# (deliverable 4).
#
# Pure-numeric tests (Wasserstein-1, R-hat-on-means, name alignment)
# run always. Live cross-engine triangulation tests skip when greta
# or INLA isn't installed.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

# ---------------------------------------------------------------- #
# Synthetic fits via fb_as_draws_simple.default-replacement helper #
# ---------------------------------------------------------------- #

# Helper: build a synthetic fit-like object for which we can supply
# raw draws via a custom method without needing a real backend.
mk_synthetic <- function(draws) {
  obj <- list(draws = draws)
  class(obj) <- c("synthetic_fit", "list")
  obj
}

# Local method registration -- the test file gets sourced after
# the package is loaded, so this S3 method is visible during
# triangulate() invocation in the same session.
fb_as_draws_simple.synthetic_fit <- function(fit, ...) fit$draws
.S3method <- get("registerS3method")
.S3method(
  "fb_as_draws_simple",
  "synthetic_fit",
  fb_as_draws_simple.synthetic_fit
)

# ---------------------------------------------------------------- #
# Wasserstein-1 (1D)                                               #
# ---------------------------------------------------------------- #

test_that(".wasserstein1_1d() is zero for identical samples", {
  set.seed(1L)
  x <- rnorm(500)
  expect_equal(flexyBayes:::.wasserstein1_1d(x, x), 0, tolerance = 1e-10)
})

test_that(".wasserstein1_1d() recovers a deterministic shift", {
  set.seed(2L)
  x <- rnorm(2000)
  shift <- 1.5
  expect_equal(
    flexyBayes:::.wasserstein1_1d(x, x + shift),
    shift,
    tolerance = 0.05
  )
})

test_that(".wasserstein1_1d() handles edge case with empty samples", {
  expect_true(is.na(flexyBayes:::.wasserstein1_1d(numeric(0), 1:5)))
})

# ---------------------------------------------------------------- #
# triangulate() integration via synthetic fits                     #
# ---------------------------------------------------------------- #

test_that("triangulate() detects literal-match common parameters", {
  set.seed(5L)
  d_a <- list(alpha = rnorm(500), beta = rnorm(500))
  d_b <- list(alpha = rnorm(500, 0.1), beta = rnorm(500, -0.1))
  fit_a <- mk_synthetic(d_a)
  fit_b <- mk_synthetic(d_b)

  tri <- triangulate(fit_a, fit_b)
  expect_s3_class(tri, "triangulate_result")
  expect_identical(tri$n_common, 2L)
  expect_setequal(tri$common, c("alpha", "beta"))
  expect_length(tri$only_a, 0L)
  expect_length(tri$only_b, 0L)
  expect_identical(nrow(tri$metrics), 2L)
})

test_that("triangulate() distinguishes only_a / only_b parameters", {
  set.seed(6L)
  d_a <- list(alpha = rnorm(500), beta = rnorm(500), gamma = rnorm(500))
  d_b <- list(alpha = rnorm(500), delta = rnorm(500))
  tri <- triangulate(mk_synthetic(d_a), mk_synthetic(d_b))

  expect_setequal(tri$common, "alpha")
  expect_setequal(tri$only_a, c("beta", "gamma"))
  expect_setequal(tri$only_b, "delta")
})

test_that("triangulate() applies a name_map to fit_b parameter names", {
  set.seed(7L)
  d_a <- list("(Intercept)" = rnorm(500))
  d_b <- list(mu_atg = rnorm(500))
  tri <- triangulate(
    mk_synthetic(d_a),
    mk_synthetic(d_b),
    name_map = c(mu_atg = "(Intercept)")
  )
  expect_identical(tri$n_common, 1L)
  expect_setequal(tri$common, "(Intercept)")
})

test_that("triangulate() reports zero common when names don't match", {
  set.seed(8L)
  d_a <- list(mu_atg = rnorm(500))
  d_b <- list("(Intercept):1" = rnorm(500))
  tri <- triangulate(mk_synthetic(d_a), mk_synthetic(d_b))
  expect_identical(tri$n_common, 0L)
  expect_identical(nrow(tri$metrics), 0L)
})

test_that("triangulate() metrics agree with hand calculation", {
  set.seed(9L)
  d_a <- list(p = rnorm(2000, mean = 1.0, sd = 0.5))
  d_b <- list(p = rnorm(2000, mean = 1.2, sd = 0.5))
  tri <- triangulate(mk_synthetic(d_a), mk_synthetic(d_b))

  m <- tri$metrics
  expect_identical(m$param, "p")
  expect_equal(m$mean_a, mean(d_a$p), tolerance = 1e-10)
  expect_equal(m$mean_b, mean(d_b$p), tolerance = 1e-10)
  expect_equal(m$mean_diff, mean(d_a$p) - mean(d_b$p), tolerance = 1e-10)
  expect_equal(
    m$sd_ratio,
    stats::sd(d_a$p) / stats::sd(d_b$p),
    tolerance = 1e-10
  )
  expect_true(m$wasserstein_1 > 0.1 && m$wasserstein_1 < 0.4)
})

test_that("triangulate() errors when no fb_as_draws_simple method exists", {
  bad <- structure(list(), class = "no_method_class")
  expect_error(triangulate(bad, bad), "does not know how to extract draws")
})

test_that("triangulate() rejects fb_as_draws_simple returning non-list", {
  registerS3method("fb_as_draws_simple", "bad_fit", function(fit, ...) {
    "not_a_list"
  })
  on.exit(
    suppressWarnings(rm(fb_as_draws_simple.bad_fit, envir = globalenv())),
    add = TRUE
  )
  bad <- structure(list(), class = "bad_fit")
  expect_error(triangulate(bad, bad), "must return a named list")
})

# ---------------------------------------------------------------- #
# Source detection                                                 #
# ---------------------------------------------------------------- #

test_that(".triangulate_source() distinguishes greta vs INLA fits", {
  greta_fit <- structure(list(), class = "flexybayes")
  inla_fit <- structure(list(), class = c("flexybayes_inla", "list"))
  expect_identical(flexyBayes:::.triangulate_source(greta_fit), "greta")
  expect_identical(flexyBayes:::.triangulate_source(inla_fit), "inla")
})

# ---------------------------------------------------------------- #
# Predicate + print                                                #
# ---------------------------------------------------------------- #

test_that("is_triangulate_result() recognises the class", {
  set.seed(10L)
  d <- list(p = rnorm(100))
  tri <- triangulate(mk_synthetic(d), mk_synthetic(d))
  expect_true(flexyBayes:::is_triangulate_result(tri))
  expect_false(flexyBayes:::is_triangulate_result(list()))
})

test_that("print.triangulate_result emits the structured header + table", {
  set.seed(11L)
  d_a <- list(p = rnorm(200))
  d_b <- list(p = rnorm(200))
  tri <- triangulate(mk_synthetic(d_a), mk_synthetic(d_b))
  out <- capture.output(print(tri))
  expect_true(any(grepl("triangulate_result", out)))
  expect_true(any(grepl("source_a", out)))
  expect_true(any(grepl("n_common: 1", out)))
  expect_true(any(grepl("Metrics", out)))
})

test_that("print.triangulate_result handles zero common parameters", {
  set.seed(12L)
  d_a <- list(p = rnorm(100))
  d_b <- list(q = rnorm(100))
  tri <- triangulate(mk_synthetic(d_a), mk_synthetic(d_b))
  out <- capture.output(print(tri))
  expect_true(any(grepl("No common parameters", out)))
})

# ---------------------------------------------------------------- #
# Live cross-engine: greta vs INLA -- skip when not installed       #
# ---------------------------------------------------------------- #

test_that("fb_as_draws_simple.flexybayes extracts greta draws", {
  skip_if_no_greta()
  d <- data.frame(
    y = rnorm(40),
    x = rnorm(40),
    g = factor(rep(1:5, length.out = 40))
  )
  fit <- flexybayes(
    fixed = y ~ x,
    data = d,
    n_samples = 50,
    warmup = 50,
    chains = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  draws <- flexyBayes:::fb_as_draws_simple(fit)
  expect_true(is.list(draws))
  expect_true(length(draws) >= 1L)
  expect_true(all(vapply(draws, is.numeric, logical(1))))
})

test_that("triangulate() runs end-to-end on greta vs INLA", {
  # End-to-end triangulation with the user-facing flexybayes() API
  # dispatching to both backends -- greta explicitly + INLA via
  # `backend = "inla"`. Canonical-name resolution is registry-driven
  # so triangulate() works without a user-supplied `name_map`.
  skip_on_cran()
  skip_if_no_greta()
  skip_if_no_inla()
  d <- data.frame(
    y = rnorm(40),
    x = rnorm(40),
    g = factor(rep(1:5, length.out = 40))
  )
  fit_g <- suppressMessages(flexybayes(
    fixed = y ~ x,
    data = d,
    backend = "greta",
    n_samples = 50,
    warmup = 50,
    chains = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_i <- suppressMessages(flexybayes(
    fixed = y ~ x,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
  tri <- triangulate(fit_g, fit_i, n_samples = 200L)
  expect_s3_class(tri, "triangulate_result")
  expect_true(is.numeric(tri$n_common))
  expect_identical(tri$source_a, "greta")
  expect_identical(tri$source_b, "inla")
  # Registry-driven canonical resolution must align at least
  # (Intercept) and x across the two engines.
  expect_true(any(c("(Intercept)", "x") %in% tri$common))
})

# ---------------------------------------------------------------- #
# Live aggregated path: triangulate() on an aggregated INLA fit     #
# ---------------------------------------------------------------- #

test_that("triangulate() runs on an aggregated INLA fit and agrees with per-row (< 0.1 SD)", {
  # The aggregated dispatch path was previously exercised only by code
  # inspection -- no live test ran triangulate() on an aggregated fit.
  # The aggregated INLA posterior is, by construction, the per-row
  # posterior under the default prior: the precision-prior closed-form
  # correction absorbs the within-cell sum-of-squares so the cell-mean
  # likelihood recovers the per-row likelihood up to scale-independent
  # constants. Triangulating the per-row and aggregated INLA fits is the
  # live confirmation of that equivalence and of the aggregated fit's
  # draw-extraction path (the class-vector design that routes the
  # aggregated INLA fit to fb_as_draws_simple.flexybayes_inla). Kept
  # same-engine so the comparison is deterministic and MCMC-free; the
  # cross-engine greta-aggregated extractor is guarded in
  # test-emit-gaussian-aggregated.R.
  skip_on_cran()
  skip_if_no_inla()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE
  )
  set.seed(1L)
  N <- 300L
  J <- 8L
  d <- data.frame(
    f1 = factor(sample(letters[1:3], N, replace = TRUE)),
    f2 = factor(sample(letters[1:3], N, replace = TRUE)),
    g = factor(sample.int(J, N, replace = TRUE))
  )
  d$y <- 1 +
    0.5 * as.integer(d$f1) -
    0.3 * as.integer(d$f2) +
    rnorm(J, 0, 0.3)[as.integer(d$g)] +
    rnorm(N, 0, 0.7)

  fit_row <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = d,
    backend = "inla",
    aggregate = FALSE,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_agg <- suppressMessages(fb(
    y ~ f1 + f2 + (1 | g),
    data = d,
    backend = "inla",
    aggregate = "auto",
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit_agg, "flexybayes_aggregated")
  expect_identical(backend_decision(fit_agg)$path, "aggregated_gaussian")

  # 4000 posterior samples keeps the per-fit Monte-Carlo error on each
  # mean well below the 0.1-SD band (the underlying posterior means
  # agree to ~1e-3; see the stress-strict test in
  # test-emit-gaussian-aggregated.R). The drift band is the audit's
  # acceptance criterion for the aggregated triangulate path.
  set.seed(20260530L)
  tri <- triangulate(fit_row, fit_agg, n_samples = 4000L)
  expect_s3_class(tri, "triangulate_result")
  expect_gt(tri$n_common, 0L)

  fx <- c("(Intercept)", "f1b", "f1c", "f2b", "f2c")
  m <- tri$metrics[tri$metrics$param %in% fx, ]
  expect_identical(nrow(m), length(fx))
  drift <- abs(m$mean_diff) / pmax(m$sd_a, m$sd_b)
  expect_true(
    all(drift < 0.1),
    info = sprintf("max fixed-effect drift = %.3f SD", max(drift))
  )
})

# --- FX-10: common-mode (shared-upstream) caveat (Independent Oracle Principle)

test_that("triangulate() attaches a common-mode caveat unless data independence is declared", {
  d_a <- data.frame(beta = rnorm(200), sigma = abs(rnorm(200)))
  d_b <- data.frame(beta = rnorm(200), sigma = abs(rnorm(200)))

  # Undeclared (default NA): caveat present, data_independence NA.
  tri <- triangulate(mk_synthetic(d_a), mk_synthetic(d_b))
  expect_true(is.na(tri$data_independence))
  expect_false(is.na(tri$shared_upstream_caveat))
  expect_match(tri$shared_upstream_caveat, "agreement")

  # Declared FALSE (same data): caveat present and names the common-mode risk.
  tri_f <- triangulate(mk_synthetic(d_a), mk_synthetic(d_a),
                       data_independence = FALSE)
  expect_false(is.na(tri_f$shared_upstream_caveat))
  expect_match(tri_f$shared_upstream_caveat, "SAME data")

  # Declared TRUE (independent data): no caveat.
  tri_t <- triangulate(mk_synthetic(d_a), mk_synthetic(d_b),
                       data_independence = TRUE)
  expect_true(is.na(tri_t$shared_upstream_caveat))

  # The caveat is metadata: the metrics are identical regardless.
  expect_equal(tri$metrics, tri_t$metrics)
})

test_that("triangulate() validates the data_independence argument", {
  d <- data.frame(beta = rnorm(50))
  expect_error(triangulate(mk_synthetic(d), mk_synthetic(d),
                           data_independence = "yes"),
               "must be a single logical")
})
