# test-aggregate-plan-joint-key.R -- the cell key is a set of variables,
# not a list of terms.
#
# `y ~ a * b` produces three terms over two variables. The planner
# multiplied the term-level counts, so a replicated 4-by-4 factorial that
# the runtime aggregator (.fb_stream_key_cols(), which has always keyed
# on the union) compresses 320 rows into 16 cells was estimated at
# 4 * 4 * 16 = 256 and refused as `compression_unproductive`. The
# 2026-08-16 adversarial review found it as P1-3.
#
# The second half of that finding is what a live probe then showed: once
# the planner admits the model, the aggregated INLA emit cannot fit it.
# It names model-matrix columns in the INLA formula, and an interaction
# column is called `a2:b2`, which INLA does not resolve even backticked
# (`object 'aa2:bb2' not found`, on both the binomial and the Poisson
# emit). The plan refuses that shape by name rather than handing INLA a
# formula it cannot read.

.joint_key_factorial <- function(n_rep = 20L, levels = c(4L, 4L),
                                 seed = 20260816L) {
  set.seed(seed)
  args <- lapply(
    seq_along(levels),
    function(i) factor(paste0(letters[i], seq_len(levels[[i]])))
  )
  names(args) <- letters[seq_along(levels)]
  args$rep <- seq_len(n_rep)
  d <- do.call(expand.grid, args)
  d$rep <- NULL
  d$y <- stats::rbinom(nrow(d), 1L, 0.4)
  d
}


# ---------------------------------------------------------------- #
# (1) The joint key                                                  #
# ---------------------------------------------------------------- #

test_that("y ~ a*b keys on two variables, K = 16 not 256", {
  d <- .joint_key_factorial()
  fb <- fb_from_brms(y ~ a * b, data = d)
  plan <- flexyBayes:::.fb_aggregation_plan(fb, flexyBayes:::.fb_dataset(d))

  expect_identical(plan$K_est, 16L)
  expect_identical(plan$N, 320L)
  # One entry per distinct key variable, not one per term.
  expect_length(plan$cell_key_terms, 2L)
  expect_setequal(
    vapply(plan$cell_key_terms, function(x) x$label, character(1L)),
    c("a", "b")
  )
})

test_that("y ~ a*b*c keys on three variables, K = 64", {
  d <- .joint_key_factorial(levels = c(4L, 4L, 4L))
  fb <- fb_from_brms(y ~ a * b * c, data = d)
  plan <- flexyBayes:::.fb_aggregation_plan(fb, flexyBayes:::.fb_dataset(d))

  expect_identical(plan$K_est, 64L)
  expect_identical(plan$N, 1280L)
  expect_length(plan$cell_key_terms, 3L)
})

test_that("a variable shared by a fixed and a random term is counted once", {
  set.seed(20260816L)
  n <- 600L
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(4L)), n, replace = TRUE))
  )
  d$y <- stats::rnorm(n)
  # `f` names a fixed effect and a grouping factor. The joint key is
  # {f, g}, so K = 20, not 5 * 5 * 4.
  fb <- fb_from_brms(y ~ f + (1 | f) + (1 | g), data = d)
  plan <- flexyBayes:::.fb_aggregation_plan(fb, flexyBayes:::.fb_dataset(d))
  expect_identical(plan$K_est, 20L)
  expect_length(plan$cell_key_terms, 2L)
})

test_that("the simple non-interaction case is unchanged", {
  set.seed(1L)
  n <- 1000L
  d <- data.frame(
    y = stats::rnorm(n),
    f = factor(sample(letters[1:5], n, replace = TRUE)),
    g = factor(sample.int(20L, n, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = d)
  plan <- flexyBayes:::.fb_aggregation_plan(fb, flexyBayes:::.fb_dataset(d))
  expect_true(plan$eligible)
  expect_identical(plan$K_est, 100L)
})


# ---------------------------------------------------------------- #
# (2) The engine limit the corrected key exposes                     #
# ---------------------------------------------------------------- #

test_that("an interaction design column refuses the aggregated route", {
  d <- .joint_key_factorial()
  fb <- fb_from_brms(y ~ a * b, data = d)
  plan <- flexyBayes:::.fb_aggregation_plan(fb, flexyBayes:::.fb_dataset(d))

  expect_false(plan$eligible)
  expect_identical(
    plan$reason_codes,
    "aggregated_interaction_column_not_representable"
  )
  # The compression the model would have had is still reported -- it is
  # what makes the refusal legible.
  expect_identical(plan$K_est, 16L)
  expect_equal(plan$compression_est, 0.05)
})

test_that("aggregate = TRUE on an interaction refuses typed and explains", {
  d <- .joint_key_factorial()
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ a * b,
      data = d,
      family = "binomial",
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_identical(
    err$reason_codes,
    "aggregated_interaction_column_not_representable"
  )
  msg <- conditionMessage(err)
  expect_match(msg, "cannot reference a column whose name contains a colon",
               fixed = TRUE)
  expect_match(msg, "aggregate = FALSE", fixed = TRUE)
})

test_that("the per-row route the refusal names does fit the model", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .joint_key_factorial(n_rep = 20L)
  fit <- suppressMessages(flexybayes(
    y ~ a * b,
    data = d,
    family = "binomial",
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_inla")
  # The interaction coefficients are present, on INLA's own naming.
  expect_true(any(grepl(":", rownames(fit$inla$summary.fixed), fixed = TRUE)))
})
