# fb_select(): glance() and print() name the out-of-fold loss and its
# direction. Oracle: the loss and the ELPD computed by hand from the per-unit
# values stored on the result.

.ll_cache <- new.env(parent = emptyenv())

.ll_data <- function(J = 24L, m = 6L) {
  set.seed(20260924L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(J)), each = m)))
  per_site <- function() stats::setNames(stats::rnorm(J), levels(d$site))
  clay <- per_site()
  rain <- per_site()
  d$clay <- as.numeric(clay[as.character(d$site)])
  d$rain <- as.numeric(rain[as.character(d$site)])
  d$plot_z <- stats::rnorm(nrow(d))
  u <- stats::setNames(stats::rnorm(J, sd = 0.6), levels(d$site))
  d$yield <- 3 + 0.8 * d$clay + 0.5 * d$rain + 0.2 * d$plot_z +
    as.numeric(u[as.character(d$site)]) + stats::rnorm(nrow(d), sd = 0.7)
  d
}

.ll_benchmark <- function() {
  if (is.null(.ll_cache$benchmark)) {
    withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
    d <- .ll_data()
    base <- flexybayes(yield ~ plot_z, random = ~ site, data = d,
                       backend = "inla", verbose = FALSE)
    .ll_cache$benchmark <- fb_select(
      base, candidates = ~ clay + rain,
      resamples = fb_resamples(unit = ~ site, v = 6L),
      workflow = "unpenalized_benchmark", seed = 20260924L
    )
  }
  .ll_cache$benchmark
}

# The mean over held-out units of a per-unit column, a unit's repeats
# averaged first.
.ll_unit_mean <- function(clusters, model, column) {
  rows <- clusters[clusters$model == model, , drop = FALSE]
  mean(tapply(rows[[column]], rows$cluster, mean))
}


test_that("glance() names the loss, and its values are that loss", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  sel <- .ll_benchmark()
  g <- glance(sel)

  expect_identical(g$loss, "nlpd")
  expect_false(any(
    c("metric", "protected_only", "full", "improvement") %in% names(g)
  ))
  clusters <- as.data.frame(sel$clusters)
  protected <- .ll_unit_mean(clusters, "protected_only", "nlpd")
  full <- .ll_unit_mean(clusters, "full", "nlpd")
  expect_equal(g$loss_protected_only, protected, tolerance = 1e-12)
  expect_equal(g$loss_full, full, tolerance = 1e-12)
  expect_equal(g$loss_reduction, protected - full, tolerance = 1e-12)
  # The candidates carry real effects, so the loss falls when they enter.
  expect_gt(g$loss_reduction, 0)
})


test_that("an elpd score is labelled as the negated elpd it holds", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  sel <- .ll_benchmark()
  # The Bayesian workflow stores the metric it was asked for, `elpd`, beside
  # scores that are -elpd per unit.
  sel$increment$metric <- "elpd"
  g <- glance(sel)
  expect_identical(g$loss, "neg_elpd")
  expect_false(identical(g$loss, "elpd"))
})


test_that("print() says the scores are a loss and lower is better", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  sel <- .ll_benchmark()
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "out-of-fold loss (nlpd, lower is better)", fixed = TRUE)
  expect_match(out, "reduction", fixed = TRUE)

  sel$increment$metric <- "elpd"
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(
    out, "out-of-fold loss (neg_elpd, lower is better)", fixed = TRUE
  )
  expect_false(grepl("out-of-fold elpd", out, fixed = TRUE))
})


test_that("the Bayesian route's loss is minus the elpd per held-out unit", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_high_missingness_warning = TRUE
  )
  d <- .ll_data(J = 12L, m = 6L)
  base <- flexybayes(yield ~ plot_z, random = ~ site, data = d,
                     backend = "inla", verbose = FALSE)
  sel <- suppressWarnings(fb_select(
    base, candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 3L),
    seed = 20260924L
  ))
  expect_identical(sel$workflow, "bayes")
  g <- glance(sel)
  expect_identical(g$loss, "neg_elpd")

  # Each stored cluster value is the summed log predictive density of that
  # site's held-out rows; the ELPD is their sum over the U sites.
  clusters <- as.data.frame(sel$clusters)
  n_units <- length(unique(clusters$cluster))
  expect_identical(g$n_units, n_units)
  for (model in c("protected_only", "full")) {
    rows <- clusters[clusters$model == model, , drop = FALSE]
    elpd <- sum(tapply(rows$elpd, rows$cluster, mean))
    loss <- g[[paste0("loss_", model)]]
    expect_equal(loss, -elpd / n_units, tolerance = 1e-12)
    # The elpd is negative and the loss positive: the sign glance() once
    # labelled `elpd`.
    expect_lt(elpd, 0)
    expect_gt(loss, 0)
  }
  expect_equal(
    g$loss_reduction, g$loss_protected_only - g$loss_full,
    tolerance = 1e-12
  )
})
