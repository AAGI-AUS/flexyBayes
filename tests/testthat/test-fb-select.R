# =============================================================================
# fb_select(penalty = NULL): the B0 deliverable.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 2, 3.1, 11.
#
# The partition principle is the thing under test. The base model, the
# extended base model and a penalised selection are one computation under
# three partitions of the design, so `penalty = NULL` is not a degenerate
# case kept for tidiness -- it is the lambda = 0 corner, it exercises every
# part of the substrate except the provider adapter, and it is what lets the
# eventual penalty be a small diff against machinery already proved.
#
# What these tests hold it to: the gate refuses what it should before
# anything runs, the increment is computed against a protected-only baseline
# from the SAME folds, and the certificate records what was done including
# where it was approximate.
# =============================================================================


fb_sel_verb_data <- function(seed = 20260912L, J = 24L, m = 6L) {
  set.seed(seed)
  d <- data.frame(
    site = factor(rep(sprintf("S%02d", seq_len(J)), each = m))
  )
  per_site <- function() {
    stats::setNames(stats::rnorm(J), levels(d$site))
  }
  clay <- per_site()
  rain <- per_site()
  noise <- per_site()
  d$clay <- as.numeric(clay[as.character(d$site)])
  d$rain <- as.numeric(rain[as.character(d$site)])
  d$noise <- as.numeric(noise[as.character(d$site)])
  # A within-site numeric fixed effect. With only factor predictors the fit
  # takes the aggregated route, which carries no parsed model and no data.
  d$plot_z <- stats::rnorm(nrow(d))
  u <- stats::setNames(stats::rnorm(J, sd = 0.6), levels(d$site))
  d$yield <- 3 + 0.8 * d$clay + 0.5 * d$rain + 0.2 * d$plot_z +
    as.numeric(u[as.character(d$site)]) + stats::rnorm(nrow(d), sd = 0.7)
  d
}

fb_sel_verb_base <- function(d) {
  skip_on_cran()
  skip_if_not_installed("INLA")
  flexybayes(yield ~ plot_z, random = ~ site, data = d,
             backend = "inla", verbose = FALSE)
}


# --- the happy path ---------------------------------------------------- #

test_that("fb_select returns an increment over a protected-only baseline", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d),
    candidates = ~ clay + rain + noise,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    seed = 20260912L
  )

  expect_s3_class(sel, "fb_select")
  expect_true(all(c("protected_only", "full") %in% sel$scores$model))
  # The baseline is scored on the SAME folds, not on a separate run, which is
  # what makes the difference an increment rather than two numbers.
  expect_identical(nrow(sel$scores), 2L)
  expect_identical(
    sel$increment$improvement,
    sel$increment$protected_only - sel$increment$full
  )
  # Three real signals against one noise column: the candidates must help.
  expect_gt(sel$increment$improvement, 0)
})


test_that("a candidate block of pure noise does not earn a large increment", {
  # The negative control. Without it the test above passes for a function
  # that always reports an improvement.
  d <- fb_sel_verb_data()
  real <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L), seed = 20260912L
  )
  noise_only <- fb_select(
    fb_sel_verb_base(d), candidates = ~ noise,
    resamples = fb_resamples(unit = ~ site, v = 6L), seed = 20260912L
  )
  expect_gt(real$increment$improvement, noise_only$increment$improvement)
})


test_that("the same seed gives the same folds and the same answer", {
  d <- fb_sel_verb_data()
  base <- fb_sel_verb_base(d)
  args <- list(
    base = base, candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L), seed = 4242L
  )
  a <- do.call(fb_select, args)
  b <- do.call(fb_select, args)
  expect_equal(a$increment$improvement, b$increment$improvement,
               tolerance = 1e-12)
})


# --- the gate ----------------------------------------------------------- #

test_that("a penalty refuses, and names what is implemented", {
  d <- fb_sel_verb_data()
  expect_error(
    fb_select(fb_sel_verb_base(d), ~ clay, penalty = list(alpha = 1),
              resamples = fb_resamples(unit = ~ site, v = 6L), seed = 1L),
    "a penalty is not implemented"
  )
})


test_that("a non-flexyBayes object refuses before anything runs", {
  expect_error(fb_select(structure(list(), class = "lm"), ~ x),
               "must be a fitted flexyBayes model")
})


test_that("an aggregated-route fit refuses by name, not by falling over", {
  # An intercept-only fixed model takes the aggregated route and carries no
  # parsed model and no data. It must say so rather than fail somewhere
  # downstream with a message about row counts.
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- fb_sel_verb_data()
  agg <- flexybayes(yield ~ 1, random = ~ site, data = d,
                    backend = "inla", verbose = FALSE)
  expect_error(
    fb_select(agg, ~ clay, resamples = fb_resamples(unit = ~ site, v = 6L),
              seed = 1L),
    "does not carry the parsed model"
  )
})


# --- the certificate ----------------------------------------------------- #

test_that("the certificate records the covariance, the estimand and the unit", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L), seed = 20260912L
  )
  cert <- sel$certificate
  expect_s3_class(cert, "fb_select_certificate")
  expect_identical(cert$covariance$estimand, "plugin_mean_covariance")
  expect_identical(cert$covariance$builder, "block_cholesky")
  expect_identical(cert$resampling$prediction, "new_cluster")
  expect_identical(cert$resampling$seed, 20260912L)
  # Convergence is the tri-state verdict, not a recomputed local rule.
  expect_identical(cert$base$convergence, .fb_convergence_verdict(sel$base)$ok)
  expect_true(is.character(cert$approximations))
})


test_that("glance reports the increment and how it was scored", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay + rain,
    resamples = fb_resamples(unit = ~ site, v = 6L), seed = 20260912L
  )
  g <- glance(sel)
  expect_identical(nrow(g), 1L)
  expect_true(all(c(
    "metric", "protected_only", "full", "improvement", "n_units",
    "se_unit", "prediction", "covariance", "builder"
  ) %in% names(g)))
  expect_identical(g$n_units, 24L)
})


test_that("print states the estimand rather than only the number", {
  d <- fb_sel_verb_data()
  sel <- fb_select(
    fb_sel_verb_base(d), candidates = ~ clay,
    resamples = fb_resamples(unit = ~ site, v = 6L), seed = 20260912L
  )
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "plugin_mean_covariance", fixed = TRUE)
  expect_match(out, "new_cluster", fixed = TRUE)
  expect_match(out, "independent", fixed = TRUE)
})
