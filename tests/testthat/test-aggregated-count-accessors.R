# =============================================================================
# Accessors on aggregated and streamed count fits.
#
# An aggregated binomial or Poisson fit is one INLA model over cell totals, so
# INLA's fitted values are one per cell. fitted() and residuals() dispatched to
# the per-row INLA method, which read those cell values and then errored on the
# data the aggregated object does not carry ("missing value where TRUE/FALSE
# needed"). The row-level values the object did store were the fixed part of
# the linear predictor only, without the random effects, so returning them
# would have been wrong by up to 3.5 counts on a Poisson fixture.
#
# The contract tested here:
# - fitted() and residuals() return one value per row, in row order, equal to
#   the per-row fit of the same model (one INLA thread, tolerance 1e-6 against
#   a measured difference below 1e-8, and 5e-7 with a within-cell offset);
# - a streamed fit holds only the first rows it read, so where that head
#   sample is shorter than the data it refuses by name rather than return the
#   sample as if it were every row;
# - nobs() and glance()$nobs are the integer summary()$n_design reports;
# - every generic that reaches an aggregated or streamed fit returns a value
#   or refuses with a registered reason code, and the test counts the pairs
#   it checked.
# =============================================================================

.acc_one_thread <- function(env = parent.frame()) {
  threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  withr::defer(INLA::inla.setOption(num.threads = threads), envir = env)
}

.acc_counts <- function(n_rows, seed) {
  set.seed(seed)
  env <- sample(1:4, n_rows, replace = TRUE)
  geno <- sample(1:12, n_rows, replace = TRUE)
  eta <- c(0, 1.2, -0.7, 0.5)[env] + stats::rnorm(12, 0, 0.6)[geno]
  data.frame(
    env = factor(env),
    geno = factor(geno),
    y = stats::rpois(n_rows, exp(eta - 0.5)),
    b = stats::rbinom(n_rows, 1, stats::plogis(eta - 1)),
    m = sample(1:5, n_rows, replace = TRUE)
  )
}

.acc_fit <- function(route, fixed, data, family, ...) {
  prior <- fb_prior(sd(group = "geno") ~ uniform(lower = 0, upper = 10))
  suppressMessages(switch(
    route,
    per_row = flexybayes(
      fixed, random = ~geno, data = data, family = family, prior = prior,
      backend = "inla", aggregate = FALSE, verbose = FALSE, ...
    ),
    aggregated = flexybayes(
      fixed, random = ~geno, data = data, family = family, prior = prior,
      backend = "inla", aggregate = TRUE, verbose = FALSE, ...
    ),
    stream = flexybayes_stream(
      fixed, random = ~geno, source = data, family = family, prior = prior,
      backend = "inla", chunk_rows = 50000, verbose = FALSE, ...
    )
  ))
}

test_that("fitted() and residuals() on aggregated count fits match the per-row fit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  .acc_one_thread()
  d <- .acc_counts(1200L, 314L)

  for (case in list(list(fixed = y ~ env, family = "poisson"),
                    list(fixed = b ~ env, family = "binomial"))) {
    agg <- .acc_fit("aggregated", case$fixed, d, case$family)
    row <- .acc_fit("per_row", case$fixed, d, case$family)
    expect_s3_class(agg, "flexybayes_aggregated")
    y <- d[[all.vars(case$fixed)[1L]]]

    f_agg <- fitted(agg)
    expect_length(f_agg, nrow(d))
    expect_true(all(is.finite(f_agg)))
    expect_equal(unname(f_agg), unname(fitted(row)), tolerance = 1e-6)

    r_agg <- residuals(agg)
    expect_length(r_agg, nrow(d))
    expect_equal(unname(r_agg), unname(y - f_agg), tolerance = 1e-12)
  }
})

test_that("an offset inside each cell reaches the aggregated fitted values", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  .acc_one_thread()
  d <- .acc_counts(1200L, 315L)
  d$E <- stats::runif(nrow(d), 1, 10)
  d$y <- stats::rpois(nrow(d), d$E * exp(0.2 * as.integer(d$env)))

  agg <- .acc_fit("aggregated", y ~ env + offset(log(E)), d, "poisson")
  row <- .acc_fit("per_row", y ~ env + offset(log(E)), d, "poisson")
  expect_equal(unname(fitted(agg)), unname(fitted(row)), tolerance = 1e-6)
  # Residuals sit near zero, where a relative tolerance against the per-row
  # fit measures nothing useful; they are held to their definition instead.
  expect_equal(unname(residuals(agg)), unname(d$y - fitted(agg)), tolerance = 1e-12)
})

test_that("aggregated fitted values follow the rows, not the cells", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  .acc_one_thread()
  d <- .acc_counts(1200L, 316L)
  set.seed(2L)
  perm <- sample.int(nrow(d))

  f_orig <- fitted(.acc_fit("aggregated", y ~ env, d, "poisson"))
  f_perm <- fitted(.acc_fit("aggregated", y ~ env, d[perm, ], "poisson"))
  expect_equal(unname(f_perm), unname(f_orig[perm]), tolerance = 1e-6)
})

test_that("a missing response keeps the aggregated count route refused by name", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .acc_counts(400L, 317L)
  d$y[3L] <- NA
  expect_error(
    .acc_fit("aggregated", y ~ env, d, "poisson"),
    class = "flexybayes_refusal_aggregation_response_incomplete"
  )
})

test_that("a streamed fit returns row values only when it holds every row", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  .acc_one_thread()

  # 1,200 rows: the head sample is every row, so the values exist and equal
  # the in-memory aggregated fit's; an exposure column scales each row.
  d <- .acc_counts(1200L, 318L)
  agg <- .acc_fit("aggregated", y ~ env, d, "poisson")
  small <- .acc_fit("stream", y ~ env, d, "poisson")
  expect_equal(unname(fitted(small)), unname(fitted(agg)), tolerance = 1e-6)
  expect_equal(unname(residuals(small)), unname(residuals(agg)), tolerance = 1e-6)

  exposed <- .acc_fit("stream", y ~ env, d, "poisson", exposure = "m")
  offset_fit <- .acc_fit("aggregated", y ~ env + offset(log(m)), d, "poisson")
  expect_equal(unname(fitted(exposed)), unname(fitted(offset_fit)), tolerance = 1e-6)

  # 110,000 rows: the head sample is the first 100,000 only.
  big <- .acc_fit("stream", y ~ env, .acc_counts(110000L, 319L), "poisson")
  for (accessor in list(fitted, residuals, stats::model.matrix)) {
    expect_error(
      accessor(big),
      class = "flexybayes_refusal_streamed_rows_not_retained"
    )
  }
})

test_that("nobs() and glance() report the integer the summary reports", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  .acc_one_thread()
  d <- .acc_counts(1200L, 320L)

  for (fit in list(.acc_fit("aggregated", y ~ env, d, "poisson"),
                   .acc_fit("stream", y ~ env, d, "poisson"))) {
    expect_identical(stats::nobs(fit), 1200L)
    expect_identical(stats::nobs(fit), summary(fit)$n_design)
    expect_identical(suppressMessages(glance(fit))$nobs, 1200L)
  }
})

test_that("every generic reaching an aggregated or streamed fit returns or refuses by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  .acc_one_thread()

  # The generics come from the package's own S3 registrations, so a method
  # added for any class an aggregated fit inherits joins the enumeration and
  # changes the count asserted at the end.
  ns <- readLines(system.file("NAMESPACE", package = "flexyBayes"))
  reg <- regmatches(ns, regexec("^S3method\\(([^,]+),([^)]+)\\)$", ns))
  reg <- do.call(rbind, reg[lengths(reg) == 3L])
  classes <- c("flexybayes", "flexybayes_inla", "flexybayes_aggregated")
  generics <- sort(unique(reg[reg[, 3L] %in% classes, 2L]))

  calls <- list(
    "anova" = function(fit) anova(fit),
    "augment" = function(fit) augment(fit),
    "bayesplot::pp_check" = function(fit) bayesplot::pp_check(fit),
    "canonical_names" = function(fit) canonical_names(fit),
    "coef" = function(fit) coef(fit),
    "confint" = function(fit) confint(fit),
    # ref_grid() reaches emm_basis() after recover_data(), and wraps a
    # recover_data() error in its own untyped "Perhaps a 'data' or 'params'
    # argument is needed"; calling recover_data() first keeps the refusal's
    # own condition where there is one.
    "emmeans::emm_basis" = function(fit) {
      emmeans::recover_data(fit)
      emmeans::ref_grid(fit)
    },
    "emmeans::recover_data" = function(fit) emmeans::recover_data(fit),
    "family" = function(fit) family(fit),
    "fb_as_draws_simple" = function(fit) fb_as_draws_simple(fit, n_samples = 20L),
    "fb_log_posterior" = function(fit) fb_log_posterior(fit),
    "fitted" = function(fit) fitted(fit),
    "formula" = function(fit) formula(fit),
    "glance" = function(fit) glance(fit),
    "insight::get_data" = function(fit) insight::get_data(fit),
    "logLik" = function(fit) logLik(fit),
    "loo::loo" = function(fit) loo::loo(fit),
    "marginaleffects::get_coef" = function(fit) marginaleffects::get_coef(fit),
    "marginaleffects::get_predict" = function(fit) marginaleffects::get_predict(fit),
    "marginaleffects::get_vcov" = function(fit) marginaleffects::get_vcov(fit),
    "marginaleffects::set_coef" = function(fit) marginaleffects::set_coef(fit, coef(fit)),
    "model.matrix" = function(fit) stats::model.matrix(fit),
    "nlme::ranef" = function(fit) nlme::ranef(fit),
    "nobs" = function(fit) stats::nobs(fit),
    "plot" = function(fit) {
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      plot(fit)
    },
    "posterior::as_draws" = function(fit) posterior::as_draws(fit, n_draws = 20L),
    "posterior::as_draws_df" = function(fit) posterior::as_draws_df(fit, n_draws = 20L),
    "posterior::as_draws_matrix" = function(fit) posterior::as_draws_matrix(fit, n_draws = 20L),
    "predict" = function(fit) predict(fit),
    "print" = function(fit) utils::capture.output(print(fit)),
    "prior_summary" = function(fit) utils::capture.output(print(prior_summary(fit))),
    "ranef" = function(fit) ranef(fit),
    "residuals" = function(fit) residuals(fit),
    "rstantools::prior_summary" = function(fit) {
      utils::capture.output(print(rstantools::prior_summary(fit)))
    },
    "summary" = function(fit) utils::capture.output(print(summary(fit))),
    "tidy" = function(fit) tidy(fit),
    "update" = function(fit) update(fit),
    "vcov" = function(fit) vcov(fit),
    "vpredict" = function(fit) vpredict(fit, ratio ~ sd_geno^2),
    "sommer::vpredict" = function(fit) {
      if (requireNamespace("sommer", quietly = TRUE)) {
        sommer::vpredict(fit, ratio ~ sd_geno^2)
      }
    }
  )
  # A generic with no call written here fails the test by name.
  expect_identical(setdiff(generics, names(calls)), character(0))

  fits <- list(
    aggregated = .acc_fit("aggregated", y ~ env, .acc_counts(1200L, 321L), "poisson"),
    streamed_head_sample = .acc_fit("stream", y ~ env, .acc_counts(110000L, 322L), "poisson")
  )
  for (pkg in c(
    "bayesplot", "emmeans", "insight", "loo", "marginaleffects", "nlme",
    "posterior", "rstantools"
  )) {
    skip_if_not_installed(pkg)
  }

  checked <- 0L
  for (g in generics) {
    for (nm in names(fits)) {
      out <- tryCatch(
        suppressWarnings(suppressMessages(calls[[g]](fits[[nm]]))),
        error = function(e) e
      )
      if (inherits(out, "error")) {
        expect_true(
          inherits(out, "flexybayes_refusal") && !is.null(out$reason_code),
          label = paste0(g, " on the ", nm, " fit: ", conditionMessage(out))
        )
      } else {
        succeed()
      }
      checked <- checked + 1L
    }
  }
  # 39 generics on two fits, and sommer's vpredict() when sommer is installed.
  expect_identical(checked, 78L + 2L * ("sommer::vpredict" %in% generics))
})
