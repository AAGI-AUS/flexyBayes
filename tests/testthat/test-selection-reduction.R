# =============================================================================
# The random-term reduction: builder B.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 4.3.
#
# Builder A needs V block-diagonal. The deposit's design is crossed -- site
# with year, variety across both -- so it is not, and a dense V there would
# be about 43 GB. This is the route for that case: where a candidate is
# constant within a random term's levels, the problem reduces to weighted
# regression on that term's estimated effects, at 3,427 rows instead of
# 73,341.
#
# The admissibility rule carries the scientific weight. A static soil layer
# is constant within every year, so it does not VARY across the year term's
# levels and cannot enter that regression at all. The concern that a static
# layer might be charged for the year component is not something to check
# after the fact -- there is no column for it to occupy.
#
# These tests are deterministic and use mocks where a real fit is not needed.
# The agreement between the two builders (oracle 10) is a separate matter and
# is measured before it is asserted.
# =============================================================================


# --- admissibility, the load-bearing rule ----------------------------- #

test_that("a covariate constant within a term's levels and varying across is admissible", {
  d <- data.frame(
    site = rep(c("a", "b", "c"), each = 4L),
    yr = rep(c("1", "2"), times = 6L),
    clay = rep(c(10, 20, 30), each = 4L),
    stringsAsFactors = FALSE
  )
  adm <- .fb_sel_reduction_admissible(d, "clay", d$site)
  expect_true(adm$constant_within)
  expect_true(adm$varies_across)
  expect_true(adm$admissible)
})


test_that("a site-constant layer is NOT identified for the year term", {
  # The whole point. `clay` varies by site and is constant within each site.
  # Across YEARS it takes every value, so it is not constant within a year;
  # and even where it were, it would not distinguish years. Either way it
  # cannot enter the year regression, which is what makes the site / year
  # asymmetry structural rather than a thing to audit afterwards.
  d <- data.frame(
    site = rep(c("a", "b", "c"), each = 4L),
    yr = rep(c("1", "2"), times = 6L),
    clay = rep(c(10, 20, 30), each = 4L),
    stringsAsFactors = FALSE
  )
  adm <- .fb_sel_reduction_admissible(d, "clay", d$yr)
  expect_false(adm$constant_within)
  expect_false(adm$admissible)
})


test_that("a year-varying covariate is not admissible for the site term", {
  d <- data.frame(
    site = rep(c("a", "b", "c"), each = 4L),
    yr = rep(c("1", "2"), times = 6L),
    rain = rep(c(100, 200), times = 6L),
    stringsAsFactors = FALSE
  )
  expect_false(.fb_sel_reduction_admissible(d, "rain", d$site)$admissible)
  expect_true(.fb_sel_reduction_admissible(d, "rain", d$yr)$admissible)
})


test_that("a constant column is refused even though it is constant within", {
  # Constant within every level, trivially. It says nothing about the term,
  # and admitting it would put a zero-variance column into the design.
  d <- data.frame(
    site = rep(c("a", "b", "c"), each = 4L),
    flat = 7,
    stringsAsFactors = FALSE
  )
  adm <- .fb_sel_reduction_admissible(d, "flat", d$site)
  expect_true(adm$constant_within)
  expect_false(adm$varies_across)
  expect_false(adm$admissible)
})


test_that("admissibility is decided per column, not for the set", {
  d <- data.frame(
    site = rep(c("a", "b", "c"), each = 4L),
    yr = rep(c("1", "2"), times = 6L),
    clay = rep(c(10, 20, 30), each = 4L),
    rain = rep(c(100, 200), times = 6L),
    stringsAsFactors = FALSE
  )
  adm <- .fb_sel_reduction_admissible(d, c("clay", "rain"), d$site)
  expect_identical(adm$admissible, c(TRUE, FALSE))
  expect_identical(adm$column, c("clay", "rain"))
})


# --- the reduction's mechanics ----------------------------------------- #

fb_red_fixture <- function(m_per_trial = 4L) {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(20260912L)
  J <- 12L
  d <- data.frame(
    trial = factor(rep(paste0("T", sprintf("%02d", seq_len(J))),
                       each = m_per_trial)),
    stringsAsFactors = FALSE
  )
  x <- stats::setNames(stats::rnorm(J), levels(d$trial))
  d$x <- as.numeric(x[as.character(d$trial)])
  u <- stats::setNames(stats::rnorm(J, sd = 0.7), levels(d$trial))
  # A within-trial fixed effect. An intercept-only fixed model routes to the
  # aggregated path, which builds an object with no parsed model and no data
  # and is a different shape entirely.
  d$z <- stats::rnorm(nrow(d))
  d$y <- 2 + 1.5 * d$x + 0.3 * d$z +
    as.numeric(u[as.character(d$trial)]) + stats::rnorm(nrow(d), sd = 0.5)
  d
}


test_that("the reduction returns one row per level, matched by name", {
  d <- fb_red_fixture()
  fit <- flexybayes(y ~ z, random = ~ trial, data = d, backend = "inla",
                    verbose = FALSE)
  red <- .fb_sel_reduction(fit, "trial", "x")

  expect_identical(length(red$u), length(red$levels))
  expect_identical(nrow(red$T), length(red$levels))
  expect_identical(rownames(red$T), red$levels)

  # Matched by NAME. The covariate value on each row must be the one that
  # belongs to that level, which a positional match would silently break.
  expected <- vapply(
    red$levels,
    function(l) d$x[as.character(d$trial) == l][[1L]],
    numeric(1L)
  )
  expect_equal(unname(red$T[, "x"]), unname(expected), tolerance = 1e-12)
})


test_that("a balanced design carries no approximation; an unbalanced one does", {
  d <- fb_red_fixture()
  fit <- flexybayes(y ~ z, random = ~ trial, data = d, backend = "inla",
                    verbose = FALSE)
  red <- .fb_sel_reduction(fit, "trial", "x")
  expect_true(red$balanced)
  expect_length(red$approximations, 0L)

  # Drop two rows from one trial. Omega's off-diagonal is then not zero and
  # the diagonal weights become an approximation, which must be declared
  # rather than adopted quietly.
  d2 <- d[-c(1L, 2L), ]
  fit2 <- flexybayes(y ~ z, random = ~ trial, data = d2, backend = "inla",
                     verbose = FALSE)
  red2 <- .fb_sel_reduction(fit2, "trial", "x")
  expect_false(red2$balanced)
  expect_identical(red2$approximations, "omega_diagonal_approximation")
})


test_that("the weights are the posterior spread, and they whiten to identity", {
  d <- fb_red_fixture()
  fit <- flexybayes(y ~ z, random = ~ trial, data = d, backend = "inla",
                    verbose = FALSE)
  red <- .fb_sel_reduction(fit, "trial", "x")
  w <- .fb_sel_reduction_whitener(red)

  expect_identical(w$builder, "random_term_reduction")
  expect_identical(w$n, length(red$levels))

  Omega <- diag(red$omega_diag, nrow = length(red$omega_diag))
  WO <- .fb_sel_whiten(w, Omega)
  WOWt <- t(.fb_sel_whiten(w, t(WO)))
  expect_equal(WOWt, diag(1, length(red$levels)), tolerance = 1e-9)
})


# --- refusals ----------------------------------------------------------- #

test_that("an inadmissible candidate refuses, naming the column and why", {
  d <- fb_red_fixture()
  d$within <- stats::rnorm(nrow(d))  # varies WITHIN each trial
  fit <- flexybayes(y ~ z, random = ~ trial, data = d, backend = "inla",
                    verbose = FALSE)
  expect_error(
    .fb_sel_reduction(fit, "trial", "within"),
    "not constant within levels"
  )
})


test_that("a term the fit does not carry refuses by name", {
  d <- fb_red_fixture()
  fit <- flexybayes(y ~ z, random = ~ trial, data = d, backend = "inla",
                    verbose = FALSE)
  expect_error(.fb_sel_reduction(fit, "nosuchterm", "x"), "no random term")
})


test_that("every refusal this builder can throw is registered", {
  codes <- c(
    "selection_reduction_term_not_in_fit",
    "selection_reduction_level_mismatch",
    "selection_candidates_not_constant_within_term"
  )
  expect_true(all(codes %in% fb_refusals()$reason_code))
})


# --- the stage-1 estimator, and why the default is what it is ---------- #

test_that("absorbed stage-1 effects agree with the one-stage route; shrunk ones do not", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  # The reduction cites the two-stage proposition, which says the term is
  # "absorbed into X0" -- absorbed means FIXED. Reading ranef() instead reads
  # posterior means of a RANDOM term, shrunk by
  # lambda = sigma_u^2 / (sigma_u^2 + sigma_e^2 / m), and a regression of
  # shrunk effects on a covariate returns lambda times the coefficient.
  #
  # Measured over 10 replicates at J = 30, m = 4 on 2026-09-12: the absorbed
  # route sits 0.19 per cent from the one-stage coefficient and the shrunk
  # route 7.42 per cent low. The bounds below carry headroom on both and are
  # a regression gate, not a restatement of those two numbers.
  J <- 30L
  m <- 4L
  gaps <- t(vapply(seq_len(4L), function(s) {
    set.seed(s)
    d <- data.frame(
      trial = factor(rep(sprintf("T%02d", seq_len(J)), each = m))
    )
    x <- stats::setNames(stats::rnorm(J), levels(d$trial))
    d$x <- as.numeric(x[as.character(d$trial)])
    d$z <- stats::rnorm(nrow(d))
    u <- stats::setNames(stats::rnorm(J, sd = 0.7), levels(d$trial))
    d$y <- 2 + 1.5 * d$x + 0.3 * d$z +
      as.numeric(u[as.character(d$trial)]) + stats::rnorm(nrow(d), sd = 0.9)

    fit_a <- flexybayes(y ~ x + z, random = ~ trial, data = d,
                        backend = "inla", verbose = FALSE)
    fx <- summary(fit_a)$fixed
    g_a <- fx$estimate[fx$term == "x"]

    base <- flexybayes(y ~ z, random = ~ trial, data = d,
                       backend = "inla", verbose = FALSE)
    coef_of <- function(mode) {
      r <- .fb_sel_reduction(base, "trial", "x", absorb = mode)
      stats::coef(stats::lm(
        r$u ~ r$T[, "x"], weights = 1 / r$omega_diag
      ))[[2L]]
    }
    c(fixed = abs(coef_of("fixed") / g_a - 1),
      random = abs(coef_of("random") / g_a - 1))
  }, numeric(2L)))

  expect_lt(mean(gaps[, "fixed"]), 0.02)
  # The shrunk route must remain visibly worse. If this ever stops being
  # true the attenuation has gone away and the default should be revisited,
  # which is a change worth reading rather than one to pass silently.
  expect_gt(mean(gaps[, "random"]), 3 * mean(gaps[, "fixed"]))
})


test_that("the shrunk route declares itself on the object", {
  d <- fb_red_fixture()
  fit <- flexybayes(y ~ z, random = ~ trial, data = d, backend = "inla",
                    verbose = FALSE)
  expect_identical(.fb_sel_reduction(fit, "trial", "x")$absorb, "fixed")
  shrunk <- .fb_sel_reduction(fit, "trial", "x", absorb = "random")
  # A coefficient from the shrunk route is on a different scale from builder
  # A's. It must not be pairable with one without the object saying so.
  expect_true("shrunk_stage_one_effects" %in% shrunk$approximations)
})
