# =============================================================================
# The autoregressive field that vanished while the fit said it converged.
#
# On a holed grid the latent field can lose its signal entirely: the field
# standard deviation runs to a floor around 0.01, both correlations sit at
# approximately zero with credible intervals spanning almost the whole of
# [-1, 1], and the variance the field should have carried is absorbed into
# the nugget. Measured over ten hole patterns on one 12 x 10 grid at 10 per
# cent missing responses, five lost the field -- while the convergence block
# read "Mode status: 0 (0 = converged)" and "Numerical confirm: PASS" on
# every one of them, and the augmentation record reported that it had
# completed the design.
#
# The intervals were the only honest tell and they are easy to read past.
# The fit now says it.
#
# The warning is scoped to the FIELD. A grouping factor whose variance
# component collapses is a different fact about a different structure --
# frequently the thing a hierarchical model exists to show -- and it keeps
# the quieter `collapsed` cell in the variance-component table with no
# warning attached. Pinned below.
# =============================================================================

.sfc_ar1_mat <- function(n, rho) {
  rho^abs(outer(seq_len(n), seq_len(n), "-"))
}

# A field of the requested strength on a complete grid. `sd_s` near zero
# with both correlations at zero is a design with no spatial signal at
# all, which is what a collapsed fit looks like from the inside.
.sfc_sim <- function(n_row, n_col, rho_r, rho_c, sd_s, sd_e, seed) {
  set.seed(seed)
  sigma <- sd_s^2 * kronecker(
    .sfc_ar1_mat(n_row, rho_r),
    .sfc_ar1_mat(n_col, rho_c)
  )
  z <- as.numeric(
    t(chol(sigma + diag(1e-8, nrow(sigma)))) %*% stats::rnorm(nrow(sigma))
  )
  d <- expand.grid(
    col = factor(seq_len(n_col)),
    row = factor(seq_len(n_row)),
    KEEP.OUT.ATTRS = FALSE
  )
  d$y <- 5 + z + stats::rnorm(nrow(d), 0, sd_e)
  d[, c("y", "row", "col")]
}

.sfc_field_fit <- function(d) {
  suppressMessages(flexybayes(
    y ~ 1,
    random = ~ ar1(row):ar1(col),
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
}

# ---------------------------------------------------------------- #
# 1. A collapsed field warns, and says what to do about it          #
# ---------------------------------------------------------------- #

test_that("a field with no spatial signal warns at fit time", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = FALSE
  )
  d <- .sfc_sim(10L, 8L, 0, 0, sd_s = 1e-3, sd_e = 1, seed = 20260819L)

  w <- tryCatch(.sfc_field_fit(d), warning = function(w) w)
  expect_s3_class(w, "warning")

  msg <- conditionMessage(w)
  # What happened, in the fit's own words.
  expect_match(msg, "did not identify", fixed = TRUE)
  expect_match(msg, "carries no spatial signal", fixed = TRUE)
  expect_match(msg, "independent-errors model", fixed = TRUE)
  # Why it converged anyway.
  expect_match(msg, "converged mode regardless", fixed = TRUE)
  # The measured trigger.
  expect_match(msg, "incomplete grid raises the risk", fixed = TRUE)
  # All three remedies.
  expect_match(msg, "fb_complete_grid()", fixed = TRUE)
  expect_match(msg, "brms engine", fixed = TRUE)
  expect_match(msg, "fb_prior()", fixed = TRUE)
  expect_match(msg, "silence_spatial_collapse_warning", fixed = TRUE)
})

test_that("the collapsed field is flagged in the varcomp table too", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = TRUE
  )
  d <- .sfc_sim(10L, 8L, 0, 0, sd_s = 1e-3, sd_e = 1, seed = 20260819L)
  fit <- .sfc_field_fit(d)

  vc <- summary(fit)$varcomp
  field <- vc[vc$component == "sd_spatial", , drop = FALSE]
  expect_identical(nrow(field), 1L)
  expect_identical(field$note, "collapsed")

  # The warning and the cell are two readings of one table, so a fit that
  # warns must also flag, and the reasons must name every parameter that
  # lost its identification: the field SD and both correlations.
  reasons <- flexyBayes:::.fb_spatial_collapse_reasons(fit)
  expect_gte(length(reasons), 2L)
  expect_true(any(grepl("field standard deviation", reasons, fixed = TRUE)))
  expect_true(any(grepl("rho_row is unidentified", reasons, fixed = TRUE)))
})

# ---------------------------------------------------------------- #
# 2. It stays quiet where it should                                 #
# ---------------------------------------------------------------- #

test_that("a field that carries real signal draws no warning", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = FALSE
  )
  d <- .sfc_sim(10L, 8L, 0.7, 0.5, sd_s = 1.5, sd_e = 0.6, seed = 20260818L)

  expect_no_warning(fit <- .sfc_field_fit(d))
  expect_identical(flexyBayes:::.fb_spatial_collapse_reasons(fit), character(0))

  vc <- summary(fit)$varcomp
  expect_identical(vc$note[vc$component == "sd_spatial"], "")
})

test_that("a non-spatial fit is outside the warning's scope entirely", {
  # A small grouping variance is a different fact about a different
  # structure. It keeps the quiet `collapsed` cell and raises nothing --
  # the teaching in the getting-started tutorial depends on that.
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = FALSE
  )
  set.seed(7L)
  d <- data.frame(g = factor(rep(seq_len(8L), each = 12L)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 3 + 0.5 * d$x + stats::rnorm(nrow(d), 0, 1)

  expect_no_warning(
    fit <- suppressMessages(flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      backend = "inla",
      verbose = FALSE
    ))
  )
  expect_identical(flexyBayes:::.fb_spatial_collapse_reasons(fit), character(0))
})

# ---------------------------------------------------------------- #
# 3. The rule itself, without paying for a fit                      #
# ---------------------------------------------------------------- #

.sfc_shell <- function(vc, medians, field = TRUE) {
  attr(vc, "posterior_median") <- medians
  random_terms <- if (field) {
    list(list(type = "ar1_spatial", row_var = "row", col_var = "col"))
  } else {
    list(list(type = "simple", var = "g"))
  }
  structure(
    list(
      fb = list(random_terms = random_terms),
      extras = list(variance_comps = vc)
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

.sfc_vc <- function(field_upper, rho_lower, rho_upper) {
  data.frame(
    component = c("sigma", "sd_spatial", "rho_row"),
    estimate = c(1, 0.01, 0),
    sd = c(0.1, 0.006, 0.7),
    q2.5 = c(0.8, 0.003, rho_lower),
    q97.5 = c(1.2, field_upper, rho_upper),
    stringsAsFactors = FALSE
  )
}

test_that("the field SD rule is a comparison against the nugget scale", {
  medians <- c(sigma = 1)

  # Upper bound at 2.7% of the nugget: the measured collapse, which the
  # 1% display cut used for ordinary components does not reach.
  expect_length(
    flexyBayes:::.fb_spatial_collapse_reasons(
      .sfc_shell(.sfc_vc(0.027, -0.1, 0.1), medians)
    ),
    1L
  )
  # A field an order of magnitude larger is not at any boundary.
  expect_identical(
    flexyBayes:::.fb_spatial_collapse_reasons(
      .sfc_shell(.sfc_vc(0.9, -0.1, 0.1), medians)
    ),
    character(0)
  )
})

test_that("a correlation is called unidentified only when it spans the space", {
  medians <- c(sigma = 1)

  # [-0.99, 0.99] is the whole parameter space: the data constrained
  # nothing.
  reasons <- flexyBayes:::.fb_spatial_collapse_reasons(
    .sfc_shell(.sfc_vc(0.9, -0.99, 0.99), medians)
  )
  expect_length(reasons, 1L)
  expect_match(reasons, "rho_row is unidentified", fixed = TRUE)

  # A correlation that is merely small, with a bounded interval, is a
  # finding rather than a failure.
  expect_identical(
    flexyBayes:::.fb_spatial_collapse_reasons(
      .sfc_shell(.sfc_vc(0.9, -0.2, 0.25), medians)
    ),
    character(0)
  )
})

test_that("a fit carrying no field is never examined", {
  # Same numbers, no field term: the rule must not read a grouping
  # factor's variance component as a lost spatial field.
  vc <- data.frame(
    component = c("sigma", "sd_g"),
    estimate = c(1, 0.001),
    sd = c(0.1, 0.0006),
    q2.5 = c(0.8, 0.0003),
    q97.5 = c(1.2, 0.002),
    stringsAsFactors = FALSE
  )
  expect_identical(
    flexyBayes:::.fb_spatial_collapse_reasons(
      .sfc_shell(vc, c(sigma = 1), field = FALSE)
    ),
    character(0)
  )
})
