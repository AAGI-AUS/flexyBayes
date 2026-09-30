# =============================================================================
# An AR1 field the data do not support is reported absent.
#
# At correlations near zero the field is the same model as the nugget, and
# under the half-normal default its SD stays well away from zero with no
# signal in the data. The fit is compared with the same model without the
# field; a log Bayes factor below -log(100) reports the field absent, with
# a warning at fit time and `collapsed` in summary(fit)$varcomp.
#
# The warning is scoped to the FIELD. A grouping factor whose variance
# component collapses keeps the quieter `collapsed` cell with no warning.
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
  # All four remedies.
  expect_match(msg, "fit the model without the field", fixed = TRUE)
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

  # The warning and the cell read one comparison, so a fit that warns
  # must also flag, and the reason states the comparison.
  reasons <- flexyBayes:::.fb_spatial_collapse_reasons(fit)
  expect_gte(length(reasons), 1L)
  expect_true(any(grepl("log Bayes factor", reasons, fixed = TRUE)))
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

test_that("a field with no spatial signal is flagged on every seed", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = TRUE
  )
  for (seed in c(20260819L, 11L, 12L)) {
    d <- .sfc_sim(10L, 8L, 0, 0, sd_s = 1e-3, sd_e = 1, seed = seed)
    fit <- .sfc_field_fit(d)
    vc <- summary(fit)$varcomp
    expect_identical(
      vc$note[vc$component == "sd_spatial"],
      "collapsed",
      info = paste("seed", seed)
    )
  }
})

test_that("fields of SD 1 and 1.5 are not flagged on any seed", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = FALSE
  )
  fields <- list(
    moderate = list(rho_r = 0.6, rho_c = 0.4, sd_s = 1, sd_e = 0.5),
    strong = list(rho_r = 0.8, rho_c = 0.6, sd_s = 1.5, sd_e = 0.3)
  )
  for (nm in names(fields)) {
    for (seed in c(20260819L, 11L, 12L)) {
      a <- fields[[nm]]
      d <- .sfc_sim(10L, 8L, a$rho_r, a$rho_c, a$sd_s, a$sd_e, seed = seed)
      label <- paste(nm, "field, seed", seed)
      expect_no_warning(fit <- .sfc_field_fit(d), message = "did not identify")
      expect_identical(
        flexyBayes:::.fb_spatial_collapse_reasons(fit),
        character(0),
        info = label
      )
      vc <- summary(fit)$varcomp
      expect_identical(vc$note[vc$component == "sd_spatial"], "", info = label)
    }
  }
})

test_that("a weak field is flagged only where a REML test also finds none", {
  # SD 0.4 beside a nugget of 1. A REML likelihood-ratio test (5%, critical
  # value 7.05) detects the field on seeds 11 and 12 (statistic 7.24 and
  # 10.60) and not on seed 20260819 (3.93).
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = TRUE
  )
  detected <- c(`11` = TRUE, `12` = TRUE, `20260819` = FALSE)
  for (seed in names(detected)) {
    d <- .sfc_sim(10L, 8L, 0.6, 0.4, 0.4, 1, seed = as.integer(seed))
    vc <- summary(.sfc_field_fit(d))$varcomp
    expect_identical(
      vc$note[vc$component == "sd_spatial"],
      if (detected[[seed]]) "" else "collapsed",
      info = paste("weak field, seed", seed)
    )
  }
})

test_that("the fit records the comparison with the model without the field", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = TRUE
  )
  d <- .sfc_sim(10L, 8L, 0, 0, sd_s = 1e-3, sd_e = 1, seed = 20260819L)
  fit <- .sfc_field_fit(d)

  test <- fit$extras$field_test
  expect_type(test, "list")
  expect_match(test$statistic, "log Bayes factor", fixed = TRUE)
  expect_true(test$refit)
  expect_true(is.finite(test$value))
  expect_equal(test$threshold, -log(100))
  expect_equal(test$mlik_field, fit$inla$mlik[1L, 1L])
  expect_equal(
    test$value,
    test$mlik_field - test$mlik_no_field,
    tolerance = 1e-8
  )
  expect_true(test$value < test$threshold)
  expect_false(test$field_supported)
  expect_output(print(fit), "log Bayes factor", fixed = TRUE)
})

test_that("a correlation away from zero is a field, with no refit", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = TRUE
  )
  d <- .sfc_sim(10L, 8L, 0.8, 0.6, sd_s = 1.5, sd_e = 0.3, seed = 11L)
  test <- .sfc_field_fit(d)$extras$field_test
  expect_false(test$refit)
  expect_true(test$field_supported)
  expect_match(test$note, "excludes zero", fixed = TRUE)
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

test_that("a recorded comparison decides, whatever the SD and correlations", {
  medians <- c(sigma = 1)
  shell <- .sfc_shell(.sfc_vc(0.027, -0.99, 0.99), medians)

  shell$extras$field_test <- list(
    value = 2,
    threshold = -log(100),
    field_supported = TRUE
  )
  expect_identical(
    flexyBayes:::.fb_spatial_collapse_reasons(shell),
    character(0)
  )

  shell <- .sfc_shell(.sfc_vc(0.9, -0.2, 0.25), medians)
  shell$extras$field_test <- list(
    value = -7.4,
    threshold = -log(100),
    field_supported = FALSE
  )
  reasons <- flexyBayes:::.fb_spatial_collapse_reasons(shell)
  expect_length(reasons, 1L)
  expect_match(reasons, "log Bayes factor -7.40 for the field", fixed = TRUE)
})
