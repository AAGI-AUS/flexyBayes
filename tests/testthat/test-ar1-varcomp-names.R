# The correlations of a separable AR1 field are named after the factor they
# run along, rho_<factor>, whichever factor the term writes first. The oracle
# is ASReml-R 4.2 on the simulated trials of data-raw/at_ar1_reference.R,
# which names each correlation by its factor (ROW:RANGE!ROW!cor).

.vn_ref <- function() {
  path <- system.file("extdata", "at_ar1_reference.rds", package = "flexyBayes")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "extdata", "at_ar1_reference.rds"
    )
  }
  readRDS(path)
}

.vn_fit <- function(...) {
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE
  )
  suppressMessages(flexybayes(..., backend = "inla", verbose = FALSE))
}

# Posterior median minus the ASReml estimate, in posterior SDs, for each
# correlation ASReml reports; NA where the fit has no component of that name.
.vn_gap <- function(fit, asreml_varcomp) {
  cor_rows <- asreml_varcomp[grepl("!cor$", asreml_varcomp$term), ]
  factor_name <- sub("^.*!([A-Za-z]+)!cor$", "\\1", cor_rows$term)
  section <- ifelse(
    grepl("^TRIAL_", cor_rows$term),
    sub("^TRIAL_([^!]+)!.*$", "__TRIAL\\1", cor_rows$term),
    ""
  )
  name <- paste0("rho_", factor_name, section)
  vc <- fit$extras$variance_comps
  med <- attr(vc, "posterior_median")
  i <- match(name, as.character(vc$component))
  stats::setNames(
    (unname(med[name]) - cor_rows$component) / vc$sd[i],
    name
  )
}

test_that("DEF-202: a field's correlations keep their factor's name in either order", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  ref <- .vn_ref()
  for (residual in list(~ ar1(ROW):ar1(RANGE), ~ ar1(RANGE):ar1(ROW))) {
    fit <- .vn_fit(y ~ VARIETY, residual = residual, data = ref$data$one_trial)
    gap <- .vn_gap(fit, ref$fits$one_trial$varcomp)
    expect_named(gap, c("rho_ROW", "rho_RANGE"))
    expect_false(anyNA(gap), label = deparse(residual))
    expect_true(
      all(abs(gap) <= 0.3, na.rm = TRUE),
      label = paste(deparse(residual), sprintf("%s %.2f", names(gap), gap))
    )
  }
})

test_that("DEF-202: one field per section names its correlations by factor", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  ref <- .vn_ref()
  fit <- .vn_fit(
    y ~ TRIAL + VARIETY,
    residual = ~ dsum(~ ar1(RANGE):ar1(ROW) | TRIAL),
    data = ref$data$dsum_field
  )
  gap <- .vn_gap(fit, ref$fits$dsum_field$varcomp)
  expect_setequal(
    names(gap),
    paste0(rep(c("rho_ROW", "rho_RANGE"), 3L), "__TRIALT", rep(1:3, each = 2L))
  )
  expect_false(anyNA(gap))
  expect_true(
    all(abs(gap) <= 0.3, na.rm = TRUE),
    label = paste(sprintf("%s %.2f", names(gap), gap), collapse = "; ")
  )
})

test_that("DEF-202: summary, the field panel and the draws agree on the names", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  d <- agridat::gilmour.serpentine
  d$rowf <- factor(d$row)
  d$colf <- factor(d$col)
  fits <- list(
    row_first = .vn_fit(yield ~ rep, random = ~ gen + ar1(rowf):ar1(colf),
                        data = d),
    col_first = .vn_fit(yield ~ rep, random = ~ gen + ar1(colf):ar1(rowf),
                        data = d)
  )
  med <- lapply(fits, function(f) {
    attr(f$extras$variance_comps, "posterior_median")
  })
  for (f in fits) {
    s <- summary(f)
    m <- attr(f$extras$variance_comps, "posterior_median")
    for (fac in c("rowf", "colf")) {
      nm <- paste0("rho_", fac)
      expect_true(nm %in% s$varcomp$component)
      panel <- s$spatial_field$median[
        s$spatial_field$parameter ==
          paste0("correlation along ", fac, " (", nm, ")")
      ]
      expect_equal(unname(m[[nm]]), panel, tolerance = 1e-8)
    }
    draws <- fb_variance_draws(f, 50L, seed = 1L)
    expect_true(all(c("rho_rowf", "rho_colf") %in% names(draws)))
    vp <- vpredict(f, r ~ rho_rowf, n_draws = 50L)
    expect_identical(nrow(vp), 1L)
  }
  # The rows of this trial are far more correlated than its columns
  # (0.99 against 0.81 on the row-first fit), so a swap would move the
  # medians by about 0.2.
  expect_equal(med$row_first[["rho_rowf"]], med$col_first[["rho_rowf"]],
               tolerance = 0.02)
  expect_equal(med$row_first[["rho_colf"]], med$col_first[["rho_colf"]],
               tolerance = 0.02)
})

.vn_one_factor_fit <- function() {
  set.seed(31L)
  d <- expand.grid(t = factor(1:8), g = factor(1:12))
  d$y <- stats::rnorm(nrow(d))
  .vn_fit(y ~ 1, random = ~ ar1(t):g, data = d)
}

test_that("DEF-225: a field with one AR1 factor names its correlation by factor", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fit <- .vn_one_factor_fit()
  s <- summary(fit)
  expect_true("rho_t" %in% s$varcomp$component)
  expect_false("rho_row" %in% s$varcomp$component)
  panel <- s$spatial_field$parameter
  expect_true(any(grepl("(rho_t)", panel, fixed = TRUE)))
  draws <- fb_variance_draws(fit, 20L, seed = 1L)
  expect_true("rho_t" %in% names(draws))
})

test_that("DEF-225: the old name on a new one-factor fit refuses and names the new one", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fit <- .vn_one_factor_fit()
  expect_error(
    vpredict(fit, r ~ rho_row, n_draws = 20L),
    "rho_t",
    class = "flexybayes_refusal_vpredict_component_not_in_posterior"
  )
  expect_error(
    fb_marginal_cov(fit),
    class = "flexybayes_refusal_marginal_cov_not_buildable"
  )
})

test_that("DEF-225: the names a field is given are rho_<factor> for one factor, rho_row/rho_col for two", {
  one <- flexyBayes:::.fb_ar1_rho_names(list(type = "ar1_spatial", row_var = "t",
                                             col_var = "g", col_ar1 = FALSE))
  expect_identical(one[["row"]], "rho_t")
  res <- flexyBayes:::.fb_ar1_rho_names(list(type = "ar1", var = "t"))
  expect_identical(res[["row"]], "rho_t")
  two <- flexyBayes:::.fb_ar1_rho_names(list(type = "ar1_spatial", row_var = "row",
                                             col_var = "col", col_ar1 = TRUE))
  expect_identical(unname(two), c("rho_row", "rho_col"))
})

test_that("DEF-225: a fit saved before 0.18.0 keeps rho_row in its stored table", {
  vc <- data.frame(
    component = c("sigma", "sd_g", "sd_spatial", "rho_row"),
    estimate = c(1, 0.5, 0.7, 0.4),
    sd = c(0.1, 0.1, 0.1, 0.1),
    q2.5 = c(0.8, 0.3, 0.5, 0.1),
    q97.5 = c(1.2, 0.7, 0.9, 0.7),
    stringsAsFactors = FALSE
  )
  saved <- structure(
    list(extras = list(variance_comps = vc)),
    class = c("flexybayes_inla", "flexybayes")
  )
  expect_true("rho_row" %in% flexyBayes:::.fb_summary_varcomp(saved)$component)
  expect_false(
    "rho_t" %in% flexyBayes:::.fb_summary_varcomp(saved)$component
  )
})
