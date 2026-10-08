# The correlations of a separable AR1 field, as summary()$varcomp names
# them and as the field panel names them.

test_that("DEF-202: varcomp names an AR1 field's correlations by factor", {
  skip(paste(
    "DEF-202 open: rho_row and rho_col name the first and second factor",
    "of the term, so ar1(col):ar1(row) reports the column correlation",
    "as rho_row"
  ))
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  d <- agridat::gilmour.serpentine
  d$row <- factor(d$row)
  d$col <- factor(d$col)
  fit <- withr::with_options(
    list(flexyBayes.inla_threads = "1:1",
         flexyBayes.silence_default_prior_note = TRUE),
    flexybayes(yield ~ rep, random = ~ gen + ar1(col):ar1(row),
               data = d, verbose = FALSE)
  )
  s <- summary(fit)
  along_row <- s$spatial_field$median[
    s$spatial_field$parameter == "correlation along row (rho_row)"
  ]
  med <- attr(fit$extras$variance_comps, "posterior_median")
  expect_equal(unname(med[["rho_row"]]), along_row, tolerance = 1e-8)
})
