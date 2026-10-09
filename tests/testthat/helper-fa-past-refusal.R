# Fit fa(f, k):g on INLA past the dispatch refusal
# (factor_analytic_not_repeatable), for the tests that check the fit itself
# against a known truth or an ASReml reference.
.fa_fit_past_refusal <- function(fixed, random, data, residual = NULL) {
  testthat::local_mocked_bindings(
    .refuse_unrepresentable_structures = function(fb) invisible(NULL),
    .package = "flexyBayes"
  )
  args <- list(
    fixed = fixed,
    random = random,
    data = data,
    backend = "inla",
    verbose = FALSE
  )
  if (!is.null(residual)) {
    args$residual <- residual
  }
  suppressMessages(do.call(flexybayes, args))
}
