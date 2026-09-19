# helper-test-tiers.R -- the long test tier.
#
# A block whose measured time would dominate a local `devtools::test()` or
# `devtools::check()` calls skip_if_not_long_tier() as its first line after
# any skip_on_cran() / skip_if_not_installed() guards, and changes nothing
# else. It then runs only when the environment variable FLEXYBAYES_RUN_LONG
# is "true". Its assertions are unchanged, so the long tier is a subset of
# the suite, not a variant of it.
#
# The release suite sets the variable (tools/run_test_suite.R), and criterion
# D11 of tools/definition_of_done.R fails when any block was skipped here or
# when a block tagged here is missing from the suite record. Tag top-level
# test_that() blocks with a literal description, so D11 can find them.
# CONTRIBUTING.md ("Test tiers") says how to choose a block for this tier.

skip_if_not_long_tier <- function() {
  if (!identical(Sys.getenv("FLEXYBAYES_RUN_LONG"), "true")) {
    testthat::skip("FLEXYBAYES_RUN_LONG != \"true\"; long tier")
  }
  invisible(TRUE)
}
