# smoke.R -- fast local sanity profile for the flexyBayes test suite.
#
# Pure-R subset only: preflight + dispatch + IR shape + canonical names
# + backend-decision metadata.  No greta MCMC, no INLA, no tarball
# build.  Target wall-time: < 30 s on a laptop.
#
# Run from the workspace root:
#   Rscript flexyBayes/tools/profiles/smoke.R

suppressMessages(devtools::load_all("flexyBayes", quiet = TRUE))

t0 <- Sys.time()
res <- testthat::test_dir(
  path = "flexyBayes/tests/testthat",
  filter = "preflight|dispatch|backend-decision|backend-routing|canonical-names|ir-shape|aggregate-plan|prior-summary|review-code|fb-from-asreml",
  reporter = "summary",
  stop_on_failure = FALSE
)
df <- as.data.frame(res)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("\n=== smoke profile ===\n")
cat(sprintf(
  "PASS: %d  FAIL: %d  SKIP: %d  WALL: %.1fs\n",
  sum(df$passed),
  sum(df$failed),
  sum(df$skipped),
  elapsed
))

if (sum(df$failed) > 0L) {
  cat("\nSmoke profile FAILED.  Inspect output above.\n")
  quit(save = "no", status = 1L)
}
