# stress.R -- env-gated full-suite probe.
#
# Sets FLEXYBAYES_RUN_STRESS=true so the four currently-gated stress
# tests run (greta MCMC backend-fidelity + known-covariance vm/ped
# Phase B-greta MCMC + 1e8-row preflight).  Wall-time depends on the
# greta MCMC paths -- typically 5 -- 15 minutes on a laptop.
#
# Run from the workspace root:
#   Rscript flexyBayes/tools/profiles/stress.R

Sys.setenv(FLEXYBAYES_RUN_STRESS = "true", NOT_CRAN = "true")
suppressMessages(devtools::load_all("flexyBayes", quiet = TRUE))

t0 <- Sys.time()
res <- testthat::test_dir(
  path = "flexyBayes/tests/testthat",
  reporter = "progress",
  stop_on_failure = FALSE
)
df <- as.data.frame(res)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("\n=== stress profile ===\n")
cat(sprintf(
  "PASS: %d  FAIL: %d  SKIP: %d  WALL: %.1fs\n",
  sum(df$passed),
  sum(df$failed),
  sum(df$skipped),
  elapsed
))

if (sum(df$failed) > 0L) {
  cat("\nStress profile FAILED.  Inspect output above.\n")
  quit(save = "no", status = 1L)
}
