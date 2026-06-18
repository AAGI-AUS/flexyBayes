# Benchmark / probe profiles

Four named profiles drive the flexyBayes verification surface (v040-plan
section 11). Each is a standalone `Rscript` run from the workspace root.

| Profile | Run | Gate criterion |
|---|---|---|
| `smoke` | `Rscript flexyBayes/tools/profiles/smoke.R` | Fast local sanity: the test tally is green (PASS > 0, FAIL == 0). Tally-only; no backend probes. |
| `release` | `Rscript flexyBayes/tools/profiles/release.R` | Release readiness: full tally green **and** the source tree builds to a tarball that passes `R CMD check --as-cran` at the structural-NOTE baseline. |
| `stress` | `FLEXYBAYES_RUN_STRESS=true Rscript flexyBayes/tools/profiles/stress.R` | Backend fidelity: the greta / INLA optional full-MCMC probes (stress-gated, skipped on CRAN/CI) reproduce their reference posteriors within Monte-Carlo tolerance. |
| `huge_n_probe` | `FLEXYBAYES_RUN_HUGE_N=true Rscript flexyBayes/tools/profiles/huge_n_probe.R` | Huge-N trust: the preflight, aggregation-plan, and prediction-kernel file-output paths all handle the `N = 5e6` data shape without fault (a fault is a release blocker). |

The `smoke` and `release` profiles run unconditionally; `stress` and
`huge_n_probe` are env-gated so they never run on CRAN / CI. The smoke
counterpart of `huge_n_probe` (the format-resolution contract at
`N = 5e6`) runs in the ordinary test suite via
`tests/testthat/test-huge-n-probe.R`.
