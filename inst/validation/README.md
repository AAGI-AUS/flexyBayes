# Evidence snapshot

Four text artefacts recording the numerical checks the package's public
claims rest on. They are copies, not the working trees: the scripts, the
fitted objects and the rendered reports live outside this package, and each
entry below names the script that produced its artefact so a reader can
reproduce it rather than take the numbers on trust.

Everything here is plain text or CSV. Nothing in this directory is read by
package code, and nothing is regenerated at build or check time.

Locate the directory from R with

``` r
system.file("validation", package = "flexyBayes")
```

## `oracle_het.log`

**What.** Two probes of the heterogeneous-variance emit against
independent implementations, on simulated data where the truth is known.
Probe 1 puts `diag(site):gen` against ASReml and `lme4`'s `dummy()`
expansion; probe 2 puts the `dsum(~ units | site)` residual against ASReml
and `nlme::varIdent`, and against a direct `brms::brm()` call under brms's
own priors, so the package's injected prior is measured rather than
assumed. Every probe asserts the parameter *count* before comparing a
value: three free variances and zero covariances, and no scalar `sigma` on
either Bayesian arm.

**Generating script.** `design-preserving-missingness/oracle_heterogeneous.R`
(workspace, outside this package).

**Data.** Simulated. Probe 1: 3 sites x 40 genotypes x 4 replicates = 480
plots, genotype SD by site `c(0.5, 1.2, 2.2)`, residual SD 0.6, seed
`20260814`. Probe 2: 3 sites x 30 genotypes x 4 replicates = 360 plots,
genotype SD 0.8, residual SD by site `c(0.4, 1.0, 2.0)`, seed `7`.

**Run.** 2026-08-15 17:45 ACST. asreml 4.2.0.392, lme4 2.0.1, nlme 3.1.168,
brms 2.23.0, rstan 2.32.7, R 4.5.2.

**Reading it.** The three per-site residual variances the package returns
through `flexybayes()` are the canonical set quoted in `NEWS.md`, in the
comment at `R/emit_brms.R`, and in the review brief. They are one run's
realisation at bulk ESS above 8,000; the percentage agreement with ASReml,
not the digits, is what the prose leads with.

## `report_sweep.txt`

**What.** The missing-response sweep: five designs (randomised complete
blocks, incomplete blocks, row-column, MET, spatial) crossed with
missingness fractions, each at 25 seeds, comparing the augmentation path
against ASReml and against deletion. The report covers which cells are
banked, which methods refused, the recovered variance components, and where
the estimates degrade.

**Generating scripts.** `design-preserving-missingness/sweep_full.R` runs
the cells and banks them; `design-preserving-missingness/report_sweep.R`
renders this text from the bank;
`design-preserving-missingness/topup_asreml_spatial.R` supplied the
superseding ASReml spatial rows the header notes.

**Data.** Simulated across the design catalogue, 25 seeds per cell, 170
banked cells, 66,388 rows in total.

**Run.** 2026-08-14. The sweep takes about six hours and is not re-run
casually.

## `head_to_head_summary.csv`

**What.** One row per fitted variance component across five models, with
the simulating truth, the ASReml estimate, the flexyBayes estimate, each
engine's signed error against the truth, and the distance between the two
fits. Accuracy, agreement and the difference between them are kept as
separate columns on purpose: agreement says whether two engines are doing
the same thing, and only the error columns say which answer is closer to
the truth.

**Generating script.** `head-to-head/run_comparison.R` (workspace).

**Data.** Simulated, five models -- randomised complete blocks at two
sizes, incomplete blocks, a MET with genotype variance by site, and a MET
with residual variance by site. Seed `20260815`, 4 chains x 2000
post-warmup draws on the sampled arms.

**Run.** 2026-08-15. asreml 4.2.0.392, brms 2.23.0, INLA 25.10.19,
flexyBayes 0.9.0, R 4.5.2.

## `head_to_head_models.csv`

**What.** The per-model companion to the row above: observation count,
wall-clock seconds for each engine, the parameter count each engine fitted,
which flexyBayes backend the model routed to, and the sampling diagnostics
where the route was a sampler. The parameter counts are the structural
check -- a timing comparison between engines fitting different numbers of
parameters would be measuring the wrong thing.

**Generating script, data, run.** As `head_to_head_summary.csv` above; both
files are projections of the same `comparison.rds`.
