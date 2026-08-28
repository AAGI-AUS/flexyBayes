# Evidence snapshot

The numerical checks the package's public claims rest on: the evidence
artefacts, the scenario registry the validation ladder reads, and the
execution-grid ledger. They are copies, not the working trees: the scripts, the
fitted objects and the rendered reports live outside this package, and each
entry below names the script that produced its artefact so a reader can
reproduce it rather than take the numbers on trust.

Everything here is plain text, CSV or YAML. Nothing in this directory is read by
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

## `benchmark_scaling.csv` and `benchmark_scaling.md`

**What.** The evidence behind the scaling claim in `DESCRIPTION` and
`README.md`. `benchmark_scaling.md` says what was measured and where the
claim stops; the CSV is the slim table, eighteen rows across three studies:
the boundary between the per-row and the streamed-aggregated fit path
(where the per-row path stops, and what the streamed one costs at the same
size), the extreme end (partitioned `.fst` shards to five billion rows),
and a re-measurement of the fast boundary cells against this release. Times
and peak memory are the operating system's numbers for a fresh subprocess
per cell (`/usr/bin/time -l`), not estimates taken inside a long-lived
session.

**Generating scripts.** `tools/benchmark_bigdata_boundary.R` and
`tools/stress_bigdata_extreme.R`, in this repository and excluded from the
build. The `boundary-confirm` rows were produced by the same boundary
design run against an installed tarball rather than a loaded source tree.

**Data.** Simulated. Six environments crossed with 60 genotypes (360 cells)
for the boundary rows, six by 200 (1,200 cells) for the extreme rows;
gaussian random intercept, INLA backend, seed `42`.

**Run.** `boundary` and `extreme` on 2026-05-31 against the development
tree released the same day as 0.6.0; `boundary-confirm` on 2026-08-19
against an installed tarball of 0.9.2. Apple silicon macOS, 32 GB RAM, 10
cores, R 4.5.2, INLA 25.10.19, fst 0.9.8.

**Reading it.** The claim rests on the `peak_mb` column, not on the wall
clock: the per-row path's memory grows with the row count and stops between
one and five million rows, while the streamed path holds one chunk plus the
cell accumulator, so 708 MB at ten million rows becomes 846 MB at five
billion. The `intercept` column is the cross-check on the other half --
the aggregated path recovers the per-row estimate to five decimal places at
every size where both were run, so the scaling is not bought with an
approximation. One run per cell on one machine: this is an
order-of-magnitude envelope, not a distribution, and it compares the
package's two paths against each other rather than against another package.

## `scenarios.yaml`

**What.** The numerical-validation registry: one row per validation
scenario, in the schema the validation ladder reads
(`Config/rpkg/validationTier` in `DESCRIPTION`, `V2` for this package).
Forty-one study rows across five studies, plus three stability cells and
one calibration study: forty-five registered scenarios in all. Every row is an existing study. Nothing was simulated
to reach a floor, and where a floor is not reached the shortfall is
written down below rather than filled with a row that would not survive
being read.

**Oracle classes present.** M-A (a declared expected behaviour at the
edge of the domain), M-B (an independently written REML that shares no
code with the package), M-D (ASReml, `lme4`, `glmmTMB`, `nlme`,
`MCMCglmm`), M-E (simulation with known truth).

**`result_cache` and the R gate.** Each row names the artefact its
numbers come from. Rows whose study ships a copy here point at that copy
(`report_sweep.txt`, `oracle_het.log`, `head_to_head_summary.csv`,
`execution_grid/ledger.csv`); rows whose evidence is a workspace
artefact point at it by its workspace-relative path. The consequence is
worth stating plainly: the R-side gate `check_scenarios.R` reads each
`result_cache` with `readRDS()` and fails on anything that is not a
serialised R object, so **that gate does not pass on this registry
today** -- the shipped artefacts are text and CSV projections of fitted
objects that are deliberately kept outside the package, and the fitted
objects themselves are not shipped. The structural gates (F20, F21, F22
in `rpkg_audit.py`) read the registry itself and do pass. Making the R
gate pass would mean either shipping the fitted objects or teaching the
gate to read a text artefact, and that is a decision about the gate
rather than about the evidence.

## `execution_grid/`

**What.** The F28 execution-grid evidence (recipe 62). `ledger.csv` has
one row per grid cell -- the claim-derived cross product of the
capability matrix, the prior routes, the malformed-prior corpus, the
family allowlist, the structure and grammar surfaces, and the
prior-translation table -- each carrying the outcome a live call
produced against an installed build. `capability_matrix.csv` is the
machine-readable export of the generated capability table, so ledger
coverage of every claimed cell can be checked without parsing R.
`roster_diff.csv` records every response family an installed engine
carries and the entry allowlist refuses, with the layer that owns each
boundary. `misses.csv` is the M1-M6 register: one row per open miss from
the last run. `sessioninfo.txt` names the engines the cells ran against.

**Generating script.** `tools/execution_grid.R`, in this repository and
excluded from the build. It runs against an installed build or a built
tarball, never `pkgload::load_all()`, with one `callr::r()` child per
cell.

**Reading it.** `class` is the outcome vocabulary -- `fit`,
`refuse_typed`, `error_untyped`, `construct_accept`, `crash`, `timeout`.
`expected` is derived from a claim surface before the run and
`expect_src` in the wide record names the surface line by line.
`verdict` is the comparison. The two numbers that matter are the count
of `error_untyped` and the count of `DIVERGENT`, and both are meant to
be zero.

## Validation tier V2: what is declared, and what V3 would add

`DESCRIPTION` declares `Config/rpkg/validationTier: V2`. V2 is what the
studies registered above actually support.

An earlier line of this package declared V3 -- the tier the ladder
assigns when outputs feed agronomic decisions -- and waived four of its
floors. That waiver was written for 0.9.2 and expired at the next
release. It is not re-issued: a second consecutive waiver on the same
four floors is the pattern a waiver exists to prevent, and declaring the
tier the evidence supports is better than declaring one above it and
carrying an exception. The four floors below are the next validation
arc, and they are recorded here rather than left to be discovered.

- Reason (F21, stochastic calibration): the registry carries one
  calibration scenario, `calibration-interval-coverage-recovery`. It is
  a real interval-coverage study run through `flexybayes()` with known
  truth, and it is short of the tier on three counts: 10 replicates per
  cell against the V3 floor of 1,000, a `0.7.0.9000` candidate rather
  than this release, and a backend set that still included the native
  engine withdrawn entirely in 0.9.3 (see `NEWS.md`). The package has no
  simulation-based calibration of its own estimator at all: the
  2026-06-02 multi-backend SBC study calibrated the ENGINES through
  direct `cmdstanr`, `INLA::inla()`, and that withdrawn engine's own
  sampler calls, so it is evidence about Stan, INLA, and that engine
  rather than about this package's emit, and registering it here would
  have claimed something it does not show. Writing a calibration entry
  that the study does not support was the alternative, and it is worse
  than a waiver.
- Reason (F20, named cells): the V3 breadth table asks for a
  null-recovery scenario at two or more sample sizes. The registered
  studies simulate non-null truths throughout, so no row claims one. The
  near-boundary signal cells the tier also asks for are covered
  (`stability-spatial-field-collapse` at a variance running to its
  floor, `stability-funnel-divergence-brms` at a between-cluster SD of
  0.02).
- Reason (registry currency): the five registered studies were run
  between 2026-06-08 and 2026-08-19 against builds from `0.7.0.9000` to
  `0.9.2`. The two head-to-head studies and the alternatives hunt are on
  the 0.9.x line; the missing-response sweep and the coverage study
  predate it.
- Reason (F24 and F25, the Mathematical Foundations Document): a
  declared V3 tier on a numerical-method package asks for `math/` --
  a foundations document carrying the derivation of every terminal
  identity the code implements, and a per-equation grounding register
  tying each one to a function, an anchor, its tests, its scenarios and
  an independent grounding reference. The package has neither file. The
  work is scheduled as **its own validation arc**, paired with the
  simulation-based-calibration study that lifts the F21 waiver above
  and sharing its derivation content with the planned methods paper,
  because the two need the same algebra written down once: the
  aggregated sufficient-statistic identities, the prior-translation
  closed forms on the standard-deviation scale, and the structured
  covariance parameterisations. Assembling a foundations document from
  the existing sources inside this release would meet the file check
  and not the floor -- the derivations would be sketched rather than
  complete, and the grounding column would carry the package's own
  tests as the authority for identities the package itself asserts,
  which is the self-oracle failure the register exists to prevent. A
  document that fails its own gate is worse than a recorded absence.

- Owner: Max Moldovan

These four floors are the next validation arc, not a permanent state:
a simulation-based calibration and an interval-coverage study of at
least 1,000 replicates, run through `flexybayes()` on the two active
engines and including a null-recovery cell at two sample sizes, is what
lifts the F20 and F21 waivers, and `math/foundations.qmd` plus a
`math/math_map.yaml` whose every equation carries a function, an
anchor, its tests, its scenarios and a grounded independent reference
is what lifts F24 and F25. Both are scheduled into the same arc.
