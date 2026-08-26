# Scaling benchmark: what "scales to large agricultural datasets" rests on

The evidence behind the scaling claim in `DESCRIPTION` and `README.md`.
The numbers are in `benchmark_scaling.csv` beside this file; this page says
what was measured, on what, and where the claim stops.

The mechanism is not a faster solver. It is that the models this package
targets are fitted to data that is **massively replicated** -- the same
handful of environments crossed with the same genotypes, over and over --
and for the exponential-family models in scope the likelihood of a *cell*
(a set of rows sharing a linear predictor) depends on those rows only
through a small vector of sufficient statistics that is **additive across
any partition of the rows**. So the statistics can be accumulated one chunk
at a time and the chunks discarded. The fit that follows is a fit to `K`
cells, not to `N` rows, and it is algebraically the same fit.

## What was measured

Two studies, one design family. Six environments crossed with genotypes,
gaussian random intercept, INLA backend, `y ~ env` with `random = ~ geno`.
Every cell runs in a **fresh subprocess under `/usr/bin/time -l`**, so
wall-clock and peak resident memory are the operating system's numbers for
that process and not an estimate taken inside a long-lived session.

| study | what it contrasts | K | sizes |
|---|---|---:|---|
| `boundary` | the per-row path against the streamed-aggregated path, over the range where both are feasible, up to the point where the per-row path is not | 360 | 1e5 to 5e7 rows |
| `extreme` | the streamed path reading partitioned `.fst` shards from disk | 1200 | 1e7 to 5e9 rows |
| `boundary-confirm` | the `boundary` cells that finish quickly, re-measured against this release | 360 | 1e5 to 1e7 rows |

`boundary` and `extreme` were run on 2026-05-31 against the development tree
released the same day as 0.6.0. `boundary-confirm` was run on 2026-08-19
against an installed tarball build of **0.9.2**, this release, so the claim
carries a current execution record rather than a historical one.

**Hardware and stack.** Apple silicon macOS, 32 GB RAM, 10 cores, R 4.5.2,
INLA 25.10.19, `fst` 0.9.8. Both runs used the same machine class; the
2026-08-19 run is on an Apple M5.

## The boundary: where the per-row path stops

| method | N | wall (s) | peak (MB) | compression |
|---|---:|---:|---:|---:|
| per-row | 100,000 | 5.17 | 1,346 | -- |
| per-row | 500,000 | 15.26 | 6,386 | -- |
| per-row | 1,000,000 | 18.96 | 12,352 | -- |
| per-row | 5,000,000 | **infeasible** -- exceeded the memory cap | | |
| streamed | 100,000 | 2.95 | 399 | 3.6e-03 |
| streamed | 1,000,000 | 2.80 | 470 | 3.6e-04 |
| streamed | 5,000,000 | 2.97 | 779 | 7.2e-05 |
| streamed | 10,000,000 | 2.92 | 1,234 | 3.6e-05 |
| streamed | 50,000,000 | 4.05 | 3,712 | 7.2e-06 |

The per-row path's memory grows with the row count and it stops between one
and five million rows on a 32 GB machine. The streamed path's wall-clock is
roughly flat across a 500-fold increase in `N`, because the fit is always a
fit to 360 cells. The streamed rows here read an **in-memory** frame, so
their `peak_mb` still carries the whole frame; the bounded-memory result is
the `extreme` study below, where the data is read from disk.

## The extreme end: a billion rows through a flat memory envelope

| N | wall (s) | peak (MB) | disk (GB) | recovered intercept |
|---:|---:|---:|---:|---:|
| 10,000,000 | 3.31 | 708 | 0.09 | 0.98542 |
| 100,000,000 | 5.67 | 796 | 0.94 | 0.98552 |
| 1,000,000,000 | 31.90 | 844 | 9.40 | 0.98565 |
| 5,000,000,000 | 150.87 | 846 | 46.98 | 0.98568 |

The column that carries the claim is `peak (MB)`: **708 MB at ten million
rows, 846 MB at five billion**, a 500-fold increase in data through a
roughly flat memory envelope, because only one chunk plus the 1,200-cell
accumulator is ever resident. Wall-clock grows roughly linearly with `N` and
is dominated by the one-pass read of the shards. The recovered intercept is
stable to four decimal places across all four scales.

## The 0.9.2 confirmation

| method | N | wall (s) | peak (MB) | intercept (0.9.2) | intercept (2026-05-31) |
|---|---:|---:|---:|---:|---:|
| per-row | 100,000 | 5.10 | 1,346 | 1.36265 | 1.36265 |
| per-row | 1,000,000 | 21.99 | 9,583 | 1.34918 | 1.34918 |
| streamed | 100,000 | 2.23 | 337 | 1.36265 | 1.36265 |
| streamed | 1,000,000 | 2.30 | 410 | 1.34918 | 1.34918 |
| streamed | 10,000,000 | 2.70 | 1,162 | 1.35327 | 1.35327 |

Same seed and same design, so the estimates are comparable directly: the
release reproduces every recovered intercept **to five decimal places**, and
the time and memory envelope is unchanged. That is the point of the
confirmation -- not that the release is faster, but that the route the
original study measured is still the route the release takes.

## The 2026-08-22 ceilings study: a realistic multi-term MET design

The studies above use a single random intercept (`random = ~ geno`) -- close
to the cheapest latent field INLA can be asked to carry. A real
multi-environment trial (MET) carries several crossed and nested random
terms at once, and that changes where the *per-row* path stops by more than
an order of magnitude. This study measures that, on
`agridat::barrero.maize` grown by appending copies of the trial network with
new years -- new environments and new genotype-by-year and
genotype-by-location-by-year levels; the genotype list itself does not grow
-- fitting the reduced, one-residual-variance variant of the Barrero Table 6
model:

```
fixed:    yield ~ loc + yearf + loc:yearf
random:   ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf
residual: ~ units
```

on INLA against ASReml, each fit in its own capped subprocess, same 32 GB /
10-core Apple silicon machine as above. Two rungs were run and both are
logged in this workspace as raw console output, not extrapolated:

| rung | N | random effects | flexyBayes / INLA | ASReml |
|---|---:|---:|---|---|
| k = 64 | 911,808 | 322,598 | **preflight refused** before any fit ran -- the design exceeds the memory ceiling (`run_scaling_k64.log`) | subprocess crashed (`could not start R ... has crashed or was killed`) |
| k = 128 | 1,823,616 | 1,031,084 | ran past preflight, then the INLA subprocess **segmentation-faulted after 2,494.4 s** (about 41.6 minutes) -- `"The inla-program exited with an error"`, a raw, untyped death (`run_inla_k128.log`) | not run at this rung |

Neither rung is a completed fit on either engine. Refusing the 911,808-row
rung outright is this package's own contract working as intended -- nothing
ran, nothing was wasted -- while the 1,823,616-row rung passes preflight and
then dies mid-solve after most of an hour, which is exactly the raw-engine-
death class the C2 fix (`inla_program_failed`) now catches, so a future run
at this size gets a typed refusal in seconds rather than 41 minutes of
silence followed by an opaque INLA string. **The flexyBayes/INLA ceiling on
this design family is therefore bracketed between 911,808 and 1,823,616
rows, and has not been measured or bisected**: there is no logged rung in
this range, on this design, where the fit actually completes.

A 911,808-row / 515,756-random-effect / 340.9-second flexyBayes *success* is
quoted in `SCALE_STRATEGY_2026-08-22.md` and `EXEC_SPEC_v0.9.3_scale_2026-08-23.md`
(workspace planning documents, not part of the package) and was carried from
there into an earlier draft of this record's CSV as a `success` row. It does
not check out against the run this workspace actually logged at that row
count: `run_scaling_k64.log` and
`dir_outcomes_barrero_scaling_22Aug2026/tab_scaling_ladder_k64.csv` both
record a **preflight refusal** with **322,598** random effects, not 515,756
-- a figure this write-up reproduced independently by re-running the ladder
script's own level-counting logic and matching it exactly. No artefact
anywhere in this workspace carries a 515,756-effect, 340.9-second, success
record for 911,808 rows. The CSV beside this file has been corrected rather
than carrying the unverified number forward a fourth time.

ASReml's own boundary on this same design family -- 58.4 s at 683,856 rows,
70.7 s at 797,832 rows, then insufficient workspace from 854,820 rows -- is
`SCALE_STRATEGY_2026-08-22.md` §1's own figure and was not independently
re-run in this session; it is at least consistent in direction with the
crash this study observed for ASReml at 911,808 rows.

## Where the claim stops

1. **Model scope.** The exactness holds for gaussian-identity,
   binomial-logit and poisson-log GLMMs with **factor** fixed effects and
   **random-intercept** grouping. Other links, continuous fixed effects,
   random slopes, structured covariance and smooths are outside the
   aggregated path and are refused before any fit, not silently approximated.
2. **Compression, not row count, is the binding quantity.** Aggregation buys
   nothing when `K` approaches `N`. The planner reports the compression so a
   user sees this before fitting, and `aggregate = TRUE` on a design with no
   compression to find is refused rather than run.
3. **Residual prior on the aggregated Gaussian path.** That path composes
   the within-cell sum-of-squares correction with INLA's residual precision
   prior; a user-supplied residual prior is not threaded through the
   composition. The count families carry no such caveat.
4. **Numerical envelope.** Counts and sums accumulate as doubles, exact for
   integer counts to `2^53`. At extreme effective sample size the custom
   log-precision prior's log-density is recentred to approximately zero at
   its mode so INLA's hyperparameter integration stays stable; the
   recentring is an exact reparametrisation and leaves small-`N` fits
   unchanged.
5. **One machine.** Every number here is one hardware class and one run per
   cell. They are an order-of-magnitude record of the memory and time
   envelope, not a distribution, and nothing here is a comparison against
   another package.

## Equivalence: the scaling is not bought with an approximation

The claim would be worth little if the aggregated path answered a different
question. Three checks, all in the shipped test suite:

- `tests/testthat/test-stream-aggregate.R` -- the streamed sufficient
  statistics equal a single-pass in-memory aggregation to a relative 1e-8,
  are invariant to where the chunk boundary falls, and the count-family sums
  are exact.
- `tests/testthat/test-aggregation-equivalence-backend.R` -- at N = 6e4 the
  per-row and streamed-aggregated INLA fits agree on the fixed effects to
  1e-4 and on the hyperparameters to 5e-3, for gaussian, binomial and
  poisson.
- The `boundary` study's own `intercept` column: at each `N` where both
  paths were run (1e5 and 1e6), they recover the same intercept to five
  decimal places, on the 2026-05-31 study and on the 0.9.2 confirmation
  alike (see the CSV).

## Reproducing it

```
Rscript tools/benchmark_bigdata_boundary.R    # the boundary study
Rscript tools/stress_bigdata_extreme.R        # the extreme study
```

Both scripts are in the repository and excluded from the build. The
`extreme` study writes tens of gigabytes of `.fst` shards and is not run
casually. The `boundary-confirm` rows were produced by the same boundary
design run against an installed tarball rather than a loaded source tree.
