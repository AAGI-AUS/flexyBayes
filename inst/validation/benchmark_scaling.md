# Scaling benchmark: the streamed fit of a count model

This page records what the streamed, aggregated fit of a binomial or Poisson
model costs in time and memory, where the per-row fit stops, and where the
evidence stops. It supports the "Exact sufficient-statistic aggregation
(binomial, Poisson)" row of the capability table in `README.md`.
`DESCRIPTION` makes no scaling claim. The numbers are in
`benchmark_scaling.csv` beside this file.

The mechanism is not a faster solver. The models this package targets are
fitted to data that is heavily replicated: the same handful of environments
crossed with the same genotypes, over and over. For a binomial or Poisson
model, the likelihood of a *cell* (the rows that share a linear predictor)
depends on those rows only through a few sufficient statistics that add up
across any partition of the rows. The statistics can therefore be
accumulated one chunk at a time and each chunk discarded, and the fit that
follows is a fit to `K` cells rather than `N` rows with the same posterior.

In 0.12.0 this holds for the count families only. A Gaussian model is fitted
per row: `aggregate = TRUE` and `flexybayes_stream()` refuse it by name
(`aggregation_residual_prior_not_composed`), because the aggregated Gaussian
likelihood does not carry a stated or default residual prior. The
section "Earlier record: the streamed Gaussian route" below explains what
that means for the measurements published before this release.

## What was measured

Two studies on one design family: a Poisson model with a log link, six
environments crossed with genotypes, `y ~ env` with `random = ~ geno`, a
stated prior `sd(group = "geno") ~ uniform(lower = 0, upper = 3)`, INLA
engine. Every cell runs in a **fresh subprocess under `/usr/bin/time -l`**,
so wall-clock time and peak resident memory are the operating system's
numbers for that process, not an estimate taken inside a long-lived session.

| study | what it contrasts | K | sizes |
|---|---|---:|---|
| `boundary` | the per-row fit against the streamed fit, over the range where both run | 360 | 1e5 to 5e7 rows |
| `extreme` | the streamed fit reading partitioned `.fst` shards from disk | 1,200 | 1e7 to 5e9 rows |

**Machine and stack.** Apple M5, 32 GB RAM, 10 cores, macOS 26.6.2; R 4.6.1,
INLA 26.8.7, `fst` 0.9.8. Run on 2026-09-15 against the development tree of
0.12.0, with no other computation running.

## The boundary: per-row against streamed

| method | N | wall (s) | peak (MB) | intercept |
|---|---:|---:|---:|---:|
| per-row | 100,000 | 3.08 | 992 | 1.27183 |
| per-row | 500,000 | 7.87 | 3,955 | 1.26747 |
| per-row | 1,000,000 | 14.00 | 7,938 | 1.27112 |
| per-row | 5,000,000 | 108.49 | 19,542 | 1.26970 |
| streamed | 100,000 | 1.48 | 450 | 1.27183 |
| streamed | 1,000,000 | 1.43 | 500 | 1.27112 |
| streamed | 5,000,000 | 1.54 | 742 | 1.26970 |
| streamed | 10,000,000 | 1.71 | 991 | 1.26862 |
| streamed | 50,000,000 | 2.93 | 3,042 | 1.26884 |

The per-row fit's memory grows with the row count: it completed at five
million rows, using 19.5 GB of a 32 GB machine, and was not run beyond that.
The streamed fit's wall-clock time stays between 1.4 and 2.9 seconds across a
500-fold increase in `N`, because it is always a fit to 360 cells. At every
size where both were run (1e5, 1e6 and 5e6 rows) the two recover the same
intercept to five decimal places.

The streamed rows in this study read an **in-memory** data frame, so their
peak memory includes the whole frame. The bounded-memory result is the
`extreme` study, where the data is read from disk.

## The extreme end: five billion rows in flat memory

| N | wall, fit (s) | peak (MB) | disk (GB) | intercept |
|---:|---:|---:|---:|---:|
| 10,000,000 | 2.27 | 712 | 0.03 | 1.29528 |
| 100,000,000 | 5.37 | 730 | 0.28 | 1.29592 |
| 1,000,000,000 | 37.71 | 734 | 2.81 | 1.29595 |
| 5,000,000,000 | 183.27 | 734 | 14.06 | 1.29591 |

Peak memory is **712 MB at ten million rows and 734 MB at five billion**, a
500-fold increase in data with memory nearly unchanged, because only one
chunk and the 1,200-cell accumulator are ever resident. The fit time grows
roughly linearly with `N` and is dominated by the single pass over the
shards. Writing the shards is not included in the fit time: it took 0.9, 6.5,
61.3 and 302.0 seconds at the four sizes. The recovered intercept varies by
less than 0.0007 across the four sizes.

The CSV records each of these rows with `family = "poisson"`. The `family`
column was added in this revision, because the table now holds a Poisson
study and a Gaussian one.

## A realistic multi-environment design: where the per-row fit stops

The studies above use a single random intercept, close to the cheapest
latent field INLA can be asked to carry. A multi-environment trial carries
several crossed and nested random terms at once, and that moves the point
where the *per-row* fit stops by more than an order of magnitude. On
2026-08-22, against 0.9.2 on a 32 GB, 10-core Apple silicon machine, the
`agridat::barrero.maize` trial network was grown by appending copies of it
with new years (new environments and new genotype-by-year and
genotype-by-location-by-year levels; the genotype list itself does not grow)
and the reduced, one-residual-variance form of the Barrero Table 6 model was
fitted per row:

```
fixed:    yield ~ loc + yearf + loc:yearf
random:   ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf
residual: ~ units
```

| N | random effects | flexyBayes on INLA |
|---:|---:|---|
| 911,808 | 322,598 | **refused by the design-memory preflight** before any fit ran |
| 1,823,616 | 1,031,084 | passed the preflight, then the INLA subprocess **failed with a segmentation fault after 2,494.4 s** (about 41.6 minutes) |

Neither size produced a completed fit. The package now reports the second
kind of failure by name (`inla_program_failed`) instead of passing on the
engine's message. **The per-row INLA ceiling on this design family lies
between 911,808 and 1,823,616 rows and has not been measured more closely.**
ASReml, run beside it, failed at the smaller size and was not run at the
larger one.

An earlier revision of this record listed a completed flexyBayes fit at
911,808 rows. No run log supports it: the run recorded at that size is the
preflight refusal above, and the CSV row was corrected accordingly.

Gaussian models are fitted per row in 0.12.0 as they were in 0.9.2, so this
study describes the route 0.12.0 takes. Its two rows carry
`family = "gaussian"`.

## Where the evidence stops

1. **Model scope.** Aggregation is exact for binomial models with a logit
   link and Poisson models with a log link, with **factor** fixed effects and
   **random-intercept** grouping. Other links, continuous fixed effects,
   random slopes, structured covariance and smooths are outside the
   aggregated route and are refused before any fit, not approximated.
   Gaussian models are refused on this route in 0.12.0, as described above.
2. **Compression, not row count, is the binding quantity.** Aggregation buys
   nothing when `K` approaches `N`. The planner reports the compression, so it
   is visible before fitting.
3. **Evidence.** The marginal likelihood of an aggregated fit is not reported
   (`aggregated_evidence_unavailable`): the engine's value omits a constant
   that depends on the data.
4. **Numerical range.** Counts and sums accumulate as doubles, which are
   exact for integer counts up to `2^53`.
5. **One machine, one run per cell.** These numbers record the order of
   magnitude of the time and memory envelope on one hardware class. They are
   not a distribution, and nothing here compares against another package's
   aggregation. A machine with less memory stops the per-row fit earlier:
   five million rows used 19.5 GB here.

## Equivalence: the scaling does not rest on an approximation

The scaling would be worth little if the aggregated fit answered a different
question. Three checks:

- `tests/testthat/test-stream-aggregate.R`: the streamed sufficient
  statistics equal a single-pass in-memory aggregation to a relative 1e-8,
  do not depend on where the chunk boundaries fall, and the count-family sums
  are exact.
- `tests/testthat/test-aggregation-equivalence-backend.R`: for binomial and
  Poisson models, the per-row and streamed INLA fits at N = 6e4, and the
  per-row and in-memory aggregated fits at N = 3e4, agree on the fixed
  effects to 1e-4 and on the hyperparameters to 5e-3.
- The `boundary` study's `intercept` column: at 1e5, 1e6 and 5e6 rows the
  per-row and streamed fits recover the same intercept to five decimal
  places.

## Earlier record: the streamed Gaussian route (2026-05-31 and 2026-08-19)

Before 0.12.0 this page measured the same two studies on a **Gaussian**
random-intercept model. On 2026-05-31, against the development tree released
as 0.6.0 (R 4.5.2, INLA 25.10.19), the per-row fit stopped between one and
five million rows, exceeding the memory cap at five million, and the streamed
fit from `.fst` shards held peak memory between 708 MB at ten million rows
and 846 MB at five billion. On 2026-08-19 an installed build of 0.9.2
reproduced the fast cells of that study to five decimal places in the
recovered intercept.

The streamed Gaussian route those measurements describe is **refused in
0.12.0**. Its aggregated likelihood combined the within-cell sums of squares
with a fixed prior on the residual precision, so a stated or default residual
prior was not the prior fitted. The route returns only once the aggregated
Gaussian likelihood carries the residual prior exactly, and it will be
measured again then. The Gaussian rows were removed from the CSV, and the
figures above are the ones that record carried.

## Reproducing it

Both studies draw the same data-generating model:

```r
set.seed(2024L)
env_effect  <- rnorm(6, 0, 0.3)
geno_effect <- rnorm(n_geno, 0, 0.25)        # 60 genotypes (boundary), 200 (extreme)
# each row: environment and genotype sampled uniformly,
# y ~ Poisson(exp(1 + env_effect[env] + geno_effect[geno]))
prior <- fb_prior(sd(group = "geno") ~ uniform(lower = 0, upper = 3))

# boundary, rows drawn with set.seed(42L)
flexybayes(y ~ env, random = ~geno, data = tab, family = "poisson",
  backend = "inla", aggregate = FALSE, prior = prior)
flexybayes_stream(y ~ env, random = ~geno, source = tab, family = "poisson",
  backend = "inla", chunk_rows = 2e6, prior = prior)

# extreme: shards of 5e7 rows written with fst::write_fst(compress = 60),
# shard i drawn with set.seed(1000L + i)
flexybayes_stream(y ~ env, random = ~geno, source = shard_paths,
  family = "poisson", chunk_rows = 5e6, prior = prior, fit = TRUE)
```

Each fit ran as its own `Rscript` process under `/usr/bin/time -l`, and the
five-billion-row size was run only after the one-billion-row size projected
it under three hours. The extreme study writes 14 GB of shards at its largest
size and is not run casually. `tools/benchmark_bigdata_boundary.R` and
`tools/stress_bigdata_extreme.R` generate the earlier Gaussian design.
