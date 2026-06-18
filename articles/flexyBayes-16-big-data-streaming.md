# Big data: streaming exact aggregation to a billion rows

A collaborator shares a dataset with hundreds of millions to a billion
rows, already written to disk in an efficient columnar format – a
multi-season multi-environment trial, a genomic panel, a sensor archive.
You want to fit a Bayesian mixed model, but the data will not fit in
memory, and pushing every row through a sampler would be hopeless. This
vignette shows how
[`flexybayes_stream()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes_stream.md)
fits such a model exactly, with memory that does not grow with the row
count. Fitting code is shown for reference – set `eval = TRUE` locally
with the relevant backends installed to run it.

## The idea: aggregate exactly, then fit

The data is almost always massively replicated – the same handful of
environments crossed with the same genotypes, over and over. For the
exponential-family models flexyBayes targets, the likelihood of a *cell*
(a set of rows that share the same linear predictor) depends on the rows
only through a few *sufficient statistics*, and those statistics are
additive across any partition of the rows. They can therefore be
accumulated one chunk at a time and the chunks discarded.

For a cell $`k`$ with $`n_k`$ rows sharing the linear predictor
$`\eta_k`$:

- gaussian (identity link) needs $`\sum y_i`$ and $`\sum y_i^2`$. The
  cell log-likelihood is algebraically identical to the per-row form, so
  this is exact, not approximate.
- poisson (log link) needs $`\sum y_i`$ and $`\sum E_i`$ (exposure),
  because independent Poissons with a shared rate add.
- binomial (logit link) needs $`\sum y_i`$ (successes) and $`\sum m_i`$
  (trials), because independent Binomials with a shared probability add.

In every case the cell statistics are sums, and a sum over the whole
dataset is the sum of the sums over its chunks. The fitted posterior is
the full-data posterior. The cost of the fit is set by the number of
distinct cells $`K`$, not the number of rows $`N`$.

## Fitting from an `.fst` file

The realistic delivery format is a columnar file on disk.
[`flexybayes_stream()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes_stream.md)
reads it in row-range chunks, accumulates the sufficient statistics, and
fits the $`K`$-cell model through INLA (or greta). Peak memory is one
chunk plus the small cell accumulator, never the full table.

``` r

library(fst)

# A dataset already on disk (here written for illustration).
path <- "trial.fst"

fit <- flexybayes_stream(
  yield ~ environment,         # fixed factor(s)
  random = ~ genotype,         # random-intercept grouping factor
  source = path,               # an .fst path on disk
  family = "gaussian",
  backend = "inla",
  chunk_rows = 5e6             # rows read per chunk
)

summary(fit)
```

For a dataset shared as several partitioned shards – the usual shape for
genuinely large data – pass the vector of shard paths and they are
streamed in turn with the same bounded footprint.

``` r

shards <- sort(list.files("trial_parts", pattern = "\\.fst$",
                          full.names = TRUE))
fit <- flexybayes_stream(
  yield ~ environment, random = ~ genotype,
  source = shards, family = "gaussian", backend = "inla"
)
```

Binomial and poisson models fit the same way. A 0/1 response is treated
as Bernoulli (one trial per row). For poisson, name an exposure column
with `exposure =` if the data carries one.

``` r

fit_bin <- flexybayes_stream(
  infected ~ region, random = ~ paddock,
  source = "survey.fst", family = "binomial", backend = "inla"
)

fit_pois <- flexybayes_stream(
  count ~ treatment, random = ~ site,
  source = "counts.fst", family = "poisson", backend = "inla"
)
```

To inspect the compression a design will achieve before committing to a
fit, pass `fit = FALSE` and the function returns the aggregation
summary.

``` r

agg <- flexybayes_stream(
  yield ~ environment, random = ~ genotype,
  source = path, fit = FALSE
)
agg          # prints N, K, and the compression ratio
```

## The in-memory entry aggregates too

When the data does fit in memory, the ordinary entry
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
reaches the same exact aggregation through `aggregate = "auto"` (the
default) on `backend = "inla"`, for the gaussian, binomial, and poisson
families. The two routes return the same posterior.

``` r

fit_a <- flexybayes(yield ~ environment, random = ~ genotype,
                    data = trial_df, backend = "inla")   # auto-aggregates
fit_b <- flexybayes(yield ~ environment, random = ~ genotype,
                    data = trial_df, backend = "inla",
                    aggregate = FALSE)                   # per row
```

## What it costs: per-row versus aggregated

The table below compares the two paths at sample sizes where both are
feasible, on a 6-environment by 60-genotype design (360 cells) with a
gaussian random intercept and the INLA backend. Memory is the operating
system peak resident size of a fresh process.

| method | N | wall (s) | peak (MB) | intercept |
|----|---:|---:|---:|---:|
| per row | 100,000 | 5.2 | 1,346 | 1.36265 |
| per row | 500,000 | 15.3 | 6,386 | 1.35696 |
| per row | 1,000,000 | 19.0 | 12,352 | 1.34918 |
| per row | 5,000,000 | infeasible – exceeded the memory budget |  |  |
| streamed | 100,000 | 2.9 | 399 | 1.36265 |
| streamed | 1,000,000 | 2.8 | 470 | 1.34918 |
| streamed | 50,000,000 | 4.0 | 3,712 | 1.35482 |

The per-row path’s time and memory grow with $`N`$ – 12 GB at one
million rows – and it became infeasible at five million on a 32 GB
machine. The streamed path returns the same intercept (1.34918 against
1.34918 at one million rows) in about three seconds.

## How far it scales

The next table streams a partitioned `.fst` dataset on a 6-environment
by 200-genotype design (1200 cells). The dataset is written as
50-million-row shards and streamed in five-million-row chunks.

| N             |    K | wall (s) | peak (MB) | on disk (GB) | intercept |
|---------------|-----:|---------:|----------:|-------------:|----------:|
| 10,000,000    | 1200 |      3.3 |       708 |         0.09 |   0.98542 |
| 100,000,000   | 1200 |      5.7 |       796 |         0.94 |   0.98552 |
| 1,000,000,000 | 1200 |     31.9 |       844 |         9.40 |   0.98565 |
| 5,000,000,000 | 1200 |    150.9 |       846 |        46.98 |   0.98568 |

Peak memory is essentially flat – 708 MB at ten million rows, 846 MB at
five billion – because only one chunk plus the 1200-cell accumulator is
ever resident. A five-billion-row gaussian mixed model fits exactly in
about two and a half minutes of streaming, reading 47 GB of shards once
from disk. The recovered intercept is stable across all four scales: the
data volume grows five-hundred-fold while the memory footprint and the
fitted estimate do not move.

## When aggregation does not apply

The method compresses only when the design is low-cardinality. A
continuous fixed covariate gives one cell per distinct value, so there
is no compression, and the model is refused before any data is read.
Random slopes, structured covariance, and smooth terms make the linear
predictor vary within a cell, which breaks the cell-constant property,
and they are refused with a named reason code rather than silently
approximated. A non-Bernoulli binomial response needs per-row trial
counts that the cell sums cannot recover, so it too falls back to the
per-row path. The guiding principle is that flexyBayes scales by
preserving the model’s meaning, and refuses loudly where exact
aggregation does not hold, rather than hiding an approximation behind
the same interface. For any of these models, fit the data per row with
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).
