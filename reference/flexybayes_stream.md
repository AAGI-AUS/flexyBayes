# Fit a mixed model to an out-of-core dataset by streaming aggregation

`flexybayes_stream()` is the large-data entry point: it reads the data a
chunk at a time from `source`, accumulates additive sufficient
statistics per design cell, and fits the resulting aggregated model
through the same INLA emit path as
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).
The aggregated likelihood is algebraically identical to the per-row
likelihood, so the posterior is the full-data posterior, not an
approximation.

## Usage

``` r
flexybayes_stream(
  fixed,
  random = NULL,
  source,
  family = "gaussian",
  trials = NULL,
  exposure = NULL,
  backend = "inla",
  chunk_rows = 5e+06,
  prior = NULL,
  fit = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- fixed:

  A two-sided formula `response ~ fixed_terms`. The fixed terms must be
  factors or factor interactions. A continuous term forces one cell per
  row and is refused. [`offset()`](https://rdrr.io/r/stats/offset.html)
  terms are allowed on the Gaussian and Poisson families (see Details).

- random:

  A one-sided formula of random-intercept grouping factors, ASReml style
  (for example `~ geno`), or `NULL` for no random terms.

- source:

  The data source. One of: a length-1 character path to an `.fst` file,
  a length-1 character path to a delimited file readable by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html),
  an in-memory `data.frame` / `data.table` (chunked internally, mainly
  for testing), or a function of one argument `i` returning the `i`-th
  chunk as a `data.frame` and `NULL` once the chunks are exhausted.

- family:

  A length-1 character family: `"gaussian"`, `"binomial"`, or
  `"poisson"`.

- trials:

  For `family = "binomial"`, the name of the trials-count column, or
  `NULL` for Bernoulli (one trial per row).

- exposure:

  For `family = "poisson"`, the name of the exposure column, or `NULL`
  for unit exposure. An `offset(log(E))` in `fixed` is the same model as
  `exposure = "E"`.

- backend:

  The estimation backend. `"inla"` is the default and the only active
  choice: INLA is the one engine with an aggregated emit. Any other
  value, including a formerly-registered engine name this package
  withdrew entirely, refuses with `unknown_backend`.

- chunk_rows:

  The number of rows to read per chunk. Larger chunks read faster but
  use more peak memory. The default 5e6 keeps peak memory modest while
  amortising read overhead.

- prior:

  An optional
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  object. When `NULL`, each random-effect standard deviation takes the
  default
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  applies, `half_normal(scale = 3)` for Poisson and
  `half_normal(scale = 5)` for binomial, on the link scale. The fixed
  effects keep INLA's default priors.
  [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
  lists the prior on every parameter of the fit.

- fit:

  When `TRUE` (default) the aggregated model is fitted and a
  `<flexybayes>` object is returned. When `FALSE` the function returns
  the `<fb_aggregated>` sufficient-statistics object without fitting,
  which shows the compression a design will achieve.

- verbose:

  When `TRUE`, report streaming progress and the compression achieved.

- ...:

  Further arguments passed to the aggregated emit.

## Value

When `fit = TRUE`, a `<flexybayes>` fit object carrying
`extras$aggregation_meta$streamed == TRUE`. When `fit = FALSE`, an
`<fb_aggregated>` object.

## Details

Fit a cell-aggregatable mixed model to a dataset held on disk, reading
it in row-range chunks and accumulating exact per-cell sufficient
statistics instead of materialising the full table in memory.

The method scales by cell count `K` (the number of distinct factor
design combinations), not by row count `N`. A replicated factorial
design collapses billions of rows to a few thousand cells and fits in
seconds. A design that keys on a continuous covariate does not compress
and is refused before any chunk is read. Continuous fixed effects,
random slopes, structured covariance, and smooth terms break the
cell-constant linear-predictor property and are likewise refused. Pass
the data to
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
per row for those models. The model needs at least one factor cell-key
term (a fixed factor, a factor interaction, or a random grouping factor
named in `random`) to aggregate by. An intercept-only model
(`fixed = y ~ 1`, no `random`) has nothing to compress and is refused by
name (`aggregation_route_unavailable`).

Supported families for a fit are `"binomial"` and `"poisson"`. A
`"gaussian"` fit refuses by name
(`aggregation_residual_prior_not_composed`): the aggregated Gaussian
route cannot yet carry the residual prior, so fit such data per row with
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).
`fit = FALSE` still returns the Gaussian sufficient statistics. For
binomial data supply `trials` to name the column of trial counts (a 0/1
response is treated as Bernoulli with one trial per row). For poisson
data supply `exposure` to name an exposure column. A fit on the
aggregated representation reports no marginal likelihood
(`aggregated_evidence_unavailable`).

An [`offset()`](https://rdrr.io/r/stats/offset.html) in `fixed` is
carried exactly on two of the three families. On `"gaussian"` it is
subtracted from the response before the cell statistics are taken, which
reproduces the per-row likelihood however the offset varies inside a
cell. On `"poisson"`, `offset(o)` multiplies each row's exposure by
`exp(o)`, so `offset(log(E))` and `exposure = "E"` fit the same model.
Given both, they multiply. On `"binomial"` a cell total of rows whose
offsets differ is not binomial, and the call refuses with
`fixed_offset_unsupported`.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the in-memory entry point, and
[`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
for prior specification.

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("fst", quietly = TRUE)
# INLA is distributed from its own repository, not from CRAN.
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1L)
  n <- 1e6
  df <- data.frame(
    env = factor(sample(letters[1:6], n, replace = TRUE)),
    geno = factor(sample(1:50, n, replace = TRUE)),
    count = rpois(n, 3)
  )
  path <- tempfile(fileext = ".fst")
  fst::write_fst(df, path)
  fit <- flexybayes_stream(count ~ env, random = ~ geno, source = path,
                           family = "poisson", backend = "inla")
  summary(fit)
}
}
```
