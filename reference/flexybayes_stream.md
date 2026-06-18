# Fit a mixed model to an out-of-core dataset by streaming aggregation

`flexybayes_stream()` is the large-data entry point: it reads the data a
chunk at a time from `source`, accumulates additive sufficient
statistics per design cell, and fits the resulting aggregated model
through the same greta or INLA emit path as
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
  backend = c("inla", "greta"),
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
  factors or factor interactions; a continuous term forces one cell per
  row and is refused.

- random:

  A one-sided formula of random-intercept grouping factors, ASReml style
  (for example `~ geno`), or `NULL` for no random terms.

- source:

  The data source. One of: a length-1 character path to an `.fst` file;
  a length-1 character path to a delimited file readable by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html);
  an in-memory `data.frame` / `data.table` (chunked internally, mainly
  for testing); or a function of one argument `i` returning the `i`-th
  chunk as a `data.frame` and `NULL` once the chunks are exhausted.

- family:

  A length-1 character family: `"gaussian"`, `"binomial"`, or
  `"poisson"`.

- trials:

  For `family = "binomial"`, the name of the trials-count column, or
  `NULL` for Bernoulli (one trial per row).

- exposure:

  For `family = "poisson"`, the name of the exposure / offset column, or
  `NULL` for unit exposure.

- backend:

  The estimation backend, `"inla"` (default) or `"greta"`.

- chunk_rows:

  The number of rows to read per chunk. Larger chunks read faster but
  use more peak memory; the default 5e6 keeps peak memory modest while
  amortising read overhead.

- prior:

  An optional
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  object. When `NULL` the backend default priors are used.

- fit:

  When `TRUE` (default) the aggregated model is fitted and a
  `<flexybayes>` object is returned. When `FALSE` the function returns
  the `<fb_aggregated>` sufficient-statistics object without fitting –
  useful for inspecting the compression a design will achieve.

- verbose:

  When `TRUE`, report streaming progress and the compression achieved.

- ...:

  Further arguments passed to the aggregated emit (for example
  `n_samples`, `warmup`, `chains` for the greta backend).

## Value

When `fit = TRUE`, a `<flexybayes>` fit object carrying
`extras$aggregation_meta$streamed == TRUE`. When `fit = FALSE`, an
`<fb_aggregated>` object.

## Details

Fit a cell-aggregatable mixed model to a dataset held on disk, reading
it in row-range chunks and accumulating exact per-cell sufficient
statistics rather than materialising the full table in memory.

The method scales by cell count `K` (the number of distinct factor
design combinations), not by row count `N`. A replicated factorial
design collapses billions of rows to a few thousand cells and fits in
seconds; a design that keys on a continuous covariate does not compress
and is refused before any chunk is read. Continuous fixed effects,
random slopes, structured covariance, and smooth terms break the
cell-constant linear-predictor property and are likewise refused – pass
the data to
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
per row for those models.

Supported families are `"gaussian"` (identity link), `"binomial"`, and
`"poisson"`. For binomial data supply `trials` to name the column of
trial counts (a 0/1 response is treated as Bernoulli with one trial per
row). For poisson data supply `exposure` to name an exposure / offset
column.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the in-memory entry point;
[`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
for prior specification.

## Author

Max Moldovan, <max.moldovan@adelaide.edu.au>

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("fst", quietly = TRUE)
set.seed(1L)
n <- 1e6
df <- data.frame(
  env = factor(sample(letters[1:6], n, replace = TRUE)),
  geno = factor(sample(1:50, n, replace = TRUE)),
  y = rnorm(n)
)
path <- tempfile(fileext = ".fst")
fst::write_fst(df, path)
fit <- flexybayes_stream(y ~ env, random = ~ geno, source = path,
                         backend = "inla")
summary(fit)
}
```
