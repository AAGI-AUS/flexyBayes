# Predict from a flexybayes model

Predict from a flexybayes model

## Usage

``` r
# S3 method for class 'flexybayes'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  se.fit = FALSE,
  chunk_size = NULL,
  allow_new_levels = c("population", "sample", "refuse"),
  output_file = NULL,
  format = c("auto", "csv", "rds", "fst"),
  interop = FALSE,
  ...
)
```

## Arguments

- object:

  A flexybayes object

- newdata:

  Optional new data frame for prediction. If NULL, returns fitted values
  from the original data.

- type:

  `"link"` for linear predictor, `"response"` for response scale.

- se.fit:

  Logical: return standard errors?

- chunk_size:

  Optional integer. When supplied, and when `newdata` has more rows than
  `chunk_size`, prediction iterates over chunks of this size and
  concatenates the results. Default `NULL` preserves the single-pass
  behaviour. Per-chunk prediction uses the same factor-level dictionary
  so `chunk_size` does not change the numerical output – only the
  wall-time and peak memory profile. A typical setting at large
  `newdata` is `chunk_size = 10000L`.

- allow_new_levels:

  One of `"population"` (default), `"sample"`, or `"refuse"`. Policy for
  handling factor levels in `newdata` that are not in the fit-time
  dictionary. `"population"` sets unknown-level rows to NA on the
  affected column and emits a warning naming the count; downstream
  prediction returns NA for these rows. `"refuse"` raises a structured
  stop on the first unknown level. `"sample"` (active since v0.3.5)
  layers a fresh `Normal(0, tau_<group>)` random-effect realisation onto
  each unknown row per posterior draw – caller's
  [`set.seed()`](https://rdrr.io/r/base/Random.html) controls
  reproducibility. The in-memory return reports the posterior-mean
  prediction (sampled-RE contribution averages toward zero across
  draws); the file-backed return (`output_file = ...`) captures the
  proper per-row posterior interval reflecting sampled-RE uncertainty.
  Only consulted when the fit carries an `extras$fb_dataset` slot (fits
  produced under v0.3.4+); legacy fits skip dictionary resolution
  entirely.

- output_file:

  Optional character path. When supplied, prediction is written to disk
  under the format resolved by `format`. The file is a tabular structure
  with columns `point`, `lower`, `upper` (95\\ columns of `newdata`.
  Refuses if `output_file` already exists (no silent overwrite). Returns
  the path invisibly. Requires a fit with posterior draws on
  `$greta$draws` (greta-backend fits including
  `fb_brms(..., backend = "greta")`); INLA fits do not currently expose
  per-draw access and route to a structured refusal. Default `NULL`
  (in-memory return).

- format:

  One of `"auto"` (default), `"csv"`, `"rds"`, `"fst"`. Only consulted
  when `output_file` is supplied. `"auto"` resolves to: `"csv"` when
  `interop = TRUE`; `"fst"` when `nrow(newdata) >= 1e6` and `fst` is
  installed; `"rds"` otherwise. `"fst"` requested without `fst`
  installed raises a structured refusal naming the install command. The
  fst path is 30–40x faster than rds at N \>= 1e6 per the Stage-3B-shape
  benchmark (`benchmark_results/fst_stage3b_2026-05-23`).

- interop:

  Logical. When `TRUE`, the format-resolution rule under
  `format = "auto"` prefers `"csv"` (universally readable) over rds /
  fst (R-only / fst-only). Useful for handing the prediction grid to a
  non-R consumer. Default `FALSE`. Only consulted when `output_file` is
  supplied and `format = "auto"`.

- ...:

  Additional arguments (ignored)

## Value

If `output_file` is supplied: invisibly returns the path that was
written to. Otherwise: if `se.fit = FALSE`, a numeric vector. If
`se.fit = TRUE`, a list with `fit` and `se.fit`.
