# Predict from a brms-passthrough flexybayes fit

Delegates to brms's `posterior_epred()` (response-scale posterior mean)
or `posterior_linpred()` (linear-predictor scale) on the live `brmsfit`
carried at `$brms`. The parent
[`predict.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes.md)
path uses a `$glm$linear.predictors` point estimate that handles only
the original-data case; this subclass override accepts `newdata` and
returns the posterior-mean prediction (per-row mean over draws), or the
full posterior matrix when `summary = FALSE`.

## Usage

``` r
# S3 method for class 'flexybayes_brms'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  re_formula = NULL,
  se.fit = FALSE,
  summary = TRUE,
  ...
)
```

## Arguments

- object:

  A `flexybayes_brms` object.

- newdata:

  Optional data.frame at which to predict. When omitted, returns the
  in-sample posterior summary.

- type:

  `"response"` (default; posterior_epred) or `"link"`
  (posterior_linpred).

- re_formula:

  Forwarded to brms; `NULL` (default) includes all random effects, `NA`
  excludes them (population-level).

- se.fit:

  Logical: if `TRUE`, returns a list with `fit` (posterior mean) and
  `se.fit` (posterior SD).

- summary:

  Logical: if `TRUE` (default), summarise across draws to a numeric
  vector; if `FALSE`, return the `draws x rows` posterior matrix.

- ...:

  Forwarded to `brms::posterior_epred()` / `brms::posterior_linpred()`.

## Details

Population-level vs. group-level prediction follows brms's `re_formula`
convention: the default `re_formula = NULL` includes all random effects;
pass `re_formula = NA` for population-level predictions only.
