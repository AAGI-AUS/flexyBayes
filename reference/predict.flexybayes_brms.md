# Predict from a brms-passthrough flexybayes fit

Delegates to brms's `posterior_epred()` (response-scale posterior mean)
or `posterior_linpred()` (linear-predictor scale) on the live `brmsfit`
carried at `$brms`. The parent `predict.flexybayes()` path uses a
`$glm$linear.predictors` point estimate that handles only the
original-data case; this subclass override accepts `newdata` and returns
the posterior-mean prediction (per-row mean over draws), or the full
posterior matrix when `summary = FALSE`.

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
  classify = NULL,
  level = 0.95,
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

- classify:

  The factors to break a marginal-means table down by: a character value
  (`"Variety"`, `"Variety:env"`) or a one-sided formula (`~ Variety`).
  `NULL` (the default) is the historical behaviour. See
  [`predict.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_inla.md)
  for how the two active engines' prediction paths differ.

- level:

  Credible level for the classify table's interval, as a proportion.
  Default `0.95`.

- ...:

  Forwarded to `brms::posterior_epred()` / `brms::posterior_linpred()`.

## Value

One of four shapes, by argument. With `classify` set, the marginal-means
table for those factors. Otherwise with `summary = FALSE`, the
`draws x rows` posterior matrix; with `se.fit = TRUE`, a list of `fit`
(posterior mean per row) and `se.fit` (posterior SD per row); and by
default a numeric vector of posterior means, one per row of `newdata` or
of the fitted data.

## Details

Population-level vs. group-level prediction follows brms's `re_formula`
convention: the default `re_formula = NULL` includes all random effects;
pass `re_formula = NA` for population-level predictions only.
