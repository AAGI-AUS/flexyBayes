# Predict from a brms-passthrough flexybayes fit

Delegates to brms's `posterior_epred()` (response-scale posterior mean)
or `posterior_linpred()` (linear-predictor scale) on the live `brmsfit`
carried at `$brms`. The parent `predict.flexybayes()` path uses a
`$glm$linear.predictors` point estimate that handles only the
original-data case. This subclass override accepts `newdata` and returns
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
  allow_new_levels = FALSE,
  mean_by = NULL,
  n_draws = 1000L,
  ...
)
```

## Arguments

- object:

  A `flexybayes_brms` object.

- newdata:

  Optional data.frame at which to predict. When omitted, returns the
  in-sample posterior summary. On a model with an
  [`offset()`](https://rdrr.io/r/stats/offset.html) brms evaluates the
  offset on `newdata`, so `newdata` needs every column the offset reads.
  An absent one refuses with `newdata_offset_missing` instead of being
  taken as zero.

- type:

  `"response"` (default. posterior_epred) or `"link"`
  (posterior_linpred).

- re_formula:

  Forwarded to brms; `NULL` (default) includes all random effects, `NA`
  excludes them (population-level).

- se.fit:

  Logical: if `TRUE`, returns a list with `fit` (posterior mean) and
  `se.fit` (posterior SD).

- summary:

  Logical: if `TRUE` (default), summarise across draws to a numeric
  vector. If `FALSE`, return the `draws x rows` posterior matrix.

- classify:

  The factors to break a predicted-value table down by: a character
  value (`"Variety"`, `"Variety:env"`) or a one-sided formula
  (`~ Variety`). The factors may be fixed or random, and the table
  follows ASReml's default rules, described under
  [`predict.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_inla.md).
  On this engine the values come from the posterior draws. A level
  combination the data never held is drawn from its fitted normal
  distribution. `NULL` (the default) is the historical behaviour.

- level:

  Credible level for the classify table's interval, or for the
  predictive interval with `allow_new_levels = TRUE`, as a proportion.
  Default `0.95`.

- allow_new_levels:

  Logical. `TRUE` returns the posterior predictive of the rows of
  `newdata`, observation noise included, in which a level the fit did
  not see is drawn from its fitted standard deviation:
  `brms::posterior_predict(allow_new_levels = TRUE, sample_new_levels = "gaussian")`.
  The table matches the INLA method's. See
  [`predict.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_inla.md)
  for the structures it is built for and the refusals. Its `estimate`,
  `sd`, `lower` and `upper` are the mean, standard deviation and
  quantiles of the draws. `FALSE` (the default) leaves the paths above
  unchanged.

- mean_by:

  With `allow_new_levels = TRUE`, the columns of `newdata` whose
  combinations define the sets of rows to average, as a character vector
  or a one-sided formula (`~ TRIAL`). `NULL` predicts each row.

- n_draws:

  With `allow_new_levels = TRUE`, the number of posterior draws used, at
  most the number the fit holds. Default `1000`.

- ...:

  Forwarded to `brms::posterior_epred()` / `brms::posterior_linpred()`.

## Value

One of five shapes, by argument. With `classify` set, the marginal-means
table for those factors. With `allow_new_levels = TRUE`, the
`fb_predict_new` table described in
[`predict.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_inla.md).
Otherwise with `summary = FALSE`, the `draws x rows` posterior matrix.
With `se.fit = TRUE`, a list of `fit` (posterior mean per row) and
`se.fit` (posterior SD per row). And by default a numeric vector of
posterior means, one per row of `newdata` or of the fitted data.

## Details

Population-level vs. group-level prediction follows brms's `re_formula`
convention: the default `re_formula = NULL` includes all random effects.
Pass `re_formula = NA` for population-level predictions only.
