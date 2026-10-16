# Population-level predictions from a per-row INLA fit

Fixed-effect (population-level) predictions: the linear predictor is
`X beta` with random effects held at their population mean (zero). On
the identity link the response- and link-scale predictions coincide.
This is the prediction surface marginaleffects uses for average
predictions and slopes.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  se.fit = FALSE,
  classify = NULL,
  level = 0.95,
  allow_new_levels = FALSE,
  mean_by = NULL,
  n_draws = 1000L,
  re_formula = NA,
  summary = TRUE,
  ...
)

# S3 method for class 'fb_predict_new'
print(x, digits = 4L, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- newdata:

  Optional data frame. Defaults to the fit data. Ignored, with a
  warning, when `classify` is supplied. On a model with an
  [`offset()`](https://rdrr.io/r/stats/offset.html) the prediction
  includes the offset evaluated on `newdata`, so `newdata` needs every
  column the offset reads. An absent one refuses with
  `newdata_offset_missing` instead of being taken as zero.

- type:

  `"response"` or `"link"`.

- se.fit:

  Logical: also return delta-method standard errors from the
  fixed-effect covariance.

- classify:

  The factors to break a predicted-value table down by: a character
  value (`"Variety"`, `"Variety:env"`) or a one-sided formula
  (`~ Variety`). The factors may be fixed or random. `NULL` (the
  default) is the historical behaviour.

- level:

  Credible level for the classify table's interval, or for the
  predictive interval with `allow_new_levels = TRUE`, as a proportion.
  Default `0.95`.

- allow_new_levels:

  Logical. `FALSE` (the default) is the population-level prediction
  above. `TRUE` returns the posterior predictive of the rows of
  `newdata` (see Details), in which every random effect enters: a level
  the fit saw at its posterior draw, and a level it did not see drawn
  from its fitted variance component.

- mean_by:

  With `allow_new_levels = TRUE`, the columns of `newdata` whose
  combinations define the sets of rows to average: a character vector or
  a one-sided formula. `mean_by = ~ TRIAL` predicts each trial mean,
  `~ TRIAL + VARIETY + NRATE` each cell mean. `NULL` (the default)
  predicts each row.

- n_draws:

  With `allow_new_levels = TRUE`, the number of posterior draws. Default
  `1000`.

- re_formula:

  `NA` (the default): the prediction is population-level, the fixed part
  of the model, as `re_formula = NA` gives on a brms fit. Any other
  value refuses with `method_argument_not_supported`.
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
  predictions at the fitted levels of the random effects.

- summary:

  `TRUE` (the default): the prediction is the posterior mean. `FALSE`
  refuses with `method_argument_not_supported`.
  [`fb_as_draws_simple()`](https://aagi-aus.github.io/flexyBayes/reference/fb_as_draws_simple.md)
  returns the posterior draws.

- ...:

  Ignored, except that an argument of ASReml's
  [`predict()`](https://rdrr.io/r/stats/predict.html) given with
  `classify` is refused by name.

- x:

  An `fb_predict_new` table, printed with a banner naming what its
  numbers are.

- digits:

  Number of significant digits for the printed table.

## Value

With `classify`, a data frame of class `fb_predict_classify`. With
`allow_new_levels = TRUE`, a data frame of class `fb_predict_new`: one
row per set, holding the `mean_by` columns (or `row`), `n_rows`,
`estimate`, `sd`, `lower` and `upper` (and `status` on a fit with
`fixed_alias = "drop"`, where a set holding a row the fit does not
estimate is `NA`), with the joint predictive draws, draws by sets, in
`attr(, "draws")`. Otherwise a numeric vector of predictions, or a list
`fit` / `se.fit` when `se.fit = TRUE`.

## Details

With `gen` random, `classify = "gen"` returns each genotype's fixed part
plus its random effect, and on a multi-environment fit with `gen:env`
random it also adds that genotype's `gen:env` effects averaged over
environments. This matches ASReml's default
[`predict()`](https://rdrr.io/r/stats/predict.html): equal-weight
averaging over the fixed-term variables not in `classify`, a covariate
held at its mean, variables that appear only in random terms ignored,
and every random term whose variables are all in `classify` or averaged
over included. The intervals integrate over INLA's hyperparameter
configurations, so they include the uncertainty in the variance
components. A table with no random term is computed by emmeans, with
intervals from the Gaussian approximation of the fixed-effect posterior.
The arguments of ASReml's
[`predict()`](https://rdrr.io/r/stats/predict.html) that change these
rules (`present`, `average`, `ignore`, `use`, `except`, `only`,
`associate`, `levels`, `vcov`, `sed`, `parallel`, `aliased`,
`design.points`) are refused with `classify_argument_not_supported`.

With `allow_new_levels = TRUE`, the prediction for a set of rows is the
mean response over them, observation noise included, as an observed
trial or cell mean at a new locality would be. Given one posterior draw
it is normal, with mean `x'beta` averaged over the rows plus the mean of
the effects at levels the fit saw, and variance the sum over random
terms of `sigma^2 * sum_l (n_l / n)^2` over the unseen levels `l` (`n_l`
of the `n` rows at level `l`), plus `sigma_field^2` times the mean AR1 x
AR1 correlation over pairs of rows in a trial the fit did not see, plus
the residual variance over `n`. The prediction is the mixture of these
normals over the draws: `estimate` and `sd` are its exact mean and
standard deviation, and `lower` and `upper` its quantiles. The fixed
effects and the seen effects are drawn jointly from INLA's
approximation, and the variances from their continuous joint posterior
(as
[`fb_variance_draws()`](https://aagi-aus.github.io/flexyBayes/reference/fb_variance_draws.md)
draws them). Every random number comes from R's generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call makes
the result repeatable.

The new-level draw is built for random intercepts, interactions
(`TRIAL:REP`), `diag(f):g` and `id(trial):ar1(row):ar1(col)`, over an
independent residual or `dsum(~ units | f)`, for a Gaussian response on
the identity link. Any other structure is refused with
`predict_new_levels_structure_not_supported`, a new level with no fitted
variance (a new level of `f` in `diag(f):g`, a new residual section, a
row or column off the fitted lattice) with
`predict_new_levels_no_variance`, and another family with
`predict_new_levels_family_not_supported`.

On a fit with `fixed_alias = "drop"`, a row that is not a linear
combination of the rows the model was fitted to has no estimable
prediction: it is returned as `NA`, and the `status` attribute (a
`status` element with `se.fit = TRUE`, a `status` column with
`classify`) reads `"Aliased"` there and `"Estimable"` elsewhere, as in
asreml's [`predict()`](https://rdrr.io/r/stats/predict.html).
