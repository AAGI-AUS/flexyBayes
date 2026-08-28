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
  ...
)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- newdata:

  Optional data frame; defaults to the fit data. Ignored, with a
  warning, when `classify` is supplied.

- type:

  `"response"` or `"link"`.

- se.fit:

  Logical: also return delta-method standard errors from the
  fixed-effect covariance.

- classify:

  The factors to break a marginal-means table down by: a character value
  (`"Variety"`, `"Variety:env"`) or a one-sided formula (`~ Variety`).
  `NULL` (the default) is the historical behaviour.

- level:

  Credible level for the classify table's interval, as a proportion.
  Default `0.95`.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

With `classify`, a data frame of class `fb_predict_classify`. Otherwise
a numeric vector of predictions, or a list `fit` / `se.fit` when
`se.fit = TRUE`.

## Details

The `classify` path builds a marginal-means table through the emmeans
seam (the same construction
[`predict.flexybayes_brms()`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_brms.md)
uses), whose interval on this engine comes from the Gaussian
approximation of the joint fixed-effect posterior rather than from
INLA's own marginals. The printed table names that.
