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
  ...
)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- newdata:

  Optional data frame; defaults to the fit data.

- type:

  `"response"` or `"link"`.

- se.fit:

  Logical: also return delta-method standard errors from the
  fixed-effect covariance.

- ...:

  Ignored.

## Value

A numeric vector of predictions, or a list `fit` / `se.fit` when
`se.fit = TRUE`.
