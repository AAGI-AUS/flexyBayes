# marginaleffects support: population-level predictions (greta backend)

marginaleffects support: population-level predictions (greta backend)

## Usage

``` r
# S3 method for class 'flexybayes'
get_predict(model, newdata = NULL, type = "response", ...)

# S3 method for class 'flexybayes_inla'
get_predict(model, newdata = NULL, type = "response", ...)
```

## Arguments

- model:

  A `flexybayes_inla` fit.

- newdata:

  Data frame to predict on (default: fit data).

- type:

  Prediction scale (identity link only).

- ...:

  Ignored.

## Value

A data frame with `rowid` and `estimate`.
