# marginaleffects support: population-level predictions (default method)

marginaleffects support: population-level predictions (default method)

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

  Ignored. Present for compatibility with the generic.

## Value

A data frame with `rowid` and `estimate`.
