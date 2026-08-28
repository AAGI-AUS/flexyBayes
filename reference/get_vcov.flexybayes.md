# marginaleffects support: covariance (default method)

marginaleffects support: covariance (default method)

## Usage

``` r
# S3 method for class 'flexybayes'
get_vcov(model, ...)

# S3 method for class 'flexybayes_inla'
get_vcov(model, ...)
```

## Arguments

- model:

  A `flexybayes_inla` fit.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Fixed-effect covariance matrix.
