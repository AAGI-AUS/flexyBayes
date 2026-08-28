# marginaleffects support: fixed-effect coefficients (default method)

marginaleffects support: fixed-effect coefficients (default method)

## Usage

``` r
# S3 method for class 'flexybayes'
get_coef(model, ...)

# S3 method for class 'flexybayes_inla'
get_coef(model, ...)
```

## Arguments

- model:

  A `flexybayes_inla` fit.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Named numeric vector of coefficients.
