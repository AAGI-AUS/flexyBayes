# marginaleffects support: fixed-effect coefficients (greta backend)

marginaleffects support: fixed-effect coefficients (greta backend)

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

  Ignored.

## Value

Named numeric vector of coefficients.
