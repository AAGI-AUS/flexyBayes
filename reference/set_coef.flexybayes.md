# marginaleffects support: set coefficients (default method)

marginaleffects support: set coefficients (default method)

## Usage

``` r
# S3 method for class 'flexybayes'
set_coef(model, coefs, ...)

# S3 method for class 'flexybayes_inla'
set_coef(model, coefs, ...)
```

## Arguments

- model:

  A `flexybayes_inla` fit.

- coefs:

  Replacement coefficient vector.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

The fit with a coefficient override attached.
