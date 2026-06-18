# marginaleffects support: set coefficients (greta backend)

marginaleffects support: set coefficients (greta backend)

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

  Ignored.

## Value

The fit with a coefficient override attached.
