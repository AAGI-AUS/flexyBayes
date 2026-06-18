# Fixed-effect coefficients of a per-row INLA fit

Posterior means of the fixed effects, read from the INLA fit's
`summary.fixed` slot (treatment-contrast basis). These are the
coefficients consumed by `emmeans::emmeans()` and
`marginaleffects::predictions()` via the flexyBayes support methods.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
coef(object, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- ...:

  Ignored.

## Value

Named numeric vector of fixed-effect posterior means.
