# Response residuals from a per-row INLA fit

Observed response minus the posterior-mean fitted value (on the response
scale). The response is recovered from the fit's fixed-effect formula
evaluated against the stored data, so a transformed response
(`log(y) ~ ...`) residualises on the modelled scale.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
residuals(object, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- ...:

  Ignored.

## Value

A numeric vector of response residuals, one per observation.
