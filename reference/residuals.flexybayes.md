# Extract residuals

Extract residuals

## Usage

``` r
# S3 method for class 'flexybayes'
residuals(object, ...)
```

## Arguments

- object:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A numeric vector of response residuals, the observed value minus the
posterior-mean fitted value, one per row of the fitted data. A row whose
response was missing and carried as latent has no observed value and
returns `NA`.
