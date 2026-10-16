# Residuals of an aggregated count fit

The observed response minus
[`fitted.flexybayes_aggregated()`](https://aagi-aus.github.io/flexyBayes/reference/fitted.flexybayes_aggregated.md),
one value per row. A streamed fit that keeps fewer rows than it was
fitted to refuses by name (`streamed_rows_not_retained`).

## Usage

``` r
# S3 method for class 'flexybayes_aggregated'
residuals(object, ...)
```

## Arguments

- object:

  A `flexybayes_aggregated` fit.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A numeric vector of response residuals, one per row of the data, in row
order.
