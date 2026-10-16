# Fitted values of an aggregated count fit

An aggregated binomial or Poisson fit is one INLA model over cell
totals. Each row's fitted value is the posterior mean of its cell's
inverse-linked linear predictor, random effects included, times the
row's own exposure: `exp(offset)` on a Poisson model, and the `exposure`
or `trials` column of a streamed fit. The values equal those of a
per-row fit of the same model.

## Usage

``` r
# S3 method for class 'flexybayes_aggregated'
fitted(object, ...)
```

## Arguments

- object:

  A `flexybayes_aggregated` fit.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A numeric vector of fitted values on the response scale, one per row of
the data, in row order.

## Details

A streamed fit keeps only the first rows it read. When those are fewer
than the rows it was fitted to, the method refuses by name
(`streamed_rows_not_retained`), not return them as the data.
