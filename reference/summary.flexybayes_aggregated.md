# Summarise a flexybayes_aggregated object

Posterior summary read off the aggregated INLA fit's `summary.fixed` +
`summary.hyperpar` slots. Shows the compression line when N/K \>= 2.

## Usage

``` r
# S3 method for class 'flexybayes_aggregated'
summary(object, ...)
```

## Arguments

- object:

  a `<flexybayes_aggregated>` object.

- ...:

  unused.

## Value

invisibly returns the posterior summary list.
