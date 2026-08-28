# Print a flexybayes_aggregated object

Brief one-screen summary of a fit produced by
`flexybayes(..., aggregate = "auto"/TRUE)` or
`fb_brms(..., aggregate = ...)`. The header names the fit's own family
(`aggregated-gaussian`, `aggregated-binomial`, `aggregated-poisson`).
Includes the `exactness` field and the cell compression ratio (when N/K
\>= 2).

## Usage

``` r
# S3 method for class 'flexybayes_aggregated'
print(x, ...)
```

## Arguments

- x:

  A `<flexybayes_aggregated>` object, as returned by a fit run on the
  aggregated representation.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the one-screen summary it prints.
