# Print a flexybayes_aggregated object

Brief one-screen summary of an aggregated-gaussian fit produced by
`flexybayes(..., aggregate = "auto"/TRUE)` or
`fb_brms(..., aggregate = ...)`. Includes the `exactness` field and the
cell compression ratio (when N/K \>= 2).

## Usage

``` r
# S3 method for class 'flexybayes_aggregated'
print(x, ...)
```

## Arguments

- x:

  a `<flexybayes_aggregated>` object.

- ...:

  unused.

## Value

invisibly returns `x`.
