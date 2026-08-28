# Print method for an internal `<fb_aggregated>` summary

Diagnostic print of the sufficient-statistics aggregation: compression
ratio, cell count, total N, and the per-cell key columns.

## Usage

``` r
# S3 method for class 'fb_aggregated'
print(x, ...)
```

## Arguments

- x:

  An `<fb_aggregated>` object, the sufficient-statistics carrier the
  aggregation layer builds.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the compression and cell-key lines
it prints.
