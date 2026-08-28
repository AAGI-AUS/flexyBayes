# Print method for an internal `<fb_aggregation_plan>` summary

Diagnostic print of the model-level aggregation plan: eligibility,
reason codes, cell-key contributions, estimated cell count, and the
compression ratio.

## Usage

``` r
# S3 method for class 'fb_aggregation_plan'
print(x, ...)
```

## Arguments

- x:

  An `<fb_aggregation_plan>` object, the model-level plan the
  aggregation layer builds.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the eligibility and cell-count
lines it prints.
