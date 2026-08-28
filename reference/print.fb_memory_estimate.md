# Print method for an internal `<fb_memory_estimate>` carrier

Renders the per-term INLA memory breakdown introduced at v0.3.10.

## Usage

``` r
# S3 method for class 'fb_memory_estimate'
print(x, ...)
```

## Arguments

- x:

  An `<fb_memory_estimate>` object holding the per-term INLA memory
  breakdown.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the breakdown it prints.
