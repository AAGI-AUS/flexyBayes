# Coerce an `<fb_memory_estimate>` carrier to a numeric scalar

Returns `x$total` (bytes, after the INLA overhead multiplier).
Registered against the `as.double` S3 generic so that
`as.numeric(<fb_memory_estimate>)` dispatches correctly — the base-R
`as.numeric` is `.Primitive("as.double")`, so dispatch fires on
`as.double`, not on `as.numeric`.

## Usage

``` r
# S3 method for class 'fb_memory_estimate'
as.double(x, ...)
```

## Arguments

- x:

  An `<fb_memory_estimate>` carrier holding the per-term breakdown and
  its `total`.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

A numeric scalar, the estimated total in bytes, or `NA_real_` when the
carrier holds no total.
