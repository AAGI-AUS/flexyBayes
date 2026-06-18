# Print method for the brms-passthrough flexybayes subclass

Mirrors `print.flexybayes` (call info + run time + diagnostics) with a
brms-specific footer (the live `brmsfit` lives at `$brms`; the GLM shim
at `$glm`; `$extras` carries the same diagnostics as the greta path).

## Usage

``` r
# S3 method for class 'flexybayes_brms'
print(x, ...)
```

## Arguments

- x:

  A `flexybayes_brms` object.

- ...:

  Ignored.
