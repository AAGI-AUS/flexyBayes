# Print method for the brms-passthrough flexybayes subclass

Opens with the header every engine's print shares, then adds the sampler
diagnostics and a brms-specific footer (the live `brmsfit` lives at
`$brms`; the GLM shim at `$glm`; `$extras` carries the same diagnostics
as the INLA path).

## Usage

``` r
# S3 method for class 'flexybayes_brms'
print(x, ...)
```

## Arguments

- x:

  A `flexybayes_brms` object.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the description it prints.
