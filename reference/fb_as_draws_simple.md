# Extract per-parameter posterior draws from a model fit

S3 generic used by
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
to extract a named list of numeric posterior-draw vectors from each fit.
Methods exist for the `flexybayes` (greta backend) and `flexybayes_inla`
(INLA backend) classes; user-defined methods can extend the generic.

## Usage

``` r
fb_as_draws_simple(fit, ...)

# S3 method for class 'flexybayes'
fb_as_draws_simple(fit, ...)

# S3 method for class 'flexybayes_inla'
fb_as_draws_simple(fit, n_samples = 1000L, ...)

# S3 method for class 'flexybayes_brms'
fb_as_draws_simple(fit, ...)

# Default S3 method
fb_as_draws_simple(fit, ...)
```

## Arguments

- fit:

  a model fit object.

- ...:

  method-specific arguments (e.g., `n_samples` for INLA).

## Value

a named list of numeric vectors.
