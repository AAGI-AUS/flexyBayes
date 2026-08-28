# Extract per-parameter posterior draws from a model fit

S3 generic used by
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
to extract a named list of numeric posterior-draw vectors from each fit.
The methods that reach an active engine are `flexybayes_brms` and
`flexybayes_inla`. User-defined methods can extend the generic.

## Usage

``` r
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

  A model fit object carrying a posterior, dispatched on by the methods
  listed above.

- ...:

  Method-specific arguments, such as `n_samples` for the INLA method.

## Value

A named list of numeric vectors, one element per parameter, each element
holding that parameter's posterior draws.

## Examples

``` r
# A named list of draw vectors, with no posterior dependency.
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
                    verbose = FALSE)
  dr <- fb_as_draws_simple(fit)
  names(dr)
}
# }
```
