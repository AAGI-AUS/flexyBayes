# Advance a deferred-execution object into its fit

Generic for the inspect-then-fit pattern. The `<flexybayes_review>`
method restores the RNG snapshot captured at review-object construction,
runs the deferred fit via the backend driver, caches the result
in-place, and returns the fit. A second call returns the cached fit.

## Usage

``` r
proceed(x, ...)

# S3 method for class 'flexybayes_review'
proceed(x, ...)
```

## Arguments

- x:

  A `<flexybayes_review>` object holding the deferred fit.

- ...:

  Reserved for future deferred-execution classes.

## Value

The fit object the originating call would have returned, of class
`flexybayes`. A second call returns the cached fit.

## Deprecated

`proceed()` is deprecated in flexyBayes 0.13.0 and will be removed in a
later release. Call
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
with the same arguments and `review_code = FALSE`, and set `seed` for a
reproducible brms fit. A call warns with a condition of class
`flexybayes_deprecated`.

## Examples

``` r
# The direct call fits the model proceed() would have fitted. The Stan
# compile is the expensive part of running this, so the fit is given a
# small budget.
if (requireNamespace("brms", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  # \donttest{
    fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "brms",
                      chains = 1L, n_samples = 200L, warmup = 100L,
                      seed = 1L)
    class(fit)
  # }
}
```
