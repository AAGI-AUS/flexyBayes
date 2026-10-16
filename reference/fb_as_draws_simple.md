# Extract per-parameter posterior draws from a model fit

S3 generic that extracts a named list of numeric posterior-draw vectors
from a fit, keyed by the engine's own parameter names. The methods that
reach an active engine are `flexybayes_brms` and `flexybayes_inla`. The
`brmsfit` method reads a brms fit inside a flexyBayes fit, a two-part
non-linear fit from
[`fb_refit()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refit.md)
under the names of the linear model. User-defined methods can extend the
generic.

## Usage

``` r
fb_as_draws_simple(fit, ...)

# S3 method for class 'flexybayes_inla'
fb_as_draws_simple(fit, n_samples = 1000L, ...)

# S3 method for class 'flexybayes_brms'
fb_as_draws_simple(fit, ...)

# S3 method for class 'brmsfit'
fb_as_draws_simple(fit, ...)

# Default S3 method
fb_as_draws_simple(fit, ...)
```

## Arguments

- fit:

  A model fit object carrying a posterior, dispatched on by the methods
  listed above.

- ...:

  Method-specific arguments: `n_samples` for the INLA method. The brms
  method returns every stored draw and refuses `n_samples` with
  `method_argument_not_supported`.

## Value

A named list of numeric vectors, one element per parameter, each element
holding that parameter's posterior draws.

## Details

On an INLA fit the draws come from the fitted approximation: each draw
takes one of INLA's hyperparameter configurations, chosen with its
posterior weight, and draws the latent field from that configuration's
Gaussian approximation, as `INLA::inla.posterior.sample()` does. Every
random number comes from R's generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call makes
the draws repeatable.

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
