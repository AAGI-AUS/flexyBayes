# Posterior draws from a flexyBayes fit

Methods for the posterior package's `as_draws()`, `as_draws_df()` and
`as_draws_matrix()` generics, so a fitted model hands its posterior to
the Bayesian workflow ecosystem – posterior summaries, bayesplot
displays, anything reading a `draws` object – in the shape those tools
expect.

## Usage

``` r
# S3 method for class 'flexybayes'
as_draws_df(x, n_draws = 1000L, ...)

# S3 method for class 'flexybayes'
as_draws(x, n_draws = 1000L, ...)

# S3 method for class 'flexybayes'
as_draws_matrix(x, n_draws = 1000L, ...)
```

## Arguments

- x:

  A fitted `flexybayes` object from either active engine.

- n_draws:

  Number of posterior draws to return, matching the default of
  [`fb_as_draws_simple()`](https://aagi-aus.github.io/flexyBayes/reference/fb_as_draws_simple.md).
  Used on the INLA path, where the draws are sampled from the fitted
  approximation; ignored on the brms path, which returns the draws the
  sampler kept.

- ...:

  Passed to the underlying draws extraction.

## Value

A posterior draws object holding one column per canonical parameter: a
`draws_df` from `as_draws_df()` and from `as_draws()`, a `draws_matrix`
from `as_draws_matrix()`.

## Details

The parameter names are canonical and the same on every engine:
`(Intercept)`, one entry per fixed-effect term, `sigma` for the residual
standard deviation, and `sd_<group>` for the standard deviation of each
random-effect group. Variance components are on the **standard-deviation
scale** whatever the engine stored: INLA's precision hyperparameters are
transformed on the way out, so `sd_g` is a draw of a standard deviation
and not of a precision. This is the same view
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
compares, and it is built by the same internal seam, so a name in one is
a name in the other.

## Where the draws come from

On a **brms** fit the draws are the sampler's own: the kept post-warmup
iterations, renamed. `n_draws` is ignored, because the number of draws
was fixed when the model was sampled.

On an **INLA** fit there is no sampler. The fit is a nested Laplace
approximation, and the draws are sampled *from that fitted
approximation* by `INLA::inla.posterior.sample()`, reading the
configuration store every flexyBayes INLA fit is built with
(`control.compute = list(config = TRUE)`). Two consequences follow, and
they are the ones Tutorials 01 and 10 teach:

- The **fitting** is deterministic. The Laplace approximation draws no
  random numbers, so a `seed` passed to
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  changes nothing about the fit itself.

- The **sampling** step is not, and it is not reproducible from
  [`set.seed()`](https://rdrr.io/r/base/Random.html) alone. The
  hyperparameter draws follow R's random stream, but the latent field –
  the intercept and the fixed-effect terms – is drawn by INLA's own
  generator, which `INLA::inla.posterior.sample()` seeds at random
  unless its own `seed` argument is given, and flexyBayes leaves that
  argument at its default. `?INLA::inla.posterior.sample` states that
  reproducing a sample needs that seed *and* R's RNG state fixed. Treat
  the result as one sample of the posterior rather than a fixed object:
  Monte-Carlo error in anything computed from it shrinks with `n_draws`,
  so raise `n_draws` before reading a small difference between two sets
  of draws as a difference between two posteriors.

## See also

[`fb_as_draws_simple()`](https://aagi-aus.github.io/flexyBayes/reference/fb_as_draws_simple.md)
for the same draws as a plain named list, with no posterior dependency;
[`canonical_names()`](https://aagi-aus.github.io/flexyBayes/reference/canonical_names.md)
for the name map itself;
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
for the priors those draws were taken under;
[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md),
whose `$converge` slot carries the engine's own convergence diagnostics;
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
for a two-engine comparison over the same canonical view.
