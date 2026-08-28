# Approximate leave-one-out cross-validation for a flexyBayes fit

A method for the loo package's `loo()` generic. On a fit from the
**brms** engine the call passes through to `brms::loo()`, which computes
PSIS-LOO from the pointwise log-likelihood Stan stored, and returns
loo's own object – `elpd_loo`, `p_loo`, `looic` and the Pareto-k
diagnostics, unchanged.

## Usage

``` r
# S3 method for class 'flexybayes'
loo(x, ...)
```

## Arguments

- x:

  A fitted `flexybayes` object.

- ...:

  Passed to `brms::loo()` on the brms path; ignored on the refusal
  paths.

## Value

On a brms-engine fit, the `loo` object `brms::loo()` returns. Otherwise
the method does not return: it raises a refusal of class
`flexybayes_refusal_loo_requires_sampler_draws`, or
`flexybayes_refusal_fit_lacks_posterior_draws` when the fit carries no
posterior.

## Details

On a fit from the **INLA** engine the method refuses by name. INLA
returns a nested Laplace approximation to the posterior, not draws of
the log-likelihood, so there is no quantity for importance sampling to
reweight and no leave-one-out estimate to report. The refusal names the
information criteria the fit did compute – WAIC at `fit$inla$waic$waic`
and DIC at `fit$inla$dic$dic` on a per-row fit – so the alternative is
in the message rather than left to be found. Neither is a leave-one-out
estimate and neither carries Pareto-k diagnostics.

A fit that carries no posterior at all is refused as
`fit_lacks_posterior_draws` instead, so the two states are told apart by
the condition class rather than by reading the message.

## See also

[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md),
whose `$converge` slot carries the engine's own diagnostics;
[`as_draws_df.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/as_draws_df.flexybayes.md)
for the posterior itself;
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
for comparing two engines' fits of one model;
[`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md)
for the refusal vocabulary.
