# Posterior predictive check for a flexyBayes fit

A method for the bayesplot package's `pp_check()` generic. On a fit from
the **brms** engine the call delegates to `brms::pp_check()`: datasets
are simulated from the posterior predictive distribution and shown
against the observed response. The default display overlays the
densities of the replicated datasets on the density of the response
(`type = "dens_overlay"`); every other bayesplot check type, and every
argument of `brms::pp_check()`, is reachable through `...`.

## Usage

``` r
# S3 method for class 'flexybayes'
pp_check(object, ...)
```

## Arguments

- object:

  A fitted `flexybayes` object.

- ...:

  Passed to `brms::pp_check()` on the brms path – `type`, `ndraws`,
  `group` and the rest; ignored on the refusal path.

## Value

On a brms-engine fit, the ggplot2 object `brms::pp_check()` returns.
Otherwise the method does not return: it raises a refusal of class
`flexybayes_refusal_pp_check_requires_predictive_draws`.

## Details

On a fit from the **INLA** engine the method refuses by name. A nested
Laplace approximation returns marginal densities, not simulated
datasets, so there is nothing to overlay – and a display of observed
against fitted values shown under this name would be a different
diagnostic wearing the check's title. The refusal names the residual
displays the fit does answer: `plot(fit, type = "residuals")` and, on a
fit carrying a design index, `plot(fit, type = "variogram")`.

`plot(fit, type = "pp_check")` runs this same code, so the two entry
points cannot disagree about what a check is.

## See also

[`plot.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
for the other displays, including the residual diagnostics this refusal
points at;
[`as_draws_df.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/as_draws_df.flexybayes.md)
for the posterior itself;
[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md),
whose `$converge` slot carries the engine's own diagnostics.
