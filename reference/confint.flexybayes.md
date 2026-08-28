# Credible intervals for the fixed effects of a flexyBayes fit

Returns posterior quantile-based credible intervals, not frequentist
confidence intervals: the bounds would be empirical quantiles of the
fixed-effect posterior draws the fit carries.

## Usage

``` r
# S3 method for class 'flexybayes'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A flexybayes fit.

- parm:

  Character vector of parameter names to return, or `NULL` (the default)
  for every fixed effect.

- level:

  Credible level for the interval, as a proportion. The default `0.95`
  returns the 2.5th and 97.5th posterior percentiles.

- ...:

  Ignored, present for compatibility with the generic.

## Value

Does not return: raises the classed `fit_lacks_posterior_draws` refusal.

## Details

No active engine stores its posterior in the shape this bare fallback
method reads. A fit from an active engine reaches its own method –
[`confint.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes_inla.md)
for INLA,
[`confint.flexybayes_brms()`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes_brms.md)
for brms – so this method always refuses by name rather than returning
an empty interval matrix, which would read as "no fixed effects" rather
than "this fit cannot answer the question".
