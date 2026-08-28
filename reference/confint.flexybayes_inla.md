# Credible intervals for the fixed effects of a per-row INLA fit

Quantiles of INLA's own posterior marginals for the fixed effects, not
frequentist confidence intervals and not a normal approximation to them.
The bounds come from `INLA::inla.qmarginal()` applied to the marginal
densities the fit stores, so any credible level is available rather than
only the 0.95 that `summary.fixed` tabulates.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit carrying fixed-effect marginals.

- parm:

  Character vector of coefficient names to return, or `NULL` (the
  default) for every fixed effect.

- level:

  Credible level for the interval, as a proportion. The default `0.95`
  returns the 2.5th and 97.5th posterior percentiles.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A matrix with one row per fixed effect and two columns holding the lower
and upper credible bounds, named for the percentiles used.

## Details

Before 0.9.0 an INLA fit did not inherit the `flexybayes` parent class,
so [`confint()`](https://rdrr.io/r/stats/confint.html) on one reached
[`stats::confint.default`](https://rdrr.io/r/stats/confint.html) and
failed on a missing `vcov` contract. This method is the INLA sibling of
[`confint.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes.md).
