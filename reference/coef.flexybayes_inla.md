# Coefficients of a per-row INLA fit

Posterior means of the fixed effects, read from the INLA fit's
`summary.fixed` slot (treatment-contrast basis). These are the
coefficients consumed by `emmeans::emmeans()` and
`marginaleffects::predictions()` via the flexyBayes support methods.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
coef(object, what = c("fixed", "random", "missing", "all"), ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- what:

  Which part of the fit to return: `"fixed"` (the default), `"random"`,
  `"missing"` or `"all"`. See
  [`coef.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes.md)
  for the shape of each.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

For the default `what = "fixed"`, a named numeric vector of fixed-effect
posterior means; otherwise as documented for
[`coef.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes.md).

## Details

An INLA fit carries no `$glm` slot, so the fixed vector is read here
rather than inherited; every other value of `what` is resolved by the
shared body
[`coef.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes.md)
uses, so the two engines answer the same question the same way.

A factor with a non-syntactic level (for example a level containing a
space, `"low N"`) is legalised internally
([`make.names()`](https://rdrr.io/r/base/make.names.html)) before the
INLA emit, so the fit itself never dies on it; this method restores the
user's own, unlegalised label in the names it returns, as do
[`ranef()`](https://aagi-aus.github.io/flexyBayes/reference/ranef.md),
[`summary()`](https://rdrr.io/r/base/summary.html) and
`predict(classify = )` on the same fit.
