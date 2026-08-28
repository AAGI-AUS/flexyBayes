# Plug-in conditional log-likelihood of a flexyBayes fit

Evaluates the conditional log-likelihood at the posterior-mean fitted
values – a plug-in quantity, not a posterior summary, and not the
marginal likelihood. It exists so
[`AIC()`](https://rdrr.io/r/stats/AIC.html)-style comparisons of two
fits from the same engine have something to read.

## Usage

``` r
# S3 method for class 'flexybayes'
logLik(object, ...)
```

## Arguments

- object:

  A `flexybayes` fit carrying a response vector, fitted values, and a
  recorded response family.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A `logLik` object: the scalar log-likelihood with `df` and `nobs`
attributes taken from the fit's recorded model information.

## Details

The method computes what the object carries and refuses by name when it
carries too little. Three requirements are checked before any arithmetic
runs: the response vector in `$glm$y`, the fitted values from
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and a recorded
family. A Gaussian fit uses the residual root-mean-square as its plug-in
scale. Families outside the Gaussian, binomial and Poisson set have no
plug-in form here and are refused rather than reported as `NA`, because
a silent `NA` propagates into
[`anova()`](https://rdrr.io/r/stats/anova.html) and
[`AIC()`](https://rdrr.io/r/stats/AIC.html) as though the comparison had
been made.

brms fits reach
[`logLik.flexybayes_brms()`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes_brms.md)
instead (brms's own `log_lik()`), and INLA fits reach
[`logLik.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes_inla.md),
which states that INLA reports a marginal log-likelihood and returns
`NA` with that message.
