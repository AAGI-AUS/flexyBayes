# Log-likelihood of a per-row INLA fit (not computed)

INLA reports a *marginal* log-likelihood (the model evidence, available
through [`summary()`](https://rdrr.io/r/base/summary.html)), not the
*conditional* model log-likelihood that the
[`logLik()`](https://rdrr.io/r/stats/logLik.html) generic denotes and
that the greta / brms backends expose. Returning the marginal quantity
under the `logLik` name would conflate two different things, so this
method honestly returns `NA` (with the degrees of freedom and
observation count filled in) and a one-line note. This also lets
downstream summaries (for example
[`glance()`](https://generics.r-lib.org/reference/glance.html)) degrade
gracefully instead of erroring with "no applicable method".

## Usage

``` r
# S3 method for class 'flexybayes_inla'
logLik(object, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- ...:

  Ignored.

## Value

A `logLik` object whose value is `NA_real_`, carrying `df` and `nobs`
attributes.
