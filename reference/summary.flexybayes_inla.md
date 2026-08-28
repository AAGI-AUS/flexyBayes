# Summary method for flexybayes_inla

Internal S3 method. Builds and prints the same `summary.flexybayes`
object every active engine returns, so `summary(fit)$varcomp` answers on
an INLA fit as it does on a brms one. Before 0.9.1 this method had its
own dialect – a bare four-slot list of INLA's own tables, not comparable
with what the other engine returned and carrying no variance-component
table at all.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
summary(object, ...)
```

## Arguments

- object:

  A `flexybayes_inla` object, the fit an INLA run returns.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, an object of class `c("summary.flexybayes", "list")`. See
[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md)
for the slots.

## Details

The variance components reach the standard-deviation scale through
INLA's precision marginals rather than by transforming a tabulated point
estimate; see `.inla_variance_comps()`. INLA's own precision-scale
hyperparameter table still prints, unmodified, beneath them.
