# Summarise a fit run on the aggregated representation

Returns the same eleven-slot `summary.flexybayes` object every other
engine returns, so `summary(fit)$varcomp` answers on an aggregated fit
as it does on a per-row one. Aggregation is the default route for an
ordinary Gaussian call with a random term, so until 0.9.1 that slot was
`NULL` on the commonest fit the package produces: this method returned
its own list of the aggregated posterior's raw pieces and carried no
variance-component table at all.

## Usage

``` r
# S3 method for class 'flexybayes_aggregated'
summary(object, ...)
```

## Arguments

- object:

  A `<flexybayes_aggregated>` object, as returned by a fit run on the
  aggregated representation.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, an object of class `c("summary.flexybayes", "list")` with the
slots documented at
[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md).
Printed as a side effect.

## Details

Two things are aggregated-specific and both appear in the printed header
rather than as extra slots. The banner names the representation
(`aggregated-gaussian`, `aggregated-binomial`, `aggregated-poisson`)
where a per-row fit names its engine, and the `aggregation` line gives
the row-to-cell compression the fit ran under. The `model` slot carries
the same compression in its text, because a reader comparing the
fixed-effect tables of an aggregated and a per-row fit otherwise has
nothing on the object saying the two were computed over different
numbers of rows.

The variance components come from the engine's own posterior. On the
INLA route that is the hyperparameter marginal transformed to the
standard-deviation scale, the same construction the per-row route uses,
and not the reciprocal square root of a tabulated precision mean. A
route that recorded posterior means alone reports the means with the
three interval columns `NA` and the row's `note` cell reading
`"no interval recorded"`.

The raw aggregated pieces this method used to return – `beta_means`,
`beta_vcov`, `sigma_means`, `tau_means` – are unchanged at
`fit$extras$summary`.

## See also

[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md)
for the slot-by-slot contract, and
[`print.flexybayes_aggregated()`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes_aggregated.md)
for the one-screen fit description.
