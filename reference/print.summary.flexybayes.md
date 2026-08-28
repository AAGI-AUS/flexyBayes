# Print a flexyBayes summary

Prints the fixed-effect table, the variance components with the prior
each one carried, the engine's own convergence diagnostics, and whatever
engine-native panels the fit earns – INLA's hyperparameter and
spatial-field tables, brms's per-level residual table.

## Usage

``` r
# S3 method for class 'summary.flexybayes'
print(x, digits = 4L, ...)
```

## Arguments

- x:

  A `summary.flexybayes` object, as returned by
  [`summary()`](https://rdrr.io/r/base/summary.html) on a fitted model.

- digits:

  Number of significant digits for the printed tables, including the
  scales inside the `prior` cell.

- ...:

  Ignored, present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the tables it prints.

## Details

Two sentences on the output are deliberate. The banner above the
variance components names the estimator, because the table looks like an
ASReml variance-component table and is not one: every number in it is a
posterior summary. On a fit carrying a latent field, the sentence below
the table names the fourth parameter, because a separable AR1 field plus
an independent nugget is a different model from ASReml's three-parameter
nugget-free residual and the tables look alike.

The numbers inside the `prior` cell are rounded for display only. The
stored cell carries the resolved prior's own string at full precision,
which is what
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
reports and what the two are compared on.
