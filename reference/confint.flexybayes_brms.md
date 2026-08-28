# Credible intervals on the brms path

Uses the brms posterior draws directly (the parent `confint.flexybayes`
refuses unconditionally with `fit_lacks_posterior_draws`, since no
active engine reaches it without its own override). Returns
quantile-based credible bounds over the `b_<term>` rows; row names are
stripped of the brms `b_` prefix to align with
[`coef()`](https://rdrr.io/r/stats/coef.html).

## Usage

``` r
# S3 method for class 'flexybayes_brms'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `flexybayes_brms` object.

- parm:

  Subset of fixed-effect names to return (NULL = all).

- level:

  Credible level (default 0.95).

- ...:

  Ignored. Present for compatibility with the generic.

## Value

A numeric matrix with one row per fixed-effect term and two columns
holding the lower and upper credible bounds at `level`. Row names are
the term names with the brms `b_` prefix removed.
