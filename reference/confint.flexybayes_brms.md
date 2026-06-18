# Credible intervals on the brms path

Uses the brms posterior draws directly (the parent `confint.flexybayes`
reads `$greta$draws`, which is `NULL` on the brms-passthrough path).
Returns quantile-based credible bounds over the b\_ rows; row names are
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

  Ignored.
