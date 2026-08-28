# Summary for flexybayes GLM-compatible object

Similar to [`summary.glm()`](https://rdrr.io/r/stats/summary.glm.html)
but with Bayesian posterior statistics instead of p-values.

## Usage

``` r
# S3 method for class 'flexybayes_glm'
summary(object, ...)
```

## Arguments

- object:

  A `flexybayes_glm` object, reached as `fit$glm` on any fitted
  `flexybayes` object.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, the printed coefficient table as a data.frame, or `NULL` when
the object carries no fixed effects.
