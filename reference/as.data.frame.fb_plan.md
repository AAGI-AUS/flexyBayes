# Coerce an `<fb_plan>` to data.frame — one row, stable columns

Stable column ordering by the internal vector `.FB_PLAN_DF_COLS`; adding
new fields appends rather than reorders.

## Usage

``` r
# S3 method for class 'fb_plan'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `<fb_plan>` object as returned by
  [`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md).

- row.names:

  Ignored. Present for compatibility with the generic.

- optional:

  Ignored. Present for compatibility with the generic.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

A one-row data.frame whose columns follow the internal
`.FB_PLAN_DF_COLS` order, with one column per recorded plan field.
