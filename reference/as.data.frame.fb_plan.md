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

  an `<fb_plan>` object.

- row.names:

  unused.

- optional:

  unused.

- ...:

  unused.
