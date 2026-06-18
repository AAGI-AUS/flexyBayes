# Tidy a per-row INLA fit into a one-row-per-term data frame

The INLA backend returns a `flexybayes_inla` object that does not
inherit from `flexybayes`, so it needs its own
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) method. The
fixed-effect summary is read directly off INLA's `summary.fixed` table,
whose `mean`, `sd`, and `0.025quant` / `0.975quant` columns map cleanly
onto the `broom`-canonical `estimate`, `std.error`, `conf.low`, and
`conf.high`.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A `flexybayes_inla` fit.

- conf.int:

  Logical. Whether to attach the credible-interval columns. Defaults to
  `TRUE`.

- conf.level:

  Numeric in `(0, 1)`. Accepted for generic compatibility; INLA reports
  the 95% marginal bounds, so a non-0.95 request is noted.

- ...:

  Currently unused; present for generic compatibility.

## Value

A `data.frame` with one row per fixed-effect term and the columns
`term`, `estimate`, `std.error`, and (when `conf.int = TRUE`) `conf.low`
/ `conf.high`.

## Details

Because the INLA fixed-effect intervals come from the marginal
posteriors INLA has already integrated, the `conf.level` argument is
accepted for generic compatibility but only the 95% bounds INLA reports
are returned; a one-off message notes this when a different level is
requested rather than silently ignoring it.

## See also

[`tidy.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes.md)

## Examples

``` r
if (FALSE) { # requireNamespace("INLA", quietly = TRUE)
if (FALSE) { # \dontrun{
fit <- flexybayes(yield ~ env, data = dat, backend = "inla")
tidy(fit)
} # }
}
```
