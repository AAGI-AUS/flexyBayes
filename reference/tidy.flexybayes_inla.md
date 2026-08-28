# Tidy a per-row INLA fit into a one-row-per-term data frame

The INLA backend returns a `flexybayes_inla` object whose internal
layout differs from the other engines' – it has no `$glm` shim – so it
needs its own [`tidy()`](https://generics.r-lib.org/reference/tidy.html)
method even though it shares the `flexybayes` parent class. The
fixed-effect summary is read directly off INLA's `summary.fixed` table,
whose `mean`, `sd`, and `0.025quant` / `0.975quant` columns map cleanly
onto the `broom`-canonical `estimate`, `std.error`, `conf.low`, and
`conf.high`.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
tidy(
  x,
  conf.int = TRUE,
  conf.level = 0.95,
  effects = c("fixed", "random"),
  ...
)
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

- effects:

  Character. Which effects to return: `"fixed"` for the population-level
  coefficients or `"random"` for the variance-component summary.
  Defaults to `"fixed"`.

- ...:

  Currently unused; present for generic compatibility.

## Value

A `data.frame` with one row per term and the columns `term`, `estimate`,
`std.error`, and (when `conf.int = TRUE`) `conf.low` / `conf.high`. The
rows are the fixed-effect terms under `effects = "fixed"` and the
variance components under `effects = "random"`.

## Details

Because the INLA fixed-effect intervals come from the marginal
posteriors INLA has already integrated, the `conf.level` argument is
accepted for generic compatibility but only the 95% bounds INLA reports
are returned; a one-off message notes this when a different level is
requested rather than silently ignoring it.

The `effects` argument is the same one
[`tidy.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes.md)
carries, and means the same thing on both: `"fixed"` returns the
population-level coefficients, `"random"` the variance-component
summary. Before 0.9.1 this method had no such argument, so
`effects = "random"` on an INLA fit was absorbed by `...` and the
FIXED-effect table came back under a call that asked for variance
components – a wrong answer with no error. An unrecognised value now
stops.

## See also

[`tidy.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes.md)

## Examples

``` r
if (FALSE) { # requireNamespace("INLA", quietly = TRUE)
# \donttest{
set.seed(1)
dat <- data.frame(yield = rnorm(40),
                  env = factor(rep(1:4, each = 10)),
                  block = factor(rep(1:5, times = 8)))
fit <- flexybayes(yield ~ env, random = ~block, data = dat,
                  backend = "inla", verbose = FALSE)
tidy(fit)
tidy(fit, effects = "random")
# }
}
```
