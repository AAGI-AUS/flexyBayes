# Tidy a flexyBayes fit into a one-row-per-term data frame

Turns a flexyBayes fit into a flat, `broom`-style `data.frame`: one row
per model term, with stable, documented columns. The method is
registered against the
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) generic
(re-exported by `broom`), so `broom::tidy(fit)` and
`generics::tidy(fit)` both dispatch here.

## Usage

``` r
# S3 method for class 'flexybayes'
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

  A flexyBayes fit: `flexybayes` (greta backend) or `flexybayes_brms`
  (brms backend, which inherits this method).

- conf.int:

  Logical. Whether to attach credible intervals. Defaults to `TRUE`.

- conf.level:

  Numeric in `(0, 1)`. The credible level for the intervals. Defaults to
  `0.95`.

- effects:

  Character. Which effects to return: `"fixed"` for the population-level
  (fixed) coefficients or `"random"` for the variance-component summary.
  Defaults to `"fixed"`.

- ...:

  Currently unused; present for generic compatibility.

## Value

A `data.frame` with one row per term and columns:

- term:

  Character. The coefficient or variance-component name.

- estimate:

  Numeric. The posterior mean.

- std.error:

  Numeric. The posterior standard deviation.

- conf.low:

  Numeric. The lower credible bound (present when `conf.int = TRUE`).

- conf.high:

  Numeric. The upper credible bound (present when `conf.int = TRUE`).

An empty `data.frame` is returned when the requested effects are absent
(for example `effects = "random"` on a fixed-effects-only fit).

## Details

This is the supported accessor for cross-engine summaries. The hub
returns backend-specific objects – `flexybayes` (greta),
`flexybayes_brms` (brms), `flexybayes_inla` (INLA) – whose internal
layouts differ. Tidying through this generic yields the same columns
across all three, so a greta-versus-INLA triangulation table can be
assembled by `rbind`-ing two
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) outputs
rather than reaching into each backend's slots by hand.

The credible intervals are posterior quantile-based intervals, not
frequentist confidence intervals; they are reported in the
`broom`-canonical `conf.low` / `conf.high` columns. The `std.error`
column carries the posterior standard deviation of each term, again
under the `broom`-canonical (dotless to the user, dotted in the column
name) label.

## See also

[`glance.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/glance.flexybayes.md),
[`augment.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/augment.flexybayes.md),
[`tidy.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes_inla.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- flexybayes(yield ~ env, data = dat, backend = "greta")
tidy(fit)
tidy(fit, effects = "random")
} # }
```
