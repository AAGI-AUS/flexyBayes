# Posterior functionals of the variance components

The exact posterior of any function of the fitted variance components: a
heritability, a genetic correlation, a variance ratio, a repeatability.
Each functional is written as a two-sided formula whose left side names
the result and whose right side is an expression over the component
names that [`summary()`](https://rdrr.io/r/base/summary.html) reports.

## Usage

``` r
vpredict(object, ...)

# S3 method for class 'flexybayes'
vpredict(object, ..., n_draws = 1000L, prob = 0.95, draws = FALSE)
```

## Arguments

- object:

  A fitted `flexybayes` model.

- ...:

  One or more two-sided formulas, `name ~ expression`.

- n_draws:

  Number of posterior draws to evaluate over.

- prob:

  Width of the credible interval.

- draws:

  When `TRUE`, return the matrix of functional draws (one column per
  functional) instead of the summary.

## Value

A `data.frame` with one row per functional and the columns `functional`,
`estimate` (posterior mean), `median`, `std.error` (posterior standard
deviation), `conf.low` and `conf.high`. With `draws = TRUE`, a numeric
matrix of draws instead.

## Details

Components are on the **standard-deviation** scale, as everywhere else
in this package, so a variance ratio squares them. Writing a name the
fit does not carry is refused, and the refusal lists the names it does.

The interval is a quantile interval of the functional's own posterior.
Where the delta method would put a heritability's interval below zero,
this one cannot: every draw of `sd_g^2 / (sd_g^2 + sigma^2)` lies in
\[0, 1\], so the interval does too.

The method is also registered on sommer's `vpredict()` generic, so
`vpredict(fit, ...)` reaches it when sommer is attached after
flexyBayes.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- flexybayes(yield ~ 1, random = ~ gen, data = trial)
vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2))
} # }
```
