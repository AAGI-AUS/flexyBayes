# Fit a Dirichlet distribution to compositional data

Estimates the concentration vector \\\alpha\\ of a Dirichlet
distribution from a matrix of compositional rows (each row a composition
on the simplex), by maximum likelihood via base
[`optim()`](https://rdrr.io/r/stats/optim.html) (dependency-free,
deterministic).

## Usage

``` r
fb_dirichlet(x, method = "ml", conf_level = 0.95, eps = 1e-06)
```

## Arguments

- x:

  A numeric matrix or `data.frame` of compositional rows: at least two
  columns, all entries non-negative, with at least four rows.

- method:

  Character. Currently only `"ml"` (maximum likelihood).

- conf_level:

  Numeric in `(0, 1)`. The interval level for the parameter summary.
  Defaults to `0.95`.

- eps:

  Numeric. The boundary nudge applied to exact zeros / ones. Defaults to
  `1e-06`.

## Value

An object of class `c("fb_dirichlet_fit", "fb_family_fit")`: a list with
`estimates` (a `data.frame` of `term` / `estimate` / `std.error` /
`conf.low` / `conf.high`, one row per component), `mean_composition`
(numeric, summing to one), `method`, `n_obs`, `n_components`, and
`logLik`.

## Details

Rows are renormalised to sum to one, and exact zeros or ones are nudged
inside the open simplex by a small `eps` (a Dirichlet density is `-Inf`
on the simplex boundary). The concentrations are optimised on the log
scale to keep them positive, and standard errors are recovered from the
Hessian by the delta method.

The fitted mean composition (`alpha / sum(alpha)`) is reported on
`fit$mean_composition`.

## See also

[`fb_family_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_family_dirichlet.md),
[`rdirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/rdirichlet.md),
[`tidy.fb_dirichlet_fit()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.fb_dirichlet_fit.md)

## Examples

``` r
set.seed(1)
X <- flexyBayes:::rdirichlet(300L, alpha = c(2, 5, 3))
fit <- flexyBayes:::fb_dirichlet(X)
fit$estimates
#>    term estimate std.error conf.low conf.high
#> c1   c1 1.886471 0.1117041 1.667535  2.105407
#> c2   c2 4.673664 0.2811266 4.122666  5.224662
#> c3   c3 2.792916 0.1667429 2.466106  3.119727
```
