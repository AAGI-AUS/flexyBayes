# Fit a generalised extreme value distribution to block maxima

Estimates the location, scale, and shape of a GEV distribution from a
vector of block-maxima observations by maximum likelihood (base
[`optim()`](https://rdrr.io/r/stats/optim.html), dependency-free and
deterministic given the data).

## Usage

``` r
fb_gev(y, return_periods = c(10, 50, 100), conf_level = 0.95)
```

## Arguments

- y:

  Numeric vector of block-maxima observations. Must contain at least
  four finite values.

- return_periods:

  Numeric vector of return periods (in blocks) at which to report return
  levels. Defaults to `c(10, 50, 100)`.

- conf_level:

  Numeric in `(0, 1)`. The interval level for the parameter summary.
  Defaults to `0.95`.

## Value

An object of class `c("fb_gev_fit", "fb_family_fit")`: a list with
`estimates` (a `data.frame` of `term` / `estimate` / `std.error` /
`conf.low` / `conf.high`), `return_levels` (a `data.frame` of
`return_period` / `return_level`), `method`, `n_obs`, and `logLik`.

## Details

Maximum likelihood maximises the GEV log-likelihood directly, with the
scale optimised on the log scale to keep it positive and the standard
errors recovered from the observed-information (Hessian) matrix by the
delta method on the scale.

Return levels (the level exceeded once per `m` blocks on average) are
computed from the fitted parameters and reported on `fit$return_levels`
for the requested return periods.

A scalable Bayesian GEV belongs on INLA's native `gev` / `bgev` family
and is planned for a future release; the greta backend ships no GEV
distribution.

## See also

[`fb_family_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_family_gev.md),
[`rgev()`](https://aagi-aus.github.io/flexyBayes/reference/rgev.md),
[`tidy.fb_gev_fit()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.fb_gev_fit.md)

## Examples

``` r
set.seed(1)
y <- rgev(200L, location = 10, scale = 2, shape = 0.1)
fit <- fb_gev(y)
fit$estimates
#>       term    estimate std.error    conf.low  conf.high
#> 1 location 10.19181337 0.1453225  9.90698659 10.4766402
#> 2    scale  1.82586952 0.1085912  1.61303475  2.0387043
#> 3    shape  0.08191751 0.0526013 -0.02117914  0.1850142
```
