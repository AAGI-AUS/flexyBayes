# Report inference-backend readiness

Checks which inference backends flexyBayes can route to are installed
and usable in the current session, returning a small table you can
inspect before fitting. The check is read-only: it probes package
availability without building a model or starting a fit.

## Usage

``` r
fb_backend_status(deep = TRUE)
```

## Arguments

- deep:

  Logical; currently unused (reserved for a future backend whose
  readiness probe is expensive enough to need an opt-out). Default
  `TRUE`.

## Value

A data frame of class `fb_backend_status` with one row per backend and
the columns `backend`, `installed` (logical), `usable` (logical), and
`note` (a human-readable status, including the install command when a
backend is absent). A `print` method renders it as a readiness table.

## Details

`installed` records whether the backend's R package is present. `usable`
records whether the backend can actually run a fit now.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry, and the
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
/
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
single-engine pins.

## Examples

``` r
fb_backend_status()
#> flexyBayes backend readiness
#> ----------------------------------------------------------------
#>   [--] INLA    not installed: install.packages('INLA', repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'))
#>   [--] brms    not installed: install.packages('brms')
#> 
#>   No active inference backend is usable -- install at least one of the above before fitting.
```
