# Report inference-backend readiness

Checks which inference backends flexyBayes can route to are installed
and usable in the current session, returning a small table you can
inspect before fitting. The check is read-only: it probes package
availability and, for greta, the reachability of its Python / TensorFlow
stack, without building a model or starting a fit.

## Usage

``` r
fb_backend_status(deep = TRUE)
```

## Arguments

- deep:

  Logical; if `TRUE` (default) greta's readiness is probed by
  initialising its Python / TensorFlow stack (the probe's output is
  captured, never leaked). If `FALSE`, greta is reported as installed
  without touching Python – a fast, non-invasive check whose `usable`
  value is `NA` ("installed, not probed").

## Value

A data frame of class `fb_backend_status` with one row per backend and
the columns `backend`, `installed` (logical), `usable` (logical), and
`note` (a human-readable status, including the install command when a
backend is absent). A `print` method renders it as a readiness table.

## Details

`installed` records whether the backend's R package is present. `usable`
records whether the backend can actually run a fit now – for greta this
additionally requires a working Python / TensorFlow stack, so a greta
that is `installed` but not `usable` needs
`greta::install_greta_deps()`. The dormant opt-in `gretaR` engine is
reported separately by
[`gretaR_status()`](https://aagi-aus.github.io/flexyBayes/reference/gretaR_status.md).

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry, the
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
/
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
/
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
single-engine pins, and
[`gretaR_status()`](https://aagi-aus.github.io/flexyBayes/reference/gretaR_status.md)
for the dormant gretaR slot.

## Examples

``` r
# \donttest{
# Probing greta initialises its Python / TensorFlow stack, so the first
# call can be slow; wrapped in \donttest{} for that reason.
fb_backend_status()
#> flexyBayes backend readiness
#> ----------------------------------------------------------------
#>   [--] greta   not installed: install.packages('greta', repos = c('https://greta-dev.r-universe.dev', getOption('repos')))
#>   [--] INLA    not installed: install.packages('INLA', repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'))
#>   [--] brms    not installed: install.packages('brms')
#> 
#>   No active inference backend is usable -- install at least one of the above before fitting.
# }
```
