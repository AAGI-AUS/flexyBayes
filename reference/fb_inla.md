# Fit a flexyBayes model via the INLA engine

Engine pin: fits the model with INLA (integrated nested Laplace
approximation) only. This is sugar for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)`(..., backend = "inla")`
and accepts the same arguments and grammars (an ASReml `fixed` /
`random` / `residual` specification or a brms-style bar-grouped formula,
see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the full argument list). The model must be latent-Gaussian feasible.
If it is not, the shared `lgm_gate()` raises a structured refusal naming
the offending term, exactly as `flexybayes(backend = "inla")` does.

## Usage

``` r
fb_inla(...)
```

## Arguments

- ...:

  Arguments passed to
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  (e.g. `fixed`, `random`, `residual`, `data`, `family`, `prior`,
  `syntax`). The `backend` argument is pinned to `"inla"` and must not
  be supplied.

## Value

An object of class `"flexybayes"` (specifically a `flexybayes_inla`
fit). See
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the structure.

## Details

Sampling-control arguments (`n_samples`, `warmup`, `chains`, `seed`,
`control`, `mcmc_verbose`) are accepted for call-compatibility with the
other engine pins but are inert under INLA's deterministic Laplace
approximation. `seed` and `control` say so when supplied: the fit is
repeatable without a seed, because nothing in the approximation is
drawn.

## Deprecated

`fb_inla()` is deprecated in flexyBayes 0.13.0 and will be removed in a
later release. Call `flexybayes(..., backend = "inla")`, which takes the
same arguments. A call warns with a condition of class
`flexybayes_deprecated`.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry that picks a backend.
[`fb_from_asreml()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_asreml.md)
/
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
for the ingest layer.

Other flexyBayes engine pins:
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)

## Examples

``` r
set.seed(1L)
df <- data.frame(
  yield = rnorm(40),
  geno  = factor(rep(letters[1:8], 5)),
  env   = factor(rep(c("a", "b"), 20))
)
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  fit <- flexybayes(yield ~ env, random = ~ geno, data = df,
                    backend = "inla")
}
#> 
#> -- flexyBayes: INLA fit ----------------------------------------
#>   formula: fb_standardised_response ~ 1 + env + f(geno, model = "iid", <prior: half-normal-SD(5) on the standardised response = half-normal-SD(4.433) on the response>)
#>   family:  gaussian
#> ------------------------------------------------------------ 
#> 
# }
```
