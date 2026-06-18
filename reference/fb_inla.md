# Fit a flexyBayes model via the INLA engine

Engine pin: fits the model with INLA (integrated nested Laplace
approximation) only. This is sugar for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)`(..., backend = "inla")`
and accepts the same arguments and grammars (an ASReml `fixed` /
`random` / `rcov` specification or a brms-style bar-grouped formula –
see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the full argument list). The model must be latent-Gaussian feasible;
if it is not, the shared `lgm_gate()` raises a structured refusal naming
the offending term, exactly as `flexybayes(backend = "inla")` does.

## Usage

``` r
fb_inla(...)
```

## Arguments

- ...:

  Arguments passed to
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  (e.g. `fixed`, `random`, `rcov`, `data`, `family`, `prior`, `syntax`).
  The `backend` argument is pinned to `"inla"` and must not be supplied.

## Value

An object of class `"flexybayes"` (specifically a `flexybayes_inla`
fit); see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the structure.

## Details

Sampling-control arguments (`n_samples`, `warmup`, `chains`,
`mcmc_verbose`) are accepted for call-compatibility with the other
engine pins but are inert under INLA's deterministic Laplace
approximation.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry that picks a backend;
[`fb_from_asreml()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_asreml.md)
/
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
for the ingest layer.

Other flexyBayes engine pins:
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md),
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)

## Examples

``` r
df <- data.frame(
  yield = rnorm(40),
  geno  = factor(rep(letters[1:8], 5)),
  env   = factor(rep(c("a", "b"), 20))
)
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  fit <- fb_inla(yield ~ env, random = ~ geno, data = df)
}
# }
```
