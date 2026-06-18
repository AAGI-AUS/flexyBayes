# Fit a flexyBayes model via the greta engine

Engine pin: fits the model with greta (full Hamiltonian Monte Carlo)
only. This is sugar for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)`(..., backend = "greta")`
and accepts the same grammars — an ASReml `fixed` / `random` / `rcov`
specification, a brms-style bar-grouped formula, or a native
`greta_model` graph built with `greta::model()`. A formula is lowered
through the shared emit path; a native graph is fit directly by
`greta::mcmc()` and returned as a `flexybayes_direct_greta` object.

## Usage

``` r
fb_greta(...)
```

## Arguments

- ...:

  Arguments passed to
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  (e.g. `fixed`, `random`, `rcov`, `data`, `family`, `prior`, `syntax`),
  or a native `greta_model` / greta-source IR as the model-spec slot.
  The `backend` argument is pinned to `"greta"`; a conflicting `backend`
  value raises a structured refusal. The pre-v0.5.0 `model = `
  native-graph argument is remapped to the model-spec slot for
  call-compatibility, so `fb_greta(model = m)` still fits a native
  graph.

## Value

For a formula: an object of class `"flexybayes"` (a greta fit). For a
native `greta_model`: a `flexybayes_direct_greta` object (subclass of
`"flexybayes"`) carrying the MCMC draws, a GLM-compatible shim, and the
canonical-name map; see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the shared structure.

## Details

To attach a canonical-name map to a native graph (so
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
can align it with an INLA or Stan fit), build the intermediate
representation first and pass it on: `fb_greta(fb_from_greta(model,`
`canonical_names = c(...)))`.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry that picks a backend;
[`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md)
for building a greta-source IR with a canonical-name map.

Other flexyBayes engine pins:
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md),
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# live greta fit -- needs a working Python / TensorFlow stack
data(sleepstudy, package = "lme4")
fit <- fb_greta(Reaction ~ Days + (1 | Subject), data = sleepstudy,
                n_samples = 200, warmup = 200, chains = 1,
                verbose = FALSE, mcmc_verbose = FALSE)
coef(fit)
} # }
```
