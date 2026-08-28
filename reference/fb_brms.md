# Fit a flexyBayes model via the brms (Stan) engine

Engine pin: fits the model with Stan through brms only. This is sugar
for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)`(..., backend = "brms")`
and accepts the same arguments and grammars — an ASReml `fixed` /
`random` / `residual` specification or a brms-style bar-grouped formula
(see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the full argument list). flexyBayes builds the intermediate
representation, translates the prior, calls `brms::brm()`, and wraps the
result; the live `brmsfit` is available on the `$brms` slot for brms's
own posterior tooling (`loo()`, `posterior_predict()`,
`bayes_factor()`).

## Usage

``` r
fb_brms(...)
```

## Arguments

- ...:

  Arguments passed to
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  (e.g. `formula` / `fixed`, `random`, `residual`, `data`, `family`,
  `prior`, `syntax`). The `backend` argument is pinned to `"brms"`; a
  conflicting `backend` value raises a structured refusal (the redundant
  `backend = "brms"` is accepted). The pre-v0.5.0 `formula = ` argument
  is remapped to the universal entry's model-spec slot for
  call-compatibility.

## Value

An object of class `"flexybayes_brms"` (a subclass of `"flexybayes"`)
carrying the live `brmsfit` on `$brms`; see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the shared structure.

## Details

The brms / Stan engine cannot represent every ASReml structured-
covariance term (`fa`, `us`, `ar1`, or a block-diagonal / low-rank
`vm()` / `ped()` carrier) or a `low_rank` smooth approximation; such a
model raises a structured refusal naming the offending construct. When
the model is latent-Gaussian feasible, re-fit with
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md);
flexyBayes has no active engine for a structured-covariance term that is
neither brms-representable nor latent-Gaussian feasible (see `NEWS.md`,
0.9.3).

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry that picks a backend;
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
for the other active engine pin;
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
for building a brms-grammar IR.

Other flexyBayes engine pins:
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)

## Examples

``` r
# Held back from the default example run because the brms route compiles
# a Stan program before it samples, and the compile alone takes far
# longer than an example should. The sampler settings are the smallest
# that reach an effective sample size the fit does not warn about, and
# the seed makes that reproducible. On a single chain they cost about
# half a second next to the compile.
# \donttest{
if (requireNamespace("brms", quietly = TRUE) &&
    requireNamespace("lme4", quietly = TRUE)) {
  data(sleepstudy, package = "lme4")
  fit <- fb_brms(Reaction ~ Days + (1 | Subject), data = sleepstudy,
                 chains = 1, n_samples = 4000, warmup = 1000, seed = 1)
  coef(fit)
}
# }
```
