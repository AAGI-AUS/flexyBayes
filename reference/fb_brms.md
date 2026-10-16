# Fit a flexyBayes model via the brms (Stan) engine

Engine pin: fits the model with Stan through brms only. This is sugar
for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)`(..., backend = "brms")`
and accepts the same arguments and grammars, an ASReml `fixed` /
`random` / `residual` specification or a brms-style bar-grouped formula
(see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the full argument list). flexyBayes builds the intermediate
representation, translates the prior, calls `brms::brm()`, and wraps the
result. The live `brmsfit` is available on the `$brms` slot for brms's
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
  `prior`, `syntax`). The `backend` argument is pinned to `"brms"`. A
  conflicting `backend` value raises a structured refusal (the redundant
  `backend = "brms"` is accepted). The pre-v0.5.0 `formula = ` argument
  is remapped to the universal entry's model-spec slot for
  call-compatibility.

## Value

An object of class `"flexybayes_brms"` (a subclass of `"flexybayes"`)
carrying the live `brmsfit` on `$brms`. See
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the shared structure.

## Details

flexyBayes has no brms lowering for some ASReml structures: `fa()`, the
autoregressive fields, a blocks or low-rank `vm()` / `ped()` carrier,
and a `low_rank` smooth approximation. Such a model raises a structured
refusal naming the term. Several of these fit on INLA through
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md).
The capability table in the README lists what each engine fits.

## Deprecated

`fb_brms()` is deprecated as of 0.13.0. Use
`flexybayes(..., backend = "brms")`.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry that picks a backend.
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
for the other active engine pin and
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
  fit <- flexybayes(
    Reaction ~ Days + (1 | Subject), data = sleepstudy,
    backend = "brms", chains = 1, n_samples = 4000, warmup = 1000,
    seed = 1
  )
  coef(fit)
}
# }
```
