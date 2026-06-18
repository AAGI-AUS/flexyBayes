# Fit a flexyBayes model via the brms (Stan) engine

Engine pin: fits the model with Stan through brms only. This is sugar
for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)`(..., backend = "brms")`
and accepts the same arguments and grammars — an ASReml `fixed` /
`random` / `rcov` specification or a brms-style bar-grouped formula (see
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
  (e.g. `formula` / `fixed`, `random`, `rcov`, `data`, `family`,
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

The brms / Stan engine cannot represent an ASReml structured-covariance
term (`vm`, `ped`, `fa`, `us`, `ar1`) or a `low_rank` smooth
approximation; such a model raises a structured refusal naming the
offending construct. Re-fit with
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
(full MCMC) or, when the model is latent-Gaussian feasible,
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md).

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal entry that picks a backend;
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
/
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
for the other engine pins;
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
for building a brms-grammar IR.

Other flexyBayes engine pins:
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md),
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)

## Examples

``` r
# \donttest{
if (requireNamespace("brms", quietly = TRUE) &&
    requireNamespace("lme4", quietly = TRUE)) {
  data(sleepstudy, package = "lme4")
  fit <- fb_brms(Reaction ~ Days + (1 | Subject), data = sleepstudy,
                 chains = 1)
  coef(fit)
}
# }
```
