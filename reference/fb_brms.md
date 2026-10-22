# Fit a flexyBayes model via the brms (Stan) engine

Engine pin: fits the model with Stan through brms only. This is sugar
for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)`(..., backend = "brms")`
and accepts the same arguments and grammars, an ASReml `fixed` /
`random` / `residual` specification or a brms-style bar-grouped formula
(see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the full argument list). flexyBayes builds the intermediate
representation, translates the prior, calls
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html), and
wraps the result. The live `brmsfit` is available on the `$brms` slot
for brms's own posterior tooling (`loo()`, `posterior_predict()`,
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
#> 
#> -- flexyBayes: brms (Stan passthrough) fit ------------------
#>   formula: Reaction ~ Days + (1 | Subject)
#>   family : gaussian (identity link)
#>   priors : 2 row(s)
#> ------------------------------------------------------------ 
#> 
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.27 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 5000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 5000 [  4%]  (Warmup)
#> Chain 1: Iteration:  400 / 5000 [  8%]  (Warmup)
#> Chain 1: Iteration:  600 / 5000 [ 12%]  (Warmup)
#> Chain 1: Iteration:  800 / 5000 [ 16%]  (Warmup)
#> Chain 1: Iteration: 1000 / 5000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 1001 / 5000 [ 20%]  (Sampling)
#> Chain 1: Iteration: 1200 / 5000 [ 24%]  (Sampling)
#> Chain 1: Iteration: 1400 / 5000 [ 28%]  (Sampling)
#> Chain 1: Iteration: 1600 / 5000 [ 32%]  (Sampling)
#> Chain 1: Iteration: 1800 / 5000 [ 36%]  (Sampling)
#> Chain 1: Iteration: 2000 / 5000 [ 40%]  (Sampling)
#> Chain 1: Iteration: 2200 / 5000 [ 44%]  (Sampling)
#> Chain 1: Iteration: 2400 / 5000 [ 48%]  (Sampling)
#> Chain 1: Iteration: 2600 / 5000 [ 52%]  (Sampling)
#> Chain 1: Iteration: 2800 / 5000 [ 56%]  (Sampling)
#> Chain 1: Iteration: 3000 / 5000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 3200 / 5000 [ 64%]  (Sampling)
#> Chain 1: Iteration: 3400 / 5000 [ 68%]  (Sampling)
#> Chain 1: Iteration: 3600 / 5000 [ 72%]  (Sampling)
#> Chain 1: Iteration: 3800 / 5000 [ 76%]  (Sampling)
#> Chain 1: Iteration: 4000 / 5000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 4200 / 5000 [ 84%]  (Sampling)
#> Chain 1: Iteration: 4400 / 5000 [ 88%]  (Sampling)
#> Chain 1: Iteration: 4600 / 5000 [ 92%]  (Sampling)
#> Chain 1: Iteration: 4800 / 5000 [ 96%]  (Sampling)
#> Chain 1: Iteration: 5000 / 5000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.248 seconds (Warm-up)
#> Chain 1:                0.493 seconds (Sampling)
#> Chain 1:                0.741 seconds (Total)
#> Chain 1: 
#> (Intercept)        Days 
#>   251.32745    10.47849 
# }
```
