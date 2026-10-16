# Construct an inference-engine specification

Names a concrete inference engine and its tuning options. The result is
passed as the `backend` argument of the fitting verbs:
`flexybayes(..., backend = fb_engine("brms", chains = 4L))`. The bare
string form (`backend = "brms"`) remains valid and is equivalent to the
default `fb_engine()` for that engine.

## Usage

``` r
fb_engine(name, opts = list(), ...)
```

## Arguments

- name:

  Character(1): the engine, one of `"inla"`, `"brms"`.

- opts:

  Named list of tuning options. Recognised names are `n_samples`,
  `warmup`, `chains` and `adapt_delta`, described under "Recognised
  options" below. An unrecognised name is an error.

- ...:

  Tuning options given individually, merged into `opts` (e.g.
  `fb_engine("brms", chains = 4L)`).

## Value

An `fb_engine` object: a classed list with elements `name`, `paradigm`
(one of `mcmc`, `laplace`, `vb`, `map`), `toolchain_status` (one of
`ready`, `requires_install`, `unavailable`), and `opts`.

## Details

`name` and the derived `paradigm` are closed vocabularies. `"auto"` is a
routing directive, not an engine. Pass `backend = "auto"` directly for
automatic routing.

## Recognised options

All four are sampler controls, read by the routes that emit through Stan
and ignored by the rest. An unrecognised name is an error at
construction, not a silently dropped argument.

- `n_samples`:

  Post-warmup draws to keep, per chain.

- `warmup`:

  Warmup iterations, per chain.

- `chains`:

  Number of chains.

- `adapt_delta`:

  Target acceptance probability for Stan's step-size adaptation, between
  0 and 1. Raising it makes the sampler take smaller steps, which costs
  time and removes the divergent transitions a smaller step would have
  avoided. brms's own default is 0.8.

## Deprecated

`fb_engine()` is deprecated as of 0.13.0. Pass the backend and engine
settings directly to
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).

## See also

[`fb_approx()`](https://aagi-aus.github.io/flexyBayes/reference/fb_approx.md)

## Examples

``` r
plan <- flexybayes(y ~ x, data = data.frame(y = 1:4, x = 4:1),
                   backend = "brms", chains = 4L, plan = TRUE)
#> flexyBayes: by default the residual and each random-effect standard deviation have a half_normal(scale = 6.455) prior, in the units of the response, with the scale divided by sd(x) for a random slope on x. See ?fb_prior to set other priors; options(flexyBayes.silence_default_prior_note = TRUE) hides this note.
plan$backend_chosen
#> [1] "brms"
```
