# Construct an inference-engine specification

Names a concrete inference engine and its tuning options. The result is
passed as the `backend` argument of the fitting verbs:
`flexybayes(..., backend = fb_engine("greta", chains = 4L))`. The bare
string form (`backend = "greta"`) remains valid and is equivalent to the
default `fb_engine()` for that engine.

## Usage

``` r
fb_engine(name, opts = list(), ...)
```

## Arguments

- name:

  Character(1): the engine, one of `"greta"`, `"inla"`, `"brms"`.

- opts:

  Named list of tuning options. Recognised names are `n_samples`,
  `warmup`, `chains`; an unrecognised name is an error.

- ...:

  Tuning options given individually, merged into `opts` (e.g.
  `fb_engine("greta", chains = 4L)`).

## Value

An `fb_engine` object: a classed list with elements `name`, `paradigm`
(one of `mcmc`, `laplace`, `vb`, `map`), `toolchain_status` (one of
`ready`, `requires_install`, `unavailable`), and `opts`.

## Details

`name` and the derived `paradigm` are closed vocabularies. `"auto"` is a
routing directive, not an engine; pass `backend = "auto"` directly for
automatic routing.

## See also

[`fb_approx()`](https://aagi-aus.github.io/flexyBayes/reference/fb_approx.md)

## Examples

``` r
e <- fb_engine("greta", chains = 4L)
e$paradigm
#> [1] "mcmc"
e$toolchain_status
#> [1] "requires_install"
```
