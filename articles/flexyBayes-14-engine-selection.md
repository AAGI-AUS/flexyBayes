# Choosing an engine: the universal entry and the engine pins

flexyBayes fits one model specification through several backends. This
vignette explains how to let the package choose a backend for you, how
to force a particular one, and how the same call reaches every backend
the model supports. Fitting code is shown for reference – set
`eval = TRUE` locally with the relevant backends installed to run it.

## One entry, every backend

[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
(and its full-name twin
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md))
is the universal entry. It detects the input grammar from the call
shape, builds a backend-agnostic intermediate representation, and
dispatches to a backend. A no-`backend` call uses `backend = "auto"`,
which routes to INLA when the model is latent-Gaussian feasible and to
greta otherwise, so a single call reaches whichever engine the model
supports.

``` r

data(sleepstudy, package = "lme4")

# brms / lme4-style grammar (a bar-grouped formula): detected automatically.
fit <- fb(Reaction ~ Days + (1 | Subject), data = sleepstudy)

# ASReml grammar: fixed / random / rcov.
fit_asreml <- fb(
  fixed  = Reaction ~ Days,
  random = ~ Subject,
  data   = sleepstudy
)
```

The grammar is read from the call, never guessed from intent. A
bar-grouped formula such as `Reaction ~ Days + (1 | Subject)` is read as
brms grammar; a bar-free `fixed` formula paired with `random` / `rcov`
is ASReml grammar. The lone ambiguous case – a bare `y ~ x` with no
grouping, which means the same model either way – defaults to ASReml.
Force the reading with `syntax`:

``` r

# Treat a bar-free formula as brms grammar explicitly.
fit_brms_grammar <- fb(Reaction ~ Days, data = sleepstudy, syntax = "brms")
```

## Letting `auto` choose, and reading its decision

Under `backend = "auto"` the call consults the latent-Gaussian gate. If
the model is feasible and INLA is installed, it fits via INLA (a fast,
deterministic Laplace approximation). Otherwise it falls back to greta
(full Hamiltonian Monte Carlo), with a one-time note. The recorded
decision – including the engines considered and why each was or was not
chosen – is available after the fit:

``` r

fit <- fb(Reaction ~ Days + (1 | Subject), data = sleepstudy,
          backend = "auto")
backend_decision(fit)
```

To preview the routing without fitting, pass `plan = TRUE`, which
returns a planning object describing the chosen backend, the
representation plan, and the prediction plan.

``` r

fb(Reaction ~ Days + (1 | Subject), data = sleepstudy, plan = TRUE)
```

## Forcing one engine: the pins

When you want a specific engine, the pins make the intent explicit. Each
is exactly
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
with `backend` fixed, and each accepts the same grammars and arguments
as the universal entry:

| Pin | Equivalent to | Inference |
|----|----|----|
| [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md) | `fb(..., backend = "greta")` | full MCMC (Hamiltonian Monte Carlo) |
| [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md) | `fb(..., backend = "inla")` | Laplace approximation |
| [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) | `fb(..., backend = "brms")` | Stan via brms (Hamiltonian Monte Carlo) |

``` r

fit_g <- fb_greta(Reaction ~ Days + (1 | Subject), data = sleepstudy)
fit_i <- fb_inla(Reaction ~ Days + (1 | Subject),  data = sleepstudy)
fit_s <- fb_brms(Reaction ~ Days + (1 | Subject),  data = sleepstudy)
```

A pin fits via one engine only, so a conflicting `backend` argument is a
clear error rather than a silent override:

``` r

# Refused: fb_inla() pins INLA, so it cannot also take backend = "greta".
fb_inla(Reaction ~ Days + (1 | Subject), data = sleepstudy,
        backend = "greta")
```

Not every engine can represent every model. The brms / Stan engine, for
example, cannot represent an ASReml structured-covariance term (`vm`,
`ped`, `fa`, `us`, `ar1`) or a `low_rank` smooth approximation. Such a
request raises a structured refusal that names the offending construct
and points to an engine that can fit it – the package never silently
fits a different model than the one asked for.

## Native greta models

If you have already built a model with native greta primitives, pass the
`greta::model()` object straight to
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).
It is fit directly by `greta::mcmc()` and returned with the standard
post-fit surface ([`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)).

``` r

library(greta)
y     <- as_data(sleepstudy$Reaction)
x     <- as_data(sleepstudy$Days)
b0    <- normal(0, 100)
b1    <- normal(0, 100)
sigma <- uniform(0, 5 * sd(sleepstudy$Reaction))
distribution(y) <- normal(b0 + b1 * x, sigma)
m <- model(b0, b1, sigma)

fit_native <- fb(m)
```

To align a native fit with an INLA or Stan fit under
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md),
the parameters need canonical names. Attach them once by building the
intermediate representation with
[`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md),
then pass the result to
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md):

``` r

ir <- fb_from_greta(
  m,
  canonical_names = c(b0 = "(Intercept)", b1 = "Days", sigma = "sigma")
)
fit_native <- fb(ir)
```

The universal entry gains no extra argument for this: the canonical map
travels on the intermediate representation, and the `fb_from_*()`
adapters are the place to prepare a specification before fitting.

## Triangulation: the reason the choice is offered

Because the same specification fits on engines that use different
inference machinery, you can fit it twice and compare. Close agreement
across an MCMC engine and a Laplace engine is evidence that the result
does not depend on the approximation; disagreement is a signal worth
investigating before trusting the fit.

``` r

fit_inla  <- fb_inla(Reaction ~ Days + (1 | Subject), data = sleepstudy)
fit_greta <- fb_greta(Reaction ~ Days + (1 | Subject), data = sleepstudy)
triangulate(fit_inla, fit_greta)
```

## Summary

- [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  /
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  is the universal entry; `backend = "auto"` reaches whichever engine
  the model supports.
- [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
  /
  [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
  /
  [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
  pin one engine each; a conflicting `backend` is refused.
- A native greta model is fit by passing it (or its
  [`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md)
  representation) to
  [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).
- An unsupported (model, engine) pair raises a structured refusal that
  names the offending construct.
- [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  compares two fits across engines as a robustness check.
