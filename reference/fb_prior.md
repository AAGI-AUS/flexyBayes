# Specify priors via the PC-canonical hybrid DSL

The flexyBayes prior DSL (domain-specific language) lives on the
standard-deviation scale (never precision / variance) and accepts two
canonical idioms: distributional (`half_normal(scale = 1)`) and
tail-quantile / PC (penalised complexity) (`pc(upper = 1, prob = 0.01)`
meaning `Pr(sigma > 1) = 0.01`). The PC idiom is the cross-engine
interlingua, it survives every backend because it is a probability
statement, not a distributional name.

## Usage

``` r
fb_prior(...)
```

## Arguments

- ...:

  One or more two-sided formulas of the form
  `target ~ distribution(args)`, one per parameter being given a prior.

  Examples:

  - `sigma ~ pc(upper = 2, prob = 0.05)`

  - `sd(group = "subject") ~ half_normal(scale = 1)`

  - `b("treatment") ~ student_t(df = 4, scale = 2.5)`

  - `cor(group = "subject") ~ lkj(eta = 2)`

  - `sd(group = "subject") ~ uniform(lower = 0, upper = 5)`

  Supported distribution families: `pc`, `half_normal`, `half_cauchy`,
  `student_t`, `normal`, `exponential`, `lkj`, `cauchy`, `gamma`,
  `uniform`, and `pooled` (see "Partially pooled SDs" below).
  `uniform()` on a variance component sits outside the PC-canonical
  interlingua, but both active backends represent it faithfully on the
  SD scale: the INLA backend via an expression-prior on the
  log-precision, and brms as a bounded `uniform()` prior on each simple
  random-effect (and residual) SD.

## Value

An `fb_prior` object, an S3 class inheriting from list, whose `$specs`
element carries the parsed `target` and `spec` pairs.

## Details

v0.1 supports the targets and distributions listed above in the file
header. Calls outside the supported set raise a structured error with
the supported list.

Argument matching is strict. Each distribution has a fixed parameter
list, and a call is matched against it by base R's own rules, so
`half_normal(1)` and `half_normal(scale = 1)` are the same
specification. An argument name the distribution does not have, a
duplicated argument, a missing required argument, a non-scalar or
non-finite value, or a value outside its domain (a non-positive scale, a
tail probability outside `(0, 1)`) is refused at construction with a
condition carrying a `flexybayes_refusal_*` class: see
[`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md).
Both halves of a PC prior are required: `pc(upper = 1)` states no
probability, and is refused instead of completed with a default the
caller never wrote.

## What each backend can carry

`fb_prior()` is engine-neutral: it records what was asked for. Whether a
given (target, distribution) pair can be *carried* depends on the
backend, and the two differ.

- **brms** takes any of the ten distributions on `sigma` and
  `sd(group = )`: the parameter is bounded below at zero and brms
  renormalises a two-sided density over that support, which is why
  `normal(0, s)` and `half_normal(s)` are the same prior there, and
  `normal` / `student_t` on `b()`. It carries no
  [`cor()`](https://rdrr.io/r/stats/cor.html) or
  [`smooth()`](https://rdrr.io/r/stats/smooth.html) prior, because no
  model this package emits on brms has either parameter.

- **INLA** takes `half_normal`, `uniform`, `pc` and `exponential` on
  `sigma`, `sd(group = )` and
  [`smooth()`](https://rdrr.io/r/stats/smooth.html) (`exponential` is
  the PC prior written in its distributional form), `normal` on `b()`,
  which arrives through `control.fixed`, and `normal` on
  `loading(group = )`, the loadings of `fa(f, k):group` relative to the
  anchor level. With no loading prior stated, each loading takes
  `normal(mean = 1, sd = 0.316)`, which is INLA's own prior on a copy
  scaling,
  [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
  lists it.

A specification the selected backend cannot carry is refused at the fit,
naming the row and both remedies, re-express it in a distribution that
engine does carry, or switch backend. Before 0.9.2 the INLA route
discarded such a row silently while
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
printed it as applied, so a prior-sensitivity analysis could vary a
scale and compare the engine's default with itself.
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
lists each prior with the engine slot it reached, as the engine reports
it after fitting, and names each stated or default prior that reached no
slot, with the reason.

A prior naming a coefficient or a variance component the model does not
have is likewise refused at the fit, against the engine's own parameter
list and in the vocabulary the caller wrote.

## The default prior

A fit that states neither `prior` nor `prior_vc_sd` runs under the
package default: `half_normal(scale = s)` on the residual standard
deviation and on the standard deviation of every variance component the
default reaches, on both engines. `s` is `5 * sd(y)` on an identity
link, `3` on a log link and `5` on a logit link. On an uncorrelated
random slope on a numeric covariate `x` it is `s / sd(x)`. The default
has no upper bound.
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
prints the scale a fit ran under. A Gaussian identity-link fit on INLA
runs on the response divided by `sd(y)`, where the same prior is
`half_normal(scale = 5)`. The printed formula and
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
show it on both scales.

## Partially pooled SDs

`pooled()` puts one shared prior on a set of per-level standard
deviations, instead of a separate prior on each. The target names the
term, and `by` the factor whose levels are pooled:

- `sd(group = "TRIAL:VARIETY", by = "TRIAL") ~ pooled(scale = half_normal(scale = 0.5))`
  for the genotype-by-trial SDs of `diag(TRIAL):VARIETY`.

- `sigma(by = "TRIAL") ~ pooled(scale = half_normal(scale = 0.5))` for
  the residual SDs of `dsum(~ units | TRIAL)`.

Either set can be pooled without the other. The log SDs of a pooled set
are normal with an estimated mean and an estimated SD, the pooling SD.
`scale` is the prior on the pooling SD and has no default. It can be any
prior brms accepts for a standard deviation. The mean has a
`normal(log(s / 5), 2)` prior, where `s` is the scale of the package's
default SD prior. On an identity link that is `normal(log(sd(y)), 2)`.

Pooling shrinks the per-level SDs toward each other, so the estimates
differ from ASReml's, which fits each SD separately. Only brms fits it,
because INLA cannot put a shared prior on separate variance parameters:
`backend = "inla"` refuses with `pooled_prior_requires_brms`, and
`backend = "auto"` uses brms. Printed output labels the fit "per-trial
SDs partially pooled", with the `by` factor's name in place of "trial".

## Examples

``` r
p <- fb_prior(
  sigma                 ~ pc(upper = 2, prob = 0.05),
  sd(group = "subject") ~ half_normal(scale = 1),
  b("treatment")        ~ student_t(df = 4, scale = 2.5)
)
p
#> <fb_prior> 3 specifications
#>   sigma ~ pc(upper = 2, prob = 0.05)
#>   sd(group = "subject") ~ half_normal(scale = 1)
#>   b("treatment") ~ student_t(df = 4, scale = 2.5)

# Model C with both per-trial SD blocks partially pooled
fb_prior(
  sd(group = "TRIAL:VARIETY", by = "TRIAL") ~
    pooled(scale = half_normal(scale = 0.5)),
  sigma(by = "TRIAL") ~ pooled(scale = half_normal(scale = 0.5))
)
#> <fb_prior> 2 specifications
#>   sd(group = "TRIAL:VARIETY", by = "TRIAL") ~ pooled(scale = half_normal(scale = 0.5))
#>   sigma(by = "TRIAL") ~ pooled(scale = half_normal(scale = 0.5))
```
