# Specify priors via the PC-canonical hybrid DSL

The flexyBayes prior DSL (domain-specific language) lives on the
standard-deviation scale (never precision / variance) and accepts two
canonical idioms: distributional (`half_normal(scale = 1)`) and
tail-quantile / PC (penalised complexity) (`pc(upper = 1, prob = 0.01)`
meaning `Pr(sigma > 1) = 0.01`). The PC idiom is the cross-engine
interlingua – it survives every backend because it is a probability
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
  `uniform`. Note that `uniform()` on a variance component sits outside
  the PC-canonical interlingua, but both active backends represent it
  faithfully on the SD scale: the INLA backend via an expression-prior
  on the log-precision, and brms as a bounded `uniform()` prior on each
  simple random-effect (and residual) SD.

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
condition carrying a `flexybayes_refusal_*` class – see
[`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md).
Both halves of a PC prior are required: `pc(upper = 1)` states no
probability, and is refused rather than completed with a default the
caller never wrote.

## What each backend can carry

`fb_prior()` is engine-neutral: it records what was asked for. Whether a
given (target, distribution) pair can be *carried* depends on the
backend, and the two differ.

- **brms** takes any of the ten distributions on `sigma` and
  `sd(group = )` – the parameter is bounded below at zero and brms
  renormalises a two-sided density over that support, which is why
  `normal(0, s)` and `half_normal(s)` are the same prior there – and
  `normal` / `student_t` on `b()`. It carries no
  [`cor()`](https://rdrr.io/r/stats/cor.html) or
  [`smooth()`](https://rdrr.io/r/stats/smooth.html) prior, because no
  model this package emits on brms has either parameter.

- **INLA** takes `half_normal`, `uniform`, `pc` and `exponential` on
  `sigma`, `sd(group = )` and
  [`smooth()`](https://rdrr.io/r/stats/smooth.html) (`exponential` is
  the PC prior written in its distributional form), and `normal` on
  `b()`, which arrives through `control.fixed`.

A specification the selected backend cannot carry is refused at the fit,
naming the row and both remedies – re-express it in a distribution that
engine does carry, or switch backend. Before 0.9.2 the INLA route
discarded such a row silently while
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
printed it as applied, so a prior-sensitivity analysis could vary a
scale and compare the engine's default with itself. Whatever
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
shows on a fit now reached the engine.

A prior naming a coefficient or a variance component the model does not
have is likewise refused at the fit, against the engine's own parameter
list and in the vocabulary the caller wrote.

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
```
