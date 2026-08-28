# Summarise a flexyBayes fit

Prints the fixed-effect posterior summaries, the variance components
with the prior each one carried, the engine's own convergence
diagnostics, and whatever engine-native panels the fit earns. On a brms
fit with a sectioned residual it additionally prints the residual
variance and standard deviation for each level of the sectioning factor,
computed from the draws rather than by transforming a posterior mean.

## Usage

``` r
# S3 method for class 'flexybayes'
summary(object, ...)
```

## Arguments

- object:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

Invisibly, an object of class `c("summary.flexybayes", "list")` carrying
the slots below. Printed as a side effect.

- `fixed`:

  Data frame: `term`, `estimate`, `std.error`, `conf.low`, `conf.high`.
  Posterior mean, posterior standard deviation and a credible interval,
  not a sampling-theory estimate and its confidence interval.

- `varcomp`:

  Data frame: `component`, `estimate`, `std.error`, `conf.low`,
  `conf.high`, `prior`, `note`. One row per variance component (and per
  correlation on a fit carrying a latent field), on the
  standard-deviation scale, under the canonical component name. `prior`
  is a projection of
  [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md).
  `note` is `"collapsed"` when the component's 97.5% quantile sits below
  1% of the posterior-median residual standard deviation – a display
  heuristic for a posterior piled against zero, not a test.

- `random`:

  Named list, one data frame per grouping factor, with the columns
  `group`, `level`, `estimate`, `std.error`, `conf.low`, `conf.high`.
  Empty when the model carries no random terms. The same object
  [`ranef()`](https://aagi-aus.github.io/flexyBayes/reference/ranef.md)
  returns.

- `missing`:

  Data frame of the unobserved design cells – ASReml's `mv` factor –
  with the columns `row`, `estimate`, `std.error`, `conf.low`,
  `conf.high`, followed by any design index variables the fit recorded.
  `row` indexes the data the engine was handed. The posterior is the
  engine's own: INLA's fitted-value marginal at that row, brms's sampled
  `mi()` response. Zero rows – never `NULL` – when every response was
  observed. The same table [`coef()`](https://rdrr.io/r/stats/coef.html)
  returns for `what = "missing"`.

- `converge`:

  List, in the engine's own terms: R-hat, effective sample sizes and
  divergent transitions where a sampler ran; mode status, marginal
  likelihood and the largest Kullback-Leibler divergence where INLA's
  Laplace approximation did. No R-hat is invented for an approximation
  that has none.

- `n_design`:

  Rows the engine was handed.

- `n_observed`:

  How many of those carried an observed response.

- `na_action`:

  The missing-response record the fit carries, or `NULL` on a fit
  assembled without that layer.

- `model`:

  One-line description of the random and residual structure, derived
  from the model representation rather than from any engine's emitted
  formula.

- `engine`:

  `"inla"` or `"brms"`.

- `call`:

  The recorded call.

- `spatial_field`:

  **Present only on a fit carrying an autoregressive latent field**, and
  absent – not `NULL`-valued – otherwise. A data frame of `parameter`,
  `median`, `lower`, `upper`: the field's correlation and
  standard-deviation parameters and the nugget standard deviation, on
  the scale a reader thinks in rather than INLA's precision scale. The
  point estimate is the posterior median, because the standard
  deviations are read off precision marginals and only a monotone
  summary survives the reciprocal square root exactly. The same table
  the printed field panel renders, unrounded.

## Details

The returned object is the same eleven-slot `summary.flexybayes` object
on every active engine, so `summary(fit)$varcomp` answers whichever
engine ran the fit. Before 0.9.1 the two engines returned two
incomparable objects – INLA's four-slot list and brms's `brmssummary` –
and neither carried a variance-component table at all.

A fit whose model carries an autoregressive latent field gains one
further slot, `spatial_field`, which is the field's own parameters on
the correlation and standard-deviation scales. It is the one
engine-native slot the object carries, it is present only where the
model has a field, and the eleven above it are on every fit.

## See also

[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
for the full resolved prior, and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) for the design and
observed counts on their own.
