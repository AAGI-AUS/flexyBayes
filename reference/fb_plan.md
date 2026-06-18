# Plan a flexyBayes fit without firing the backend

Returns the dispatch + preflight + representation + memory decision the
routing layer would make, without running MCMC or Laplace approximation.
Useful for verifying the backend chosen, memory estimate, and any
structural refusals before paying the fit cost.

## Usage

``` r
fb_plan(
  formula,
  data,
  backend = c("auto", "greta", "inla", "brms"),
  priors = NULL,
  known_matrices = list(),
  family = "gaussian",
  link = NULL,
  weights = NULL,
  aggregate = "auto",
  memory_ceiling_gb = NULL,
  predict_plan = NULL,
  ...
)
```

## Arguments

- formula:

  a brms-style two-sided formula (`y ~ x + s(z) + (1 | g)`).

- data:

  a data.frame.

- backend:

  one of `"greta"`, `"inla"`, `"brms"`, `"auto"` (default `"auto"`).

- priors:

  optional
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  list; defaults to the v0.2 uniform-on-SD default.

- known_matrices:

  named list of structured-covariance matrices referenced by `vm()` or
  `ped()` terms.

- family, link:

  standard family/link arguments.

- weights:

  optional observation weights.

- aggregate:

  `"auto"` / `TRUE` / `FALSE` — as on
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).

- memory_ceiling_gb:

  optional override for the preflight memory ceiling (defaults to
  `flexyBayes.preflight_ceiling_gb` option, or
  `flexyBayes.preflight_ram_fraction` x available RAM).

- predict_plan:

  optional `list(newdata = ..., chunk_size = ...)` to compute a
  prediction-shape plan. Plan-only; does not fire
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- ...:

  currently unused; reserved for future plan inputs.

## Value

an `<fb_plan>` classed list. See
[`print.fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_plan.md)
for the surface;
[`summary.fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/summary.fb_plan.md)
for the verbose dump;
[`as.data.frame.fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/as.data.frame.fb_plan.md)
for the programmatic-consumer shape.
