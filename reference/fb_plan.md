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
  random = NULL,
  residual = NULL,
  backend = c("auto", "inla", "brms"),
  priors = NULL,
  known_matrices = list(),
  family = "gaussian",
  link = NULL,
  weights = NULL,
  aggregate = "auto",
  memory_ceiling_gb = NULL,
  predict_plan = NULL,
  syntax = c("auto", "asreml", "brms"),
  ...
)
```

## Arguments

- formula:

  Either a brms-style two-sided formula such as
  `y ~ x + s(z) + (1 | g)`, or (when `random` and/or `residual` are
  supplied) the ASReml-grammar fixed-effects formula, e.g.
  `yield ~ loc + yearf`. The same formula the fit would receive.

- data:

  A data.frame holding every variable the formula names. Planning reads
  its dimensions and column types, never its values.

- random:

  An optional one-sided ASReml-grammar random-effects formula (e.g.
  `~ gen + gen:loc`). `NULL` (default) for a brms-style `formula`, or
  for a fixed-effects-only ASReml model.

- residual:

  An optional one-sided ASReml-grammar residual structure formula (e.g.
  `~ dsum(~ units | env)`). `NULL` (default) for the ordinary
  independent-residual model.

- backend:

  A single string, one of `"auto"` (the default), `"inla"`, or `"brms"`.
  Chooses the engine the plan reports on, with `"auto"` planning the
  route dispatch would take.

- priors:

  An optional
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  list. Defaults to the uniform-on-SD default the fit would inject.

- known_matrices:

  A named list of structured-covariance matrices referenced by `vm()` or
  `ped()` terms. Empty when the model names none.

- family, link:

  Standard family and link arguments, given as single strings.

- weights:

  An optional numeric vector of observation weights, one per row of
  `data`.

- aggregate:

  `"auto"`, `TRUE`, or `FALSE`, read exactly as on
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).

- memory_ceiling_gb:

  An optional numeric override for the preflight memory ceiling.
  Defaults to the `flexyBayes.preflight_ceiling_gb` option, or to
  `flexyBayes.preflight_ram_fraction` of available RAM.

- predict_plan:

  An optional `list(newdata = ..., chunk_size = ...)` requesting a
  prediction-shape plan. Plan-only, so it never fires
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- syntax:

  Grammar override, one of `"auto"` (default, detected from `formula` /
  `random` / `residual` exactly as
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  does), `"asreml"`, or `"brms"`.

- ...:

  Currently unused, reserved for future plan inputs.

## Value

An `<fb_plan>` classed list. See
[`print.fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_plan.md)
for the one-screen surface,
[`summary.fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/summary.fb_plan.md)
for the verbose dump, and
[`as.data.frame.fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/as.data.frame.fb_plan.md)
for the programmatic-consumer shape.

## Details

Accepts either grammar
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
accepts: a single brms-style two-sided formula
(`y ~ x + s(z) + (1 | g)`), or the ASReml triple (`formula` as the
fixed-effects side, plus `random` / `residual`). Before 0.9.3,
`fb_plan()` was brms-formula-only and silently dropped `random` /
`residual` into its unused `...` when a caller supplied them
ASReml-style, planning a fixed-effects-only model and reporting a
backend the full model would not reach live (FS-22). Both spellings now
route through the same grammar detection
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
uses (`syntax = "auto"` by default), so
`fb_plan(y ~ x, random = ..., residual = ..., data = ...)` and
`flexybayes(fixed = y ~ x, random = ..., residual = ..., data = ..., plan = TRUE)`
plan the identical model. Note the two spellings of the fixed part:
`fb_plan()` takes it positionally as `formula`, the entry point takes it
as `fixed`.

## Examples

``` r
# Planning never fires the backend, so this runs without INLA or
# brms installed.
set.seed(1)
d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))

# brms-style grammar.
p1 <- fb_plan(y ~ x + (1 | g), data = d)
print(p1)
#> == flexyBayes plan ======================================
#>   Will fit:                no  (no active route: auto_no_active_route)
#>   Backend requested:       auto
#>   Backend chosen:          NA
#>   Path:                    auto_no_active_route
#>   Routing policy version:  stage5a_v1
#>   Gate outcome:            accept
#>   Aggregation:             not eligible (continuous_cell_key_data_dependent)
#>   Cov validation policy:   n/a
#>   Representation plan:
#>     x                        -> indexed_fixed_numeric           indexed (0.0 MB); aggregation-eligible
#>     (1 | g)                  -> indexed_random_intercept        indexed (0.0 MB); aggregation-eligible
#>   Rejected routes:
#>     inla     -> backend_not_installed
#>   Memory estimate:         ~ 0.0 MB (preflight)
#>     INLA per-term total:  ~ 0.0 MB (overhead 2.0x)
#>       - (1 | g)                indexed_random_intercept        0.0 MB
#>       - (fixed effects)        fixed_model_matrix              0.0 MB
#>   Representation:          exact
#>   Engine:                  NA
#> ==========================================================

# ASReml-style grammar plans the identical model (FS-22 parity).
p2 <- fb_plan(y ~ x, random = ~ g, data = d)
identical(p1$backend_chosen, p2$backend_chosen)
#> [1] TRUE

# The representation table and the memory estimate are both on the
# plan object without ever fitting.
as.data.frame(p1)[, c("backend_chosen", "gate_outcome")]
#>   backend_chosen gate_outcome
#> 1           <NA>       accept
```
