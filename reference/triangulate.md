# Cross-engine posterior triangulation

Compute per-parameter posterior comparison metrics across two Bayesian
fits produced via different backends (greta / INLA / brms passthrough).
The signature feature of flexyBayes.

## Usage

``` r
triangulate(
  fit_a,
  fit_b,
  name_map = NULL,
  transform_a = NULL,
  transform_b = NULL,
  n_samples = 1000L,
  data_independence = NA
)
```

## Arguments

- fit_a:

  a fit object with a `fb_as_draws_simple` method (e.g., `flexybayes`
  from greta, `flexybayes_inla` from INLA).

- fit_b:

  a second fit object, typically from a different backend.

- name_map:

  named character vector or list mapping fit_b's parameter names (left)
  to canonical names matching fit_a (right). Optional.

- transform_a, transform_b:

  named list of functions, each taking a numeric vector of posterior
  draws and returning a numeric vector of the same length. Names key
  parameters in `fit_a` / `fit_b` (using each fit's *original* parameter
  names). Common use: pass
  `transform_b = list("Precision for g" = function(x) 1 / sqrt(x))` to
  convert INLA's precision draws to standard-deviation scale so they
  line up with greta's `sigma_g`. Optional.

- n_samples:

  integer: number of posterior samples to draw for fits whose extractor
  needs sampling (e.g., INLA via `INLA::inla.posterior.sample`).

- data_independence:

  single logical declaring whether the two fits were built on
  independently-sourced data. `triangulate()` measures inter-fit
  *agreement*, and the backend-independence registry certifies code (not
  data) independence – so if both fits share the same upstream data, a
  fabricated data fact is common-mode and their agreement cannot detect
  it. `TRUE` declares the data independent (no caveat); `FALSE` (same
  data) or `NA` (the default, undeclared) attach a
  `shared_upstream_caveat` field to the result, surfaced prominently by
  the print method, so agreement is never silently mistaken for
  corroboration of a shared data fact (Independent Oracle Principle).

## Value

a `triangulate_result` S3 object (list). Key fields: `metrics`
(data.frame, one row per common parameter), `common` (character),
`only_a`, `only_b`, `n_common`, `source_a`, `source_b`.

## Details

For each parameter present in both fits (post `name_map`), the returned
table reports: posterior means, posterior SDs, Q2.5 / Q97.5 differences
(tail drift), Wasserstein-1 distance (1D empirical), the SD ratio, and
an R-hat-on-means scalar that compares between-engine vs within-engine
variance. Parameters present in only one fit are reported in `only_a` /
`only_b`.

Backends use different parameter naming conventions; supply
`name_map = c(<fit_b name> = <canonical name>, ...)` to align them.
Without alignment, only literal-match parameter names are compared.

Backends also use different *parameter scales* for variance components –
INLA reports precision (`Precision for g`), greta reports standard
deviation (`sigma_g`). Supply `transform_a` / `transform_b` – named
lists of one-argument functions keyed by parameter name – to put the two
posteriors on a common scale before comparison. Transforms are applied
first; then `name_map` aligns the (already-transformed) fit_b names to
fit_a's canonical names. Names in `transform_b` therefore refer to
fit_b's *original* parameter names, not the post-`name_map` canonical
names.

## Matched priors

When one of the inputs is an aggregated-gaussian fit (cell-level
sufficient statistics rather than the per-row likelihood), the
posteriors being compared are only directly comparable when the two fits
share priors. The aggregated path combines the cell-mean likelihood with
a precision prior carrying a closed-form correction that absorbs the
within-cell sum-of-squares; under the *default* prior this recovers the
per-row posterior to numerical precision, so the aggregated fit is
tagged `prior_parametrization = "per_row_equivalent"` (visible in the
aggregated [`print()`](https://rdrr.io/r/base/print.html) /
[`summary()`](https://rdrr.io/r/base/summary.html) and in
[`canonical_names()`](https://aagi-aus.github.io/flexyBayes/reference/canonical_names.md)).
When an explicit prior is supplied the fit is tagged `"custom"`: the
equivalence against a *default-prior* per-row fit no longer holds, and
on the aggregated INLA path the observation-precision prior is not
plumbed through, so a custom residual prior is silently not applied
there. Before reading the agreement metrics on a custom-prior aggregated
fit, confirm both inputs carry the same prior with
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md).

## Examples

``` r
# Live INLA posterior sampling can fail in restricted-process
# check environments (the `inla.posterior.sample` parallelism
# check trips); the example uses `\dontrun{}` deliberately. On
# an interactive install with greta + INLA + sn available it
# runs in a few seconds.
if (FALSE) { # \dontrun{
if (requireNamespace("greta", quietly = TRUE) &&
    requireNamespace("INLA",  quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(40), x = rnorm(40),
                  g = factor(rep(1:5, 8)))
  fit_g <- fb(y ~ x + (1 | g), data = d, backend = "greta",
              n_samples = 100, warmup = 100, chains = 1,
              verbose = FALSE)
  fit_i <- fb(y ~ x + (1 | g), data = d, backend = "inla",
              verbose = FALSE)
  prec_to_sd <- function(x) 1 / sqrt(x)
  tri <- triangulate(
    fit_g, fit_i,
    transform_b = list(
      "Precision for g"                         = prec_to_sd,
      "Precision for the Gaussian observations" = prec_to_sd
    ),
    name_map = c(
      "(Intercept):1"                           = "mu_atg",
      "x:1"                                     = "beta_x",
      "Precision for g"                         = "sigma_g",
      "Precision for the Gaussian observations" = "sigma_e_atg"
    )
  )
  print(tri)
}
} # }
```
