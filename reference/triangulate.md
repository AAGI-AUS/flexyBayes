# Cross-engine posterior triangulation

Compute per-parameter posterior comparison metrics across two Bayesian
fits produced on different backends (INLA and the brms / Stan
passthrough).

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

  A fit object carrying a `fb_as_draws_simple` method, which means
  `flexybayes_brms` from the brms path or `flexybayes_inla` from the
  INLA path.

- fit_b:

  A second fit object of the same kind, typically produced by the other
  backend on the same model and data.

- name_map:

  An optional named character vector or list mapping `fit_b`'s parameter
  names (left) to canonical names matching `fit_a` (right).

- transform_a, transform_b:

  An optional named list of functions, each taking a numeric vector of
  posterior draws and returning a numeric vector of the same length.
  Names key parameters in `fit_a` / `fit_b` (using each fit's *original*
  parameter names). Common use: pass
  `transform_b = list("Precision for g" = function(x) 1 / sqrt(x))` to
  convert INLA's precision draws to standard-deviation scale so they
  line up with brms's `sd_g__Intercept`. Optional, and normally
  unnecessary –
  [`canonical_names()`](https://aagi-aus.github.io/flexyBayes/reference/canonical_names.md)
  supplies the common transforms.

- n_samples:

  A single integer giving the number of posterior samples to draw for
  fits whose extractor needs sampling, such as INLA via
  `INLA::inla.posterior.sample`.

- data_independence:

  A single logical declaring whether the two fits were built on
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

A `triangulate_result` S3 object (a list). Key fields: `metrics` (a
data.frame with one row per common parameter and the columns `param`,
`mean_a`, `mean_b`, `mean_diff`, `sd_a`, `sd_b`, `sd_ratio`,
`q025_diff`, `q975_diff`, `wasserstein_1`, `mean_shift_sd`, `w1_sd`,
`verdict`, `reason`), `status` (one of `"concordant"`, `"discordant"`,
`"inconclusive"`), `status_reasons` (character, why a status is
inconclusive or discordant), `comparability` (the fingerprint verdict),
`diagnostics` (the per-fit diagnostic summaries), `common` (character),
`only_a`, `only_b`, `n_common`, `n_compared`, `source_a`, `source_b`.

## Details

For each parameter present in both fits (post `name_map`), the returned
table reports each fit's posterior mean and SD, the difference in means,
the ratio of SDs, the Q2.5 and Q97.5 differences (tail drift), and the
1D empirical Wasserstein-1 distance. Parameters present in only one fit
are reported in `only_a` / `only_b`.

`triangulate()` is a diagnostic on a gated overlap, not the package's
validation contract. Agreement is not validation: the parser, the
intermediate representation, the prior interlingua and the data
preparation are shared, so two engines can agree on a mistranslated
model. What agreement does support is the narrower claim that neither
engine's approximation is distorting the posterior of the compared
parameters, and it supports even that only inside the three gates
described below.

## The three gates

**Comparability.** Fits made by
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
carry a model fingerprint: the canonical formula triple, the family and
link, the data dimensions, column names and a content digest, and the
prior recorded for each variance component. Two fingerprints that
disagree raise a typed refusal
(`flexybayes_refusal_triangulate_incomparable_fits`) naming the first
element that differs. There is no override argument: comparing two
models is a different operation, and it is not this one. A fit built
outside the package carries no fingerprint, so comparability cannot be
verified and the overall status is `inconclusive`.

**Diagnostics.** A brms fit passes when every model parameter has R-hat
at or below 1.01 and bulk effective sample size at or above 400, with no
divergent transitions; an INLA fit passes its own numerical-confirm gate
(mode status and a finite marginal likelihood). A failing fit yields
status `inconclusive` and no parameter verdicts: disagreement between an
unconverged fit and a converged one is not a finding.

**Matched priors.** A variance component (`sigma`, `sd_<group>`,
`cor_<group>`, `rho_*`) is compared only when both fits record the same
prior for it – the shared default this package injects before dispatch,
or an identical explicit
[`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md).
Components outside that record are reported as `not_compared` with the
reason, most often that the term is outside the default-prior walker and
each engine chose its own hyperprior. Fixed-effect coefficients are
compared without a prior match: both engines put a vague normal on them,
and the fixed structure itself is already part of the fingerprint.

## Thresholds and what they are not

A parameter is called concordant when the absolute mean difference is at
most 0.1 posterior SD, the Wasserstein-1 distance is at most 0.1
posterior SD, and the SD ratio lies in \[0.9, 1.1\]; the reference scale
is the larger of the two posterior SDs. The overall status is
`discordant` when any compared parameter is discordant, `concordant`
when all are, and `inconclusive` when a gate blocked the comparison or
nothing was left to compare.

These cut-offs are **heuristics**, in the sense Gelman et al. (2020) use
for cross-implementation checks in the Bayesian workflow: they are
chosen to flag differences a reader should look at, not derived from a
sampling distribution, and they carry no error rate. A `concordant`
verdict is not a calibration statement. flexyBayes ships no
simulation-based calibration of its own fits (Talts et al., 2018), so
whether either engine is calibrated for a given model and dataset
remains an open question and future work.

Read `discordant` as "look at this", not as "one engine is wrong".
Wasserstein-1 compares the whole marginal, so two posteriors whose means
and SDs agree can still be flagged on shape alone – which is the
expected signature of a Laplace marginal for a variance component next
to an HMC one. All three metrics are estimated from finite draws, and
the Wasserstein-1 estimate in particular is biased upward at small
`n_samples`; raise `n_samples` and re-run before reading a borderline
flag as a finding.

Backends use different parameter naming conventions; supply
`name_map = c(<fit_b name> = <canonical name>, ...)` to align them.
Without alignment, only literal-match parameter names are compared.

Backends also use different *parameter scales* for variance components –
INLA reports precision (`Precision for g`), brms reports standard
deviation (`sd_g__Intercept`). Supply `transform_a` / `transform_b` –
named lists of one-argument functions keyed by parameter name – to put
the two posteriors on a common scale before comparison. Transforms are
applied first; then `name_map` aligns the (already-transformed) fit_b
names to fit_a's canonical names. Names in `transform_b` therefore refer
to fit_b's *original* parameter names, not the post-`name_map` canonical
names.

## Aggregated fits and matched priors

When one of the inputs is an aggregated-gaussian fit (cell-level
sufficient statistics rather than the per-row likelihood), the
posteriors being compared are only directly comparable when the two fits
share priors. The aggregated path combines the cell-mean likelihood with
a precision prior carrying a closed-form correction that absorbs the
within-cell sum-of-squares. Under the *legacy scalar* prior this
recovers the per-row posterior to numerical precision, so the aggregated
fit is tagged `prior_parametrization = "per_row_equivalent"` (visible in
the aggregated [`print()`](https://rdrr.io/r/base/print.html) /
[`summary()`](https://rdrr.io/r/base/summary.html) and in
[`canonical_names()`](https://aagi-aus.github.io/flexyBayes/reference/canonical_names.md)).
Under the package's own auto-default it is tagged `"package_default"`,
which names the prior without claiming that equivalence. When an
explicit prior is supplied the fit is tagged `"custom"`: the equivalence
against a *default-prior* per-row fit no longer holds, and on the
aggregated INLA path the observation-precision prior is not plumbed
through, so a custom residual prior is silently not applied there.
Before reading the agreement metrics on a custom-prior aggregated fit,
confirm both inputs carry the same prior with
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md).

## References

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
Yao, Y., et al. (2020). Bayesian workflow. *arXiv preprint*
arXiv:2011.01808.

Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A.
(2018). Validating Bayesian inference algorithms with simulation-based
calibration. *arXiv preprint* arXiv:1804.06788.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., et al. (2021).
Rank-normalization, folding, and localization: an improved R-hat for
assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667–718.

## Examples

``` r
# Live INLA posterior sampling can fail in restricted-process
# check environments (the `inla.posterior.sample` parallelism
# check trips); the example uses `\dontrun{}` deliberately. On
# an interactive install with INLA + brms + sn available it
# runs in a few seconds.
if (FALSE) { # \dontrun{
if (requireNamespace("brms", quietly = TRUE) &&
    requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(40), x = rnorm(40),
                  g = factor(rep(1:5, 8)))
  fit_b <- fb(y ~ x + (1 | g), data = d, backend = "brms",
              n_samples = 500, warmup = 500, chains = 2,
              verbose = FALSE, mcmc_verbose = FALSE)
  fit_i <- fb(y ~ x + (1 | g), data = d, backend = "inla",
              verbose = FALSE)
  # canonical_names() aligns the names and the precision-to-SD
  # transform, so no name_map is needed in the common cases.
  print(triangulate(fit_b, fit_i))
}
} # }
```
