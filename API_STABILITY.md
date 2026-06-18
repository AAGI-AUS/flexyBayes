# API stability — flexyBayes

`flexyBayes` is in the **0.8.x** development line. Every public export carries
a `lifecycle::badge("experimental")`, so the guarantees below are deliberately
weaker than they will be at v1.0: bug fixes and additions never break callers,
but renamings, default changes, and shape changes are permitted within 0.x
through a one-minor-release `lifecycle::deprecate_warn()` cycle.

This document describes the **lean core**: the mixed-model fitting and
cross-engine triangulation surface. The orchestra-composition layer (surrogate
emulators, ensemble sources, simulator-derived priors, the synthesised
fourth-opinion slot) lives in the companion package **flexyBayesOrchestra** and
is documented there.

## Stability ladder (lifecycle)

| Stage | Meaning for `flexyBayes` |
|---|---|
| **stable** | Behaviour and signature frozen across the 1.x series; breaking changes ride a major bump and a deprecation cycle. *No exports are stable yet.* |
| **maturing** | Behaviour and signature unlikely to change in 0.x; signature additions are non-breaking. |
| **experimental** | Behaviour and signature may change in any 0.x release through a `deprecate_warn()` cycle of at least one minor release. **All current exports.** |
| **deprecated** | Still callable; emits `lifecycle::deprecate_warn()`. Removal scheduled in `NEWS.md`. |
| **superseded** | Still callable, no warning, but a newer entry point is preferred. |

## The entry surface

flexyBayes has two kinds of fitting entry point.

**Universal entries** take a `backend` argument and route the model to an
engine:

| Export | Stage | Notes |
|---|---|---|
| `flexybayes()` | experimental | The asreml-style entry (`fixed`, `random`, `rcov`); also accepts brms-style and greta-style grammar via `syntax = "auto"`. `backend = c("auto", "greta", "inla", "brms", "gretaR")`, default `"auto"`. `"auto"` never routes to Stan/brms. `prior` accepts an `fb_prior()` object. |
| `fb()` | experimental | Literal alias for `flexybayes()`; documented but not promoted. |

**Engine pins** fix one engine and therefore take **no** `backend` argument —
passing a conflicting `backend` raises an `engine_pin_backend_conflict`
refusal:

| Export | Stage | Notes |
|---|---|---|
| `fb_greta()` | experimental | Pins the greta (Hamiltonian Monte Carlo) engine. Also accepts a user-built native greta model (optionally wrapped with `fb_from_greta()` to set `canonical_names`). |
| `fb_inla()` | experimental | Pins the INLA engine (approximate inference — see "Inference semantics"). |
| `fb_brms()` | experimental | Pins the Stan/brms engine via `brms::brm()`. |

**Native-model adapters** lift a model fitted elsewhere into the flexyBayes
object so the shared diagnostics, prediction, and interop methods apply:

| Export | Stage | Notes |
|---|---|---|
| `fb_from_greta()` | experimental | Wrap a native greta model; carries `canonical_names`. |
| `fb_from_brms()` | experimental | Wrap a fitted `brmsfit`. |
| `fb_from_asreml()` | experimental | Wrap a fitted `asreml` object. |

## Cross-engine comparison

| Export | Stage | Notes |
|---|---|---|
| `triangulate()` | experimental | Pairwise comparison of two posteriors from different engines on shared parameters. Metric set for 0.8.x: Wasserstein-1 distance, tail drift, SD ratio, mean difference. *(The earlier cross-engine R-hat-on-means metric was removed — R-hat across independent engines is not a valid convergence statistic; per-fit within-engine R-hat is reported on each fit.)* Metric additions are non-breaking; removals or renamings ride a deprecation cycle. |

## Priors

| Export | Stage | Notes |
|---|---|---|
| `fb_prior()` | experimental | Prior DSL on the standard-deviation scale. Targets: `sigma`, `sd(group)`, `b(name)`, `cor(group)`, `smooth(var)`. Families: `pc`, `normal`, `student_t`, `half_normal`, `half_cauchy`, `cauchy`, `gamma`, `exponential`, `lkj`, `uniform`. The penalised-complexity (PC) family is the cross-engine interlingua (Simpson et al. 2017) for translating a prior across greta, INLA, and brms. |
| `prior_summary()` | experimental | S3 generic returning the resolved-prior view for a fit; reports which resolution path fired (`auto-default`, `user-supplied`, `legacy-scalar`). Methods for every fit subclass. |

## Covariance, engine, and approximation helpers

| Export | Stage | Notes |
|---|---|---|
| `fb_cov()`, `is_fb_cov()` | experimental | Construct / test a structured-covariance specification carried on a random term. |
| `fb_engine()`, `is_fb_engine()` | experimental | Construct / test an engine specification (engine name plus sampler-control options) passable as `backend = fb_engine(...)`. |
| `fb_approx()`, `is_fb_approx()`, `validate_approximation()` | experimental | Construct / test / validate an approximate-scheme request. Approximate schemes are gated: an unregistered scheme is refused with a structured message naming the available exact routes. |

## Planning, streaming, and big-data

| Export | Stage | Notes |
|---|---|---|
| `fb_plan()` | experimental | Returns the dispatch / aggregation plan for a model without fitting (`plan = TRUE` on a universal entry returns the same object). Explains which backend was chosen and why. |
| `flexybayes_stream()` | experimental | Streaming Gaussian sufficient-statistic aggregation for data too large to hold in memory; `fit = FALSE` returns the `<fb_aggregated>` carrier (compression ratio, `K`, `N`) without fitting. |

## Diagnostics and introspection

| Export | Stage | Notes |
|---|---|---|
| `backend_decision()` | experimental | The captured dispatch trace: `backend`, `path`, `gate_checks`, `reason`. Shape stable across 0.8.x. |
| `canonical_names()` | experimental | S3 generic returning the canonical-name registry view for a fit. Methods for every fit subclass. |
| `fb_refusals()` | experimental | The refusal vocabulary as a table (code, message template, since-version). |
| `gretaR_status()` | experimental | Reports whether the dormant `gretaR` R-native engine is detected at run time (see "Backends"). |
| `proceed()`, `cat_code()` | experimental | Companions to the `review_code = TRUE` workflow (inspect, then run, generated engine code). |

## Interoperability contract

flexyBayes fits plug into the wider ecosystem through registered S3 methods,
whose signatures are dictated by the host package's contract and are therefore
as stable as that contract:

- **Draws**: `fb_as_draws_simple()` (a named list of named draw vectors) and the
  `posterior`-compatible path.
- **broom**: `tidy()`, `glance()`, `augment()` — column names follow broom
  convention; new columns are non-breaking.
- **emmeans**: `recover_data()`, `emm_basis()`.
- **marginaleffects**: `get_coef()`, `get_predict()`, `get_vcov()`, `set_coef()`.
- **insight**: `get_data()`.
- **base/stats**: `coef()`, `confint()`, `vcov()`, `predict()`, `fitted()`,
  `residuals()`, `family()`, `formula()`, `logLik()`, `nobs()`, `summary()`,
  `anova()`, `update()`, `plot()`.

`predict.flexybayes()` accepts a `newdata` interface mirroring
`stats::predict()`, with `output_file` / `format = c("auto", "csv", "rds", "fst")`
for chunked output and `allow_new_levels = c("population", "sample", "refuse")`.
Interval semantics are *posterior expected-response* (no residual observation
noise), not posterior-predictive.

## Backends

| Backend | Availability | Inference |
|---|---|---|
| `greta` | greta-dev R-universe (archived from CRAN); needs a working Python/TensorFlow stack at run time | Hamiltonian Monte Carlo (exact up to Monte Carlo error) |
| `INLA` | from its own repository (`Additional_repositories`) | Integrated nested Laplace approximation (**approximate**) |
| `brms` | CRAN; needs a Stan toolchain | Hamiltonian Monte Carlo via Stan |
| `gretaR` | **not a declared dependency** — install it yourself | torch-native MCMC; dormant, opt-in, detected at run time via `gretaR_status()` |

## Inference semantics — read this

flexyBayes is *formula-preserving*: the model emitted to a backend faithfully
represents the formula you wrote. That is distinct from the *inference* being
exact. The greta and brms backends draw from the posterior by MCMC (exact up to
Monte Carlo error, subject to convergence diagnostics). The **INLA backend is
approximate inference** — integrated nested Laplace approximation — even when the
emitted model is a faithful translation of the formula. Treat "formula-preserving"
and "exact inference" as different claims; flexyBayes makes the first everywhere
and the second only on the MCMC backends.

## Default-prior contract

On the auto-default path — when neither an `fb_prior()` nor a legacy
`prior_vc_sd` scalar is supplied:

| Target | Default | Basis |
|---|---|---|
| Residual `sigma` (Gaussian, identity link) | `uniform(0, 5 * sd(y))` | flexyBayes heuristic (weakly-informative bounded SD prior) |
| Residual `sigma` (log link) | `uniform(0, 5 * sd(log(y + 0.5)))` | flexyBayes heuristic |
| Residual `sigma` (logit link) | `uniform(0, 5)` | flexyBayes heuristic |
| `sd(group)` for `simple` / `ide` / `id` random terms | `uniform(0, U)`, `U` by family as above | flexyBayes heuristic |
| `sd(group)` for `vm()`, `ped()` structured-cov terms | `uniform(0, U)`, `U` by family | flexyBayes heuristic |
| `sd(group)` for `at()`, `us()`, `fa()`, `ar1()`, `spl()` terms | legacy `lognormal(0, prior_vc_sd)` (per-form uniform default deferred) | flexyBayes legacy |
| Fixed-effect coefficients | `normal(0, prior_fixed_sd)`, `prior_fixed_sd = 100` | weakly-informative on the natural data scale |

The bounded-uniform-on-SD default is a flexyBayes heuristic in the spirit of
Gelman (2006), which argues for bounded / weakly-informative variance priors over
the conjugate inverse-gamma. Note that Gelman (2006) specifically recommends a
half-Cauchy for the few-groups regime; the bounded uniform is chosen here for its
cross-engine translatability, and the PC prior remains the recommended explicit
choice (`fb_prior(sigma ~ pc(upper = U, prob = p))`) when the number of groups is
small. The legacy scalar bridge (`prior_fixed_sd`, `prior_vc_sd`) preserves the
original `lognormal(0, prior_vc_sd)` semantics verbatim when `prior_vc_sd` is
passed explicitly.

## Deprecation policy

When an experimental export is renamed, restructured, or removed:

1. The old call path emits `lifecycle::deprecate_warn()` for at least one minor
   release.
2. The next minor release moves it to `lifecycle::deprecate_stop()` — the
   function still exists but signals a hard error directing callers to the
   replacement.
3. The minor release after that removes the old export.

For a default-value change, the old default remains available behind an explicit
argument for the rest of the 0.x series.

## Pinning

Production users who need stability should pin to a specific 0.8.x patch via
`renv::snapshot()` until v1.0 lands.

## New in the 0.8.x line (all experimental)

The 0.8.x line adds the exports below. All are **experimental** under the
ladder above; none changes the stability posture of the lean-core
fitting / triangulation surface.

| Export | Added | Notes |
|---|---|---|
| `triangulate_genomic()` / `triangulate_gwas()` | 0.8.0 | Genomic / GWAS cross-engine and field-standard triangulation. |
| `fb_met_summary()` | 0.8.0 | Breeder summary of a greta factor-analytic G×E fit (greta-only). |
| `fb_gblup_cv()` | 0.8.0 | Genomic-prediction accuracy by cross-validation. |
| `fb_gwas()` | 0.8.0 | EMMAX / P3D whole-genome scan. |
| `tidy()` / `glance()` / `augment()` | 0.8.1 | broom-style accessors (re-exported from `generics`). `tidy()` covers all three backends; `glance()` / `augment()` cover greta + brms. |
| `fb_gev()` / `fb_dirichlet()` | 0.8.1 | Generalised-extreme-value and Dirichlet fitters (with `fb_family_*` descriptors). |
| `fb_log_posterior()` | 0.8.2 | Constellation C4 producer; greta is the real producer, brms / INLA honestly abstain. |
| `glance.flexybayes_inla()` / `augment.flexybayes_inla()` | 0.8.3 | Explicit, classed refusals for INLA fits, pointing users to `tidy()`, `summary()`, and `fb_structured_cov()` (an INLA fit previously raised a bare "no applicable method" error). |
