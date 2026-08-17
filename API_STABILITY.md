# API stability -- flexyBayes

`flexyBayes` is in the **0.9.0** development line. Every public export carries
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
| `flexybayes()` | experimental | The asreml-style entry (`fixed`, `random`, `residual`); also accepts brms-style and greta-source grammar via `syntax = "auto"`. `backend = c("auto", "greta", "inla", "brms", "gretaR")`, default `"auto"`. The two quarantined values stay in the vocabulary so that requesting one refuses by name rather than by `match.arg()`. `"auto"` fits on INLA when `lgm_gate()` accepts and INLA is installed, otherwise on brms when brms can represent the model, and otherwise refuses with `auto_no_active_route`. There is no silent fallback and no quarantined engine is ever selected. `prior` accepts an `fb_prior()` object. |
| `fb()` | experimental | Literal alias for `flexybayes()`; documented but not promoted. |

**Engine pins** fix one engine and therefore take **no** `backend` argument --
passing a conflicting `backend` raises an `engine_pin_backend_conflict`
refusal:

| Export | Stage | Notes |
|---|---|---|
| `fb_greta()` | experimental | Pins greta, which is quarantined: the call refuses with `backend_quarantined` before any code is emitted, and a native greta model graph passed to it refuses with `native_greta_fit_quarantined`. The entry point is kept so a re-entry restores one call path rather than a released signature. |
| `fb_inla()` | experimental | Pins the INLA engine (approximate inference -- see "Inference semantics"). |
| `fb_brms()` | experimental | Pins the Stan/brms engine via `brms::brm()`. |

**Sampler arguments.** `seed` and `control` are part of the universal
signature and reach `brms::brm()` unchanged, so a brms posterior is
reproducible from the call alone: `control = list(adapt_delta = 0.95)` is the
route to `adapt_delta` and `max_treedepth`. Left `NULL` they map onto brms's
own defaults (`seed = NA`, `control = NULL`) rather than being passed through,
and `update()` repeats whatever the fit recorded. Both are no-ops on the
deterministic INLA path, which says so once per session.

**Native-model adapters** lift a model fitted elsewhere into the flexyBayes
object so the shared diagnostics, prediction, and interop methods apply:

| Export | Stage | Notes |
|---|---|---|
| `fb_from_greta()` | experimental | Wrap a native greta model graph in the intermediate representation; carries `canonical_names`. The adapter still builds the representation, but fitting the wrapped graph is quarantined with the engine (`native_greta_fit_quarantined`). |
| `fb_from_brms()` | experimental | Wrap a fitted `brmsfit`. |
| `fb_from_asreml()` | experimental | Wrap a fitted `asreml` object. |

## Cross-engine comparison

| Export | Stage | Notes |
|---|---|---|
| `triangulate()` | experimental | Pairwise comparison of two posteriors from different engines on shared parameters, behind three gates. A model fingerprint (formula triple, family and link, data digest, recorded priors) refuses mismatched fits with `triangulate_incomparable_fits`. A per-fit diagnostics gate returns status `inconclusive` instead of a parameter verdict. A matched-prior gate marks a parameter `not_compared` when the two fits do not record the same prior for it. Metric set for 0.9.x: mean difference, SD ratio, quantile (tail) drift, Wasserstein-1 distance, and the SD-scaled forms of the last two. Overall `status` is `concordant`, `discordant`, or `inconclusive`, and the thresholds are documented heuristics on a gated overlap rather than calibrated tests. *(The earlier cross-engine R-hat-on-means metric was removed -- R-hat across independent engines is not a valid convergence statistic, and per-fit within-engine R-hat is reported on each fit.)* Metric additions are non-breaking; removals or renamings ride a deprecation cycle. |

## Priors

| Export | Stage | Notes |
|---|---|---|
| `fb_prior()` | experimental | Prior DSL on the standard-deviation scale. Targets: `sigma`, `sd(group)`, `b(name)`, `cor(group)`, `smooth(var)`. Families: `pc`, `normal`, `student_t`, `half_normal`, `half_cauchy`, `cauchy`, `gamma`, `exponential`, `lkj`, `uniform`. The penalised-complexity (PC) family is the cross-engine interlingua (Simpson et al. 2017) for translating one prior specification into each engine's own parameterisation. |
| `prior_summary()` | experimental | S3 generic returning the resolved-prior view for a fit, reporting which resolution path fired (auto-default bounded uniform on SD, user-supplied `fb_prior()`, or the legacy scalar bridge). Methods for every fit subclass. |

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
| `flexybayes_stream()` | experimental | Chunked sufficient-statistic aggregation for data too large to hold in memory, for `family = "gaussian"`, `"binomial"` or `"poisson"` (the count families take `trials` / `exposure`). INLA is the default and the only engine with an aggregated emit, and `backend = "greta"` refuses with `backend_quarantined`. `fit = FALSE` returns the `<fb_aggregated>` carrier (compression ratio, `K`, `N`) without fitting. |

## Diagnostics and introspection

| Export | Stage | Notes |
|---|---|---|
| `backend_decision()` | experimental | The captured dispatch trace, eight fields: `backend`, `path`, `gate_checks`, `reason`, `preflight_summary`, `representation_plan`, `rejected_routes`, `routing_policy_version`. Shape stable across 0.9.x. |
| `canonical_names()` | experimental | S3 generic returning the canonical-name registry view for a fit. Methods for every fit subclass. |
| `fb_refusals()` | experimental | The refusal vocabulary as a table (code, message template, since-version). The registry is the count -- no document quotes a number of refusal codes. |
| `fb_backend_status()` | experimental | Which engines are installed and usable in the current session, and why an unusable one is unusable. |
| `fb_structured_cov()` | experimental | The identified covariance of a factor-analytic structured-covariance term. |
| `gretaR_status()` | experimental | Reports whether the quarantined `gretaR` engine is present at run time, and the reason it is not dispatchable (see "Backends"). |
| `proceed()`, `cat_code()` | experimental | Companions to the `review_code = TRUE` workflow (inspect, then run, generated engine code). brms is the only engine with a code slot, so the deferred-execution token is available under `backend = "brms"` and under `"auto"`, and refuses under `backend = "inla"` with `review_code_backend_unsupported`. |

## Interoperability contract

Every fit carries its engine's class first and `"flexybayes"` second
(`c("flexybayes_brms", "flexybayes", "list")`,
`c("flexybayes_inla", "flexybayes", "list")`), so a method written for the
shared class reaches both engines and an engine-specific method still wins.

flexyBayes fits plug into the wider ecosystem through registered S3 methods,
whose signatures are dictated by the host package's contract and are therefore
as stable as that contract:

- **Draws**: `fb_as_draws_simple()` (a named list of named draw vectors) and the
  `posterior`-compatible path.
- **broom**: `tidy()`, `glance()`, `augment()` -- column names follow broom
  convention, and new columns are non-breaking. `tidy()` has a method for
  every fit subclass. `glance()` and `augment()` describe a sampled fit, so
  on an INLA fit they refuse by name and point at `tidy()`, `summary()` and
  `predict()` instead.
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

Where a method has no meaning on an engine it states that rather than
returning a number. `logLik()` on an INLA fit reports that INLA returns a
marginal log-likelihood and yields `NA` with that message, and a fit carrying
neither a response vector nor a recorded family refuses instead of letting
`AIC()` or `anova()` report a comparison that was never made.

## Backends

Two engines are active. `greta` and `gretaR` are quarantined: their registry
descriptors and emit code are retained as re-entry candidates, `backend =
"auto"` never selects them, and an explicit request refuses with
`backend_quarantined`. Installing them adds no fitting capability.

| Backend | Status | Availability | Inference |
|---|---|---|---|
| `INLA` | active | from its own repository (`Additional_repositories`) | Integrated nested Laplace approximation (**approximate**) |
| `brms` | active | CRAN; needs a Stan toolchain | Hamiltonian Monte Carlo via Stan |
| `greta` | quarantined | greta-dev R-universe (archived from CRAN); needs a working Python/TensorFlow stack at run time | Hamiltonian Monte Carlo -- refuses rather than fitting |
| `gretaR` | quarantined | **not a declared dependency** -- install it yourself | torch-native MCMC; detected at run time via `gretaR_status()`, refuses rather than fitting |

## Inference semantics -- read this

flexyBayes is *formula-preserving*: the model emitted to a backend faithfully
represents the formula you wrote. That is distinct from the *inference* being
exact. The brms backend draws from the posterior by MCMC (exact up to
Monte Carlo error, subject to convergence diagnostics). The **INLA backend is
approximate inference** -- integrated nested Laplace approximation -- even when the
emitted model is a faithful translation of the formula. Treat "formula-preserving"
and "exact inference" as different claims; flexyBayes makes the first everywhere
and the second only on the MCMC backends.

## Default-prior contract

On the auto-default path -- when neither an `fb_prior()` nor a legacy
`prior_vc_sd` scalar is supplied:

| Target | Default | Basis |
|---|---|---|
| Residual `sigma` (Gaussian, identity link) | `uniform(0, 5 * sd(y))` | flexyBayes heuristic (weakly-informative bounded SD prior) |
| Residual `sigma` (log link -- Poisson, negative binomial, gamma) | `uniform(0, 3)` on the log scale | flexyBayes heuristic |
| Residual `sigma` (logit link -- binomial, beta) | `uniform(0, 5)` on the logit scale | flexyBayes heuristic |
| `sd(group)` for `simple` / `ide` / `id` random terms | `uniform(0, U)`, `U` by family as above | flexyBayes heuristic |
| `sd(group)` and `sd(slope)` for an uncorrelated slope `(x \|\| g)` | `uniform(0, U)`, `U` by family | flexyBayes heuristic |
| `sd(group)` for `vm()`, `ped()` structured-cov terms | `uniform(0, U)`, `U` by family | flexyBayes heuristic |
| `sd(group)` for a `gen:env` interaction intercept, and the per-level SDs of `diag()` / `idh()` / `at()` and `us()` | `uniform(0, U)`, `U` by family | flexyBayes heuristic (added at 0.9.0) |
| `us(f):g` level correlations | brms's own LKJ | recorded as engine-default |
| AR1 field hyperparameters (`ar1()`, `ar1(row):ar1(col)`) | INLA's own hyperpriors | recorded as engine-default |
| Any other random-term type (a multi-way interaction, `spl()`, `fa()`) | the engine's own default | recorded as engine-default |
| Fixed-effect coefficients | `normal(0, prior_fixed_sd)`, `prior_fixed_sd = 100` | weakly-informative on the natural data scale |

"Recorded as engine-default" is a contract, not a shrug: the parameter and the
reason are written into the fit's prior provenance, so `triangulate()`'s
matched-prior gate reports it as `not_compared` rather than comparing two
engines that were never asked the same question.

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
2. The next minor release moves it to `lifecycle::deprecate_stop()` -- the
   function still exists but signals a hard error directing callers to the
   replacement.
3. The minor release after that removes the old export.

For a default-value change, the old default remains available behind an explicit
argument for the rest of the 0.x series.

## Pinning

Production users who need stability should pin to a specific 0.9.x patch via
`renv::snapshot()` until v1.0 lands.

## Added on the 0.8.x line (all experimental)

The 0.8.x line added the exports below. All are **experimental** under the
ladder above; none changes the stability posture of the lean-core
fitting / triangulation surface.

| Export | Added | Notes |
|---|---|---|
| `triangulate_genomic()` / `triangulate_gwas()` | 0.8.0 | Genomic / GWAS cross-engine and field-standard triangulation. |
| `fb_met_summary()` | 0.8.0 | Breeder summary of a factor-analytic GxE fit. It is computed from realised factor-analytic effects, which only greta produced, so with greta quarantined it abstains on an INLA or brms fit with `met_summary_not_available` and names what an active engine reports instead (`summary()` for the components, `brms::VarCorr()` for a `diag()` or `us()` covariance). |
| `fb_gblup_cv()` | 0.8.0 | Genomic-prediction accuracy by cross-validation. |
| `fb_gwas()` | 0.8.0 | EMMAX / P3D whole-genome scan. |
| `tidy()` / `glance()` / `augment()` | 0.8.1 | broom-style accessors (re-exported from `generics`). `tidy()` has a method for every fit subclass. `glance()` / `augment()` describe a sampled fit and refuse by name on an INLA fit. |
| `fb_gev()` / `fb_dirichlet()` | 0.8.1 | Generalised-extreme-value and Dirichlet fitters (with `fb_family_*` descriptors). |
| `fb_log_posterior()` | 0.8.2 | Constellation C4 producer. greta was the only engine that evaluated the log density, so on every other object class the method abstains with a typed message rather than returning a number. |
| `glance.flexybayes_inla()` / `augment.flexybayes_inla()` | 0.8.3 | Explicit, classed refusals for INLA fits, pointing users to `tidy()`, `summary()`, and `fb_structured_cov()` (an INLA fit previously raised a bare "no applicable method" error). |
