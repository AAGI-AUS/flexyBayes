# API stability -- flexyBayes

`flexyBayes` is in the **0.15.0** development line. Every public export is at
the experimental `lifecycle` stage, so the guarantees below are deliberately
weaker than they will be at v1.0: bug fixes and additions never break callers,
but renamings, default changes, and shape changes are permitted within 0.x
through a deprecation cycle of at least one minor release, in which the old
call warns with a condition of class `flexybayes_deprecated`.

This document describes the **lean core**: the mixed-model fitting
surface. The orchestra-composition layer (surrogate
emulators, ensemble sources, simulator-derived priors, the synthesised
fourth-opinion slot) lives in the companion package **flexyBayesOrchestra** and
is documented there.

## Stability ladder (lifecycle)

| Stage | Meaning for `flexyBayes` |
|---|---|
| **stable** | Behaviour and signature frozen across the 1.x series; breaking changes ride a major bump and a deprecation cycle. *No exports are stable yet.* |
| **maturing** | Behaviour and signature unlikely to change in 0.x; signature additions are non-breaking. |
| **experimental** | Behaviour and signature may change in any 0.x release through a deprecation cycle of at least one minor release. **All current exports.** |
| **deprecated** | Still callable; warns with a base R condition of class `flexybayes_deprecated` that names what to call. Removal scheduled in `NEWS.md`. |
| **superseded** | Still callable, no warning, but a newer entry point is preferred. |

## The entry surface

flexyBayes has two kinds of fitting entry point.

**Universal entries** take a `backend` argument and route the model to an
engine:

| Export | Stage | Notes |
|---|---|---|
| `flexybayes()` | experimental | The asreml-style entry (`fixed`, `random`, `residual`); also accepts brms-style grammar via `syntax = "auto"`. `backend = c("auto", "inla", "brms")`, default `"auto"`. Naming a withdrawn or otherwise unrecognised backend refuses by name with `unknown_backend` rather than by `match.arg()`. `"auto"` fits on INLA when `lgm_gate()` accepts and INLA is installed, otherwise on brms when brms can represent the model, and otherwise refuses with `auto_no_active_route`. Nothing is translated silently. `auto` does change engine in one case after starting -- an INLA program failure is caught and re-routed to brms -- and that is reported by a message, never silently; the fit's class and `backend` field record which engine produced it. `prior` accepts an `fb_prior()` object. |
| `fb()` | deprecated in 0.13.0 | Alias for `flexybayes()` with the same arguments: it warns and forwards the call. Call `flexybayes()`. |

**Engine pins** fix one engine and therefore take **no** `backend` argument --
passing a conflicting `backend` raises an `engine_pin_backend_conflict`
refusal:

| Export | Stage | Notes |
|---|---|---|
| `fb_inla()` | deprecated in 0.13.0 | Pins the INLA engine (approximate inference -- see "Inference semantics"). Call `flexybayes(..., backend = "inla")`. |
| `fb_brms()` | deprecated in 0.13.0 | Pins the Stan/brms engine via `brms::brm()`. Call `flexybayes(..., backend = "brms")`. |

**Sampler arguments.** `seed` and `control` are part of the universal
signature and reach `brms::brm()` unchanged, so a brms posterior is
reproducible from the call alone: `control = list(adapt_delta = 0.95)` is the
route to `adapt_delta` and `max_treedepth`. Left `NULL` they map onto brms's
own defaults (`seed = NA`, `control = NULL`) rather than being passed through,
and `update()` repeats whatever the fit recorded. Both are no-ops on the
deterministic INLA path, which says so once per session.

**Grammar adapters** ingest a model *specification* written in another
package's syntax and return the flexyBayes internal representation
(`fb_terms`), which the entry points then fit. They do not take a fitted
object: pass the formula and data, not a `brmsfit` or an `asreml` fit.

| Export | Stage | Notes |
|---|---|---|
| `fb_from_brms()` | deprecated in 0.13.0 | Compatibility parser for brms-style syntax. Use `fb_plan(..., syntax = "brms")` to inspect a model without fitting. |
| `fb_from_asreml()` | deprecated in 0.13.0 | Compatibility parser for ASReml-style syntax. Use `fb_plan(..., syntax = "asreml")`. |

## Cross-engine comparison

| Export | Stage | Notes |
|---|---|---|
| `triangulate()` | removed in 0.14.0 | Internal from 0.11.0 and removed in 0.14.0, together with `triangulate_genomic()`, `triangulate_gwas()` and the refusal code `triangulate_incomparable_fits`. Nothing in the package called it, and agreement between two engines fitted from one parsed model cannot detect a translation error both share. The canonical draws it compared remain available on each fit through `posterior::as_draws_df()` and `canonical_names()`. |

## Priors

| Export | Stage | Notes |
|---|---|---|
| `fb_prior()` | experimental | Prior DSL on the standard-deviation scale. Targets: `sigma`, `sd(group)`, `b(name)`, `cor(group)`, `smooth(var)`. Families: `pc`, `normal`, `student_t`, `half_normal`, `half_cauchy`, `cauchy`, `gamma`, `exponential`, `lkj`, `uniform`. The penalised-complexity (PC) family is the cross-engine interlingua (Simpson et al. 2017) for translating one prior specification into each engine's own parameterisation. |
| `prior_summary()` | experimental | S3 generic returning the resolved-prior view for a fit, reporting which resolution path fired (auto-default bounded uniform on SD, user-supplied `fb_prior()`, or the legacy scalar bridge). Methods for every fit subclass. |

## Covariance, engine, and approximation helpers

| Export | Stage | Notes |
|---|---|---|
| `fb_cov()` | experimental | Construct a structured-covariance specification carried on a random term. |
| `is_fb_cov()` | deprecated in 0.13.0 | Test for that class. Use `inherits(x, "fb_cov")`. |
| `fb_engine()` | deprecated in 0.13.0 | Compatibility constructor. Pass `backend`, `chains`, `n_samples`, `warmup` and `control` directly to `flexybayes()`; `fb_select(engine = )` accepts a named sampler-settings list. |
| `is_fb_engine()` | deprecated in 0.13.0 | Test for that class. Use `inherits(x, "fb_engine")`. |
| `is_fb_approx()` | deprecated in 0.13.0 | Test for the `fb_approx` class. Use `inherits(x, "fb_approx")`. |
| `fb_approx()`, `validate_approximation()` | experimental | Construct / validate an approximate-scheme request. Approximate schemes are gated: an unregistered scheme is refused with a structured message naming the available exact routes. |

## Planning, streaming, and big-data

| Export | Stage | Notes |
|---|---|---|
| `fb_plan()` | experimental | Returns the dispatch / aggregation plan for a model without fitting (`plan = TRUE` on a universal entry returns the same object). Explains which backend was chosen and why. |
| `fb_select()` | experimental | Selection over a fitted model with an explicit workflow axis. `workflow = "bayes"` fits a posterior reference model; `workflow = "penalized_optimizer"` is a frequentist penalised-likelihood path (`glmnet` for ridge/lasso/elastic net and the package proximal provider for grouped families); `workflow = "unpenalized_benchmark"` is generalised least squares with plugged-in covariance. Every result records workflow, provider, prior and decision separately, so no optimiser is labelled Bayesian. Fold-local preprocessing, candidate availability, support history and score-contract hashes are recorded. |
| `fb_select_fold_data()` | experimental | Each fold's whitened and protected-removed data from `fb_select(keep_fold_data = TRUE)`, with its held-out covariance and variance components. |
| `fb_select_path()` | experimental | The penalised path of a kept `fb_select()` run solved again on its fold data, at the run's elastic-net alpha or another. |
| `fb_select_refit()` | experimental | The unpenalised GLS refit of a named column set on each fold of a kept `fb_select()` run, with held-out cluster scores. |
| `fb_variance_draws()` | experimental | Joint posterior draws of an INLA fit's variance components, under the names and on the scale of `summary(fit)$varcomp`. |
| `fb_refit()` | experimental | The model of a fit refitted on new data, with added fixed columns, extra brms prior rows and Stan code, a formula hook, or a reduced INLA readout. |
| `fb_marginal_cov()` | experimental | The marginal covariance of the response of a fit over given rows: one block per variance component, with its design, structure and Kronecker factors where it is separable, the residual, and any known variance or covariance. `as.matrix()` assembles it dense. |
| `fb_penalty()` | experimental | Declare the joint shrinkage prior over the candidate block: `target ~ family(...)` formulas resolved in written order against the model matrix. Families: `ridge`, `laplace`, `elastic_net`, `regularized_horseshoe`, `r2d2`, and `none()` to exempt coordinates. Targets: `term(everything())`, `term(x1, x2)`, `term(c("x1", "x2"))`, `term(matches("^SLGA_"))`. |
| `fb_resamples()` | experimental | Declare the resampling scheme a selection is scored on: folds on original units before any transformation, frozen factor levels, and an unseen level in an assessment fold refused rather than dropped. |
| `fb_decision()` | experimental | Declare how the selection decision is made: the method (`projection`, `full`, `posterior_threshold`), the metric, the rule (`one_se`, `min`, `tolerance`) and the search. Validated at construction, so a rule missing the number it needs refuses before anything is fitted. |
| `flexybayes_stream()` | experimental | Chunked sufficient-statistic aggregation for data too large to hold in memory, fitted for `family = "binomial"` or `"poisson"` (which take `trials` / `exposure`). INLA is the default and the only engine with an aggregated emit. `fit = FALSE` returns the `<fb_aggregated>` carrier (compression ratio, `K`, `N`) without fitting, for a Gaussian family too. **Changed in 0.12.0:** a Gaussian fit refuses with `aggregation_residual_prior_not_composed` before any chunk is read, because the aggregated Gaussian likelihood did not carry the residual prior; fit a Gaussian model per row. An aggregated or streamed fit reports no marginal likelihood (`aggregated_evidence_unavailable`). |

## Diagnostics and introspection

| Export | Stage | Notes |
|---|---|---|
| `backend_decision()` | experimental | The captured dispatch trace, eight fields: `backend`, `path`, `gate_checks`, `reason`, `preflight_summary`, `representation_plan`, `rejected_routes`, `routing_policy_version`. Shape stable across 0.9.x. |
| `canonical_names()` | experimental | S3 generic returning the canonical-name registry view for a fit. Methods for every fit subclass. |
| `fb_refusals()` | experimental | The refusal vocabulary as a table (code, package, description, since-version). The registry is the count -- no document quotes a number of refusal codes. |
| `fb_register_refusals()` | experimental | For packages that extend flexyBayes: adds the package's refusal codes, from its `.onLoad()`, to its own partition of the registry. A code another package or flexyBayes already registered is refused. |
| `fb_refuse()` | experimental | Raises a registered refusal as a `flexybayes_refusal` condition carrying the code and the package that registered it. |
| `fb_isolated_lapply()` | experimental | For packages that extend flexyBayes: `lapply()` with each call in a separate R process that starts from the session's settings and random-number state; values, messages and warnings come back in job order, and a process that dies is a refusal. |
| `fb_backend_status()` | experimental | Which engines are installed and usable in the current session, and why an unusable one is unusable. |
| `fb_structured_cov()` | withdrawn at 0.10.0 | The identified covariance of a factor-analytic structured-covariance term. No active engine emits an `fa()` term, so it reported nothing on every reachable fit and is no longer exported. |
| `proceed()`, `cat_code()` | deprecated in 0.13.0 | Companions to the `review_code = TRUE` workflow (inspect, then run, generated engine code). brms is the only engine with a code slot, so the deferred-execution token is available under `backend = "brms"` and under `"auto"`, and refuses under `backend = "inla"` with `review_code_backend_unsupported`. The generated code is the token's `code` element, and `flexybayes(..., return_code = TRUE)` returns it; the fit comes from calling `flexybayes()` with the same arguments and `review_code = FALSE`. |
| `fb_complete_grid()` | experimental | Completes a design grid (`index = ~ row * col`) before fitting, so a structured covariance is built over the intended index set rather than over the observed rows only. |
| `genomic_summary()` | experimental | Turns the posterior draws of a relationship-matrix term into genomic quantities for a fitted model. |
| `ranef()` | experimental | S3 generic returning random-effect predictions -- posterior summaries of the group-level effects -- from a fit. |
| `vpredict()` | experimental | S3 generic returning a posterior functional of the variance components -- a heritability, a variance ratio, a genetic correlation -- named as an expression in the canonical component names. ASReml's own verb for the same job. Computed on the draws, so the interval is a quantile of the functional rather than a delta-method approximation and cannot leave the parameter space. Refuses by name when the functional is missing or malformed (`vpredict_no_functional`), names a component the fit does not carry (`vpredict_component_not_in_posterior`), or does not return one value per draw (`vpredict_functional_not_one_value_per_draw`). |

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
  every fit subclass. `glance()` returns its one-row summary on an INLA fit
  with the sampler-specific columns (`chains`, `samples`, `max_rhat`,
  `min_ess`) `NA`, since a deterministic Laplace fit has no sampler to
  report. `augment()` has no INLA-side answer and refuses by name, pointing
  at `tidy()` and `predict()` instead.
- **emmeans**: `recover_data()`, `emm_basis()`.
- **marginaleffects**: `get_coef()`, `get_predict()`, `get_vcov()`, `set_coef()`.
- **insight**: `get_data()`.
- **base/stats**: `coef()`, `confint()`, `vcov()`, `predict()`, `fitted()`,
  `residuals()`, `family()`, `formula()`, `logLik()`, `nobs()`, `summary()`,
  `anova()`, `update()`, `plot()`.

`predict.flexybayes_inla()` and `predict.flexybayes_brms()` accept a `newdata`
interface mirroring `stats::predict()`. Both take `type`, `se.fit`, `classify`
and `level`; the brms method additionally takes `re_formula` and `summary`.
There is no parent `predict.flexybayes()` method: it was removed at 0.9.3 so
that each engine states its own prediction contract. Chunked output is
not an argument of these methods: `output_file` belongs to `fb_plan()`. The
brms method passes further arguments on to brms's own prediction functions,
so `allow_new_levels = TRUE` reaches brms's new-level policy there; the INLA
method takes no new-level argument. An earlier revision of this page
advertised three arguments that `...` silently swallowed.
Interval semantics are *posterior expected-response* (no residual observation
noise), not posterior-predictive.

Where a method has no meaning on an engine it states that rather than
returning a number. `logLik()` on an INLA fit reports that INLA returns a
marginal log-likelihood and yields `NA` with that message, and a fit carrying
neither a response vector nor a recorded family refuses instead of letting
`AIC()` or `anova()` report a comparison that was never made.

## Backends

Two engines are active. A third native engine was withdrawn entirely in
0.9.3 (see `NEWS.md`): no code path, export, registry row, or `Suggests`
entry remains, and naming it (or any other unrecognised backend) raises
an ordinary `unknown_backend` refusal. Re-entry, should it ever be
proposed, would be a fresh implementation, not a repair of retained
code.

| Backend | Status | Availability | Inference |
|---|---|---|---|
| `INLA` | active | from its own repository (`Additional_repositories`) | Integrated nested Laplace approximation (**approximate**) |
| `brms` | active | CRAN; needs a Stan toolchain | Hamiltonian Monte Carlo via Stan |

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
| Residual `sigma` (Gaussian, identity link) | `half_normal(scale = 5 * sd(y))` | flexyBayes heuristic (weakly informative SD prior with no upper bound; 0.13.0) |
| Residual `sigma` (log link -- Poisson, negative binomial, gamma) | `half_normal(scale = 3)` on the log scale | flexyBayes heuristic |
| Residual `sigma` (logit link -- binomial, beta) | `half_normal(scale = 5)` on the logit scale | flexyBayes heuristic |
| `sd(group)` for `simple` / `ide` / `id` random terms | `half_normal(scale = s)`, `s` by family as above | flexyBayes heuristic |
| `sd(group)` for an uncorrelated slope `(x \|\| g)` | `half_normal(scale = s)`, `s` by family | flexyBayes heuristic |
| `sd(slope)` for an uncorrelated slope `(x \|\| g)` on a numeric `x` | `half_normal(scale = s / sd(x))`, on both engines | flexyBayes heuristic (response units per covariate unit) |
| `sd(group)` for `vm()`, `ped()` structured-cov terms | `half_normal(scale = s)`, `s` by family | flexyBayes heuristic |
| `sd(group)` for a `gen:env` interaction intercept, and the per-level SDs of `diag()` / `idh()` / `at()` and `us()` | `half_normal(scale = s)`, `s` by family | flexyBayes heuristic (added at 0.9.0) |
| `us(f):g` level correlations | brms's own LKJ | recorded as engine-default |
| AR1 field hyperparameters (`ar1()`, `ar1(row):ar1(col)`) | INLA's own hyperpriors | recorded as engine-default |
| Any other random-term type (a multi-way interaction, `spl()`, `fa()`) | the engine's own default | recorded as engine-default |
| Fixed-effect coefficients | each engine's own default -- brms: flat on the population-level coefficients and a response-centred `student_t` on the intercept; INLA: `prec = 0.001` on the slopes and a flat intercept | `prior_fixed_sd` is applied when it is supplied, on every backend; unsupplied, the engine's default stands and `prior_summary()` names it |

"Recorded as engine-default" is a contract, not a shrug: the parameter and the
reason are written into the fit's prior provenance, so `prior_summary()` lists
it by name with the reason.

The half-normal-on-SD default is a flexyBayes heuristic in the spirit of
Gelman (2006), which argues for weakly informative variance priors over the
conjugate inverse-gamma. It is exact on both engines, and it has no upper bound
for an approximate posterior interval to run past. Its scale is the bound of the
`uniform(0, U)` default it replaced in 0.13.0. Gelman (2006) discusses a
half-Cauchy for the few-groups regime, and the PC prior remains an explicit
choice (`fb_prior(sigma ~ pc(upper = U, prob = p))`) when the number of groups is
small. The legacy scalar bridge (`prior_fixed_sd`, `prior_vc_sd`) preserves the
original `lognormal(0, prior_vc_sd)` semantics verbatim when `prior_vc_sd` is
passed explicitly, and does so on both active engines: brms receives a
`lognormal(0, prior_vc_sd)` prior row and INLA the expression prior that writes
the same density on the standard-deviation scale, so the two carry one prior and
the fit's prior record lists the components as carried rather than as engine
defaults.

## Deprecation policy

When an experimental export is renamed, restructured, or removed:

**0.10.0 did not follow this cycle.** Eleven exports were withdrawn in one step
with no deprecation-warning release ahead of them, on the ground that the package
has never been published and so has no installed base to warn. The policy below
governs from the first public release onward, and a withdrawal without a cycle
after that is a break of it, not a precedent set by this one.

1. The old call path warns for at least one minor release, with a base R
   warning of class `flexybayes_deprecated` that names what to call.
2. The next minor release makes it an error of class `defunctError` -- the
   function still exists but signals a hard error directing callers to the
   replacement.
3. The minor release after that removes the old export.

For a default-value change, the old default remains available behind an explicit
argument for the rest of the 0.x series.

## Pinning

Production users who need stability should pin to a specific 0.12.x patch via
`renv::snapshot()` until v1.0 lands.

## Added on the 0.8.x line (all experimental)

The 0.8.x line added the entries below. All were **experimental** under the
ladder above; none changed the stability posture of the lean-core
fitting surface. **Rows marked withdrawn are no longer
exported as of 0.10.0** (see `NEWS.md`); the functions remain in the package as
internals, so a caller reaching one now gets "could not find function" rather
than the behaviour described.

| Export | Added | Notes |
|---|---|---|
| `triangulate_genomic()` / `triangulate_gwas()` | 0.8.0 | Genomic / GWAS cross-engine and field-standard comparison. `triangulate_gwas()` **withdrawn at 0.10.0**; both **removed in 0.14.0**. |
| `fb_met_summary()` | 0.8.0 | **Withdrawn at 0.10.0.** Breeder summary of a factor-analytic GxE fit. It is computed from realised factor-analytic effects, which only the engine withdrawn in 0.9.3 produced, so it now abstains unconditionally with `met_summary_not_available` and names what an active engine reports instead (`summary()` for the components, `brms::VarCorr()` for a `diag()` or `us()` covariance). |
| `fb_gblup_cv()` | 0.8.0 | Genomic-prediction accuracy by cross-validation. |
| `fb_gwas()` | 0.8.0 | EMMAX / P3D whole-genome scan. **Withdrawn at 0.10.0.** |
| `tidy()` / `glance()` / `augment()` | 0.8.1 | broom-style accessors (re-exported from `generics`). `tidy()` has a method for every fit subclass. `glance()` returns a row for an INLA fit with the sampler-specific columns `NA`; `augment()` refuses by name on one. |
| `fb_gev()` / `fb_dirichlet()` | 0.8.1 | Generalised-extreme-value and Dirichlet fitters (with `fb_family_*` descriptors). **All six withdrawn at 0.10.0.** |
| `fb_log_posterior()` | 0.8.2 | **Withdrawn at 0.10.0.** Constellation C4 producer. The engine withdrawn in 0.9.3 was the only one that evaluated the log density, so the method now abstains unconditionally with a typed message (`fb_c4_unavailable`) rather than returning a number. |
| `glance.flexybayes_inla()` / `augment.flexybayes_inla()` | 0.8.3 | Added so an INLA fit stops raising a bare "no applicable method" error. `glance()` has since returned a one-row summary with the sampler columns `NA` rather than refusing (0.9.3, D15); `augment()` still refuses by name and points at `tidy()` and `predict()`. |
