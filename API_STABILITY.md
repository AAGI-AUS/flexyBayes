# API stability

This document describes the stability of the flexyBayes 0.18.0 interface.
Every export is at the experimental lifecycle stage, so the guarantees are
weaker than they will be at 1.0. Bug fixes and additions never break
callers. Renamings, default changes and changes of shape are allowed within
0.x after a deprecation cycle of at least one minor release, in which the old
call warns with a condition of class `flexybayes_deprecated`.

The document covers the mixed-model fitting interface.

## Lifecycle stages

| Stage | Meaning for flexyBayes |
|---|---|
| **stable** | Behaviour and signature fixed across the 1.x series. A breaking change needs a major version and a deprecation cycle. *No exports are stable yet.* |
| **maturing** | Behaviour and signature unlikely to change in 0.x. Added arguments do not break callers. |
| **experimental** | Behaviour and signature may change in any 0.x release after a deprecation cycle of at least one minor release. **All current exports.** |
| **deprecated** | Still callable. Warns with a base R condition of class `flexybayes_deprecated` that names what to call instead. The removal is announced in `NEWS.md`. |
| **superseded** | Still callable without a warning, but a newer entry point is preferred. |

## Fitting entry points

flexyBayes has two kinds of fitting entry point.

**Universal entry points** take a `backend` argument and send the model to an
engine:

| Export | Stage | Notes |
|---|---|---|
| `flexybayes()` | experimental | The ASReml-style entry (`fixed`, `random`, `residual`). It also accepts brms-style formulas through `syntax = "auto"`. `backend = c("auto", "inla", "brms")`, default `"auto"`. A withdrawn or unrecognised backend is refused by name with `unknown_backend`, not by `match.arg()`. `"auto"` fits on INLA when `lgm_gate()` accepts the model and INLA is installed, otherwise on brms when brms can represent the model, and otherwise refuses with `auto_no_active_route`. Nothing is translated silently. In one case `"auto"` changes engine after starting: an INLA program failure is caught and the model is refitted on brms. A message reports this, and the fit's class and `backend` field record which engine produced it. `prior` accepts an `fb_prior()` object. |
| `fb()` | deprecated in 0.13.0 | Alias for `flexybayes()` with the same arguments. It warns and forwards the call. Call `flexybayes()`. |

**Engine pins** fix one engine and take **no** `backend` argument. Passing a
conflicting `backend` raises an `engine_pin_backend_conflict` refusal.

| Export | Stage | Notes |
|---|---|---|
| `fb_inla()` | deprecated in 0.13.0 | Pins the INLA engine (approximate inference, see "Inference semantics"). Call `flexybayes(..., backend = "inla")`. |
| `fb_brms()` | deprecated in 0.13.0 | Pins the Stan/brms engine through `brms::brm()`. Call `flexybayes(..., backend = "brms")`. |

**Sampler arguments.** `seed` and `control` are part of the universal
signature and reach `brms::brm()` unchanged, so a brms posterior can be
reproduced from the call alone. `control = list(adapt_delta = 0.95)` sets
`adapt_delta`, and `max_treedepth` is set the same way. Left `NULL`, they
fall back to brms's own defaults (`seed = NA`, `control = NULL`) and are not
passed through. `update()` repeats whatever the fit recorded. Both arguments
have no effect on the deterministic INLA path, which says so once per
session.

**Grammar adapters** read a model *specification* written in another
package's syntax and return the flexyBayes internal representation
(`fb_terms`), which the entry points then fit. They do not take a fitted
object. Pass the formula and data, not a `brmsfit` or an `asreml` fit.

| Export | Stage | Notes |
|---|---|---|
| `fb_from_brms()` | deprecated in 0.13.0 | Compatibility parser for brms-style syntax. Use `fb_plan(..., syntax = "brms")` to inspect a model without fitting. |
| `fb_from_asreml()` | deprecated in 0.13.0 | Compatibility parser for ASReml-style syntax. Use `fb_plan(..., syntax = "asreml")`. |

## Cross-engine comparison

| Export | Stage | Notes |
|---|---|---|
| `triangulate()` | removed in 0.14.0 | Internal from 0.11.0 and removed in 0.14.0, together with `triangulate_genomic()`, `triangulate_gwas()` and the refusal code `triangulate_incomparable_fits`. Nothing in the package called it, and agreement between two engines fitted from one parsed model cannot detect a translation error they share. The canonical draws it compared are still available on each fit through `posterior::as_draws_df()` and `canonical_names()`. |

## Priors

| Export | Stage | Notes |
|---|---|---|
| `fb_prior()` | experimental | Prior specification on the standard-deviation scale. Targets: `sigma`, `sd(group)`, `b(name)`, `cor(group)`, `smooth(var)`. Families: `pc`, `normal`, `student_t`, `half_normal`, `half_cauchy`, `cauchy`, `gamma`, `exponential`, `lkj`, `uniform`. The penalised-complexity (PC) family (Simpson et al. 2017) is the common form from which one prior specification is translated into each engine's own parameterisation. |
| `prior_summary()` | experimental | S3 generic returning the resolved priors of a fit, with the route each came from: the default half-normal on each standard deviation, a user-supplied `fb_prior()`, or the legacy scalar arguments. Methods for every fit subclass. |

## Covariance, engine and approximation helpers

| Export | Stage | Notes |
|---|---|---|
| `fb_cov()` | experimental | Constructs a structured-covariance specification carried on a random term. |
| `is_fb_cov()` | deprecated in 0.13.0 | Tests for that class. Use `inherits(x, "fb_cov")`. |
| `fb_engine()` | deprecated in 0.13.0 | Compatibility constructor. Pass `backend`, `chains`, `n_samples`, `warmup` and `control` directly to `flexybayes()`. |
| `is_fb_engine()` | deprecated in 0.13.0 | Tests for that class. Use `inherits(x, "fb_engine")`. |
| `is_fb_approx()` | deprecated in 0.13.0 | Tests for the `fb_approx` class. Use `inherits(x, "fb_approx")`. |
| `fb_approx()`, `validate_approximation()` | experimental | Construct and validate a request for an approximate scheme. An unregistered scheme is refused with a message naming the available exact routes. |

## Planning, streaming and large data

| Export | Stage | Notes |
|---|---|---|
| `fb_plan()` | experimental | Returns the dispatch and aggregation plan for a model without fitting it (`plan = TRUE` on a universal entry point returns the same object). It says which backend was chosen and why. |
| `fb_select()` | deprecated in 0.16.0 | Refuses with `selection_removed`: penalised selection is no longer part of flexyBayes. Removal in the release after the first public release that announces their removal. The two-stage vignette fits a penalised regression on the design whitened by `fb_marginal_cov()`. |
| `fb_select_fold_data()` | deprecated in 0.16.0 | Refuses with `selection_removed`: penalised selection is no longer part of flexyBayes. Removal in the release after the first public release that announces their removal. The two-stage vignette fits a penalised regression on the design whitened by `fb_marginal_cov()`. |
| `fb_select_path()` | deprecated in 0.16.0 | Refuses with `selection_removed`: penalised selection is no longer part of flexyBayes. Removal in the release after the first public release that announces their removal. The two-stage vignette fits a penalised regression on the design whitened by `fb_marginal_cov()`. |
| `fb_select_refit()` | deprecated in 0.16.0 | Refuses with `selection_removed`: penalised selection is no longer part of flexyBayes. Removal in the release after the first public release that announces their removal. The two-stage vignette fits a penalised regression on the design whitened by `fb_marginal_cov()`. |
| `fb_variance_draws()` | experimental | Joint posterior draws of an INLA fit's variance components, under the names and on the scale of `summary(fit)$varcomp`. |
| `fb_refit()` | experimental | Refits the model of a fit on new data, with added fixed columns, extra brms prior rows and Stan code, a formula hook, or a reduced INLA output. |
| `fb_marginal_cov()` | experimental | The marginal covariance of the response of a fit over given rows: one block per variance component, with its design, structure and Kronecker factors where it is separable, the residual, and any known variance or covariance. `as.matrix()` assembles it as a dense matrix. |
| `fb_penalty()` | deprecated in 0.16.0 | Refuses with `selection_removed`: penalised selection is no longer part of flexyBayes. Removal in the release after the first public release that announces their removal. The two-stage vignette fits a penalised regression on the design whitened by `fb_marginal_cov()`. |
| `fb_resamples()` | deprecated in 0.16.0 | Refuses with `selection_removed`: penalised selection is no longer part of flexyBayes. Removal in the release after the first public release that announces their removal. The two-stage vignette fits a penalised regression on the design whitened by `fb_marginal_cov()`. |
| `fb_decision()` | deprecated in 0.16.0 | Refuses with `selection_removed`: penalised selection is no longer part of flexyBayes. Removal in the release after the first public release that announces their removal. The two-stage vignette fits a penalised regression on the design whitened by `fb_marginal_cov()`. |
| `flexybayes_stream()` | experimental | Chunked sufficient-statistic aggregation for data too large to hold in memory, fitted for `family = "binomial"` or `"poisson"` (which take `trials` or `exposure`). INLA is the default and the only engine with an aggregated emit. `fit = FALSE` returns the `<fb_aggregated>` object (compression ratio, `K`, `N`) without fitting, for a Gaussian family too. **Changed in 0.12.0:** a Gaussian fit refuses with `aggregation_residual_prior_not_composed` before any chunk is read, because the aggregated Gaussian likelihood did not carry the residual prior. Fit a Gaussian model per row. An aggregated or streamed fit reports no marginal likelihood (`aggregated_evidence_unavailable`). |

## Diagnostics and introspection

| Export | Stage | Notes |
|---|---|---|
| `backend_decision()` | experimental | The recorded dispatch trace, with eight fields: `backend`, `path`, `gate_checks`, `reason`, `preflight_summary`, `representation_plan`, `rejected_routes`, `routing_policy_version`. Its shape has not changed since 0.9.x. |
| `canonical_names()` | experimental | S3 generic returning the canonical parameter names of a fit. Methods for every fit subclass. |
| `fb_refusals()` | experimental | The refusal codes as a table (code, package, description, version added). Count its rows for the number of codes. |
| `fb_register_refusals()` | experimental | For packages that extend flexyBayes: adds the package's refusal codes, from its `.onLoad()`, to its own part of the registry. A code that another package or flexyBayes already registered is refused. |
| `fb_refuse()` | experimental | Raises a registered refusal as a `flexybayes_refusal` condition carrying the code and the package that registered it. |
| `fb_isolated_lapply()` | experimental | For packages that extend flexyBayes: `lapply()` with each call in a separate R process that starts from the session's settings and random-number state. Values, messages and warnings come back in job order, and a process that dies is a refusal. |
| `fb_extension_api()` | experimental | The version of the set of functions through which a package that extends flexyBayes reads a fit. It is 1 in this release. |
| `fb_model_terms()` | experimental | For packages that extend flexyBayes: the model description of a fit. |
| `fb_response_values()` | experimental | For packages that extend flexyBayes: a response, a column or an expression, evaluated over a data frame. |
| `fb_offset_values()` | experimental | For packages that extend flexyBayes: the summed offsets of a model over the rows of a data frame. |
| `fb_model_weights()` | experimental | For packages that extend flexyBayes: the observation weights a model was fitted with. |
| `fb_known_variance()` | experimental | For packages that extend flexyBayes: the known variance or covariance of a model and the data column that carries it. |
| `fb_engine_settings()` | experimental | For packages that extend flexyBayes: sampler settings, checked against those flexyBayes recognises. |
| `fb_convergence()` | experimental | For packages that extend flexyBayes: the convergence verdict of a fit with its evidence. |
| `fb_inla_posterior_sample()` | experimental | For packages that extend flexyBayes: joint posterior samples of an INLA fit on the scale of the data. |
| `fb_check_refit()` | experimental | For packages that extend flexyBayes: refuses a fit that cannot be refitted as the model it was fitted as. |
| `fb_fixed_design()` | experimental | For packages that extend flexyBayes: the design matrix of a fixed formula, with each column's term and ASReml's column order. |
| `fb_backend_status()` | experimental | Which engines are installed and usable in the current session, and why an unusable one is unusable. |
| `fb_structured_cov()` | withdrawn at 0.10.0 | The identified covariance of a factor-analytic structured-covariance term. No active engine emits an `fa()` term, so it reported nothing on every reachable fit and is no longer exported. |
| `proceed()`, `cat_code()` | deprecated in 0.13.0 | Companions to the `review_code = TRUE` workflow (inspect the generated engine code, then run it). brms is the only engine with generated code to review, so the deferred-execution token is available under `backend = "brms"` and under `"auto"`, and refuses under `backend = "inla"` with `review_code_backend_unsupported`. The generated code is the token's `code` element, and `flexybayes(..., return_code = TRUE)` returns it. To fit, call `flexybayes()` with the same arguments and `review_code = FALSE`. |
| `fb_complete_grid()` | experimental | Completes a design grid (`index = ~ row * col`) before fitting, so a structured covariance is built over the intended index set and not only over the observed rows. |
| `genomic_summary()` | experimental | Turns the posterior draws of a relationship-matrix term into genomic quantities for a fitted model. |
| `ranef()` | experimental | S3 generic returning random-effect predictions (posterior summaries of the group-level effects) from a fit. |
| `vpredict()` | experimental | S3 generic returning a posterior function of the variance components, such as a heritability, a variance ratio or a genetic correlation, written as an expression in the canonical component names. ASReml uses the same verb for the same job. It is computed on the draws, so the interval is a quantile of the function, not a delta-method approximation, and it cannot leave the parameter space. Refuses by name when the function is missing or malformed (`vpredict_no_functional`), names a component the fit does not carry (`vpredict_component_not_in_posterior`), or does not return one value per draw (`vpredict_functional_not_one_value_per_draw`). |

## Interoperability

Every fit carries its engine's class first and `"flexybayes"` second
(`c("flexybayes_brms", "flexybayes", "list")`,
`c("flexybayes_inla", "flexybayes", "list")`). A method written for the
shared class reaches both engines, and an engine-specific method takes
precedence.

flexyBayes fits work with other packages through registered S3 methods.
Their signatures are set by the host package and are as stable as that
package's interface.

- **Draws**: `fb_as_draws_simple()` (a named list of named draw vectors) and
  the `posterior`-compatible path.
- **broom**: `tidy()`, `glance()`, `augment()`. Column names follow broom
  conventions, and new columns do not break callers. `tidy()` has a method
  for every fit subclass. `glance()` returns its one-row summary on an INLA
  fit with the sampler columns (`chains`, `samples`, `max_rhat`, `min_ess`)
  set to `NA`, because a deterministic Laplace fit has no sampler to report.
  `augment()` has no INLA method and refuses by name, pointing to `tidy()`
  and `predict()`.
- **emmeans**: `recover_data()`, `emm_basis()`.
- **marginaleffects**: `get_coef()`, `get_predict()`, `get_vcov()`, `set_coef()`.
- **insight**: `get_data()`.
- **base/stats**: `coef()`, `confint()`, `vcov()`, `predict()`, `fitted()`,
  `residuals()`, `family()`, `formula()`, `logLik()`, `nobs()`, `summary()`,
  `anova()`, `update()`, `plot()`.

`predict.flexybayes_inla()` and `predict.flexybayes_brms()` take `newdata`
in the same way as `stats::predict()`. Both take `type`, `se.fit`,
`classify` and `level`, and the brms method also takes `re_formula` and
`summary`. There is no parent `predict.flexybayes()` method. It was removed
at 0.9.3 so that each engine states its own prediction behaviour. Chunked
output is not an argument of these methods. `output_file` belongs to
`fb_plan()`. Both methods take `allow_new_levels` and `mean_by`. With
`allow_new_levels = TRUE`, each grouping level absent from the fit is drawn
from its fitted variance, and the intervals are posterior-predictive,
including observation noise. `mean_by` returns the mean of the rows in each
group. Without `allow_new_levels`, intervals are for the *posterior expected
response* (no residual observation noise).

Where a method has no meaning on an engine, it says so and does not return a
number. `logLik()` on an INLA fit reports that INLA returns a marginal
log-likelihood and returns `NA` with that message. A fit with neither a
response vector nor a recorded family refuses, so `AIC()` and `anova()`
cannot report a comparison that was never made.

## Backends

Two engines are active. A third native engine was withdrawn entirely in
0.9.3 (see `NEWS.md`). No code path, export, registry row or `Suggests`
entry for it remains, and naming it, like any other unrecognised backend,
raises an `unknown_backend` refusal. Bringing it back would mean a new
implementation.

| Backend | Status | Availability | Inference |
|---|---|---|---|
| `INLA` | active | from its own repository (`Additional_repositories`) | Integrated nested Laplace approximation (**approximate**) |
| `brms` | active | CRAN, needs a Stan toolchain | Hamiltonian Monte Carlo via Stan |

## Inference semantics

flexyBayes preserves the formula: the model sent to a backend represents the
formula you wrote. That is a different claim from exact inference. The brms
backend draws from the posterior by MCMC, which is exact up to Monte Carlo
error and subject to convergence diagnostics. The **INLA backend is
approximate inference** (integrated nested Laplace approximation), even when
the model it receives is a faithful translation of the formula. flexyBayes
preserves the formula on every backend and gives exact inference only on the
MCMC backend.

## Default priors

On the default path, when neither an `fb_prior()` nor a legacy
`prior_vc_sd` scalar is supplied:

| Target | Default | Basis |
|---|---|---|
| Residual `sigma` (Gaussian, identity link) | `half_normal(scale = 5 * sd(y))` | flexyBayes heuristic (weakly informative SD prior with no upper bound, 0.13.0) |
| Residual `sigma` (log link: Poisson, negative binomial, gamma) | `half_normal(scale = 3)` on the log scale | flexyBayes heuristic |
| Residual `sigma` (logit link: binomial, beta) | `half_normal(scale = 5)` on the logit scale | flexyBayes heuristic |
| `sd(group)` for `simple` / `ide` / `id` random terms | `half_normal(scale = s)`, `s` by family as above | flexyBayes heuristic |
| `sd(group)` for an uncorrelated slope `(x \|\| g)` | `half_normal(scale = s)`, `s` by family | flexyBayes heuristic |
| `sd(slope)` for an uncorrelated slope `(x \|\| g)` on a numeric `x` | `half_normal(scale = s / sd(x))`, on both engines | flexyBayes heuristic (response units per covariate unit) |
| `sd(group)` for `vm()`, `ped()` structured-covariance terms | `half_normal(scale = s)`, `s` by family | flexyBayes heuristic |
| `sd(group)` for a `gen:env` interaction intercept, and the per-level SDs of `diag()` / `idh()` / `at()` and `us()` | `half_normal(scale = s)`, `s` by family | flexyBayes heuristic (added at 0.9.0) |
| `us(f):g` level correlations | brms's own LKJ | recorded as engine default |
| AR1 field hyperparameters (`ar1()`, `ar1(row):ar1(col)`) | INLA's own hyperpriors | recorded as engine default |
| Any other random-term type (a multi-way interaction, `spl()`, `fa()`) | the engine's own default | recorded as engine default |
| Fixed-effect coefficients | each engine's own default. brms: flat on the population-level coefficients and a response-centred `student_t` on the intercept. INLA: `prec = 0.001` on the slopes and a flat intercept | `prior_fixed_sd` is applied on every backend when it is supplied. When it is not, the engine's default stands and `prior_summary()` names it |

A parameter recorded as engine default is written into the fit's prior
record with the reason, and `prior_summary()` lists it by name with that
reason.

The half-normal default on a standard deviation is a flexyBayes heuristic in
the spirit of Gelman (2006), who recommends weakly informative variance
priors over the conjugate inverse-gamma. It is exact on both engines, and it
has no upper bound for an approximate posterior interval to run past. Its
scale is the bound of the `uniform(0, U)` default it replaced in 0.13.0.
Gelman (2006) discusses a half-Cauchy for the case of few groups, and the PC
prior remains available (`fb_prior(sigma ~ pc(upper = U, prob = p))`) when
the number of groups is small. The legacy scalar arguments (`prior_fixed_sd`,
`prior_vc_sd`) keep the original `lognormal(0, prior_vc_sd)` prior when
`prior_vc_sd` is passed explicitly, on both active engines. brms receives a
`lognormal(0, prior_vc_sd)` prior row, and INLA an expression prior with the
same density on the standard-deviation scale. The two engines therefore fit
the same prior, and the fit's prior record lists the components as carried,
not as engine defaults.

## Deprecation policy

**0.10.0 did not follow this policy.** Eleven exports were withdrawn in one
step, without a release that warned first, because the package had not been
released publicly. The policy applies from the first public release onward.

When an experimental export is renamed, restructured or removed:

1. The old call warns for at least one minor release, with a base R warning
   of class `flexybayes_deprecated` that names what to call instead.
2. The next minor release makes it an error of class `defunctError`. The
   function still exists but stops with an error naming the replacement.
3. The minor release after that removes the old export.

For a change of default value, the old default remains available through an
explicit argument for the rest of the 0.x series.

## Pinning

Users who need stability should pin an exact flexyBayes version, for example
with `renv::snapshot()`, until 1.0.

## Added on the 0.8.x line (all experimental)

The 0.8.x line added the entries below. All were **experimental**, and none
changed the stability of the fitting interface. **Rows marked withdrawn are
not exported from 0.10.0** (see `NEWS.md`). The functions remain in the
package as internal functions, so calling one now gives "could not find
function".

| Export | Added | Notes |
|---|---|---|
| `triangulate_genomic()` / `triangulate_gwas()` | 0.8.0 | Genomic and GWAS comparison across engines and against field-standard methods. `triangulate_gwas()` **withdrawn at 0.10.0**. Both **removed in 0.14.0**. |
| `fb_met_summary()` | 0.8.0 | **Withdrawn at 0.10.0.** Breeder summary of a factor-analytic GxE fit. It is computed from realised factor-analytic effects, which only the engine withdrawn in 0.9.3 produced, so it now always abstains with `met_summary_not_available` and names what an active engine reports instead (`summary()` for the components, `brms::VarCorr()` for a `diag()` or `us()` covariance). |
| `fb_gblup_cv()` | 0.8.0 | Genomic-prediction accuracy by cross-validation. |
| `fb_gwas()` | 0.8.0 | EMMAX / P3D whole-genome scan. **Withdrawn at 0.10.0.** |
| `tidy()` / `glance()` / `augment()` | 0.8.1 | broom-style accessors (re-exported from `generics`). `tidy()` has a method for every fit subclass. `glance()` returns a row for an INLA fit with the sampler columns `NA`, and `augment()` refuses by name on one. |
| `fb_gev()` / `fb_dirichlet()` | 0.8.1 | Generalised extreme value and Dirichlet fitters (with `fb_family_*` descriptors). **All six withdrawn at 0.10.0.** |
| `fb_log_posterior()` | 0.8.2 | **Withdrawn at 0.10.0.** Returned the log posterior of a fit as a function, for use by downstream tools. The engine withdrawn in 0.9.3 was the only one that evaluated the log density, so the method now always abstains with a typed message (`fb_c4_unavailable`) and returns no number. |
| `glance.flexybayes_inla()` / `augment.flexybayes_inla()` | 0.8.3 | Added so that an INLA fit no longer raises a bare "no applicable method" error. Since 0.9.3, `glance()` returns a one-row summary with the sampler columns `NA`. `augment()` still refuses by name and points to `tidy()` and `predict()`. |
