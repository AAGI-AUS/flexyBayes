# flexyBayes

Flexible Bayesian Mixed Models with ASReml and brms-Style Syntax

Licence: MIT. Version 0.10.0 is an experimental release: every export
is at the experimental `lifecycle` stage and the API may change within
the 0.x series.

`flexyBayes` is a multi-backend Bayesian mixed-model framework. It
routes one model specification to INLA (integrated nested Laplace
approximation) or brms (Stan passthrough), refuses by name what neither
engine can represent, and reports how far the two posteriors sit apart
when the same model is fitted on both -- a diagnostic, `triangulate()`,
not a proof of correctness. The parser, the intermediate representation
and the prior interlingua are shared, so agreement between the engines
is evidence about the samplers, not about the model.

> **Development release.** All exports are at the experimental
> `lifecycle` stage and the API may change within the 0.x series. Not on CRAN.
> **flexyBayes fits on two active engines: brms and INLA.** A third native
> engine was withdrawn entirely in 0.9.3 (see `NEWS.md`); naming it, or any
> other unrecognised backend, now raises an ordinary unknown-backend refusal.
> See `system.file("KNOWN_ISSUES.md", package = "flexyBayes")` for the
> current per-backend capability boundaries before relying on results.

## Which entry point do I use?

| If you have ... | Use | Notes |
|---|---|---|
| ASReml syntax (`fixed` / `random` / `residual`) | `fb()` / `flexybayes()` | Variance-component / agricultural workflows. |
| brms or lme4 syntax (`y ~ x + (1 \| g)`) | `fb()` / `flexybayes()` | The grammar is detected from the call, and `syntax =` forces it. |
| To force one engine | `fb_inla()` / `fb_brms()` | Single-engine pins. A conflicting `backend` is refused. |
| Two fits you want to compare across engines | `triangulate()` | Auto-resolves canonical parameter names. |

`fb()` is the short alias for `flexybayes()`, and either name is the
universal entry that spans every backend.

## Which backend will I get?

| Verb | `inla` | `brms` (Stan) | `auto` |
|---|:-:|:-:|:-:|
| `fb()` / `flexybayes()` | ✓ | ✓ | ✓ (INLA or brms via `lgm_gate()`) |
| `fb_inla()` | ✓ | – | – |
| `fb_brms()` | – | ✓ | – |

The universal entry reaches any active backend: name one with
`backend =`, or let `backend = "auto"` choose. Each `fb_<engine>()` pin
fits exactly one engine and refuses a conflicting `backend`.
`backend = "auto"` runs the LGM feasibility gate and routes to INLA on
acceptance, otherwise to brms. A model neither can represent refuses with
`auto_no_active_route` rather than fitting something else, and nothing is
translated silently. `auto` does change engine in one case after starting:
an INLA program failure is caught and re-routed to brms, reported by a
message rather than silently, and the fit's class and `backend` field
always say which engine produced it. Reach Stan explicitly with `fb_brms()` or
`fb(..., backend = "brms")`. Naming a withdrawn or otherwise
unrecognised backend raises `unknown_backend` -- see *Backend support*
below.

## Backend support

flexyBayes is standalone-functional with any one backend installed, and the
planner (see *Quick start*) needs none at all. The backends differ in install
burden and in what they offer.

| Backend | On CRAN? | Install burden | Inference | In flexyBayes |
|---|---|---|---|---|
| INLA | No (own repository) | Moderate -- binary, no compiler | Approximate (nested Laplace) | Supported |
| brms (Stan) | Yes | Heavy -- first-call Stan compile (~30--60 s) | MCMC (sampling error only) | Supported |

A third native engine was withdrawn entirely in 0.9.3 -- see `NEWS.md`.

All exports are at the **experimental** `lifecycle` stage. See
`API_STABILITY.md` in the source repository for what that guarantees. The
planner runs with no inference backend installed at all, so it is the least
you need in place to explore the package. For a worked fit with production
sampling settings and its diagnostics reported in full, follow the
*getting started* vignette.

### Backend support by model class

What each active engine does, by model class. The table below is generated
from a single R-level source and every verdict in it is re-derived from the
gate and emit code by `tests/testthat/test-capability-matrix.R`. Editing it
by hand fails that test. `fits` means the structure emits and a test
exercises it -- it does not promise that every fit converges at small
budgets, which is model-specific and always reported, so treat a high R-hat
badge as a diagnostic rather than a result. Only the two active engines
are columns; see the callout above.

<!-- capability-matrix:begin -->
| Model class | Spelling | INLA | brms | Notes |
|---|---|:-:|:-:|---|
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | fits | fits | The certified overlap class, which both engines emit and `triangulate()` compares. |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | fits | fits | INLA's likelihood allowlist is read from `INLA::inla.models()` when INLA is installed. |
| Hurdle gamma (zero mass plus a positive gamma part) | `family = "hurdle_gamma"` | refuses | fits | brms-native (`dpars` mu, shape, hu); the zero-mass probability `hu` keeps brms's own prior. INLA's likelihood roster carries no counterpart, so the family gate refuses it there and `auto` routes to brms. |
| Uncorrelated random slope | `(x \|\| g)` | refuses | fits | The three-arbitrator verification named a since-withdrawn engine as one arbitrator, so the INLA mapping stays deferred until the criterion is rebuilt around the active engines. The deferral is host-independent -- no local artefact lifts it. `auto` routes to brms. |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | refuses | fits | The indexed-slope INLA mapping shares the deferred three-arbitrator verification with the uncorrelated random slope, and refuses on every host. `auto` routes to brms. |
| Correlated random slope | `(x \| g)` | refuses | refuses | Refused at ingest, before any engine is chosen. Fit `(x \|\| g)` when the correlation is not of inferential interest. |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | refuses | fits | INLA collapses the finest strata, so it refuses rather than reporting a zero. brms emits `(1 \| a:b)`. |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | refuses | fits | One variance per level of `f`, no covariance between levels. All three spellings emit identical code. |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | refuses | fits | The correlated sibling of the diagonal structure -- `k(k+1)/2` parameters against `diag()`'s `k`. At one observation per cell the residual variance is confounded with the diagonal of the covariance: the covariance block converges and `sigma` does not, and a longer chain does not help. Replicate within cell, or put an informative prior on the residual. |
| Heterogeneous variances with one shared correlation | `~ corh(f):g` | refuses | refuses | No active engine has an equicorrelation group-level structure. Use `diag(f):g` or `us(f):g`. |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | refuses | fits | Lowered to distributional regression on log sigma, `sigma ~ 0 + f`. Refused for families with no residual scale. |
| Combined interaction random effects and heterogeneous residual (full MET) | `random = ~ gen + gen:env` with the `dsum` residual | refuses | fits | The emit carries both the group-level term and the `sigma` predictor, and a live fit samples cleanly on simulated multi-environment data. `auto` reaches brms for this class. |
| Factor-analytic genotype-by-environment covariance | `~ fa(env, k):gen` | refuses | refuses | Parsed for the formula catalogue and refused at dispatch -- no active engine emits a factor-analytic covariance. |
| Multi-trait covariance | `~ us(trait):vm(gen)` | refuses | refuses | No active engine represents a trait-by-genotype unstructured covariance. |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits | INLA takes the sparse-precision, pedigree-precision and block carriers, and brms additionally takes dense and Cholesky. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses | A latent AR1 field plus the Gaussian observation nugget -- four hyperparameters, one observation per grid node. This is not ASReml's three-parameter nugget-free residual, so the residual spelling refuses and names this one. |
| Per-trial separable AR1 field | `random = ~ at(trial):ar1(row):ar1(col)` | fits | refuses | One field realisation per level of `trial`, via INLA's `replicate =` mechanism, but the row correlation, column correlation and field SD are shared across every level -- not estimated per trial. `at(trial, level):ar1(row):ar1(col)` (a level argument, asking for a single conditioned trial or for per-trial hyperparameters) refuses by name (`at_field_per_level_hyper_not_representable`). brms has no lowering for either spelling. |
| Univariate P-spline | `~ spl(x)` | fits | refuses | Mapped to INLA's second-order random walk. brms has no lowering for the smooth basis. |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits | Precision weighting, Var(y_i) = sigma^2 / w_i (the ASReml / lme4 / glm(weights=) sense): INLA's `scale = w`; on brms a known offset on the sigma distributional parameter, NOT brms's own `weights()` addition term (a different, likelihood-power quantity per brms's own documentation). Both engines match lme4::lmer(weights=) closely on a shared simulated fixture. Any other family, or a non-identity link on Gaussian, refuses by name (`weights_requires_gaussian`); `aggregate = TRUE` alongside weights also refuses by name (`weights_not_aggregatable`). |
| Exact sufficient-statistic aggregation | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a | Exact cell-likelihood aggregation for iid exponential-family models with small cell count. The brms path has no aggregated emit. |

`fits` -- the engine emits the structure and a test exercises it. `emits` -- the engine generates the structure and no live fit has yet confirmed it samples acceptably. `refuses` -- the request raises rather than fitting something else. `n/a` -- the class does not apply to that engine's interface.

This block is generated from `.fb_capability_matrix()` by `tools/generate_capability_matrix.R`. Edit the R table, re-run the generator, and let `tests/testthat/test-capability-matrix.R` check that every verdict still matches the gate and emit code. Do not edit the rows here by hand.
<!-- capability-matrix:end -->

> **What the aggregation row is worth.** The sufficient-statistic route is
> what carries the package to dataset sizes the per-row path cannot reach:
> the per-row path runs out of memory between one and five million rows on
> a 32 GB machine, while the streamed path fits five billion rows through a
> roughly flat memory envelope, because it is always fitting the same small
> number of cells. The measured record -- sizes, wall-clock, peak memory,
> the hardware, and the model scope outside which the route refuses rather
> than approximates -- is banked with the package at
> `system.file("validation/benchmark_scaling.md", package = "flexyBayes")`.

> **MET capability, stated currently.** A multi-environment-trial model fits
> on brms, and its pieces are checked against ASReml: nested
> genotype-by-environment random effects recover every variance component
> against the REML reference on `agridat::besag.met`. `diag()` / `idh()` /
> `at()` and `us()` genotype variances emit and are validated on the
> parameter count before the values, and the heteroscedastic residual
> `dsum(~ units | env)` returns per-site posterior-mean variances of
> 0.1146 / 1.1516 / 4.6981 where ASReml gives 0.1093 / 1.1248 / 4.6079 --
> within 4.8%, largest on the smallest variance. The *combination* now fits
> too: interaction random effects and a sectioned residual in one model
> generate both blocks in the emitted Stan program and sample with every
> R-hat below 1.05 and no divergent transitions on simulated data, which is
> what moved that row from `emits` to `fits`. Two limits are worth knowing
> before planning around it. The verified fit is 120 rows, so the scale a
> national trial series needs is untested -- `seed` and `control` reach
> `brms::brm()`, so a larger series that needs a raised `adapt_delta` has a
> route to it, but nobody has run one. And
> `us(env):gen` on unreplicated data confounds its covariance diagonal with
> the residual -- the covariance converges, the residual does not, and more
> iterations make it worse. INLA refuses every heterogeneous-variance
> structure and every interaction random effect by name. See
> `system.file("KNOWN_ISSUES.md", package = "flexyBayes")` for the
> per-engine reasons.

**Breeder MET summaries.** Overall performance, stability, GxE BLUPs, factor
loadings and environment genetic correlations were computed from a
factor-analytic (`fa(env, k):gen`) fit's identified *realised* effects on the
engine withdrawn in 0.9.3 (see `NEWS.md`). No active engine produces that fit
shape, so that summary is unavailable in this release and its entry point is
no longer exported. The INLA MET route gives variance components via
`summary()`.

## Installation

```r
# INLA (approximate-inference backend) -- not on CRAN
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"))

# brms (Stan passthrough) -- on CRAN
install.packages("brms")

# flexyBayes itself (not yet on CRAN) -- install from the repository:
# install.packages("remotes")
remotes::install_github("AAGI-AUS/flexyBayes")
```

`flexyBayes` degrades gracefully when an optional engine is missing:
each backend is detected at run time, and a model sent to an
unavailable engine is refused with a clear message naming what to
install rather than failing obscurely.

## Quick start

The planner needs no inference backend to show what flexyBayes will do with a
model: it builds the intermediate representation, chooses a backend, and
reports the plan without fitting.

```r
library(flexyBayes)
data(sleepstudy, package = "lme4")

# Inspect the routing decision and representation plan -- no backend needed
plan <- flexybayes(
  fixed  = Reaction ~ Days,
  random = ~ Subject,
  data   = sleepstudy,
  plan   = TRUE
)
plan
```

To fit, install at least one backend (see *Backend support* above). The
following uses production sampling settings. The *getting started* vignette
walks through the same fit with its convergence diagnostics.

```r
fit <- flexybayes(
  fixed  = Reaction ~ Days,
  random = ~ Subject,
  data   = sleepstudy,
  n_samples = 2000, warmup = 5000, chains = 4
)

# Standard R output
summary(fit)
coef(fit)
confint(fit)

# emmeans + marginaleffects
emmeans::emmeans(fit, ~ Days, at = list(Days = c(0, 5)))
marginaleffects::avg_slopes(fit)
```

## Cross-engine triangulation

Fit the same model on two backends and compare:

```r
fit_brms <- flexybayes(Reaction ~ Days, random = ~ Subject,
                       data = sleepstudy, backend = "brms")
fit_inla <- flexybayes(Reaction ~ Days, random = ~ Subject,
                       data = sleepstudy, backend = "inla")

triangulate(fit_brms, fit_inla)
```

The result carries an overall `status` -- `concordant`, `discordant` or
`inconclusive` -- and a per-parameter table: each fit's posterior mean and
SD, the mean shift and Wasserstein-1 distance in posterior-SD units, the
SD ratio, and a per-parameter verdict.

Three gates run before any of those numbers is reported.

- **Comparability.** Each fit carries a fingerprint of the question it
  answered: canonical formula triple, family and link, data dimensions,
  column names and content digest, and the recorded prior per variance
  component. Two fits that disagree are refused by name
  (`triangulate_incomparable_fits`). There is no override argument.
- **Diagnostics.** A fit that failed its own convergence checks makes the
  status `inconclusive` and yields no parameter verdicts. Disagreement
  between an unconverged fit and a converged one is not a finding.
- **Matched priors.** A variance component is compared only when both fits
  record the same prior for it. The rest are reported as `not_compared`
  with the reason, most often that the term sits outside the shared
  default and each engine chose its own hyperprior.

The thresholds behind the verdicts (0.1 posterior SD on the mean shift and
on Wasserstein-1, an SD ratio inside [0.9, 1.1]) are heuristics for
directing attention, in the sense Gelman et al. (2020) use for
cross-implementation checks. They carry no error rate, and `concordant` is
not a calibration claim: the package ships no simulation-based calibration
of its own fits.

brms (full HMC) and INLA (a Laplace approximation) sit on different
inference paradigms, so the comparison is independent evidence about
the *inference*. It is not independent evidence about the model: both
fits come from the same parsed representation, so a mistranslation is
common-mode and triangulates perfectly. The *cross-engine triangulation*
vignette works through the disagreement patterns, and
the package keeps an internal registry of which code paths the two engines
have been certified not to share.

`canonical_names()` does the work of aligning backend-native parameter
names (brms's `sd_g__Intercept`, INLA's `Precision for g` on the
precision scale) to a single canonical name with the correct scale
transform -- no `name_map` argument needed in standard cases.

## Companion accessors

| Accessor | Returns |
|---|---|
| `backend_decision(fit)` | The captured dispatch trace: backend, path, `lgm_gate` checks, reason. |
| `prior_summary(fit)` | The resolved prior -- auto-default (weakly-informative bounded uniform on SD, with half-Cauchy advised for small group counts), user-supplied `fb_prior()`, or legacy scalar bridge. |
| `canonical_names(fit)` | The backend-native ↔ canonical-name table with per-row scale transforms. |
| `review_code = TRUE` on `flexybayes()` / `fb_brms()` | Inspect-before-fit workflow. `cat_code(rev)` prints the generated backend code, and `proceed(rev)` advances into the fit. Supported on the formula-entry verbs only. |

A fitted object exposes everything a downstream tool needs through these
exported accessors, so a pipeline can consume a fit -- its variance
components, its dispatch decision, its seed -- without reading flexyBayes
internals. The AAGI ORCHESTRA workspace, which coordinates several
analysis packages, builds its provenance records from exactly this
surface; nothing in flexyBayes depends on it.

## ASReml-hands accessors

If you arrive with an `asreml()` call in hand, *Getting started*'s
section 6, "Reading a fit the way a REML user does", walks the
accessors below one at a time. They are views over an ordinary Bayesian
fit -- the estimator is named on every table, because a posterior mean
is not a REML component.

| You type | You get |
|---|---|
| `summary(fit)` | One eleven-slot object on every engine: `$fixed`, `$varcomp`, `$random`, `$missing`, `$converge`, `$n_design`, `$n_observed`, `$na_action`, `$model`, `$engine`, `$call`. A fit carrying an autoregressive latent field adds `$spatial_field`. |
| `summary(fit)$varcomp` | Variance components on the standard-deviation scale with a credible interval, the prior each component ran under, and a note column for boundary collapse. |
| `coef(fit, what = "random")` / `ranef(fit)` | Random-effect predictions, one data frame per grouping factor. `what = "missing"` returns the unobserved design cells, `what = "all"` returns all three tables. |
| `predict(fit, classify = "Variety", level = 0.95)` | The marginal-means table, built on `emmeans`. Means and credible intervals only -- no pairwise standard-error block. |
| `nobs(fit, type = "observed")` | Observed responses, beside `nobs(fit)`, which stays the design row count the engine saw. |
| `na_action = list(y = "include", x = "fail")` | ASReml's `na.method()` vocabulary, detected by shape. An `asreml::na.method()` object is accepted directly, with no asreml dependency. |
| `fb_complete_grid(data, ~ row * col, response = "yield")` | Absent design cells reinstated with an `NA` response. A varying design factor refuses unless `unused_level =` names the level to write. |
| `plot(fit, type = "variogram")` | The empirical residual semivariance over the design array, computed on the observed rows. |
| `update(fit, random = ~ Block + Variety)` | A re-fit on either engine, carrying every recorded argument forward, `na_action` and the resolved prior included. |

There is no `wald()` method, no pairwise standard-error table, and no
covariate zero-fill. Each is a deliberate absence rather than a gap, and
the vignette says why.

### Two doors, one object

The same fit answers to the generics a Bayesian reaches for. Neither
idiom is a wrapper around the other.

| ASReml hands | Bayesian hands |
|---|---|
| `summary(fit)$varcomp` | `posterior::as_draws_df(fit)` |
| `coef(fit, what = "random")`, `ranef(fit)` | `posterior::as_draws()`, `posterior::as_draws_matrix()`, `fb_as_draws_simple(fit)` |
| `predict(fit, classify = "Variety")` | `emmeans` or `marginaleffects` on the same fit |
| `plot(fit, type = "variogram")` | `pp_check(fit)`, `plot(fit, type = "diagnostics")` |
| `summary(fit)$converge` | The same slot, reporting the engine's own diagnostics |
| No ASReml counterpart | `prior_summary(fit)`, `loo(fit)`, `triangulate(fit_a, fit_b)` |

`loo()` and `pp_check()` pass through to brms on a sampled fit and refuse
by name on a Laplace fit, naming the WAIC and DIC that fit does carry.

## Output structure

Every fit carries three top-level slots:

```r
fit$glm         # GLM-compatible shim -- works with summary(), emmeans,
                # marginaleffects, broom

fit$inla        # native INLA output (when backend = "inla")
fit$brms        # live brmsfit (when backend = "brms")

fit$extras      # BLUPs, variance components, convergence diagnostics,
                # generated code, parsed IR, run time, captured call
```

There is no third-engine slot on a fit object: naming a withdrawn or
otherwise unrecognised `backend` refuses before a fit object exists (see
*Backend support*).

## Supported ASReml syntax (reference)

```r
# Fixed effects
yield ~ env                  # fixed factor
yield ~ env + x_cov          # factor + covariate
yield ~ 0 + env              # means model (no intercept)
yield ~ env + I(x^2)         # expression terms

# Random effects
random = ~ geno                       # simple iid
random = ~ block:rep                  # nested
random = ~ vm(geno, Gmat)             # GBLUP (dense V; brms)
random = ~ vm(geno, chol = L)         # user-supplied Cholesky (brms; v0.3.7+)
random = ~ vm(geno, precision = Q)    # user-supplied sparse precision (INLA; v0.3.7+)
random = ~ ped(animal, Amat)          # pedigree (animal model)
random = ~ ped(animal, A_inv,
               use_sparse_precision = TRUE)  # sparse pedigree precision (v0.3.7+)
random = ~ diag(env):geno             # diagonal GxE (brms; idh() and at() are synonyms)
random = ~ us(env):geno               # unstructured GxE (brms)
random = ~ corh(env):geno             # equicorrelated GxE -- refused by name
random = ~ fa(env, 2):id(geno)        # factor-analytic GxE -- refused (no active engine)
random = ~ ar1(row):ar1(col)          # separable AR1 field + nugget (INLA)
random = ~ ar1(t)                     # one-dimensional AR1 field + nugget (INLA)
random = ~ spl(x_cov)                 # P-spline (INLA)

# Residual
residual = ~ units                # iid residuals (default)
residual = ~ dsum(~ units | env)  # one residual variance per environment (brms)
residual = ~ at(env):units        # the same structure, ASReml's other spelling
residual = ~ ar1(row):ar1(col)    # ASReml's nugget-free residual -- refused by
                                  #   name; write the field on the random side

# Families
family = "gaussian" | "binomial" | "poisson" | "negative_binomial" |
         "gamma" | "beta"
```

## Vignettes

Eleven vignettes ship with the package:

| # | Vignette |
|---|---|
| 01 | Getting started |
| 02 | The formula surface: ASReml-shaped terms and structured covariance |
| 03 | Foundational regression |
| 04 | Hierarchical models |
| 05 | Default priors and the `fb_prior()` DSL |
| 06 | METs and genomic selection |
| 07 | Downstream analysis |
| 08 | Spatio-temporal models |
| 09 | Cross-engine triangulation |
| 10 | Dispatch, refusals, and the backend registry |
| 11 | Streaming exact aggregation |

Vignette 10 is the technical/internals reference, and the rest target a
general audience. The deck is contiguous 01--11 as of the 0.9.3 renumber
(see `NEWS.md`): the pre-0.9.3 deck ran 00 and 06--16 with gaps at 05 and
12--15 recording earlier merges, and *From an ASReml call* (page 00) was
folded into *Getting started*'s new accessor section rather than kept as
its own page. Every number in the table above is the page's current,
shipped filename suffix.

Heavy MCMC vignettes use a `.Rmd.orig` precompile pattern. The `.Rmd`
that ships in the package tarball is the pre-evaluated static output.
Browse them with `browseVignettes("flexyBayes")` **after a full install** --
`R CMD build` then `R CMD INSTALL` the tarball, or
`devtools::install(build_vignettes = TRUE)`. A plain `install_github()` or
source-directory install does **not** build the vignettes into `inst/doc`.

Some vignettes fit at small sampling budgets for speed and say so where it
matters, and each flags any diagnostic that falls short of production
convergence thresholds rather than presenting it as a clean result. For a
convergence-clean workflow see the *getting started* and *cross-engine
triangulation* vignettes.

## Design principles

- **Indexing, not model matrices**: parameters are sized by number of
  levels (`p`), not observations (`N`), giving O(p) memory for random
  effects.
- **Non-centred parameterisation**: all random effects use NCP for
  efficient MCMC sampling.
- **Cholesky decomposition**: known covariance matrices (G, A) are
  decomposed once.
- **Triangulation as first-class output**: the package is built
  around the claim that *evidence of cross-engine agreement* is the
  signature flexyBayes deliverable, not "one more Bayesian mixed-
  model frontend".

## Correctness

flexyBayes ships an extensive `testthat` suite (`devtools::test()`) covering
the paths most likely to hide errors: fixed-effect and factor models against
`lm()` / `lme4::lmer()` references, random intercepts and the asreml-route
random slopes, structured-covariance terms, streaming-aggregation equivalence
to the per-row fit, weights, offsets, missing-response handling, backend
routing and the structured refusal taxonomy, prior translation, and
cross-engine `triangulate()` agreement. The *cross-engine triangulation* and
per-family vignettes show clean reproducible checks.

## Testing & CI

Continuous integration validates the INLA, brms, and engine-independent
surface (the ASReml / brms parsers, the intermediate representation,
`lgm_gate()`, the dispatch policy table, the refusal registry, the prior DSL,
and the `triangulate()` metrics). A third native engine was withdrawn
entirely in 0.9.3 (see `NEWS.md`); its tests were deleted along with the
engine rather than skipped, so the suite carries no dormant coverage for a
capability the package no longer offers. Run the full suite locally with
`devtools::test()`.

## Known limitations

flexyBayes refuses what it cannot yet fit rather than fitting it
silently. The current release does not cover the following. Each is a
roadmap deferral, and a request that needs one is met with a
structured refusal naming the gap, not a quiet wrong answer.

- **Scale ceiling on the per-row path (realistic multi-term design)**: on a
  crossed/nested multi-environment-trial design the flexyBayes/INLA ceiling
  is bracketed between 911,808 rows (preflight refuses) and 1,823,616 rows
  (the engine dies after 41.6 minutes), not measured or bisected -- see
  `inst/validation/benchmark_scaling.md` for the two logged rungs.
- **Multi-environment-trial scale**: the combined model -- interaction
  random effects (`gen:loc`, `gen:loc:yearf`) *together with* a
  heteroscedastic per-environment residual (`dsum(~ units | env)`) --
  fits on brms and is verified by a live fit, but on 120 simulated rows.
  A national trial series is untested. The sampler controls such a run
  would need are available (`seed` and `control` are forwarded to
  `brms::brm()`), so what is missing is the run, not the route. INLA
  refuses both halves. See `inst/KNOWN_ISSUES.md` for the status.
- **Structured GxE beyond `diag` / `us`**: `corh(env):gen` (heterogeneous
  variances with one shared correlation) and `fa(env, k):gen`
  (factor-analytic) have no active emit and refuse by name.
- **Spatial structure**: only the separable AR1 field
  (`random = ~ ar1(row):ar1(col)`) is supported, on INLA, over a complete
  grid with one observation per node. It is emitted as an AR1-by-AR1 latent
  field plus the Gaussian observation nugget -- four hyperparameters, all
  four reported by `summary()` -- and is therefore not ASReml's
  three-parameter nugget-free residual. Writing it on the residual side
  refuses by name and points here. Intrinsic
  CAR and BYM2 areal models are not implemented. You can express a custom
  spatial precision by passing your own matrix to `vm(g, precision = Q)`,
  but there is no BYM2 helper.
- **Smooth terms**: univariate penalised splines (`s(x)`, `spl(x)`) are
  supported on INLA. Multivariate and tensor-product smooths (`te()`,
  `ti()`, `t2()`) are refused and deferred to a later release.
- **Observation weights**: fit for the Gaussian family on an identity
  link, on both engines, in the ASReml / lme4 / `glm(weights=)` precision
  sense (`Var(y_i) = sigma^2 / w_i`) -- both engines match
  `lme4::lmer(weights=)` closely on a shared simulated fixture. Any other
  family, or a non-identity link on Gaussian, refuses by name
  (`weights_requires_gaussian`) rather than returning the unweighted
  posterior under a weighted call; `aggregate = TRUE` alongside weights
  also refuses by name (`weights_not_aggregatable`).
- **Hidden-Markov, multi-state, and survival models**: not supported.
  Survival / time-to-event families are refused at the family gate.
  A NIMBLE backend covering these is on the roadmap with no fixed
  release target.
- **Missing data**: flexyBayes does not impute covariates. A missing
  *predictor* is refused by default, as it is in ASReml
  (`na.method(x = "fail")`). Setting `na_action = list(y = "include", x = "omit")`
  drops the affected rows with a warning naming the count and the columns,
  and ASReml's third setting -- treat the missing covariate as zero -- is
  refused by name rather than reproduced. A missing *response* is
  retained by default (`na_action = "augment"`), carried as a latent
  quantity the engine marginalises, so the design index set a structured
  covariance is built over survives a lost plot. That preserves the
  representation, not information: under ignorability the parameter
  posterior is the same either way, and where missingness depends on the
  unobserved response both augmenting and omitting are biased. That
  identity is a statement about the posterior and not a promise about what
  an optimiser returns -- the two settings hand the engine an intact design
  and a ragged one, and on a weakly identified mode they can land in
  different places. See `?flexybayes`. Non-Gaussian missing responses on
  brms are refused.

## Requirements

- R ≥ 4.1.0
- INLA (optional, via Additional_repositories)
- brms (optional, for the Stan passthrough)

## Contributing

Contributions are welcome. See `CONTRIBUTING.md` for the development workflow
(fork, `usethis::pr_init()`, `devtools::check()` must pass, and a `NEWS.md`
bullet for any user-facing change) and `CODE_OF_CONDUCT.md` (Contributor
Covenant). Please report bugs and request features at
<https://github.com/AAGI-AUS/flexyBayes/issues>. The package's architecture
decisions are indexed in `DESIGN_DECISIONS.md`.

## Citation

If you use `flexyBayes` in your research, please cite:

```
@software{flexyBayes,
  title  = {flexyBayes: Flexible Bayesian Mixed Models with ASReml
            and brms-Style Syntax},
  author = {Moldovan, Max and Tanaka, Emi and Hui, Francis K.C. and
            Forte Deltell, Anabel},
  year   = {2026},
  url    = {https://github.com/AAGI-AUS/flexyBayes}
}
```

## License

MIT. See `LICENSE` and `LICENSE.md`.
