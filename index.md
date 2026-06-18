# flexyBayes

Flexible Bayesian Mixed Models with ASReml and brms-Style Syntax

`flexyBayes` is a multi-backend Bayesian mixed-model framework with
cross-engine posterior **triangulation** exposed as a first-class
robustness diagnostic. The package routes one model specification to
greta (MCMC), INLA (integrated nested Laplace approximation), or brms
(Stan passthrough), and tells you — with evidence — when the engines
agree, when they differ, and why a backend was refused.

> **Development release (v0.8.3).** All exports are at the experimental
> `lifecycle` stage and the API may change within the 0.x series. Not on
> CRAN. The full multi-environment-trial model (genotype-by-environment
> interaction random effects together with a heteroscedastic
> per-environment residual) is **not yet supported** by any backend.
> Read `system.file("KNOWN_ISSUES.md", package = "flexyBayes")` before
> relying on results.

## Which entry point do I use?

| If you have… | Use | Notes |
|----|----|----|
| ASReml syntax (`fixed` / `random` / `rcov`) | [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) / [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) | Variance-component / agricultural workflows. |
| brms or lme4 syntax (`y ~ x + (1 \| g)`) | [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) / [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) | The grammar is detected from the call; `syntax =` forces it. |
| A pre-built greta model | [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) | Pass the `greta::model()` object (or its [`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md) representation) straight in. |
| To force one engine | [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md) / [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md) / [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) | Single-engine pins; a conflicting `backend` is refused. |
| Two fits you want to compare across engines | [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md) | Auto-resolves canonical parameter names. |

[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
is the short alias for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md);
either name is the universal entry that spans every backend.

## Which backend will I get?

| Verb | `greta` | `inla` | `brms` (Stan) | `auto` |
|----|:--:|:--:|:--:|:--:|
| [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) / [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) | ✓ | ✓ | ✓ | ✓ (greta or INLA via `lgm_gate()`) |
| [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md) | ✓ | – | – | – |
| [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md) | – | ✓ | – | – |
| [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) | – | – | ✓ | – |

The universal entry reaches any backend: name one with `backend =`, or
let `backend = "auto"` choose. Each `fb_<engine>()` pin fits exactly one
engine and refuses a conflicting `backend`. `backend = "auto"` runs the
LGM feasibility gate and routes to INLA on acceptance, otherwise to
greta; it never routes to Stan, because brms’s first-call Stan compile
(typically 30–60 s) is incompatible with the auto-fast-path promise.
Reach Stan explicitly with
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
or `fb(..., backend = "brms")`.

## Backend support

flexyBayes is standalone-functional with any one backend installed, and
the planner (see *Quick start*) needs none at all. The backends differ
in install burden and in what they offer.

| Backend | On CRAN? | Install burden | Inference | In flexyBayes |
|----|----|----|----|----|
| INLA | No (own repository) | Moderate — binary, no compiler | Approximate (nested Laplace) | Supported |
| brms (Stan) | Yes | Heavy — first-call Stan compile (~30–60 s) | MCMC (sampling error only) | Supported |
| greta | No (greta-dev R-universe) | Heavy — Python + TensorFlow stack | MCMC (sampling error only) | Supported |
| gretaR | No (opt-in) | Heavy — torch | MCMC (out-of-process NUTS) | Dormant, opt-in |

All exports are at the **experimental** `lifecycle` stage; see
`API_STABILITY.md` in the source repository for what that guarantees.
The fastest way to explore the package without any backend is the
planner; for a worked fit with production sampling settings and honest
diagnostics, follow the *getting started* vignette.

### Backend support by model class

What each backend does, by model class, reflecting the 2026-06 empirical
validation. “Supported” means the syntax emits and the model is tested;
it does not promise that every fit converges at small budgets –
convergence is model-specific and is always reported, so treat a high
R-hat badge as a diagnostic, not a result.

| Model class | greta | INLA | brms |
|----|----|----|----|
| Gaussian LMM, simple random intercepts | validated | validated | validated |
| GLMM (binomial / Poisson / NB), simple RE | validated | validated | validated |
| Random slopes, structured covariance, GBLUP, pedigree, separable spatial, splines | supported | supported | supported |
| **Interaction / nested random effects (GxE)** | does **not** converge | **refused** | **refused** |
| **Heteroscedastic / per-stratum residual** (`dsum`) | does **not** converge | **refused** | **refused** |
| Interaction *fixed* effects, binomial path | supported | **known bug** | supported |

> **Capability boundary (read before relying on flexyBayes for METs).**
> The full multi-environment-trial model – genotype-by-environment
> *interaction random effects* together with a *heteroscedastic
> per-environment residual* – is **not fittable by any flexyBayes
> backend in this release**. ASReml and `lme4` fit it; flexyBayes does
> not yet. This is the package’s central open problem. See
> `inst/KNOWN_ISSUES.md`
> (`system.file("KNOWN_ISSUES.md", package = "flexyBayes")`) for the
> full status, the per-backend reasons, and the intended route to a fix.
> Until then, flexyBayes is an honest multi-backend orchestration and
> triangulation tool for the model classes above – not yet a stand-alone
> Bayesian replacement for full ASReml-style MET analysis.

**Breeder MET summaries.**
[`fb_met_summary()`](https://aagi-aus.github.io/flexyBayes/reference/fb_met_summary.md)
(overall performance, stability, GxE BLUPs, factor loadings, environment
genetic correlations) requires a **greta** factor-analytic
(`fa(env, k):gen`) fit — it is computed from the identified *realised*
effects, which fit on greta. On an INLA or brms fit it refuses with a
pointer to the right path; the scalable INLA MET route gives variance
components via [`summary()`](https://rdrr.io/r/base/summary.html) /
[`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md).

## Installation

``` r

# greta (R-native MCMC backend) — archived from CRAN, served via the
# greta-dev R-universe
install.packages("greta",
  repos = c("https://greta-dev.r-universe.dev", getOption("repos")))
greta::install_greta_deps()

# INLA (approximate-inference backend) — not on CRAN
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"))

# brms (Stan passthrough) — on CRAN
install.packages("brms")

# flexyBayes itself (not yet on CRAN) -- install from the repository:
# install.packages("remotes")
remotes::install_github("AAGI-AUS/flexyBayes")
```

`flexyBayes` degrades gracefully when an optional engine is missing:
each backend is detected at run time, and a model sent to an unavailable
engine is refused with a clear message naming what to install rather
than failing obscurely.

## Quick start

The planner needs no inference backend and is the fastest way to see
what flexyBayes will do with a model: it builds the intermediate
representation, chooses a backend, and reports the plan without fitting.

``` r

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
following uses production sampling settings; the *getting started*
vignette walks through the same fit with its convergence diagnostics.

``` r

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

``` r

fit_g <- flexybayes(Reaction ~ Days, random = ~ Subject,
                    data = sleepstudy, backend = "greta")
fit_i <- flexybayes(Reaction ~ Days, random = ~ Subject,
                    data = sleepstudy, backend = "inla")

triangulate(fit_g, fit_i)
#> Common parameters: 3
#> # A tibble: 3 × 5
#>   parameter   mean_diff sd_ratio tail_drift wasserstein_1
#>   <chr>           <dbl>    <dbl>      <dbl>         <dbl>
#> 1 (Intercept)    -0.04     0.997     0.012          0.04
#> 2 Days            0.001    1.001     0.003          0.001
#> 3 sigma           0.12     1.004     0.020          0.12
```

Add a third engine via brms / Stan:

``` r

fit_s <- fb_brms(Reaction ~ Days + (1 | Subject), data = sleepstudy)

triangulate(fit_g, fit_s)
triangulate(fit_i, fit_s)
```

[`canonical_names()`](https://aagi-aus.github.io/flexyBayes/reference/canonical_names.md)
does the work of aligning backend-native parameter names (greta’s
`sigma_g`, INLA’s `Precision for g` on the precision scale, brms’s
`sd_g__Intercept`) to a single canonical name with the correct scale
transform — no `name_map` argument needed in standard cases.

## Companion accessors

| Accessor | Returns |
|----|----|
| `backend_decision(fit)` | The captured dispatch trace: backend, path, `lgm_gate` checks, reason. |
| `prior_summary(fit)` | The resolved prior — auto-default (weakly-informative bounded uniform on SD; half-Cauchy advised for small group counts), user-supplied [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md), or legacy scalar bridge. |
| `canonical_names(fit)` | The backend-native ↔︎ canonical-name table with per-row scale transforms. |
| `review_code = TRUE` on [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) / [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) | Inspect-before-fit workflow; `cat_code(rev)` prints the generated backend code; `proceed(rev)` advances into the fit. Supported on the formula-entry verbs only. [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md) raises a documented refusal: the user’s model is already greta-side, so there is no flexyBayes-generated code string to inspect (model-graph pretty-printer reserved for a future release). |

## Output structure

Every fit carries three (or four) top-level slots:

``` r

fit$glm         # GLM-compatible shim — works with summary(), emmeans,
                # marginaleffects, effectsize, broom

fit$greta       # native greta output (when backend = "greta")
fit$inla        # native INLA output (when backend = "inla")
fit$brms        # live brmsfit (when backend = "brms")

fit$extras      # BLUPs, variance components, convergence diagnostics,
                # generated code, parsed IR, run time, captured call
```

## Supported ASReml syntax (reference)

``` r

# Fixed effects
yield ~ env                  # fixed factor
yield ~ env + x_cov          # factor + covariate
yield ~ 0 + env              # means model (no intercept)
yield ~ env + I(x^2)         # expression terms

# Random effects
random = ~ geno                       # simple iid
random = ~ block:rep                  # nested
random = ~ vm(geno, Gmat)             # GBLUP (dense V; both backends)
random = ~ vm(geno, chol = L)         # user-supplied Cholesky (greta only; v0.3.7+)
random = ~ vm(geno, precision = Q)    # user-supplied sparse precision (INLA preferred; v0.3.7+)
random = ~ ped(animal, Amat)          # pedigree (animal model)
random = ~ ped(animal, A_inv,
               use_sparse_precision = TRUE)  # sparse pedigree precision (v0.3.7+)
random = ~ at(env):geno               # diagonal GxE
random = ~ us(env):id(geno)           # unstructured GxE
random = ~ fa(env, 2):id(geno)        # factor-analytic GxE
random = ~ ar1(row):id(col)           # spatial AR1
random = ~ spl(x_cov)                 # P-spline

# Residual
rcov = ~ units                # iid residuals (default)
rcov = ~ at(env):units        # heterogeneous by environment

# Families
family = "gaussian" | "binomial" | "poisson" | "negative_binomial" |
         "gamma" | "beta"
```

## Vignettes

Sixteen numbered vignettes ship with the package:

1.  Getting started
2.  ASReml-shaped formulas
3.  Foundational regression
4.  Hierarchical models
5.  Structured covariance
6.  Priors and regularisation
7.  METs and genomic selection
8.  Downstream analysis
9.  Spatio-temporal models
10. Cross-engine triangulation
11. LGM feasibility filter
12. Backend internals
13. LGM feasibility and memory
14. Choosing an engine: the universal entry and the engine pins
15. Architecture: the backend registry and how a new engine joins
16. Big-data streaming

Heavy MCMC vignettes use a `.Rmd.orig` precompile pattern; the `.Rmd`
that ships in the package tarball is the pre-evaluated static output.
Browse them with `browseVignettes("flexyBayes")` **after a full
install** — `R CMD build` then `R CMD INSTALL` the tarball, or
`devtools::install(build_vignettes = TRUE)`. A plain `install_github()`
or source-directory install does **not** build the vignettes into
`inst/doc`.

Several reference vignettes (02, 04, 05, 07) fit at small sampling
budgets and print diagnostics that do not meet the convergence
thresholds; each opens with an “illustration of output *shape*, not
inference” disclaimer. For a convergence-clean workflow see the *getting
started* and *cross-engine triangulation* vignettes.

## Design principles

- **Indexing, not model matrices**: parameters are sized by number of
  levels (`p`), not observations (`N`), giving O(p) memory for random
  effects.
- **Non-centred parameterisation**: all random effects use NCP for
  efficient MCMC sampling.
- **Cholesky decomposition**: known covariance matrices (G, A) are
  decomposed once.
- **Triangulation as first-class output**: the package is built around
  the claim that *evidence of cross-engine agreement* is the signature
  flexyBayes deliverable, not “one more Bayesian mixed- model frontend”.

## Correctness

flexyBayes ships an extensive `testthat` suite (`devtools::test()`)
covering the paths most likely to hide errors: fixed-effect and factor
models against [`lm()`](https://rdrr.io/r/stats/lm.html) /
`lme4::lmer()` references, random intercepts and the asreml-route random
slopes, structured-covariance terms, streaming-aggregation equivalence
to the per-row fit, weights, offsets, missing-response handling, backend
routing and the structured refusal taxonomy, prior translation, and
cross-engine
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
agreement. The *cross-engine triangulation* and per-family vignettes
show clean reproducible checks; the structured-covariance reference
vignettes are explicitly marked non-inferential where a model does not
mix at vignette-scale budgets.

## Testing & CI

Continuous integration validates the INLA, brms, and engine-independent
surface (the ASReml / brms parsers, the intermediate representation,
`lgm_gate()`, the dispatch policy table, the refusal registry, the prior
DSL, and the
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
metrics). The **greta** integration path (Python + TensorFlow) is
**not** installed on the standard CI matrix, so a green CI run does not
include live greta fits; greta-dependent tests gate via
`skip_if_no_greta()` and are run locally / via a manual
greta-integration check. Run the full suite locally with
`devtools::test()` after `greta::install_greta_deps()`.

## Known limitations

flexyBayes refuses what it cannot yet fit rather than fitting it
silently. The current release does not cover the following; each is a
roadmap deferral, and a request that needs one is met with a structured
refusal naming the gap, not a quiet wrong answer.

- **Full multi-environment-trial (MET) models**: genotype-by-environment
  *interaction random effects* (`gen:loc`, `gen:loc:yearf`) and
  *heteroscedastic per-environment residuals* (`dsum(~ units | env)`)
  are **not fittable by any trustworthy backend in this release** – INLA
  and brms refuse them, greta expresses but does not converge on them.
  This is the package’s central open problem; see `inst/KNOWN_ISSUES.md`
  for the full status and the intended route to a fix. Use ASReml or
  `lme4` for the full MET model today.
- **Spatial structure**: only separable AR1 (`ar1(row):id(col)`) is
  supported. Intrinsic CAR and BYM2 areal models are not implemented.
  You can express a custom spatial precision by passing your own matrix
  to `vm(g, precision = Q)`, but there is no BYM2 helper.
- **Smooth terms**: univariate penalised splines (`s(x)`, `spl(x)`) are
  supported. Multivariate and tensor-product smooths (`te()`, `ti()`,
  `t2()`) are refused and deferred to a later release.
- **Hidden-Markov, multi-state, and survival models**: not supported.
  Survival / time-to-event families are refused at the family gate. A
  NIMBLE backend covering these is on the roadmap with no fixed release
  target.
- **Missing data**: flexyBayes does not impute. Greta and the brms
  (Stan) passthrough require complete cases for the model variables;
  INLA treats an `NA` response as a prediction target. Resolve
  missingness (drop or impute) before fitting, or use the prediction
  path deliberately.

## Requirements

- R ≥ 4.1.0
- greta ≥ 0.4.0 (Python + TensorFlow + TensorFlow Probability)
- INLA (optional, via Additional_repositories)
- brms (optional, for the Stan passthrough)

## Contributing

Contributions are welcome. See `CONTRIBUTING.md` for the development
workflow (fork, `usethis::pr_init()`, `devtools::check()` must pass, and
a `NEWS.md` bullet for any user-facing change) and `CODE_OF_CONDUCT.md`
(Contributor Covenant). Please report bugs and request features at
<https://github.com/AAGI-AUS/flexyBayes/issues>. The package’s
architecture decisions are indexed in `DESIGN_DECISIONS.md`.

## Citation

If you use `flexyBayes` in your research, please cite:

    @software{flexyBayes,
      title  = {flexyBayes: Flexible Bayesian Mixed Models with ASReml
                and brms-Style Syntax},
      author = {Moldovan, Max and Tanaka, Emi and Hui, Francis K.C. and
                Forte Deltell, Anabel},
      year   = {2026},
      url    = {https://github.com/AAGI-AUS/flexyBayes}
    }

## License

GPL (≥ 3)
