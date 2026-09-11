# flexyBayes

Fit linear and generalised linear mixed models from ASReml-style or
brms-style formulas, on two engines: INLA directly, and Stan through brms.

The ASReml grammar reaches INLA natively rather than through a sampler, so
the models designed experiments actually use -- multi-stratum random effects,
a residual variance per environment, unstructured and factor-analytic
genotype-by-environment covariance -- fit by approximation in seconds rather
than by MCMC in minutes. Where a model is outside what an engine can
represent, the package says so by name instead of fitting something else.

> **Development version 0.10.0.9000.** Not on CRAN. Every export is at the
> experimental `lifecycle` stage and the API may change within the 0.x
> series.

## Installation

``` r
# INLA (approximate-inference engine) -- not on CRAN
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"))

# brms (Stan) -- on CRAN
install.packages("brms")

# flexyBayes -- not yet on CRAN
# install.packages("remotes")
remotes::install_github("AAGI-AUS/flexyBayes", build_vignettes = TRUE)
```

`build_vignettes = TRUE` matters: `remotes::install_github()` does not build
them by default, and without it `vignette(package = "flexyBayes")` comes back
empty.

Either engine is optional. A model sent to an engine that is not installed is
refused with a message naming what to install.

## A worked example

An alpha-design field trial: 24 genotypes in three replicates, blocked within
replicate. The `fixed` / `random` split and the `rep:block` nesting are
written as they are in ASReml.

``` r
library(flexyBayes)
data(john.alpha, package = "agridat")

fit <- flexybayes(
  fixed  = yield ~ rep,
  random = ~ rep:block + gen,
  data   = john.alpha,
  backend = "inla"
)

summary(fit)$varcomp[, c("component", "estimate", "conf.low", "conf.high")]
#>     component  estimate  conf.low conf.high
#>         sigma 0.2926140 0.2289825 0.3715715
#>  sd_rep:block 0.2880199 0.1632093 0.4677316
#>        sd_gen 0.3973759 0.2736240 0.5610492
```

The three standard deviations sit on `lme4::lmer()`'s REML estimates for the
same model (0.2857, 0.2650, 0.3780), and each arrives with an interval rather
than a standard error.

Because the variance components are draws, any function of them has an exact
posterior -- ASReml's `vpredict()` question, answered without the delta
method:

``` r
vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2))
#>  functional  estimate    median std.error  conf.low conf.high
#>          h2 0.6328936 0.6462664 0.1020417 0.4596090 0.8002571
```

That interval cannot leave `[0, 1]`, because no draw of it can. A
delta-method interval on the same quantity has no such guarantee.

## What it fits

One table, generated from the package's own capability object and checked
against the gate and emit code on every test run. `fits` means an engine
emits the structure and a test exercises it; `refuses` means the request
raises rather than fitting something else.

<!-- capability-matrix:begin -->
| Model class | Spelling | INLA | brms | Notes |
|---|---|:-:|:-:|---|
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | fits | fits | The class both engines emit, and the one their agreement is checked on in the test suite. |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | fits | fits | INLA's likelihood allowlist is read from `INLA::inla.models()` when INLA is installed. |
| Hurdle gamma (zero mass plus a positive gamma part) | `family = "hurdle_gamma"` | refuses | fits | brms-native (`dpars` mu, shape, hu); the zero-mass probability `hu` keeps brms's own prior. INLA's likelihood roster carries no counterpart, so the family gate refuses it there and `auto` routes to brms. |
| Uncorrelated random slope | `(x \|\| g)` | fits | fits | Two independent iid fields on INLA, one for the intercept and one for the slope, checked against lme4's REML fit of the same model. Reported as `sd_<g>` and `sd_<x>_<g>`. |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | fits | fits | An ordinary treatment-coded design matrix on both engines. The INLA posterior means sit on the least-squares estimates of the same model under a weak prior. |
| Correlated random slope | `(x \| g)` | fits | fits | The uncorrelated slope with the correlation freed: one intercept SD, one slope SD and the correlation between them. INLA holds all three in a single iid2d field over a two-block index, with the slope block entering as a fixed-scaling copy weighted by `x`; brms writes `(x \| g)` verbatim. Both report `sd_<g>`, `sd_<x>_<g>` and `cor_<g>__Intercept__<x>`. The iid2d hyperparameters are one joint Wishart with no standard-deviation slot, so on INLA an explicit SD prior refuses by name and the package default is reported as an engine default rather than as a prior that reached the fit. The joint Wishart's scale matrix is set from the response's own variance, and the slope block from the response variance per unit of the slope variable. Leaving INLA's default there -- a scale matrix of ones, so a prior centred on unit variances -- makes the mode non-unique on any response that is not near unit scale: on `lme4::sleepstudy` eight identical calls returned intercept SDs from 0.5 to 54.8 and correlations from -1.00 to +0.75. Scaled, the same fit returns 25.7 to 26.9 against lme4's 24.7, with the correlation between -0.07 and +0.15 around lme4's 0.066. |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | fits | fits | An iid field over the interaction index on INLA, `(1 \| a:b)` on brms -- the same model, and the two agree with lme4's REML fit of it on `agridat::john.alpha`. Only observed cells get a level, as in lme4 and brms. |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | fits | fits | One variance per level of `f`, no covariance between levels. All three spellings emit identical code. INLA writes k independent iid fields, one per level; brms writes `(0 + f \|\| g)`. Both report the component as `sd_<g>__<f><level>` and agree with lme4's REML fit of the same model. `at(f, level)` conditions on selected levels -- a different model -- and still refuses. |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | fits | fits | One iidKd field on INLA (k = 2..5; wider refuses and names fa()), `(0 + f \| g)` on brms. INLA's Wishart hyperprior shrinks the correlations toward each other, and an explicit SD prior refuses there because iidKd has no slot for one. The correlated sibling of the diagonal structure -- `k(k+1)/2` parameters against `diag()`'s `k`. At one observation per cell the residual variance is confounded with the diagonal of the covariance: the covariance block converges and `sigma` does not, and a longer chain does not help. Replicate within cell, or put an informative prior on the residual. |
| Heterogeneous variances with one shared correlation | `~ corh(f):g` | refuses | refuses | No active engine has an equicorrelation group-level structure. Use `diag(f):g` or `us(f):g`. |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | fits | fits | One residual variance per section. INLA carries it through the multiple-likelihood interface, a matrix response with one Gaussian family per level; brms lowers it to distributional regression on log sigma, `sigma ~ 0 + f`. Both report `sigma_<level>` on the standard-deviation scale. Refused for families with no residual scale; `at(f, level):units` conditions on selected levels, which is a different model, and refuses. |
| Combined interaction random effects and heterogeneous residual (full MET) | `random = ~ gen + gen:env` with the `dsum` residual | fits | fits | The full multi-environment-trial model. Both halves refused on INLA before 0.10.0.9000, so the combination was unreachable there. It now recovers the genotype, genotype-by-environment and all three per-section residual standard deviations on simulated data. |
| Factor-analytic genotype-by-environment covariance | `~ fa(env, k):gen` | fits | refuses | Each common factor is one latent field entered per level of `f` as a SCALED COPY, and INLA's copy scaling IS the loading, so the loadings carry a posterior rather than a point estimate. Factor j is anchored at level j, the zero upper triangle ASReml uses against rotation, so an anchored loading is absent from the summary rather than reported as 1. Verified at k = 1 and k = 2. brms group-level effects are either uncorrelated or fully unstructured, with no factor-analytic form between. |
| Multi-trait covariance | `~ us(trait):vm(gen)` | refuses | refuses | No active engine represents a trait-by-genotype unstructured covariance. |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits | INLA takes the sparse-precision, pedigree-precision and block carriers, and brms additionally takes dense and Cholesky. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses | A latent AR1 field plus the Gaussian observation nugget -- four hyperparameters, one observation per grid node. This is not ASReml's three-parameter nugget-free residual, so the residual spelling refuses and names this one. |
| Per-trial separable AR1 field | `random = ~ at(trial):ar1(row):ar1(col)` | fits | refuses | One field realisation per level of `trial`, via INLA's `replicate =` mechanism, but the row correlation, column correlation and field SD are shared across every level -- not estimated per trial. `at(trial, level):ar1(row):ar1(col)` (a level argument, asking for a single conditioned trial or for per-trial hyperparameters) refuses by name (`at_field_per_level_hyper_not_representable`). brms has no lowering for either spelling. |
| Grouped temporal AR1 field | `random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)` | fits | refuses | AR1 along `t` within each level of `g` and independent across levels: one correlation, one field SD and the Gaussian nugget, three hyperparameters shared by every level. The units-by-occasions repeated-measures covariance of Verbyla and Cullis (1992), written on the random side; an occasion a unit missed stays as a design row with the response NA under `na_action = "augment"`. Either factor order is accepted. brms has no lowering for the grouped autoregressive precision and refuses. |
| Univariate P-spline | `~ spl(x)` | fits | refuses | Mapped to INLA's second-order random walk. brms has no lowering for the smooth basis. An analogue of ASReml's spline rather than a translation of it -- basis, scaling and prediction rule differ, and no equivalence has been established (see TRANSLATION.md). |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits | Precision weighting, Var(y_i) = sigma^2 / w_i (the ASReml / lme4 / glm(weights=) sense): INLA's `scale = w`; on brms a known offset on the sigma distributional parameter, NOT brms's own `weights()` addition term (a different, likelihood-power quantity per brms's own documentation). Both engines match lme4::lmer(weights=) closely on a shared simulated fixture. Any other family, or a non-identity link on Gaussian, refuses by name (`weights_requires_gaussian`); `aggregate = TRUE` alongside weights also refuses by name (`weights_not_aggregatable`). |
| Exact sufficient-statistic aggregation | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a | Exact cell-likelihood aggregation for iid exponential-family models with small cell count. The brms path has no aggregated emit. |

`fits` -- the engine emits the structure and a test exercises it. `emits` -- the engine generates the structure and no live fit has yet confirmed it samples acceptably. `refuses` -- the request raises rather than fitting something else. `n/a` -- the class does not apply to that engine's interface.

This block is generated from `.fb_capability_matrix()` by `tools/generate_capability_matrix.R`. Edit the R table, re-run the generator, and let `tests/testthat/test-capability-matrix.R` check that every verdict still matches the gate and emit code. Do not edit the rows here by hand.
<!-- capability-matrix:end -->

The full ASReml surface -- every spelling the parser accepts, and what each
one is for -- is in `vignette("flexyBayes-03-syntax")`. What no engine fits, and
why, is in `system.file("KNOWN_ISSUES.md", package = "flexyBayes")`.

The table above says whether an engine fits a structure. It does not say what
the lowering promises, and the two are different questions: some structures are
exact structural translations of the ASReml model, some share a covariance
family but not the identification or the priors, and one is an analogue built
from a different basis. Which is which, and what independently grounds each
row, is in `system.file("TRANSLATION.md", package = "flexyBayes")`.

## Documentation

Four vignettes, meant to be read in order and to be enough to start work:

| Vignette | What it covers |
|---|---|
| `flexyBayes-01-mixed-models` | Mixed models in brief: strata, variance components, what is random and why |
| `flexyBayes-02-bayesian` | Bayesian estimation in brief: priors, posteriors, and what changes in the output |
| `flexyBayes-03-syntax` | Every ASReml spelling the package covers, what it is for, and what it does not cover |
| `flexyBayes-04-illustration` | A worked analysis with ASReml and flexyBayes side by side, and what the posterior gives that REML cannot |

## Reporting problems

Bugs and feature requests: <https://github.com/AAGI-AUS/flexyBayes/issues>.
A reproducible example and `sessionInfo()` help most. `CONTRIBUTING.md`
covers the development workflow.

## Citation

``` 
@software{flexyBayes,
  title  = {flexyBayes: Flexible Bayesian Mixed Models with ASReml
            and brms-Style Syntax},
  author = {Moldovan, Max},
  year   = {2026},
  url    = {https://github.com/AAGI-AUS/flexyBayes}
}
```

Licence: MIT.
