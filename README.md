# flexyBayes

Fit linear and generalised linear mixed models from ASReml-style or
brms-style formulas, on two engines: INLA directly, and Stan through brms.

The ASReml grammar reaches INLA natively rather than through a sampler, so
the models designed experiments actually use -- multi-stratum random effects,
a residual variance per environment, unstructured and factor-analytic
genotype-by-environment covariance -- fit by approximation in seconds rather
than by MCMC in minutes. Where a model is outside what an engine can
represent, the package says so by name instead of fitting something else.

> **Development version 0.12.0.** Not on CRAN. Every export is at the
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

## Selecting over a candidate block

`fb_select()` protects the model you already have and asks what a block of
candidate fixed effects adds to it. The candidates carry a joint shrinkage
prior, the covariance structure is integrated over rather than plugged in,
and the submodel is chosen by projection predictive inference on folds cut at
the cluster.

``` r
sel <- fb_select(
  fit,
  candidates = ~ rain + clay + elevation + ph,
  penalty    = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
  resamples  = fb_resamples(unit = ~ site, v = 10L),
  decision   = fb_decision(rule = "one_se"),
  seed       = 20260912L
)
```

Posterior draws under a continuous shrinkage prior contain no exact zeros, so
"selected" is a named decision and never a property of the draws. Every
result records four independent axes -- what was fitted, what fitted it,
under what prior, and how the answer was read off -- and `benchmark = TRUE`
adds a `glmnet` arm and a generalised-least-squares arm on the same folds,
each carrying one of three declared equivalence levels. The words "same as
glmnet" appear nowhere unqualified.

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
| Uncorrelated random slope | `(x \|\| g)` | fits | fits | Two independent iid fields on INLA, one for the intercept and one for the slope, checked against lme4's REML fit of the same model. Reported as `sd_<g>` and `sd_<x>_<g>` on INLA and as `sd_<g>` and `sd_<g>__<x>` on brms; triangulate() maps both to the INLA spelling. |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | fits | fits | An ordinary treatment-coded design matrix on both engines. The INLA posterior means sit on the least-squares estimates of the same model under a weak prior. |
| Correlated random slope | `(x \| g)` | refuses | fits | The uncorrelated slope with the correlation freed: one intercept SD, one slope SD and the correlation between them. brms writes `(x \| g)` verbatim and reports `sd_<g>`, `sd_<g>__<x>` and `cor_<g>__Intercept__<x>`. Estimated correlation matrices are fitted on brms only: INLA refuses by name (`correlated_covariance_requires_brms`) and `backend = "auto"` routes to brms. The INLA route was one iid2d field under a joint Wishart prior, and paired against ASReml 4.2.0.392 and lme4 on `lme4::sleepstudy` it returned a slope variance of 46.4 at the posterior mode against 35.08 and 35.07, while rescaling the response by ten moved the intercept-slope correlation from 0.02 to 0.969 where REML holds at 0.066 (oracle row 19, defect D1 of finding S01). |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | fits | fits | An iid field over the interaction index on INLA, `(1 \| a:b)` on brms -- the same model, and the two agree with lme4's REML fit of it on `agridat::john.alpha`. Only observed cells get a level, as in lme4 and brms. |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | fits | fits | One variance per level of `f`, no covariance between levels. All three spellings emit identical code. INLA writes k independent iid fields, one per level; brms writes `(0 + f \|\| g)`. Both report the component as `sd_<g>__<f><level>` and agree with lme4's REML fit of the same model. `at(f, level)` conditions on selected levels -- a different model -- and still refuses. |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | refuses | fits | `(0 + f \| g)` on brms. Estimated correlation matrices are fitted on brms only: INLA refuses by name (`correlated_covariance_requires_brms`) and `backend = "auto"` routes to brms. The INLA route was one iidKd field under a joint Wishart prior, and paired against ASReml 4.2.0.392 on `agridat::besag.met` it returned a genetic correlation of 0.142 where ASReml returns 0.930, while eight identical calls on the response rescaled by ten none converged and returned correlations from -0.894 to +0.686 (finding O04). The correlated sibling of the diagonal structure -- `k(k+1)/2` parameters against `diag()`'s `k`. At one observation per cell the residual variance is confounded with the diagonal of the covariance: the covariance block converges and `sigma` does not, and a longer chain does not help. Replicate within cell, or put an informative prior on the residual. |
| Heterogeneous variances with one shared correlation | `~ corh(f):g` | refuses | refuses | No active engine has an equicorrelation group-level structure. Use `diag(f):g` or `us(f):g`. |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | fits | fits | One residual variance per section. INLA carries it through the multiple-likelihood interface, a matrix response with one Gaussian family per level; brms lowers it to distributional regression on log sigma, `sigma ~ 0 + f`. Both report `sigma_<level>` on the standard-deviation scale. Refused for families with no residual scale; `at(f, level):units` conditions on selected levels, which is a different model, and refuses. |
| Combined interaction random effects and heterogeneous residual (full MET) | `random = ~ gen + gen:env` with the `dsum` residual | fits | fits | The full multi-environment-trial model. Both halves refused on INLA before 0.10.0.9000, so the combination was unreachable there. It now recovers the genotype, genotype-by-environment and all three per-section residual standard deviations on simulated data. |
| Factor-analytic genotype-by-environment covariance | `~ fa(env, k):gen` | fits | refuses | Each common factor is one latent field entered per level of `f` as a SCALED COPY, and INLA's copy scaling IS the loading, so the loadings carry a posterior rather than a point estimate. Factor j is anchored at level j, the zero upper triangle ASReml uses against rotation, so an anchored loading is absent from the summary rather than reported as 1. Verified at k = 1 and k = 2. brms group-level effects are either uncorrelated or fully unstructured, with no factor-analytic form between. |
| Multi-trait covariance | `~ us(trait):vm(gen)` | refuses | refuses | No active engine represents a trait-by-genotype unstructured covariance. |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits | INLA takes the sparse-precision, pedigree-precision and block carriers, and brms additionally takes dense and Cholesky. Under `fb_select()`: shrinkage yes, projection no. A joint shrinkage prior over candidate fixed effects fits on a base carrying one of these terms, but the term lowers to `(1 \| gr(g, cov = K))` and projpred refuses a function call on the right-hand side of a group-term `\|` (probe C5), so the projection decision is red. `decision = "auto"`, the default, resolves to `full` there; naming `"projection"` refuses before the reference model is sampled. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses | A latent AR1 field plus the Gaussian observation nugget -- four hyperparameters, one observation per grid node. This is not ASReml's three-parameter nugget-free residual, so the residual spelling refuses and names this one. |
| Per-trial separable AR1 field | `random = ~ at(trial):ar1(row):ar1(col)` | fits | refuses | One field realisation per level of `trial`, via INLA's `replicate =` mechanism, but the row correlation, column correlation and field SD are shared across every level -- not estimated per trial. `at(trial, level):ar1(row):ar1(col)` (a level argument, asking for a single conditioned trial or for per-trial hyperparameters) refuses by name (`at_field_per_level_hyper_not_representable`). brms has no lowering for either spelling. |
| Grouped temporal AR1 field | `random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)` | fits | refuses | AR1 along `t` within each level of `g` and independent across levels: one correlation, one field SD and the Gaussian nugget, three hyperparameters shared by every level. The units-by-occasions repeated-measures covariance of Verbyla and Cullis (1992), written on the random side; an occasion a unit missed stays as a design row with the response NA under `na_action = "augment"`. Either factor order is accepted. brms has no lowering for the grouped autoregressive precision and refuses. |
| Univariate P-spline | `~ spl(x)` | fits | refuses | Mapped to INLA's second-order random walk. brms has no lowering for the smooth basis. An analogue of ASReml's spline rather than a translation of it -- basis, scaling and prediction rule differ, and no equivalence has been established (see TRANSLATION.md). |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits | Precision weighting, Var(y_i) = sigma^2 / w_i (the ASReml / lme4 / glm(weights=) sense): INLA's `scale = w`; on brms a known offset on the sigma distributional parameter, NOT brms's own `weights()` addition term (a different, likelihood-power quantity per brms's own documentation). Both engines match lme4::lmer(weights=) closely on a shared simulated fixture. Any other family, or a non-identity link on Gaussian, refuses by name (`weights_requires_gaussian`); `aggregate = TRUE` alongside weights also refuses by name (`weights_not_aggregatable`). |
| Exact sufficient-statistic aggregation (binomial, Poisson) | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a | Exact cell-likelihood aggregation for binomial (logit) and Poisson (log) models with a small cell count; the posterior is the per-row posterior. The marginal likelihood is not reported (`aggregated_evidence_unavailable`): the engine's value omits a constant that depends on the data. The brms path has no aggregated emit. |
| Exact sufficient-statistic aggregation (Gaussian) | `aggregate = TRUE`, `flexybayes_stream(family = "gaussian")` | refuses | n/a | Refused by name (`aggregation_residual_prior_not_composed`). The aggregated Gaussian route composes the within-cell sum of squares with a fixed prior on the residual precision, so a stated or default residual prior would not be the prior fitted. Gaussian models are fitted per row, which is what `aggregate = "auto"` does; `flexybayes_stream(fit = FALSE)` still returns the sufficient statistics. |
| Ridge shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ ridge(scale))` | n/a | fits | A joint shrinkage prior over a declared block of candidate fixed effects, compiled to one coefficient-specific `normal` row per candidate so the block is targeted exactly and no protected coefficient is touched. Only `fb_select()` takes a shrinkage prior, and it fits its reference model on brms whatever engine fitted the base, carrying the base's structure and priors through and recording the re-emission on the certificate. No call asks INLA to fit one, and INLA has no implementation of the block. |
| Laplace shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ laplace(scale))` | n/a | fits | One coefficient-specific `double_exponential` row per candidate. The posterior mean is not a lasso estimate and the draws carry no exact zeros; the alias `lasso()` says so once per session. Fitted on brms, as for ridge. A Laplace prior is Gaussian given one local scale per coefficient, so the model can be written as a latent Gaussian model, but it then carries a hyperparameter for every coefficient, and INLA's numerical integration is built for models with a few hyperparameters. |
| Elastic-net shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ elastic_net(lambda, alpha))` | n/a | fits | The canonical elastic-net density written as a `target +=` contribution over the candidate indices through `stanvars` and `check = FALSE`, so the block is targeted exactly while a protected coefficient keeps its own prior row. Fixed hyperparameters only in this release, because the normalising constant is then irrelevant. Fitted on brms, as for ridge. |
| Regularised horseshoe shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ regularized_horseshoe(p0))` | n/a | fits | brms's native `horseshoe()`, with `par_ratio` set from `p0` and calibrated against the number of resampling units rather than the number of rows. Exact block targeting ONLY when the base's fixed part is intercept-only: it is a class-wide prior, brms refuses a coefficient-specific row beside a class-wide special prior, so a base carrying protected fixed effects beyond the intercept refuses with `penalty_block_targeting_unsupported`. The two ways forward are `laplace()`, `ridge()` or `elastic_net()`, which target the block exactly, or `scope = "all_fixed"`, an opt-in under which the protected coefficients also sit under the global-local prior, the projection search forces them into every submodel, and the certificate records the named approximation `protected_under_global_shrinkage`. A candidate block containing a design column the re-emitted reference cannot standardise -- a factor contrast, an interaction product -- refuses with `penalty_candidate_columns_unscaled`. Fitted on brms, as for ridge; like the Laplace prior, it is Gaussian only given a local scale per coefficient. |
| R2D2 shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ r2d2(mean_R2, prec_R2, cons_D2))` | n/a | fits | brms's native `R2D2()`, whose scales are stated on the explained-variance grammar rather than on a coefficient scale, so `p0` is refused for it. Class-wide, with the same targeting rule as the regularised horseshoe: exact only when the base's fixed part is intercept-only, otherwise `laplace()`, `ridge()` or `elastic_net()`, or `scope = "all_fixed"` with the named approximation. Unstandardised candidate columns refuse as for the horseshoe. Fitted on brms, as for ridge. |

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

Five vignettes, meant to be read in order and to be enough to start work:

| Vignette | What it covers |
|---|---|
| `flexyBayes-01-mixed-models` | Mixed models in brief: strata, variance components, what is random and why |
| `flexyBayes-02-bayesian` | Bayesian estimation in brief: priors, posteriors, and what changes in the output |
| `flexyBayes-03-syntax` | Every ASReml spelling the package covers, what it is for, and what it does not cover |
| `flexyBayes-04-illustration` | A worked analysis with ASReml and flexyBayes side by side, and what the posterior gives that REML cannot |
| `flexyBayes-05-selection` | Selecting over a candidate block: the joint shrinkage prior, the projection decision, and what a benchmark against `glmnet` does and does not establish |

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
