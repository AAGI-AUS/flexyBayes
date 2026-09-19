# flexyBayes

Fit linear and generalised linear mixed models from ASReml-style or
brms-style formulas, on two engines: INLA directly, and Stan through brms.
The table under "What it fits" lists the models each engine fits. A model the
chosen engine cannot fit is refused with a message saying why and what to use
instead.

> **Version 0.14.0.** Not on CRAN. Every export is at the
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

The three standard deviations are close to `lme4::lmer()`'s REML estimates
for the same model (0.2857, 0.2650, 0.3780), and each comes with an interval
rather than a standard error.

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

`fb_select()` takes a fitted model, `base` below, and asks which candidate
fixed effects add to its predictions. The candidates share a joint shrinkage
prior and are compared on folds that hold out whole clusters.
`vignette("flexyBayes-05-selection")` explains the choices, and what a
comparison with `glmnet` does and does not show.

``` r
sel <- fb_select(
  base,
  candidates = ~ rain + clay + elevation + ph,
  penalty    = fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
  resamples  = fb_resamples(unit = ~ site, v = 10L),
  decision   = fb_decision(rule = "one_se"),
  seed       = 20260912L
)
```

## What it fits

The table lists the model classes at least one engine fits. The notes on each
row, and the structures neither engine fits, are in
`system.file("KNOWN_ISSUES.md", package = "flexyBayes")`.

<!-- capability-matrix:begin -->
<!-- Generated from .fb_capability_matrix() in R/capability_matrix.R by tools/generate_capability_matrix.R. Edit the R table and re-run the generator; tests/testthat/test-capability-matrix.R fails on a hand edit here. -->
| Model class | Spelling | INLA | brms |
|---|---|:-:|:-:|
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | fits | fits |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | fits | fits |
| Hurdle gamma (zero mass plus a positive gamma part) | `family = "hurdle_gamma"` | refuses | fits |
| Uncorrelated random slope | `(x \|\| g)` | fits | fits |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | fits | fits |
| Correlated random slope | `(x \| g)` | refuses | fits |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | fits | fits |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | fits | fits |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | refuses | fits |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | fits | fits |
| Baseline heteroscedastic MET: iid genotype, iid genotype-by-environment and a sectioned residual | `random = ~ gen + gen:env` with the `dsum` residual | fits | fits |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses |
| Per-trial separable AR1 field | `random = ~ at(trial):ar1(row):ar1(col)` | fits | refuses |
| Grouped temporal AR1 field | `random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)` | fits | refuses |
| Univariate P-spline | `~ spl(x)` | fits | refuses |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits |
| Exact sufficient-statistic aggregation (binomial, Poisson) | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a |
| Ridge shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ ridge(scale))` | n/a | fits |
| Laplace shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ laplace(scale))` | n/a | fits |
| Elastic-net shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ elastic_net(lambda, alpha))` | n/a | fits |
| Regularised horseshoe shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ regularized_horseshoe(p0))` | n/a | fits |
| R2D2 shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ r2d2(mean_R2, prec_R2, cons_D2))` | n/a | fits |

`fits`: the engine fits it, and a test checks the fit. `refuses`: flexyBayes stops with a message saying why and what to use instead. `n/a`: the request does not apply to that engine.
<!-- capability-matrix:end -->

Every ASReml spelling the package accepts, what each is for, and how closely
each follows ASReml are in `vignette("flexyBayes-03-syntax")`.

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
