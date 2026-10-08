<!-- README.md is generated from README.Rmd with knitr::knit(). Edit this file. -->



# flexyBayes

flexyBayes fits linear and generalised linear mixed models written in ASReml
or brms syntax. It fits them with INLA, or with Stan through brms. If the
chosen engine cannot fit a model, flexyBayes stops and says why.

Version 0.16.0. Not on CRAN, and the functions may still change.

## Installation

``` r
install.packages(c("remotes", "knitr", "rmarkdown"))
remotes::install_github("AAGI-AUS/flexyBayes", build_vignettes = TRUE)
```

Keep `build_vignettes = TRUE`. Without it, `vignette(package = "flexyBayes")`
lists nothing. Building the vignettes needs knitr, rmarkdown and pandoc
(RStudio includes pandoc).

You need at least one engine. brms is on CRAN. INLA is not:

``` r
install.packages("brms")
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"))
```

## Example

An alpha design from agridat: 24 genotypes in three replicates, with blocks
nested in replicates. The call is written as it would be in ASReml.


``` r
library(flexyBayes)
data(john.alpha, package = "agridat")

fit <- flexybayes(
  fixed   = yield ~ rep,
  random  = ~ rep:block + gen,
  data    = john.alpha,
  backend = "inla",
  verbose = FALSE
)

summary(fit)$varcomp[, c("component", "estimate", "conf.low", "conf.high")]
#>      component estimate conf.low conf.high
#> 1        sigma    0.293    0.229     0.371
#> 2 sd_rep:block    0.288    0.163     0.467
#> 3       sd_gen    0.397    0.274     0.561
```



The estimates are posterior means of the three standard deviations, with 95%
credible intervals. REML estimates from `lme4::lmer()` for the same model are
0.286, 0.265 and 0.378.

Heritability is computed on each posterior draw, so its interval stays
between 0 and 1:


``` r
vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2))
#>   functional estimate median std.error conf.low conf.high
#> 1         h2    0.636  0.639       0.1    0.429     0.821
```

## What it fits

The table lists the models at least one engine fits. Notes on each row, and
the models neither engine fits, are in
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
| Per-trial separable AR1 field | `random = ~ id(trial):ar1(row):ar1(col)` | fits | refuses |
| Grouped temporal AR1 field | `random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)` | fits | refuses |
| Univariate P-spline | `~ spl(x)` | fits | refuses |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits |
| Exact sufficient-statistic aggregation (binomial, Poisson) | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a |

`fits`: the engine fits it, and a test checks the fit. `refuses`: flexyBayes stops with a message saying why and what to use instead. `n/a`: the request does not apply to that engine.
<!-- capability-matrix:end -->

## Documentation

| Vignette | Content |
|---|---|
| `flexyBayes-01-mixed-models` | Mixed models in brief |
| `flexyBayes-02-bayesian` | Bayesian estimation in brief |
| `flexyBayes-03-syntax` | The ASReml syntax flexyBayes covers |
| `flexyBayes-04-illustration` | ASReml and flexyBayes on the same trials |
| `flexyBayes-06-two-stage` | A two-stage analysis of a series of trials |
| `flexyBayes-07-from-asreml` | Moving from ASReml-R |

## Reporting problems

Please report bugs at <https://github.com/AAGI-AUS/flexyBayes/issues>, with a
reproducible example and the output of `sessionInfo()`.

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
