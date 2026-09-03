# flexyBayes

Flexible Bayesian Mixed Models with `ASReml` and `brms`-style Syntax

Licence: MIT. Version 0.10.0 is an experimental release: every export is
at the experimental `lifecycle` stage and the API may change within the
0.x series.

`flexyBayes` is a multi-backend Bayesian mixed model framework,
effectively acting as a wrapper to estimate a range of linear and
generalized linear mixed models via INLA (integrated nested Laplace
approximation) or `brms` (wrapping around Stan).

> **Development release.** All exports are at the experimental
> `lifecycle` stage and the API may change within the 0.x series. Not on
> CRAN. **flexyBayes fits on two active engines: brms and INLA.** A
> third native engine was withdrawn entirely in 0.9.3 (see `NEWS.md`);
> naming it, or any other unrecognised backend, now raises an ordinary
> unknown-backend refusal. See
> `system.file("KNOWN_ISSUES.md", package = "flexyBayes")` for the
> current per-backend capability boundaries before relying on results.

## Installation

``` r

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
each backend is detected at run time, and a model sent to an unavailable
engine is refused with a clear message naming what to install rather
than failing obscurely.

## Quick start

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

## Which entry point do I use?

| If you are used… | Then | Notes |
|----|----|----|
| ASReml syntax (`fixed` / `random` / `residual`) | [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) / [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) | Variance-component / Aimed more at analysis of agricultural data. |
| brms/lme4 syntax (`y ~ x + (1 \| g)`) | [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) / [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md) | At more any more other mixed model analysis |
| I don’t mind! | [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md) / [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) | If you want to use a specific backend |

[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
is the short alias for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md),
and either name is the universal entry that spans every backend.

## Backend support

The backends differ in install burden and in what they offer.

| Backend | On CRAN? | Install burden | Inference | In flexyBayes |
|----|----|----|----|----|
| INLA | No (own repository) | Moderate – binary, no compiler | Approximate (nested Laplace) | Supported |
| brms (Stan) | Yes | Heavy – first-call Stan compile (~30–60 s) | MCMC (sampling error only) | Supported |

All exports are at the **experimental** `lifecycle` stage. See
`API_STABILITY.md` in the source repository for what that guarantees.

### Some model types

| Model class | Spelling | INLA | brms | Notes |
|----|----|:--:|:--:|----|
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | ✓ | ✓ |  |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | ✓ | ✓ |  |
| Hurdle gamma | `family = "hurdle_gamma"` | x | ✓ |  |
| Uncorrelated random slope | `(x \|\| g)` | x | ✓ |  |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | x | ✓ |  |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | x | ✓ | brms basically turns this to `(1 \| a:b)`. |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | x | ✓ | One variance per level of `f`, no covariance between levels. |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | x | ✓ |  |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | x | ✓ | Refused for families with no residual scale. |
| Combined interaction random effects and heterogeneous residual (full MET) | `random = ~ gen + gen:env` with the `dsum` residual | x | ✓ |  |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | ✓ | ✓ | INLA takes the sparse-precision, pedigree-precision and block carriers, and brms additionally takes dense and Cholesky forms. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | ✓ | x |  |
| Per-trial separable AR1 field | `random = ~ at(trial):ar1(row):ar1(col)` | ✓ | x | One field realisation per level of `trial`, via INLA’s `replicate =` mechanism, but the row correlation, column correlation and field SD are shared across every level. |
| Univariate P-spline | `~ spl(x)` | ✓ | x |  |
| Observation weights (Gaussian, identity link) | `weights = w` | ✓ | ✓ | Precision weighting, `Var(y_i) = sigma\^2 / w_i` (the ASReml / lme4 / glm(weights=) sense). |

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

| \#  | Vignette                                            |
|-----|-----------------------------------------------------|
| 01  | Getting started: what changes when you go Bayesian  |
| 02  | The formula surface, and what is not supported      |
| 03  | Regression and hierarchical models                  |
| 04  | Priors, and what they do to your answer             |
| 05  | Multi-environment trials and genomics               |
| 06  | Spatial and temporal structure                      |
| 07  | After the fit: summaries, comparison, triangulation |
| 08  | Big data: fitting without holding the data          |

Each vignette presents the equivalent the ASReml and flexyBayes fit,
what the posterior adds, and what it costs.

## Known limitations

- **Structured GxE beyond `diag` / `us`**: `corh(env):gen`
  (heterogeneous variances with one shared correlation) and
  `fa(env, k):gen` (factor-analytic) are current not implemented.
- **Spatial structure**: only the separable AR1 field
  (`random = ~ ar1(row):ar1(col)`) is currently supported, on INLA, over
  a complete grid with one observation per node. Intrinsic CAR and BYM2
  areal models are not implemented. You can express a custom spatial
  precision by passing your own matrix to `vm(g, precision = Q)`, but
  there is no BYM2 helper.
- **Smooth terms**: univariate penalised splines (`s(x)`, `spl(x)`) are
  supported on INLA. Multivariate and tensor-product smooths (`te()`,
  `ti()`, `t2()`) are currently not implemented.
- **Observation weights**: fit for the Gaussian family on an identity
  link, on both engines, in the ASReml / lme4 / `glm(weights=)`
  precision sense (`Var(y_i) = sigma^2 / w_i`) -Any other family, or a
  non-identity link on Gaussian is currently not possible. \<!–
- **Hidden-Markov, multi-state, and survival models**: not supported.
  Survival / time-to-event families are refused at the family gate. A
  NIMBLE backend covering these is on the roadmap with no fixed release
  target. –\>
- **Missing data**: flexyBayes does not impute covariates. A missing
  *predictor* is refused by default, as it is in ASReml
  (`na.method(x = "fail")`). Setting
  `na_action = list(y = "include", x = "omit")` drops the affected rows
  with a warning naming the count and the columns. A missing *response*
  is retained by default (`na_action = "augment"`), carried as a latent
  quantity the package marginalises, so the design index set a
  structured covariance is built over survives a lost plot. That
  preserves the representation, not information: under ignorability the
  parameter posterior is the same either way, and where missingness
  depends on the unobserved response both augmenting and omitting are
  biased. Non-Gaussian missing responses on brms are refused.

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
