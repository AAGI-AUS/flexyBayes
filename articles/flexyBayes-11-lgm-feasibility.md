# Tutorial 11: LGM feasibility: what \`lgm_gate()\` accepts and refuses

## 1. Why a feasibility filter?

INLA — integrated nested Laplace approximation (Rue, Martino, & Chopin,
2009) — is an extraordinarily fast Bayesian backend. The cost of that
speed is a structural assumption: the model must be a *latent Gaussian
model* (LGM). When the assumption holds, INLA returns posterior
marginals accurate to several decimal places. When it does not hold,
INLA still returns *something* — but that something is an approximation
whose error is hard to bound and easy to miss.

A *latent Gaussian model* has the form

``` math
\begin{aligned}
y_i \mid \boldsymbol{\eta}, \boldsymbol{\theta}
  &\sim \pi(y_i \mid \eta_i, \boldsymbol{\theta}) \\
\boldsymbol{\eta} = (\eta_1, \dots, \eta_n)
  &= \mathbf{A}\, \mathbf{x}, \\
\mathbf{x} \mid \boldsymbol{\theta}
  &\sim \mathcal{N}(\mathbf{0},\, \mathbf{Q}(\boldsymbol{\theta})^{-1})
\end{aligned}
```

where the latent field $`\mathbf{x}`$ is conditionally Gaussian given a
small ($`\le 15`$, in practice) hyperparameter vector
$`\boldsymbol{\theta}`$, and the linear predictor $`\boldsymbol{\eta}`$
is a linear combination of $`\mathbf{x}`$. The likelihood
$`\pi(y_i \mid \eta_i, \boldsymbol{\theta})`$ may be non-Gaussian
(Poisson, binomial, gamma, …) but $`y_i`$ depends on $`\mathbf{x}`$ only
through $`\eta_i`$.

Models that violate the LGM structure include:

- finite mixtures (multiple latent classes);
- distributional regression on a dispersion or shape parameter;
- non-linear predictors;
- random-effect priors that are not Gaussian (horseshoe, Cauchy,
  Student-$`t`$ on the *latent field* itself).

`lgm_gate()` is the structural-first filter that decides, *before*
calling INLA, whether your model is in the LGM class. It runs six checks
plus a seventh post-fit numerical check. Refusals are *non-silent by
construction*: every refusal carries the rule that fired, a one-line
gloss, a re-route suggestion (typically: use the greta backend), and a
documented override path.

## 2. The seven checks

| \# | Check | What it inspects |
|----|----|----|
| 1 | `family_allowlist` | the response family must appear in INLA’s likelihood roster |
| 2 | `predictor_linearity` | no `nl = TRUE` term; the linear predictor is linear in the latent field |
| 3 | `distributional_regression` | no auxiliary parameter (sigma, phi, nu, zi, hu) carries a non-intercept right-hand side |
| 4 | `re_gaussian_prior` | every random-effect prior has a Gaussian / multivariate normal / Gauss–Markov-random-field family |
| 5 | `latent_class` | no `mixture`, `hmm`, or `multistate` term; no `family = "mixture(...)"` |
| 6 | `hyperparam_budget` | the total number of hyperparameters is at most 15 (hard) — soft warning at 10 |
| 7 | `numerical_confirm` (post-fit) | INLA’s `mode$mode.status == 0`; `mlik` finite; no hyperparameters on the boundary |

Checks 1–6 fire structurally (without running INLA) on the `fb_terms`
intermediate representation (IR). Check 7 lives inside `emit_inla()` and
is non-overridable: a numerically failed fit is discarded, never quietly
returned.

## 3. The pass case

The simplest LGM-feasible model is a Gaussian random intercept on a
small dataset. We use `lme4::sleepstudy`.

``` r

library(flexyBayes)
data(sleepstudy, package = "lme4")
```

[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
builds the IR; `lgm_gate()` runs the six structural checks and either
returns the augmented IR (with the `lgm_compatible` capability) or a
structured refusal.

``` r

fb_t <- flexyBayes:::fb_from_brms(
  Reaction ~ Days + (1 | Subject),
  data = sleepstudy
)
res <- flexyBayes:::lgm_gate(fb_t)
class(res)
#> [1] "fb_terms" "list"
res$capabilities
#> [1] "lgm_compatible"      "gretaR_slot_dormant"
```

The IR carries `lgm_compatible`. Internally,
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
(the INLA engine pin, sugar for `flexybayes(..., backend = "inla")`)
uses this flag to decide whether to call `emit_inla()` (it does), or
refuse with a structured error (it does not).

``` r

fit <- fb_inla(Reaction ~ Days + (1 | Subject), data = sleepstudy,
               verbose = FALSE)
fit$num_check       # post-fit check 7
#> $pass
#> [1] TRUE
#> 
#> $reasons
#> character(0)
```

`fit$num_check` carries the verdict of check 7: `mode_status`,
`mlik_finite`, and `boundary_proximity`. All three pass for a
well-conditioned Gaussian random-intercept fit.

## 4. Refusal patterns

We now demonstrate four representative refusals. Each one resolves to a
structured `lgm_refusal` object whose print method names the rule, the
gloss, the re-route, and the override path.

### 4.1 Refusal pattern A — random slopes (caught by `fb_from_brms`, before the gate)

A brms-formula random slope is rejected at *parse* time, before the
fb_terms IR is ever built. The error message names the limitation
directly:

``` r

fb_brms(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#> Error:
#> ! Correlated random slopes (x | g) are not yet supported.
#> Uncorrelated random slopes (x || g) are supported -- they fit the
#> marginal slope and intercept variances independently and are equivalent
#> to the correlated form when the correlation is small.
#>
#> If your model needs the correlation parameter, defer to a
#> future release (structured-covariance representation).
#> 
#> Workaround: re-fit as (x || g) if the correlation is not of inferential
#> interest, or use backend = "greta" via fb_brms() with a hand-rolled
#> covariance prior.
#> 
#> Got: (Days | Subject)
```

The asreml-format entry `flexybayes(fixed = ..., random = ...)` accepts
random-slope-equivalent constructions and routes them through greta —
see the *hierarchical models* and *asreml-shaped formulas* vignettes.

### 4.2 Refusal pattern B — family outside the INLA roster (check 1)

INLA carries an internal roster of supported likelihoods (around 100
families when `INLA::inla.models()` is queried live). A family not in
that roster cannot be Laplace-approximated by INLA. We construct a
fictitious family name to make the refusal visible:

``` r

fb_t_bad <- fb_t
fb_t_bad$family <- "imaginary_invented_family"
ref <- flexyBayes:::lgm_gate(fb_t_bad)
class(ref)
#> [1] "lgm_refusal" "list"
ref
#> flexyBayes: INLA backend refused for this model.
#> Reasons (1 structural failure):
#>   - [family_allowlist] family "imaginary_invented_family" is not in the INLA likelihood allowlist (no built-in Laplace machinery).
#>     Diagnostic: fb$family = "imaginary_invented_family"
#> Re-route: try backend = "brms" via fb_brms() for the Stan passthrough
#>   (model must be in the brms corpus), or backend = "greta" for the broader
#>   flexyBayes formula path. NIMBLE is planned future work, not yet implemented.
#> For mixtures with <= 2 components, the planned INLA-within-MCMC escape
#>   (Gomez-Rubio & Rue 2018) is on the roadmap.
#> Override (not recommended): lgm_gate(fb, force = "inla",
#>   acknowledge_silent_bias_risk = TRUE, reason = "<your reason>")
#> Docs: vignette("flexyBayes-11-lgm-feasibility").
```

The structured refusal names the rule that fired (`family_allowlist`),
the diagnostic (`fb$family = "imaginary_invented_family"`), the
suggested re-route
([`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
for the Stan passthrough;
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
– equivalently `flexybayes(..., backend = "greta")` – for the broader
flexyBayes formula path), and the override path. Note that *the package
never silently approximates*: the user always sees why INLA refused.

### 4.3 Refusal pattern C — non-Gaussian RE prior (check 4)

The latent field of an LGM is Gaussian by construction. A horseshoe,
Cauchy, or Student-$`t`$ prior on the *random-effect field* (not on the
random-effect *standard deviation*, which is always a hyperparameter)
violates the assumption. The
[`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
DSL does not expose a target for the latent field directly, so this
refusal is constructed by setting an explicit non-Gaussian RE-prior
entry on the IR:

``` r

fb_t_horseshoe <- fb_t
fb_t_horseshoe$priors <- structure(
  list(re = list(list(family = "horseshoe", target = "u_Subject"))),
  class = "fb_prior"
)
ref_h <- flexyBayes:::lgm_gate(fb_t_horseshoe)
ref_h
#> flexyBayes: INLA backend refused for this model.
#> Reasons (1 structural failure):
#>   - [re_gaussian_prior] non-Gaussian random-effect prior detected (horseshoe). INLA's latent field is Gaussian by construction.
#>     Diagnostic: fb$priors$re entries with non-Gaussian family
#> Re-route: try backend = "brms" via fb_brms() for the Stan passthrough
#>   (model must be in the brms corpus), or backend = "greta" for the broader
#>   flexyBayes formula path. NIMBLE is planned future work, not yet implemented.
#> For mixtures with <= 2 components, the planned INLA-within-MCMC escape
#>   (Gomez-Rubio & Rue 2018) is on the roadmap.
#> Override (not recommended): lgm_gate(fb, force = "inla",
#>   acknowledge_silent_bias_risk = TRUE, reason = "<your reason>")
#> Docs: vignette("flexyBayes-11-lgm-feasibility").
```

Real horseshoe-on-latent-field shrinkage is best handled with the greta
backend; a future release will add a horseshoe target to the DSL on the
greta side.

### 4.4 Refusal pattern D — hyperparameter-budget exhaustion (check 6)

INLA integrates over the hyperparameter vector $`\boldsymbol{\theta}`$
on a numerical grid (or central-composite design). The grid scales
exponentially in $`\dim(\boldsymbol{\theta})`$. INLA’s authors recommend
$`\dim(\boldsymbol{\theta}) \le 15`$ for tractable Laplace
approximation. An unstructured (`us`) genotype-by-environment covariance
matrix across nine environments has $`9 \cdot 10 / 2 = 45`$ free entries
— well over the budget.

``` r

data(yan.winterwheat, package = "agridat")
fb_us <- flexyBayes:::fb_from_asreml(
  fixed = yield ~ env, random = ~ us(env):gen, rcov = NULL,
  data = yan.winterwheat
)
ref_us <- flexyBayes:::lgm_gate(fb_us)
ref_us
#> flexyBayes: INLA backend refused for this model.
#> Reasons (2 structural failures):
#>   - [hyperparam_budget] hyperparameter count = 46 exceeds hard limit 15 (CCD/grid integration intractable; INLA Laplace approximation deteriorates).
#>     Diagnostic: count = likelihood (1) + random + rcov contributions
#>   - [random_term_type_inla] random term type "us_gxe" (unstructured G×E) is outside the INLA emit allowlist. INLA's f() machinery does not currently represent this structured-covariance class without an SPDE / kronecker expansion that the current INLA emit does not produce. Re-route via backend = "greta".
#>     Diagnostic: fb$random_terms type(s): us_gxe
#> Re-route: try backend = "brms" via fb_brms() for the Stan passthrough
#>   (model must be in the brms corpus), or backend = "greta" for the broader
#>   flexyBayes formula path. NIMBLE is planned future work, not yet implemented.
#> For mixtures with <= 2 components, the planned INLA-within-MCMC escape
#>   (Gomez-Rubio & Rue 2018) is on the roadmap.
#> Override (not recommended): lgm_gate(fb, force = "inla",
#>   acknowledge_silent_bias_risk = TRUE, reason = "<your reason>")
#> Docs: vignette("flexyBayes-11-lgm-feasibility").
```

The refusal cites `hyperparam_budget` and reports the count. The
re-route is unambiguous: use the greta backend (which handles
unstructured covariances natively via Cholesky parameterisation), or
reduce to a factor-analytic approximation (`fa(env, k):gen`) with small
`k` — the latter is the standard agricultural-statistics move when
$`n_{\text{env}}`$ exceeds ten.

## 5. The override path

Sometimes a user wants to force INLA dispatch despite a refusal —
typically because they have done their own homework and know the
silent-bias risk. The override is *two-key armed*: it requires
`force = "inla"`, `acknowledge_silent_bias_risk = TRUE`, and a non-empty
`reason`:

``` r

ovr <- flexyBayes:::lgm_gate(
  fb_t_bad,
  force                        = "inla",
  acknowledge_silent_bias_risk = TRUE,
  reason                       = "exploratory: comparing with greta-fit posterior"
)
class(ovr)            # fb_terms (overrideed)
#> [1] "fb_terms" "list"
grep("lgm_force", ovr$capabilities, value = TRUE)
#> [1] "lgm_force_overridden"                                            
#> [2] "lgm_force_reason:exploratory: comparing with greta-fit posterior"
#> [3] "lgm_force_bypassed:family_allowlist"
```

The override leaves an audit trail in `$capabilities`: the literal
`reason` string and the rule(s) that were bypassed. Downstream tooling
(the
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
print method, in particular) flags forced-INLA fits so the override is
visible in any cross-engine comparison.

**The override does not apply to check 7.** A numerical failure inside
INLA (mode finder did not converge, marginal likelihood is not finite,
hyperparameter on the boundary) is non-overridable. A forced-INLA fit
that survives check 7 is *numerically* OK; a forced- INLA fit that fails
check 7 is rejected even with the override active.

## 6. Capability flags beyond the gate verdict

`fb_terms$capabilities` is a character vector that flows through the
dispatch infrastructure as a structured trail of evidence. Beyond the
gate-verdict flags (`lgm_compatible`, `lgm_force_overridden`, the
soft-warning markers, and the override audit), every LGM-compatible IR
also carries `gretaR_slot_dormant` — an informational capability flag
advertising that the package has a provisioned slot for the gretaR
backend but the slot has not yet been activated.

``` r

gs <- gretaR_status()
gs$dormancy_reason
#> [1] "slot_provisioned_not_activated"
cat(gs$activation_procedure, sep = "\n")
#> Install the gretaR package when it goes public on CRAN.
#> Verify the local gretaR build against the audit checklist.
#> Run options(flexyBayes.gretaR_activated = TRUE) in the session.
```

The slot is dormant because the gretaR package is not yet publicly
available and the package’s audit process for new backends has not yet
been wired into the registry. Activation will land in a subsequent
release when both conditions are met. Until then,
[`gretaR_status()`](https://aagi-aus.github.io/flexyBayes/reference/gretaR_status.md)
is the canonical discoverability surface for users who see the flag in
`fb_terms$capabilities` or in the print method of a fitted model.

The user-facing `backend` argument on
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
/
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
does *not* advertise `"gretaR"` while the slot is dormant; the selection
set `c("auto", "greta", "inla", "brms")` keeps the surface honest. Once
activated, the slot will widen the set to include `"gretaR"`, swap the
canonical-name mapper stub for the production mapper, and replace the
dispatch-helper branch with the real emit step.

## 7. Common pitfalls

**Mistaking the prior on $`\sigma`$ for the prior on the latent field.**
A Student-$`t`$ prior on `sigma` (the random-effect standard deviation)
is fine and does not trigger refusal — `sigma` is a hyperparameter, not
part of the latent field. A Student-$`t`$ prior on the random-effect
*coefficient vector* (not exposed in the DSL, constructed manually here)
does trigger refusal, correctly.

**Reading “lgm_compatible” as “the model is correct”.** The flag means
*only* that INLA’s structural assumptions are not violated; it does not
mean the model fits the data, that the priors are appropriate, or that
the posterior is close to the truth. The *priors and regularisation* and
*cross-engine triangulation* vignettes cover these orthogonal concerns.

**Treating a refusal as a bug.** Refusals are how the package prevents
silent approximation. The right response is one of: use the greta
backend (`backend = "greta"`), reduce model complexity (e.g.,
factor-analytic instead of unstructured), or — if you understand the
risk — invoke the override.

## 8. Active prompts

1.  Construct an `fb_terms` object with `family = "tweedie"` (a real
    distribution, supported by INLA) and run `lgm_gate()`. Does it pass?
    Why?
2.  Construct an `fb_terms` for `fa(env, 2):gen` on `yan.winterwheat`.
    Does it pass the budget check? What is the hyperparameter count?
3.  Run a forced-INLA fit on a refused model and compare its posterior
    with the greta-backend fit via
    [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md).
    The override is intended to surface bias; how visible is it?

## 9. Session information

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.5.1
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
#> 
#> time zone: Australia/Adelaide
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] bayesplot_1.15.0 emmeans_2.0.2    ggplot2_4.0.3    flexyBayes_0.8.3
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6           DBI_1.3.0              mnormt_2.1.2          
#>   [4] tfautograph_0.3.2      sandwich_3.1-1         rlang_1.2.0           
#>   [7] magrittr_2.0.5         multcomp_1.4-29        otel_0.2.0            
#>  [10] ggridges_0.5.7         matrixStats_1.5.0      e1071_1.7-17          
#>  [13] compiler_4.5.2         mgcv_1.9-3             loo_2.9.0             
#>  [16] png_0.1-9              callr_3.7.6            vctrs_0.7.3           
#>  [19] reshape2_1.4.5         stringr_1.6.0          pkgconfig_2.0.3       
#>  [22] crayon_1.5.3           backports_1.5.1        labeling_0.4.3        
#>  [25] effectsize_1.0.2       nloptr_2.2.1           MatrixModels_0.5-4    
#>  [28] torch_0.17.0           bit_4.6.0              INLA_25.10.19         
#>  [31] xfun_0.57              jsonlite_2.0.0         progress_1.2.3        
#>  [34] parallel_4.5.2         prettyunits_1.2.0      tensorflow_2.20.0     
#>  [37] R6_2.6.1               stringi_1.8.7          RColorBrewer_1.1-3    
#>  [40] reticulate_1.45.0      boot_1.3-32            parallelly_1.47.0     
#>  [43] numDeriv_2016.8-1.1    estimability_1.5.1     Rcpp_1.1.1-1.1        
#>  [46] knitr_1.51             zoo_1.8-15             parameters_0.28.3     
#>  [49] base64enc_0.1-6        Matrix_1.7-4           splines_4.5.2         
#>  [52] tidyselect_1.2.1       dichromat_2.0-0.1      abind_1.4-8           
#>  [55] agridat_1.26           codetools_0.2-20       processx_3.9.0        
#>  [58] listenv_0.10.1         gretaR_0.2.0           lattice_0.22-7        
#>  [61] tibble_3.3.1           plyr_1.8.9             bayestestR_0.17.0     
#>  [64] withr_3.0.2            bridgesampling_1.2-1   S7_0.2.2              
#>  [67] posterior_1.7.0        coda_0.19-4.1          evaluate_1.0.5        
#>  [70] marginaleffects_0.32.0 future_1.70.0          survival_3.8-3        
#>  [73] sf_1.1-0               units_1.0-1            proxy_0.4-29          
#>  [76] RcppParallel_5.1.11-2  pillar_1.11.1          tensorA_0.36.2.1      
#>  [79] whisker_0.4.1          KernSmooth_2.23-26     checkmate_2.3.4       
#>  [82] stats4_4.5.2           insight_1.4.6          reformulas_0.4.4      
#>  [85] sn_2.1.3               distributional_0.7.0   generics_0.1.4        
#>  [88] hms_1.1.4              rstantools_2.6.0       scales_1.4.0          
#>  [91] minqa_1.2.8            coro_1.1.0             globals_0.19.1        
#>  [94] xtable_1.8-8           class_7.3-23           glue_1.8.1            
#>  [97] tools_4.5.2            data.table_1.18.2.1    lme4_2.0-1            
#> [100] mvtnorm_1.3-6          grid_4.5.2             rbibutils_2.4.1       
#> [103] datawizard_1.3.0       nlme_3.1-168           cli_3.6.6             
#> [106] tfruns_1.5.4           viridisLite_0.4.3      fmesher_0.7.0         
#> [109] Brobdingnag_1.2-9      dplyr_1.2.1            gtable_0.3.6          
#> [112] greta_0.5.1            digest_0.6.39          classInt_0.4-11       
#> [115] TH.data_1.1-5          brms_2.23.0            farver_2.1.2          
#> [118] lifecycle_1.0.5        bit64_4.8.0            MASS_7.3-65
```

## References

Gómez-Rubio, V., & Rue, H. (2018). Markov chain Monte Carlo with the
integrated nested Laplace approximation. *Statistics and Computing*,
28(5), 1033–1051.

Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
inference for latent Gaussian models by using integrated nested Laplace
approximations. *Journal of the Royal Statistical Society: Series B*,
71(2), 319–392.

Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
(2017). Penalising model component complexity: A principled, practical
approach to constructing priors. *Statistical Science*, 32(1), 1–28.
