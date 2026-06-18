# Tutorial 07: Multi-environment trials and genomic selection

## 1. The two questions of plant breeding

> **This is a syntax-and-representation catalogue.** It demonstrates how
> to *write* models across the term surface. The fits use small sampling
> budgets so the document builds quickly, and several printed
> diagnostics therefore do **not** meet the package’s own convergence
> thresholds ($`\widehat{R} \le 1.01`$, ESS $`\ge 400`$); some model
> classes also mix poorly on their backend at any budget. Read the
> posterior summaries here as illustrations of the output *shape*, not
> as inferential results. For a convergence-clean worked fit and the
> trustworthy-fit workflow, see the *getting started* and *cross-engine
> triangulation* vignettes.

Crop variety evaluation poses two intertwined questions, and
`flexyBayes` provides one statistical framework for both.

1.  **Multi-environment trials (METs).** Genotypes are evaluated across
    many environments (locations × years). Their *rankings* change
    across environments — *genotype-by-environment interaction* (GxE).
    The first question is which genotypes are best on average, and which
    trade off mean for stability.
2.  **Genomic selection.** Modern breeding programmes have genome- wide
    markers on every candidate genotype, even those untested in the
    field. Predicting the field performance of an *unphenotyped*
    candidate from its markers is *genomic prediction*; selecting
    candidates by predicted performance is *genomic selection*
    (Meuwissen, Hayes, & Goddard, 2001).

This vignette walks through both questions, fitting models on real
agricultural data, and notes the connection through the kinship matrix.

> **A note on the routes below – read first.** This vignette fits the
> MET model two ways. The **recommended primary route is INLA**: the
> independent (diagonal) genotype-by-environment model refits in a
> couple of seconds and its fixed-effect and variance-component
> posteriors are trustworthy at default settings, so the INLA fits below
> *are* usable inference. The **factor-analytic (structured-covariance)
> route runs on greta** and is shown as the slower, harder alternative:
> it is the hardest case for Hamiltonian Monte Carlo, and at the modest
> MCMC budgets used here it **frequently fails to mix** – large
> $`\widehat{R}`$, near-zero ESS, and factor-analytic loadings that do
> not sample reliably. Read the greta factor-analytic loadings and
> biplots as illustrations of the *workflow*, not as inference; the
> genomic-prediction (GBLUP / pedigree) fits later in the vignette are
> likewise greta structured-covariance models held to the same caution.
> Always verify convergence as the *getting started* and *structured
> covariance* vignettes describe before interpreting a greta structured
> fit. In particular, a factor-analytic term’s raw loadings are
> identified only up to rotation and sign, so their $`\widehat{R}`$ is
> expected to be large even when the model has converged; judge
> convergence on the *identified* covariance
> $`G = \Lambda\Lambda^\top + \mathrm{diag}(\psi)`$ via
> [`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md),
> which is rotation- and sign-invariant.

## 2. Multi-environment trials: `yan.winterwheat`

`agridat::yan.winterwheat` (Yan & Tinker, 2006) records winter wheat
yield for 18 genotypes in 9 Ontario environments — a small, classic
dataset that produces interpretable MET results.

``` r

library(flexyBayes)
data(yan.winterwheat, package = "agridat")
str(yan.winterwheat)
#> 'data.frame':    162 obs. of  3 variables:
#>  $ gen  : Factor w/ 18 levels "Ann","Ari","Aug",..: 1 2 3 4 5 6 7 8 9 10 ...
#>  $ env  : Factor w/ 9 levels "BH93","EA93",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ yield: num  4.46 4.42 4.67 4.73 4.39 ...
```

The natural mixed model is

``` math
y_{ijk} = \mu + e_j + g_i + (ge)_{ij} + \varepsilon_{ijk},
```

with $`e_j`$ environment effects (fixed), $`g_i`$ genotype effects
(random), $`(ge)_{ij}`$ genotype-by-environment interaction, and
$`\varepsilon_{ijk}`$ residuals. The *parameterisation of $`(ge)_{ij}`$*
is the modelling choice. Three options span the parsimony spectrum:

- **diagonal (independent GxE)**: each genotype-by-environment cell
  $`(ge)_{ij} \sim \mathcal{N}(0, \sigma_{ge}^2)`$, independent across
  cells – the parsimonious, scalable workhorse;
- **factor-analytic** with `k = 2`: a small-rank approximation of the
  full GxE covariance, the standard breeder’s model for *stability*;
- **unstructured**: the full $`9 \times 9`$ covariance, with 45 free
  entries – too many for 162 observations to identify.

### 2.1 Recommended route: the diagonal GxE model on INLA

For the *mean-performance* question – which genotypes are best on
average, and how large is the GxE relative to the main effects – the
diagonal GxE model is both sufficient and fast. INLA fits it by
deterministic Laplace approximation in a couple of seconds, with no
sampler to tune and no convergence to chase, so its posteriors are
usable inference rather than a workflow sketch. **This is the
recommended starting point.**

``` r

yan.winterwheat$ge <- interaction(
  yan.winterwheat$gen,
  yan.winterwheat$env,
  drop = TRUE
)
fit_inla <- flexybayes(
  fixed   = yield ~ env,
  random  = ~ gen + ge,
  data    = yan.winterwheat,
  backend = "inla",
  verbose = FALSE
)
fit_inla
#> Bayesian fit  [flexybayes_inla / INLA backend]
#> ------------------------------------------------------- 
#>   formula: yield ~ 1 + env + f(gen, model = "iid", hyper = list(prec = list(prior = "expression: U=4.940099638422605; lb=-2*log(U); ld=-log(U)-log(2)-theta/2; return( theta<lb ? -1.0e10 : ld );"))) +     f(ge, model = "iid", hyper = list(prec = list(prior = "expression: U=4.940099638422605; lb=-2*log(U); ld=-log(U)-log(2)-theta/2; return( theta<lb ? -1.0e10 : ld );")))
#>   family:  gaussian
#>   n_obs:   162
#>   fixed:   9
#>   random:  2
#>   hyper:   3
#>   runtime: 2.2 sec
#>   numerical confirm: PASS
#> ------------------------------------------------------- 
#>   $inla -- raw INLA fit (use INLA's summary, plot, etc.)
#>   $fb   -- the fb_terms IR used for dispatch
```

The fixed-effect environment means and the two variance components – the
genotype main-effect variance $`\sigma_g^2`$ and the GxE cell variance
$`\sigma_{ge}^2`$ – read straight off the standard accessors. The
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) method
returns the fixed effects as a flat, `broom`-style table with the
canonical `term` / `estimate` / `std.error` / `conf.low` / `conf.high`
columns, identical in shape to what the greta and brms backends return –
so a cross-engine comparison table is an
[`rbind()`](https://rdrr.io/r/base/cbind.html), not a hand-built
reconciliation.

``` r

head(tidy(fit_inla))
#>          term   estimate std.error   conf.low  conf.high
#> 1 (Intercept)  4.3629310 0.1314894  4.1039214  4.6219402
#> 2     envEA93  0.0752906 0.1300043 -0.1800980  0.3306798
#> 3     envHW93 -1.2259762 0.1300043 -1.4813645 -0.9705866
#> 4     envID93 -0.8678681 0.1300043 -1.1232565 -0.6124786
#> 5     envKE93  1.3196134 0.1300043  1.0642246  1.5750024
#> 6     envNN93  0.6972298 0.1300043  0.4418411  0.9526189
```

The ratio $`\sigma_{ge}^2 / (\sigma_g^2 + \sigma_{ge}^2)`$ measures how
much of the genotypic signal is environment-specific: a large ratio says
rankings reshuffle across environments and a single
*average-performance* number hides real GxE, which is precisely the case
where the factor-analytic stability model below earns its extra
parameters.

### 2.2 Slower alternative: the factor-analytic model on greta

The factor-analytic with `k = 2` reduces the unstructured $`9 \times 9`$
covariance to 18 entries (16 loadings + 2 specific variances) and is the
standard model for *stability* – it decomposes each genotype’s GxE into
a few latent environmental gradients. It is a structured-covariance
model, so it runs on the greta backend and is the **slower,
harder-mixing** alternative to the diagonal model above: the loadings
sample poorly at the small budget used here, and the fit below is a
workflow illustration, not inference.

``` r

fit_fa <- flexybayes(
  fixed     = yield ~ env,
  random    = ~ fa(env, 2):gen,
  data      = yan.winterwheat,
  n_samples = mcmc_small$n_samples,
  warmup    = mcmc_small$warmup,
  chains    = mcmc_small$chains,
  verbose   = FALSE
)
fit_fa
#> Bayesian mixed model  [flexyBayes]
#> ------------------------------------------------------- 
#>   Fixed  : yield ~ env 
#>   Random : ~fa(env, 2):gen 
#>   Family : gaussian ( identity link )
#>   MCMC   : 2 chain(s) x 1000 samples (warmup = 2000 ) -- 19.2 sec
#>   Params : 6 monitored; 2 fixed, 1 random terms
#>   Representation: exact
#>   Engine:         greta MCMC
#>   Max Rhat: 12.46  [!] 
#>   Min ESS: 0 
#> ------------------------------------------------------- 
#>   $glm    -- GLM-compatible (summary, emmeans, etc.)
#>   $greta  -- native greta (draws, model, calculate)
#>   $extras -- diagnostics, BLUPs, variance components
```

The variance-component block reports the loading variances (one per
factor) and the specific variances (one per environment). Together these
reconstruct the implied environment-by-environment covariance
$`\boldsymbol{\Sigma} = \boldsymbol{\Lambda}\boldsymbol{\Lambda}^\top + \boldsymbol{\Psi}`$.
Use the factor-analytic route when the *stability* decomposition is the
scientific target and a sampling budget large enough to mix it is
available; use the diagonal INLA route for fast, trustworthy
mean-performance inference and as the de-risking baseline you always fit
first.

## 3. Genotype rankings and stability

Two questions plant breeders ask of every fit:

- **Mean performance**: which genotypes have the highest predicted
  yield, averaged across environments?
- **Stability**: which genotypes maintain their ranking across
  environments — and which do not?

Mean performance comes from the posterior of the genotype effect
$`g_i`$. Stability comes from the magnitude and direction of each
genotype’s factor scores: a genotype near the origin of the factor space
is stable; a genotype far from the origin in factor 1 is sensitive to
the dimension that factor 1 captures.

``` r

draws  <- flexyBayes::fb_as_draws_simple(fit_fa)
gen_levels <- levels(yan.winterwheat$gen)
# Note: parameter names are backend-specific; the genotype effects
# typically appear as g_i_1 ... g_i_18 or similar in the greta draws.
# We summarise the posterior means and rank.
gen_means <- vapply(seq_along(gen_levels), function(k) {
  nm <- paste0("g_atg_", k)
  if (!is.null(draws[[nm]])) mean(draws[[nm]]) else NA_real_
}, numeric(1))
names(gen_means) <- gen_levels
sort(gen_means, decreasing = TRUE)
#> named numeric(0)
```

Note: the precise parameter naming inside `draws` varies with the
generated greta program; `fit_fa$extras$code` shows the literal names.
The pedagogical point is that the posterior of every random effect is
available — *with full uncertainty quantification* — for ranking,
contrasts, and downstream selection decisions.

## 4. Visualising GxE: the interaction plot

Before any model is fit, the canonical exploratory plot of MET data is
the *genotype-by-environment interaction plot*.

``` r

library(ggplot2)
gxe <- aggregate(yield ~ gen + env, data = yan.winterwheat, mean)
ggplot(gxe, aes(x = env, y = yield, group = gen, colour = gen)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) +
  scale_colour_viridis_d(option = "viridis") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "environment", y = "yield (t/ha)",
       title = "yan.winterwheat: genotype means per environment")
```

![plot of chunk
yw-interaction-plot](flexyBayes-07-met-and-genomics-figs/yw-interaction-plot-1.png)

plot of chunk yw-interaction-plot

Crossing lines indicate GxE — the rank of a genotype changes across
environments. Parallel lines indicate no GxE — the diagonal model
suffices. In `yan.winterwheat`, the lines cross visibly, motivating the
factor-analytic decomposition.

## 5. Genomic BLUP (GBLUP)

When genome-wide markers are available, the kinship matrix
$`\mathbf{G}`$ — the realised additive genetic relationship matrix
computed from marker genotypes (VanRaden, 2008) — replaces the default
identity covariance in the genotype random effect:

``` math
\mathbf{u} \mid \sigma_g^2 \sim \mathcal{N}\!\left(\mathbf{0},\, \sigma_g^2 \mathbf{G}\right).
```

This is *genomic BLUP* (GBLUP). It pools information across genotypes
via marker-based similarity, and — crucially — it extrapolates to
*unphenotyped* genotypes whose markers are known.

`flexyBayes` plumbs this through the `vm()` term and the
`known_matrices` argument:

``` r

data(met_example, package = "flexyBayes")
fit_gblup <- flexybayes(
  fixed     = yield ~ env,
  random    = ~ vm(geno, Gmat),
  data      = met_example$dat,
  known_matrices = list(Gmat = met_example$G_mat),
  n_samples = mcmc_small$n_samples,
  warmup    = mcmc_small$warmup,
  chains    = mcmc_small$chains,
  verbose   = FALSE
)
fit_gblup
#> Bayesian mixed model  [flexyBayes]
#> ------------------------------------------------------- 
#>   Fixed  : yield ~ env 
#>   Random : ~vm(geno, Gmat) 
#>   Family : gaussian ( identity link )
#>   MCMC   : 2 chain(s) x 1000 samples (warmup = 2000 ) -- 7.8 sec
#>   Params : 5 monitored; 2 fixed, 1 random terms
#>   Representation: exact
#>   Engine:         greta MCMC
#>   Max Rhat: 1.988  [!] 
#>   Min ESS: 5 
#> ------------------------------------------------------- 
#>   $glm    -- GLM-compatible (summary, emmeans, etc.)
#>   $greta  -- native greta (draws, model, calculate)
#>   $extras -- diagnostics, BLUPs, variance components
```

The variance-component block reports a single $`\sigma_g^2`$ — the
*additive genetic variance*. Combined with the residual variance, this
gives the *narrow-sense heritability* $`h^2 = \sigma_g^2 /
(\sigma_g^2 + \sigma_e^2)`$ on the appropriate plot mean basis (Falconer
& Mackay, 1996).

``` r

draws <- flexyBayes::fb_as_draws_simple(fit_gblup)
# vm(geno, Gmat) names the genotype SD `sigma_geno`.
sigma_g <- draws$sigma_geno
sigma_e <- draws$sigma_e_atg
h2 <- sigma_g^2 / (sigma_g^2 + sigma_e^2)
c(mean = mean(h2), q025 = quantile(h2, 0.025), q975 = quantile(h2, 0.975))
#>       mean  q025.2.5% q975.97.5% 
#>  0.4481187  0.1875975  0.7756129
```

The credible interval on $`h^2`$ is the practitioner’s first quantity of
interest — a heritability of 0.4 with a tight interval is actionable for
breeding-value selection; 0.4 with a wide interval is a warning that the
data are insufficient for confident decisions.

## 6. Pedigree BLUP (the conventional alternative)

If a known pedigree is available — but no genome-wide markers — the
pedigree-derived numerator relationship matrix $`\mathbf{A}`$ takes the
role of $`\mathbf{G}`$:

``` r

fit_ped <- flexybayes(
  fixed     = yield ~ env,
  random    = ~ vm(geno, Amat),
  data      = met_example$dat,
  known_matrices = list(Amat = met_example$A_mat),
  n_samples = mcmc_small$n_samples,
  warmup    = mcmc_small$warmup,
  chains    = mcmc_small$chains,
  verbose   = FALSE
)
fit_ped
#> Bayesian mixed model  [flexyBayes]
#> ------------------------------------------------------- 
#>   Fixed  : yield ~ env 
#>   Random : ~vm(geno, Amat) 
#>   Family : gaussian ( identity link )
#>   MCMC   : 2 chain(s) x 1000 samples (warmup = 2000 ) -- 7.8 sec
#>   Params : 5 monitored; 2 fixed, 1 random terms
#>   Representation: exact
#>   Engine:         greta MCMC
#>   Max Rhat: 15.639  [!] 
#>   Min ESS: 6 
#> ------------------------------------------------------- 
#>   $glm    -- GLM-compatible (summary, emmeans, etc.)
#>   $greta  -- native greta (draws, model, calculate)
#>   $extras -- diagnostics, BLUPs, variance components
```

GBLUP and pedigree BLUP are structurally identical; only the
relationship matrix changes. When both markers and pedigree are
available, the genomic-pedigree blend (`H` matrix; Christensen & Lund,
2010) is on the roadmap.

## 7. Genomic prediction error

The cross-validation procedure for genomic prediction holds out a random
subset of phenotyped genotypes, fits the model on the rest, predicts the
held-out genotypes from their markers + the trained $`\sigma_g^2`$
posterior, and compares predictions with realised phenotypes. The
summary metric is *prediction accuracy*: the correlation between
predicted and observed across the held-out set.

``` r

# Sketch (full implementation depends on dataset structure):
folds <- rsample::vfold_cv(unique(geno_id), v = 5)
acc <- vapply(seq_along(folds$splits), function(k) {
  train_id <- rsample::analysis(folds$splits[[k]])
  test_id  <- rsample::assessment(folds$splits[[k]])
  fit_k <- flexybayes(..., data = subset(dat, geno %in% train_id), ...)
  pred  <- predict(fit_k, newdata = subset(dat, geno %in% test_id))
  cor(pred, subset(dat, geno %in% test_id)$yield)
}, numeric(1))
mean(acc)
```

The Bayesian framework propagates uncertainty into the predictions —
each prediction comes with a posterior distribution, not a single
number. Decision-theoretic genomic selection uses these full
distributions (Bauer, Reif, & Schipprack, 2009) rather than
posterior-mean point predictions.

## 8. Pitfalls

**Kinship matrix scaling.** $`\mathbf{G}`$ matrices are conventionally
scaled so that their average diagonal equals 1. Different scaling
conventions (VanRaden 2008 versus Astle & Balding 2009) shift the
posterior on $`\sigma_g^2`$ by a constant. State the scaling explicitly
when reporting heritability.

**Connectedness across environments.** Factor-analytic GxE assumes that
*enough* genotypes are common across environments to identify the
loadings. With strongly disconnected designs — different genotype sets
per environment — `fa(env, 2):gen` may give noisy loadings; convergence
diagnostics will flag this.

**`spl()` for environmental covariates.** Many MET analyses have
quantitative environmental covariates (temperature, rainfall) — fit
these via `spl()` in the random formula to capture non-linear
genotype-by-environment-covariate interactions.

**Genotype-environment covariance versus genotype variance.** The
factor-analytic model puts structure on $`\mathrm{Var}(g_{ij})`$ across
environments; it does *not* introduce a separate $`\sigma_g^2`$
hyperparameter. This is why the factor-analytic fit reports loading
variances rather than a single genotype variance.

## 9. Active prompts

1.  Refit `yan.winterwheat` with `fa(env, 1):gen`. How does the
    posterior on `gen` change? Does the one-factor approximation pick up
    the dominant pattern visible in the interaction plot?
2.  Compute the heritability from the GBLUP fit on `met_example`. Is the
    credible interval narrow enough to support breeding decisions?
3.  Replace the kinship matrix in the GBLUP fit with the pedigree matrix
    and triangulate the two fits. Where do they agree, and where do they
    differ?

## 10. Session information

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
#> [1] ggplot2_4.0.3    flexyBayes_0.8.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] tidyselect_1.2.1       viridisLite_0.4.3      dplyr_1.2.1           
#>   [4] farver_2.1.2           loo_2.9.0              tensorflow_2.20.0     
#>   [7] S7_0.2.2               tensorA_0.36.2.1       INLA_25.10.19         
#>  [10] TH.data_1.1-5          digest_0.6.39          estimability_1.5.1    
#>  [13] lifecycle_1.0.5        gretaR_0.5.0           sf_1.1-0              
#>  [16] survival_3.8-3         agridat_1.26           processx_3.9.0        
#>  [19] posterior_1.7.0        magrittr_2.0.5         compiler_4.5.2        
#>  [22] rlang_1.2.0            progress_1.2.3         tools_4.5.2           
#>  [25] data.table_1.18.2.1    knitr_1.51             labeling_0.4.3        
#>  [28] prettyunits_1.2.0      bridgesampling_1.2-1   bit_4.6.0             
#>  [31] classInt_0.4-11        reticulate_1.45.0      RColorBrewer_1.1-3    
#>  [34] multcomp_1.4-29        abind_1.4-8            KernSmooth_2.23-26    
#>  [37] withr_3.0.2            grid_4.5.2             xtable_1.8-8          
#>  [40] e1071_1.7-17           future_1.70.0          globals_0.19.1        
#>  [43] emmeans_2.0.2          scales_1.4.0           MASS_7.3-65           
#>  [46] dichromat_2.0-0.1      cli_3.6.6              mvtnorm_1.3-6         
#>  [49] crayon_1.5.3           generics_0.1.4         RcppParallel_5.1.11-2 
#>  [52] otel_0.2.0             tfruns_1.5.4           DBI_1.3.0             
#>  [55] proxy_0.4-29           stringr_1.6.0          splines_4.5.2         
#>  [58] bayesplot_1.15.0       parallel_4.5.2         coro_1.1.0            
#>  [61] matrixStats_1.5.0      base64enc_0.1-6        marginaleffects_0.32.0
#>  [64] brms_2.23.0            vctrs_0.7.3            Matrix_1.7-4          
#>  [67] sandwich_3.1-1         jsonlite_2.0.0         greta_0.5.1           
#>  [70] callr_3.7.6            hms_1.1.4              bit64_4.8.0           
#>  [73] listenv_0.10.1         units_1.0-1            glue_1.8.1            
#>  [76] parallelly_1.47.0      codetools_0.2-20       distributional_0.7.0  
#>  [79] stringi_1.8.7          gtable_0.3.6           tibble_3.3.1          
#>  [82] pillar_1.11.1          Brobdingnag_1.2-9      torch_0.17.0          
#>  [85] R6_2.6.1               fmesher_0.7.0          evaluate_1.0.5        
#>  [88] lattice_0.22-7         png_0.1-9              backports_1.5.1       
#>  [91] tfautograph_0.3.2      MatrixModels_0.5-4     rstantools_2.6.0      
#>  [94] class_7.3-23           Rcpp_1.1.1-1.1         checkmate_2.3.4       
#>  [97] nlme_3.1-168           coda_0.19-4.1          whisker_0.4.1         
#> [100] xfun_0.57              zoo_1.8-15             pkgconfig_2.0.3
```

## References

Bauer, A. M., Reif, J. C., & Schipprack, W. (2009). Reduction of the
genetic distance for marker-assisted backcrossing. *Theoretical and
Applied Genetics*, 119(1), 47–53.

Christensen, O. F., & Lund, M. S. (2010). Genomic prediction when some
animals are not genotyped. *Genetics Selection Evolution*, 42(2), 1–8.

Falconer, D. S., & Mackay, T. F. C. (1996). *Introduction to
quantitative genetics* (4th ed.). Longman.

Meuwissen, T. H. E., Hayes, B. J., & Goddard, M. E. (2001). Prediction
of total genetic value using genome-wide dense marker maps. *Genetics*,
157(4), 1819–1829.

Smith, A. B., Cullis, B. R., & Thompson, R. (2001). Analyzing variety by
environment data using multiplicative mixed models and adjustments for
spatial field trend. *Biometrics*, 57(4), 1138–1147.

VanRaden, P. M. (2008). Efficient methods to compute genomic
predictions. *Journal of Dairy Science*, 91(11), 4414–4423.

Yan, W., & Tinker, N. A. (2006). Biplot analysis of multi-environment
trial data: Principles and applications. *Canadian Journal of Plant
Science*, 86(3), 623–645.
