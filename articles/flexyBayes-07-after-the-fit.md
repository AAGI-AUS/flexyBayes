# 7. After the fit: summaries, comparison, triangulation

This page in one panel

- The problem : The fit is done. Getting the numbers out, comparing two
  models, and checking that the answer is a property of the data rather
  than of the engine.

- In ASReml :

  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html), `wald()`, and a
  change in REML log-likelihood.

- In flexyBayes :

  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html), `loo()`, and
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md).

- What the posterior adds : Any function of the parameters is computed
  on the draws. And with two engines, a check no single fit can perform
  on itself.

- What it costs : Model comparison is predictive rather than a test, and
  it needs a sampled fit.

## 1. Getting the numbers out

``` r

library(flexyBayes)
yw  <- agridat::yan.winterwheat
fit <- flexybayes(yield ~ env, random = ~ gen, data = yw,
                  backend = "inla", verbose = FALSE)
```

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) follows broom
convention: one row per term, columns that do not change name between
engines.

``` r

head(tidy(fit), 3)
#>          term    estimate std.error   conf.low  conf.high
#> 1 (Intercept)  4.36293122 0.1305793  4.1059742  4.6198880
#> 2     envEA93  0.07529038 0.1288959 -0.1777923  0.3283735
#> 3     envHW93 -1.22597659 0.1288959 -1.4790591 -0.9728932
tidy(fit, effects = "random")
#>     term  estimate  std.error  conf.low conf.high
#> 1  sigma 0.3861117 0.02352645 0.3426347 0.4348663
#> 2 sd_gen 0.3910171 0.07876814 0.2625308 0.5705438
```

[`ranef()`](https://aagi-aus.github.io/flexyBayes/reference/ranef.md)
gives the random-effect predictions – REML’s BLUPs, with intervals.

``` r

head(ranef(fit)$gen, 3)
#>   group level    estimate std.error   conf.low conf.high
#> 1   gen   Ann -0.16656260 0.1509254 -0.4646959 0.1289700
#> 2   gen   Ari  0.10539930 0.1508288 -0.1904270 0.4028784
#> 3   gen   Aug  0.04727761 0.1507772 -0.2489074 0.3442058
```

[`glance()`](https://generics.r-lib.org/reference/glance.html) gives one
row of fit-level information. On an INLA fit the sampler-specific
columns are `NA`, because a deterministic Laplace approximation has no
chains to report:

``` r

glance(fit)
#>   nobs npar logLik   family     link chains samples max_rhat min_ess run_time
#> 1  162    9     NA gaussian identity     NA      NA       NA      NA        2
```

## 2. Anything that is a function of the parameters

This is where a posterior is simply more convenient than a point
estimate.
[`fb_as_draws_simple()`](https://aagi-aus.github.io/flexyBayes/reference/fb_as_draws_simple.md)
returns the draws as a named list, and any derived quantity is computed
on them directly.

``` r

d <- fb_as_draws_simple(fit)
names(d)[1:4]
#> [1] "Predictor:1" "Predictor:2" "Predictor:3" "Predictor:4"
```

A ratio, a difference, a probability – each is one line, and each
arrives with its own uncertainty rather than needing a delta-method
approximation.

## 3. Comparing two models

Model comparison under a posterior is predictive: which model would
better predict a new observation. `loo::loo()` estimates that by
approximate leave-one-out cross-validation. It needs pointwise
likelihoods, so it needs a sampled fit. flexyBayes registers a method on
`loo`’s generic rather than exporting one of its own, so the call is
qualified or [`library(loo)`](https://mc-stan.org/loo/) is attached
first.

``` r

m1 <- flexybayes(yield ~ env, random = ~ gen, data = yw,
                 backend = "brms", chains = 2L, n_samples = 1500L,
                 warmup = 750L, seed = 1L, verbose = FALSE,
                 mcmc_verbose = FALSE)
m2 <- flexybayes(yield ~ env, random = ~ gen + gen:env, data = yw,
                 backend = "brms", chains = 2L, n_samples = 1500L,
                 warmup = 750L, seed = 1L, verbose = FALSE,
                 mcmc_verbose = FALSE)
#> Warning: There were 163 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 8 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.64, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Parts of the model have not converged (some Rhats are > 1.05). Be
#> careful when analysing the results! We recommend running more iterations and/or
#> setting stronger priors.
#> Warning: There were 163 divergent transitions after warmup. Increasing
#> adapt_delta above 0.8 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: flexyBayes: the sampler may not have converged -- 82 parameters with
#> Rhat >= 1.10 (max 1.64); min effective sample size 3. Treat the posterior with
#> caution: increase `warmup` / `n_samples`, raise `control = list(adapt_delta =
#> 0.95)` where the run reported divergent transitions, simplify the model, or
#> supply a more informative prior. Inspect the full diagnostics with summary().
#> Silence this warning via options(flexyBayes.silence_convergence_warning =
#> TRUE).
loo::loo_compare(loo::loo(m1), loo::loo(m2))
#> Warning: Found 59 observations with a pareto_k > 0.7 in model 'x$brms'. We
#> recommend to set 'moment_match = TRUE' in order to perform moment matching for
#> problematic observations.
#>        elpd_diff se_diff
#> x$brms   0.0       0.0  
#> x$brms -40.3       1.8
```

The first row is the better-predicting model. `elpd_diff` is the
expected log predictive density difference and `se_diff` its standard
error. A difference smaller than about twice its standard error is not a
distinction worth acting on.

This is not a hypothesis test and there is no p-value. The question it
answers is “which predicts better”, not “is the extra term significant”.

## 4. Triangulation

A single fit cannot tell you whether its answer is a property of the
data or of the algorithm that produced it. Two engines can.

``` r

tri <- triangulate(fit, m1)
tri$status
#> [1] "discordant"
table(tri$metrics$verdict)
#> 
#> concordant discordant 
#>          9          2
```

Three gates run before any comparison. A **fingerprint** gate refuses to
compare fits of different models, data or priors. A **diagnostics** gate
returns `inconclusive` rather than a verdict when either fit has not
converged. A **matched-prior** gate marks a parameter `not_compared`
when the two fits did not use the same prior for it.

``` r

tri$metrics[tri$metrics$verdict != "concordant",
            c("param", "mean_a", "mean_b", "sd_ratio", "verdict")]
#>     param    mean_a    mean_b  sd_ratio    verdict
#> 10  sigma 0.3860015 0.3866297 0.9606597 discordant
#> 11 sd_gen 0.3880369 0.3953745 0.9158931 discordant
```

Where the two disagree here it is on `sd_ratio` – the ratio of posterior
widths – not on location. INLA’s approximation reports a narrower
posterior for the variance components than the sampler does. Both fits
report clean convergence on their own; only the comparison shows it.

The reading is not that one engine is wrong. It is that for a quoted
interval on a variance component the sampler’s is the more conservative,
and that is worth knowing before the number goes into a paper.

## 5. Downstream packages

A fit works with the standard ecosystem: `emmeans` for estimated
marginal means, `marginaleffects` for contrasts and predictions,
`bayesplot` and `posterior` for draws, `loo` for comparison. Where a
method has no meaning on an engine it says so rather than returning
something misleading –
[`augment()`](https://generics.r-lib.org/reference/augment.html) on an
INLA fit refuses by name, because there is no retained per-row model
frame to append fitted values to.

## 6. Pitfalls

**`loo::loo()` needs a sampled fit.** An INLA fit has no pointwise
log-likelihood draws, so comparison runs on the brms side.

**Triangulation compares parameters, not conclusions.** Two fits can
agree on every component and still support different decisions if the
decision depends on a tail.

**`inconclusive` is not `discordant`.** It means a gate stopped the
comparison, usually because one fit had not converged. Fix the fit, then
compare.

**Interval semantics.**
[`predict()`](https://rdrr.io/r/stats/predict.html) returns posterior
expected-response intervals, not posterior-predictive ones: they exclude
residual observation noise.

## 7. Session information

    #> R version 4.5.2 (2025-10-31)
    #> Platform: aarch64-apple-darwin20
    #> Running under: macOS Tahoe 26.5.2
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
    #> [1] flexyBayes_0.10.0
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] tidyselect_1.2.1       dplyr_1.2.1            farver_2.1.2          
    #>   [4] loo_2.9.0              S7_0.2.2               INLA_25.10.19         
    #>   [7] TH.data_1.1-5          tensorA_0.36.2.1       digest_0.6.39         
    #>  [10] estimability_1.5.1     lifecycle_1.0.5        sf_1.1-0              
    #>  [13] StanHeaders_2.32.10    processx_3.9.0         survival_3.8-3        
    #>  [16] agridat_1.26           magrittr_2.0.5         posterior_1.7.1       
    #>  [19] compiler_4.5.2         rlang_1.3.0            tools_4.5.2           
    #>  [22] data.table_1.18.2.1    sn_2.1.3               knitr_1.51            
    #>  [25] bridgesampling_1.2-1   mnormt_2.1.2           curl_7.0.0            
    #>  [28] pkgbuild_1.4.8         classInt_0.4-11        plyr_1.8.9            
    #>  [31] RColorBrewer_1.1-3     multcomp_1.4-29        abind_1.4-8           
    #>  [34] KernSmooth_2.23-26     numDeriv_2016.8-1.1    withr_3.0.3           
    #>  [37] grid_4.5.2             stats4_4.5.2           colorspace_2.1-2      
    #>  [40] xtable_1.8-8           e1071_1.7-17           future_1.70.0         
    #>  [43] inline_0.3.21          ggplot2_4.0.3          globals_0.19.1        
    #>  [46] emmeans_2.0.2          scales_1.4.0           MASS_7.3-65           
    #>  [49] dichromat_2.0-0.1      cli_3.6.6              mvtnorm_1.3-6         
    #>  [52] reformulas_0.4.4       generics_0.1.4         otel_0.2.0            
    #>  [55] RcppParallel_5.1.11-2  reshape2_1.4.5         minqa_1.2.8           
    #>  [58] DBI_1.3.0              proxy_0.4-29           rstan_2.32.7          
    #>  [61] stringr_1.6.0          splines_4.5.2          bayesplot_1.15.0      
    #>  [64] parallel_4.5.2         matrixStats_1.5.0      marginaleffects_0.32.0
    #>  [67] brms_2.23.0            vctrs_0.7.3            boot_1.3-32           
    #>  [70] V8_8.0.1               Matrix_1.7-4           sandwich_3.1-1        
    #>  [73] jsonlite_2.0.0         callr_3.8.0            listenv_0.10.1        
    #>  [76] units_1.0-1            glue_1.8.1             parallelly_1.47.0     
    #>  [79] nloptr_2.2.1           codetools_0.2-20       distributional_0.8.1  
    #>  [82] stringi_1.8.7          gtable_0.3.6           QuickJSR_1.9.0        
    #>  [85] lme4_2.0-1             tibble_3.3.1           pillar_1.11.1         
    #>  [88] Brobdingnag_1.2-9      R6_2.6.1               fmesher_0.7.0         
    #>  [91] Rdpack_2.6.6           evaluate_1.0.5         lattice_0.22-7        
    #>  [94] rbibutils_2.4.1        backports_1.5.1        rstantools_2.6.0      
    #>  [97] class_7.3-23           MatrixModels_0.5-4     Rcpp_1.1.2            
    #> [100] coda_0.19-4.1          gridExtra_2.3          nlme_3.1-168          
    #> [103] checkmate_2.3.4        xfun_0.57              zoo_1.8-15            
    #> [106] pkgconfig_2.0.3
