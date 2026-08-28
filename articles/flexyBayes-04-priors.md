# 4. Priors, and what they do to your answer

This page in one panel

- The problem : A Bayesian fit has a prior whether or not you chose one.
  Knowing what it is, and when it is doing work, is part of reporting
  the analysis.

- In ASReml : No equivalent. REML has no prior, which is exactly why a
  variance component at the boundary has no well-behaved interval.

- In flexyBayes :

  `prior_summary(fit)` to see it, `fb_prior(…)` to change it.

- What the posterior adds : Somewhere to put what you already know. With
  few levels that is the difference between an interval and a shrug.

- What it costs : An assumption you have to state.

## 1. There is always a prior

flexyBayes does not ask you to pick one before you can fit anything. It
picks a weakly informative default, tells you at fit time, and records
it on the fit.

The default is a bounded uniform on each standard deviation, with the
upper bound scaled to the response. Working on the standard-deviation
scale rather than the variance or precision scale is deliberate: it is
the scale the quantity is reported on, so a bound is a statement you can
reason about.

``` r

library(flexyBayes)
data(sleepstudy, package = "lme4")
fit <- flexybayes(Reaction ~ Days, random = ~ Subject,
                  data = sleepstudy, backend = "inla", verbose = FALSE)
#> Warning: flexyBayes: a variance component in this fit is at the boundary --
#> sd_Subject (upper credible bound 1.514e-10 against a residual SD of 47.78). The
#> term carries no signal the data could separate from residual noise, and its
#> variance has been absorbed into the residual. The convergence block reports a
#> converged mode regardless, because the optimiser did converge, to a solution
#> with that term at zero. This is not always a property of the data: the same
#> model on a random subset of the same trial can return a well-identified
#> component. Three routes: give the term an informative prior with fb_prior() so
#> it is not left to run to a boundary; fit the same model on the other engine for
#> a second reading; or report the component as unidentified rather than as zero.
#> Silence via options(flexyBayes.silence_boundary_collapse_warning = TRUE).
prior_summary(fit)
#> <prior_summary>  backend = inla
#>   Source: auto-default bounded uniform on SD (weakly-informative; half-Cauchy advised for small J)
#>     Upper bound U = 281.6438
#>     Scale basis    = identity_sd_uniform
#>     prior_fixed_sd:  not supplied -- the fixed effects carry INLA's own control.fixed defaults -- prec = 0.001 on the slopes (a standard deviation near 32), and prec.intercept = 0, which is flat
#> 
#> <fb_prior> 2 specifications
#>   sigma ~ uniform(lower = 0, upper = 281.6438)
#>   sd(group = "Subject") ~ uniform(lower = 0, upper = 281.6438)
```

Every component is listed with the prior that was actually used, not the
prior that was requested. On a fit where an engine substituted its own
default, that is what appears.

## 2. When the prior is doing work

With plenty of levels the data dominate and the default is invisible.
With few levels it is not. The comparison below fits the same model to
all eighteen subjects and then to four.

``` r

small <- sleepstudy[sleepstudy$Subject %in% levels(sleepstudy$Subject)[1:4], ]
small$Subject <- droplevels(small$Subject)

fit_small <- flexybayes(Reaction ~ Days, random = ~ Subject,
                        data = small, backend = "inla", verbose = FALSE)
#> Warning: flexyBayes: a variance component in this fit is at the boundary --
#> sd_Subject (upper credible bound 7.892e-06 against a residual SD of 60.59). The
#> term carries no signal the data could separate from residual noise, and its
#> variance has been absorbed into the residual. The convergence block reports a
#> converged mode regardless, because the optimiser did converge, to a solution
#> with that term at zero. This is not always a property of the data: the same
#> model on a random subset of the same trial can return a well-identified
#> component. Three routes: give the term an informative prior with fb_prior() so
#> it is not left to run to a boundary; fit the same model on the other engine for
#> a second reading; or report the component as unidentified rather than as zero.
#> Silence via options(flexyBayes.silence_boundary_collapse_warning = TRUE).

rbind(
  `18 subjects` = tidy(fit, effects = "random")[2, 2:5],
  `4 subjects`  = tidy(fit_small, effects = "random")[2, 2:5]
)
#>                 estimate    std.error     conf.low    conf.high
#> 18 subjects 1.240777e-10 1.459947e-11 9.489472e-11 1.514309e-10
#> 4 subjects  4.645489e-06 1.416672e-06 2.350807e-06 7.892142e-06
```

The point estimate moves a little. The interval widens a great deal.
That width is the honest report: four levels carry little information
about the spread of a population, and nothing in the arithmetic can
manufacture more.

## 3. Changing it

[`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
takes one-sided formulas, each naming a target on the left and a
distribution on the right.

``` r

p <- fb_prior(
  sigma                 ~ pc(upper = 2, prob = 0.05),
  sd(group = "Subject") ~ half_normal(scale = 30)
)
p
#> <fb_prior> 2 specifications
#>   sigma ~ pc(upper = 2, prob = 0.05)
#>   sd(group = "Subject") ~ half_normal(scale = 30)
```

Three targets are addressable: `sigma` for the residual, `sd(group =)`
for a named random-effect group, and `b()` for a fixed-effect
coefficient.

`pc(upper, prob)` is a penalised-complexity prior (Simpson and others,
2017), read as “the probability that this standard deviation exceeds
`upper` is `prob`”. It shrinks toward zero unless the data pull it away,
which is a defensible default for a component that may not be there at
all.

``` r

fit_p <- flexybayes(Reaction ~ Days, random = ~ Subject,
                    data = small, prior = p,
                    backend = "inla", verbose = FALSE)
#> Warning: flexyBayes: a variance component in this fit is at the boundary --
#> sd_Subject (upper credible bound 3.232e-06 against a residual SD of 39.81). The
#> term carries no signal the data could separate from residual noise, and its
#> variance has been absorbed into the residual. The convergence block reports a
#> converged mode regardless, because the optimiser did converge, to a solution
#> with that term at zero. This is not always a property of the data: the same
#> model on a random subset of the same trial can return a well-identified
#> component. Three routes: give the term an informative prior with fb_prior() so
#> it is not left to run to a boundary; fit the same model on the other engine for
#> a second reading; or report the component as unidentified rather than as zero.
#> Silence via options(flexyBayes.silence_boundary_collapse_warning = TRUE).
tidy(fit_p, effects = "random")
#>         term     estimate    std.error     conf.low    conf.high
#> 1      sigma 3.984534e+01 2.603752e+00 3.489204e+01 4.509775e+01
#> 2 sd_Subject 2.870273e-06 2.181719e-07 2.400414e-06 3.231656e-06
```

## 4. Reporting

Two rules make a Bayesian analysis reproducible by someone else.

State the prior.
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
produces it in a form that can be pasted into a methods section, and it
reports what ran rather than what was asked for.

Say whether it mattered. Refitting under a different reasonable prior
and reporting both is the standard check. Where the two agree, the data
are doing the work. Where they do not, the prior is part of the result
and belongs in the write-up.

## 5. Pitfalls

**A prior on the wrong scale.** The DSL is on the standard-deviation
scale throughout. A bound of 2 means two units of the response, not two
units of variance.

**An engine default is not silence.** Rows marked `engine default` in
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
carry whatever INLA or Stan chose. Those are priors too, and they are
not always weak.

**Bounded uniform has a bound.** The default’s upper limit is scaled to
the response, and a component genuinely larger than the bound cannot be
recovered.
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
shows the number.

**A prior cannot rescue an unidentified component.** If a term carries
no signal, an informative prior returns the prior. That is the intended
behaviour, and the boundary warning is the flag that it happened.

## 6. Session information

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
