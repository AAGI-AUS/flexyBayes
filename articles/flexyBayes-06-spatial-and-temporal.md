# 6. Spatial and temporal structure

This page in one panel

- The problem : Plots near each other in a field yield alike for reasons
  the design did not plan. Ignoring that inflates the apparent precision
  of every variety comparison.

- In ASReml :

  `asreml(yield ~ 1, random = ~ gen, residual = ~ ar1(row):ar1(col))`

- In flexyBayes :

  `flexybayes(yield ~ 1, random = ~ gen + ar1(row):ar1(col))` – the
  field moves to the random side, for the reason in Section 2.

- What the posterior adds : The two correlations and the field standard
  deviation arrive with intervals, so a field that has not identified is
  visible rather than inferred from a converged flag.

- What it costs : INLA only, a complete grid, and a four-parameter model
  where ASReml fits three.

## 1. Design first, correlation second

A randomised block design already removes planned variation. Spatial
correlation is what remains: a fertility gradient, a drainage line, a
sprayer pass. Modelling it is not a substitute for a good design, and a
badly designed trial cannot be rescued here.

The practical consequence of ignoring it is not bias in the variety
means so much as over-confidence in them. Neighbouring plots that share
an unmodelled gradient carry less independent information than the
degrees-of-freedom count suggests.

## 2. Three parameters or four

ASReml writes the field as a *residual*:
`residual = ~ ar1(row):ar1(col)` is a single correlated residual with no
independent plot error – three parameters.

INLA represents a separable autoregressive process as a latent field
*plus* observation noise – four parameters. The two coincide only as the
nugget goes to zero. flexyBayes will not fit the four-parameter model
under the three-parameter name, so the residual spelling is refused and
the field is written on the random side instead.

``` r

library(flexyBayes)
sn <- agridat::stroup.nin
sn$row <- factor(sn$row)
sn$col <- factor(sn$col)
sn$yield_z <- as.numeric(scale(sn$yield))

refused <- tryCatch(
  flexybayes(yield_z ~ 1, random = ~ gen,
             residual = ~ ar1(row):ar1(col), data = sn,
             backend = "inla", verbose = FALSE),
  error = function(e) e
)
cat(strwrap(conditionMessage(refused), 72)[1:4], sep = "\n")
#> residual = ~ ar1(row):ar1(col) asks for ASReml's separable
#> autoregressive residual: a single correlated residual with no
#> independent plot error, three parameters on a field grid (row
#> correlation, column correlation, one variance). INLA represents a
```

## 3. The field on a real trial

`agridat::stroup.nin` is a Nebraska wheat trial: 56 varieties on an 11
by 22 plot array, with 18 plots unharvested. The array has 242 nodes and
the data carry 224 rows.

A separable field needs every node of the array, which is the same
requirement ASReml states for `residual = ~ ar1:ar1`. The unharvested
plots are kept as design cells with the response missing, which is the
default `na_action = "augment"` and ASReml’s `na.method(y = "include")`.

``` r

fit <- flexybayes(yield_z ~ 1, random = ~ ar1(row):ar1(col),
                  data = sn, backend = "inla", verbose = FALSE)
vc_field <- tidy(fit, effects = "random")
vc_field
#>         term  estimate  std.error  conf.low conf.high
#> 1      sigma 0.5282562 0.03436625 0.4649726 0.5996954
#> 2 sd_spatial 0.8058180 0.13904051 0.5759922 1.1197668
#> 3    rho_row 0.8321471 0.06403936 0.6810619 0.9294357
#> 4    rho_col 0.8884133 0.04463320 0.7821983 0.9550196
```

Four numbers describe the field. `rho_row` and `rho_col` are the
correlations between adjacent rows and adjacent columns. `sd_spatial` is
the field’s standard deviation and `sigma` the independent plot noise.
Here the field carries more variation than the noise does, and both
correlations are high: this trial has strong spatial structure, which is
why it is the textbook example.

Printing all four is how a *failed* field shows up. When a field does
not identify, `sd_spatial` collapses toward zero and the correlations
revert to their prior, spanning almost the whole interval. The package
raises a warning when that happens rather than leaving it to be noticed.

## 4. Adding a smooth trend: does it earn its place?

Field trials often show *both* short-range correlation and a large-scale
gradient. The textbook decomposition puts a smooth trend on the mean and
an autoregressive field on the residual (Gilmour, Cullis, & Verbyla,
1997). Both go on the random side here, with `spl()` for the trend.

`spl()` is a penalised spline, fitted on INLA as a second-order random
walk. brms refuses it, so this section is INLA only.

``` r

sn$rown <- as.numeric(as.character(sn$row))
sn$coln <- as.numeric(as.character(sn$col))

fit_trend <- flexybayes(
  yield_z ~ 1,
  random = ~ spl(rown) + spl(coln) + ar1(row):ar1(col),
  data   = sn, backend = "inla", verbose = FALSE
)
vc_trend <- tidy(fit_trend, effects = "random")
vc_trend
#>         term    estimate   std.error    conf.low  conf.high
#> 1      sigma 0.525881003 0.036681753 0.458447688 0.60224850
#> 2    sd_rown 0.009980914 0.006196669 0.003397621 0.02772970
#> 3    sd_coln 0.009642915 0.005766860 0.003341367 0.02602769
#> 4 sd_spatial 0.641975610 0.118328937 0.446169200 0.90905686
#> 5    rho_row 0.760851783 0.094199636 0.538902725 0.90391146
#> 6    rho_col 0.838559437 0.070866395 0.667958793 0.94164888
```

The two splines collapse, and the field absorbs what they were offered:

``` r

pick <- function(v, nm) v$estimate[v$term == nm]
data.frame(
  component = c("sd_spatial", "sd_rown", "sd_coln"),
  field_only  = c(pick(vc_field, "sd_spatial"), NA, NA),
  with_trend  = c(pick(vc_trend, "sd_spatial"),
                  pick(vc_trend, "sd_rown"), pick(vc_trend, "sd_coln"))
)
#>    component field_only  with_trend
#> 1 sd_spatial   0.805818 0.641975610
#> 2    sd_rown         NA 0.009980914
#> 3    sd_coln         NA 0.009642915
```

Both spline standard deviations sit near zero against a residual an
order of magnitude larger, and the field’s standard deviation falls once
the splines are present. The reading is that the autoregressive field
was already carrying the large-scale gradient, so once the splines are
offered the same variation there is nothing left for them to explain.

That is a result about `stroup.nin`, not about splines. On a trial with
a strong monotone gradient and weak short-range correlation the split
goes the other way. The point is that the fit says which, rather than
leaving the composite as an untested recommendation – and a collapsed
component here raises the boundary warning, so it cannot be read as “the
trend is zero” when it means “the trend is not separable from the
field”.

## 5. Correlation in time

The same machinery covers a regular time series: `ar1(t)` on a factor
whose levels are in time order.

``` r

set.seed(7)
n <- 120
d <- data.frame(t = factor(seq_len(n)), x = rnorm(n))
d$y <- 2 + 0.5 * d$x + as.numeric(arima.sim(list(ar = 0.7), n = n)) * 0.8

fit_t <- flexybayes(y ~ x, random = ~ ar1(t), data = d,
                    backend = "inla", verbose = FALSE)
vc_time <- tidy(fit_t, effects = "random")
vc_time
#>         term  estimate  std.error   conf.low conf.high
#> 1      sigma 0.2322977 0.12357981 0.05708243 0.7900921
#> 2 sd_spatial 1.1558414 0.14548427 0.91382879 1.4821413
#> 3    rho_row 0.7373798 0.06816154 0.58389637 0.8503652
```

The series was generated with an autoregressive coefficient of 0.7 and
the posterior recovers it. The component names come from the shared
spatial machinery, so the temporal correlation is reported as `rho_row`.

## 6. What is not covered

| Wanted | Status |
|----|----|
| ASReml’s nugget-free separable residual | refused by name (Section 2) |
| Irregular geography (SPDE / Matern) | not implemented |
| Areal data (BYM2) | not implemented |
| Continuous-time correlation (Ornstein-Uhlenbeck) | not implemented |
| Separable space-time interaction | not implemented |

For these, `R-INLA` directly or
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) with a spatial
basis is the route. A fit from either can still be compared against a
flexyBayes fit through the
[`fb_as_draws_simple()`](https://aagi-aus.github.io/flexyBayes/reference/fb_as_draws_simple.md)
interface.

## 7. Pitfalls

**Row and column must be factors, in field order.** Correlation is read
along the stored level order. Labels that sort lexically rather than
numerically give the wrong neighbour structure. Set the levels
explicitly when in doubt.

**A complete grid is required.** A missing plot, a duplicated
combination, or a non-rectangular array raises
`ar1_spatial_requires_complete_grid`. Supply the field-book rows for
every sown plot with the response set to `NA`, or pad with
[`fb_complete_grid()`](https://aagi-aus.github.io/flexyBayes/reference/fb_complete_grid.md).

**Too few rows or columns.** A 4 by 3 array cannot identify two
correlations, a field variance and a plot variance at once. The symptom
is a collapsed field with correlations at their prior width.

**`backend = "auto"` refuses rather than substituting.** brms has no
lowering for a Kronecker autoregressive precision, so a field model
under `auto` raises `auto_no_active_route` instead of quietly fitting
something else.

## 8. Session information

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
