# Extract per-parameter posterior draws from a model fit

S3 generic that extracts a named list of numeric posterior-draw vectors
from a fit, keyed by the engine's own parameter names. The methods that
reach an active engine are `flexybayes_brms` and `flexybayes_inla`. The
`brmsfit` method reads a brms fit inside a flexyBayes fit, a two-part
non-linear fit from
[`fb_refit()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refit.md)
under the names of the linear model. User-defined methods can extend the
generic.

## Usage

``` r
fb_as_draws_simple(fit, ...)

# S3 method for class 'flexybayes_inla'
fb_as_draws_simple(fit, n_samples = 1000L, ...)

# S3 method for class 'flexybayes_brms'
fb_as_draws_simple(fit, ...)

# S3 method for class 'brmsfit'
fb_as_draws_simple(fit, ...)

# Default S3 method
fb_as_draws_simple(fit, ...)
```

## Arguments

- fit:

  A model fit object carrying a posterior, dispatched on by the methods
  listed above.

- ...:

  Method-specific arguments: `n_samples` for the INLA method. The brms
  method returns every stored draw and refuses `n_samples` with
  `method_argument_not_supported`.

## Value

A named list of numeric vectors, one element per parameter, each element
holding that parameter's posterior draws.

## Details

On an INLA fit the draws come from the fitted approximation: each draw
takes one of INLA's hyperparameter configurations, chosen with its
posterior weight, and draws the latent field from that configuration's
Gaussian approximation, as
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
does. Every random number comes from R's generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call makes
the draws repeatable.

## Examples

``` r
# A named list of draw vectors, with no posterior dependency.
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
                    verbose = FALSE)
  dr <- fb_as_draws_simple(fit)
  names(dr)
}
#>  [1] "Predictor:1"                            
#>  [2] "Predictor:2"                            
#>  [3] "Predictor:3"                            
#>  [4] "Predictor:4"                            
#>  [5] "Predictor:5"                            
#>  [6] "Predictor:6"                            
#>  [7] "Predictor:7"                            
#>  [8] "Predictor:8"                            
#>  [9] "Predictor:9"                            
#> [10] "Predictor:10"                           
#> [11] "Predictor:11"                           
#> [12] "Predictor:12"                           
#> [13] "Predictor:13"                           
#> [14] "Predictor:14"                           
#> [15] "Predictor:15"                           
#> [16] "Predictor:16"                           
#> [17] "Predictor:17"                           
#> [18] "Predictor:18"                           
#> [19] "Predictor:19"                           
#> [20] "Predictor:20"                           
#> [21] "Predictor:21"                           
#> [22] "Predictor:22"                           
#> [23] "Predictor:23"                           
#> [24] "Predictor:24"                           
#> [25] "Predictor:25"                           
#> [26] "Predictor:26"                           
#> [27] "Predictor:27"                           
#> [28] "Predictor:28"                           
#> [29] "Predictor:29"                           
#> [30] "Predictor:30"                           
#> [31] "Predictor:31"                           
#> [32] "Predictor:32"                           
#> [33] "Predictor:33"                           
#> [34] "Predictor:34"                           
#> [35] "Predictor:35"                           
#> [36] "Predictor:36"                           
#> [37] "Predictor:37"                           
#> [38] "Predictor:38"                           
#> [39] "Predictor:39"                           
#> [40] "Predictor:40"                           
#> [41] "Predictor:41"                           
#> [42] "Predictor:42"                           
#> [43] "Predictor:43"                           
#> [44] "Predictor:44"                           
#> [45] "Predictor:45"                           
#> [46] "Predictor:46"                           
#> [47] "Predictor:47"                           
#> [48] "Predictor:48"                           
#> [49] "Predictor:49"                           
#> [50] "Predictor:50"                           
#> [51] "Predictor:51"                           
#> [52] "Predictor:52"                           
#> [53] "Predictor:53"                           
#> [54] "Predictor:54"                           
#> [55] "Predictor:55"                           
#> [56] "Predictor:56"                           
#> [57] "Predictor:57"                           
#> [58] "Predictor:58"                           
#> [59] "Predictor:59"                           
#> [60] "Predictor:60"                           
#> [61] "g:1"                                    
#> [62] "g:2"                                    
#> [63] "g:3"                                    
#> [64] "g:4"                                    
#> [65] "g:5"                                    
#> [66] "g:6"                                    
#> [67] "(Intercept):1"                          
#> [68] "x:1"                                    
#> [69] "Precision for the Gaussian observations"
#> [70] "Precision for g"                        
# }
```
