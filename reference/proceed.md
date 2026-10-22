# Advance a deferred-execution object into its fit

Generic for the inspect-then-fit pattern. The `<flexybayes_review>`
method restores the RNG snapshot captured at review-object construction,
runs the deferred fit via the backend driver, caches the result
in-place, and returns the fit. A second call returns the cached fit.

## Usage

``` r
proceed(x, ...)

# S3 method for class 'flexybayes_review'
proceed(x, ...)
```

## Arguments

- x:

  A `<flexybayes_review>` object holding the deferred fit.

- ...:

  Reserved for future deferred-execution classes.

## Value

The fit object the originating call would have returned, of class
`flexybayes`. A second call returns the cached fit.

## Deprecated

`proceed()` is deprecated in flexyBayes 0.13.0 and will be removed in a
later release. Call
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
with the same arguments and `review_code = FALSE`, and set `seed` for a
reproducible brms fit. A call warns with a condition of class
`flexybayes_deprecated`.

## Examples

``` r
# The direct call fits the model proceed() would have fitted. The Stan
# compile is the expensive part of running this, so the fit is given a
# small budget.
if (requireNamespace("brms", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  # \donttest{
    fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "brms",
                      chains = 1L, n_samples = 200L, warmup = 100L,
                      seed = 1L)
    class(fit)
  # }
}
#> 
#> -- flexyBayes: brms (Stan passthrough) fit ------------------
#>   formula: y ~ x + (1 | g)
#>   family : gaussian (identity link)
#>   priors : 2 row(s)
#> ------------------------------------------------------------ 
#> 
#> Compiling Stan program...
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.23 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 1: Iteration: 101 / 300 [ 33%]  (Sampling)
#> Chain 1: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 1:                0.021 seconds (Sampling)
#> Chain 1:                0.034 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.08, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Parts of the model have not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations and/or setting stronger priors.
#> Warning: flexyBayes: the sampler may not have converged: min effective sample size 15. Treat the posterior with caution: increase `warmup` / `n_samples`, raise `control = list(adapt_delta = 0.95)` where the run reported divergent transitions, simplify the model, or supply a more informative prior. Inspect the full diagnostics with summary(). Silence this warning via options(flexyBayes.silence_convergence_warning = TRUE).
#> [1] "flexybayes_brms" "flexybayes"      "list"           
```
