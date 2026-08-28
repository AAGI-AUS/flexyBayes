# Glance at a flexyBayes fit

Returns a one-row `data.frame` of model-level statistics: the response
family and link, the observation and parameter counts, the
log-likelihood, the worst convergence diagnostics across monitored
parameters, and the wall-clock run time.

## Usage

``` r
# S3 method for class 'flexybayes'
glance(x, ...)

# S3 method for class 'flexybayes_inla'
glance(x, ...)
```

## Arguments

- x:

  A flexyBayes fit (`flexybayes` or `flexybayes_brms`).

- ...:

  Currently unused; present for generic compatibility.

## Value

A one-row `data.frame` with `nobs`, `npar`, `logLik`, `family`, `link`,
`chains`, `samples`, `max_rhat`, `min_ess`, and `run_time`.

## See also

[`tidy.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes.md)

## Examples

``` r
# Not run: glance() reports sampler diagnostics, so the example needs
# a sampled fit and therefore a Stan compile. The fragment is
# complete and runs as written where a C++ toolchain is present.
if (FALSE) { # \dontrun{
set.seed(1)
dat <- data.frame(yield = rnorm(40), env = factor(rep(1:4, each = 10)))
fit <- flexybayes(yield ~ env, data = dat, backend = "brms",
                  chains = 1L, n_samples = 200L, warmup = 100L)
glance(fit)
} # }
```
