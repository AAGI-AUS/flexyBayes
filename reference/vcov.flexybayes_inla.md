# Posterior covariance of a per-row INLA fit's fixed effects

The joint posterior covariance of the fixed effects, the matrix that
contrast and marginal-mean standard errors need. It is computed exactly
from the Gaussian approximations INLA stores for the fit, averaged over
INLA's hyperparameter configurations, so it involves no random draws:
every call on the same fit returns the same matrix and the session's
random number stream is left untouched. It is the covariance of the
draws
[`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
would give, without their Monte Carlo error. The standard deviations on
its diagonal agree closely with `summary.fixed$sd`, which comes from
INLA's marginal densities.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
vcov(object, n_samples = 2000L, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- n_samples:

  Not used, as the covariance is exact. Accepted so that calls written
  for earlier versions, which estimated the matrix from `n_samples`
  posterior draws, still run.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Posterior covariance matrix of the fixed effects, with the coefficient
names of [`coef()`](https://rdrr.io/r/stats/coef.html) as dimnames.
