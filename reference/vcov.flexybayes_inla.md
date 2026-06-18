# Posterior covariance of a per-row INLA fit's fixed effects

Monte-Carlo estimate of the joint posterior covariance of the fixed
effects, computed from `inla.posterior.sample()`. The marginal standard
deviations match `summary.fixed$sd`; the off-diagonals carry the joint
dependence that contrast / marginal-mean standard errors require.
Because the estimate is sampling-based it varies slightly between calls;
raise `n_samples` for a tighter estimate.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
vcov(object, n_samples = 2000L, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- n_samples:

  Posterior sample size for the covariance estimate (default 2000).

- ...:

  Ignored.

## Value

Posterior covariance matrix of the fixed effects, with `summary.fixed`
rownames as dimnames.
