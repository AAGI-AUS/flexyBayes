# Credible intervals for flexybayes_glm

Posterior quantile credible intervals, computed from the fixed-effect
draws.

## Usage

``` r
# S3 method for class 'flexybayes_glm'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `flexybayes_glm` object, reached as `fit$glm` on any fitted
  `flexybayes` object.

- parm:

  A character vector naming the parameters to report, or `NULL` (the
  default) for every fixed effect.

- level:

  A single numeric giving the credible level, defaulting to 0.95.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

A two-column matrix of interval bounds, carrying an `interval_basis`
attribute of `"posterior_quantile"` or `"normal_approximation"`.

## Details

Where the draws are unavailable the method falls back to the
normal-approximation interval \\\hat\beta \pm z\_{\alpha/2}
\mathrm{sd}\\, and marks the returned matrix with
`attr(, "interval_basis") == "normal_approximation"`. The two agree only
for a symmetric posterior; on a skewed one they can differ materially,
which is why the basis is reported rather than assumed. Earlier versions
documented this method as quantile-based while always returning the
approximation.
