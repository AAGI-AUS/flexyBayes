# Log-likelihood on the brms path

Delegates to `brms::log_lik()` then sums pointwise log-likelihood across
observations and averages across draws. The `df` attribute carries the
parameter count from `$extras$model_info`; `nobs` carries the
observation count.

## Usage

``` r
# S3 method for class 'flexybayes_brms'
logLik(object, ...)
```

## Arguments

- object:

  A `flexybayes_brms` object.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

An object of class `logLik`: the pointwise log-likelihood summed over
observations and averaged across draws, carrying `df` (parameter count)
and `nobs` (observation count) attributes.
