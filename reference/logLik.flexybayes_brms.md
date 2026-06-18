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

  Ignored.
