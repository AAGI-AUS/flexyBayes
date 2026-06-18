# In-sample fitted values from a per-row INLA fit

Returns INLA's posterior-mean fitted values (response scale) for the
observed rows, taken from `summary.fitted.values`. Without this method
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) dispatched to
[`stats::fitted.default`](https://rdrr.io/r/stats/fitted.values.html),
which silently returned `NULL` for an INLA fit because the object does
not populate a `$glm` slot.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
fitted(object, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- ...:

  Ignored.

## Value

A numeric vector of posterior-mean fitted values, one per observation.
