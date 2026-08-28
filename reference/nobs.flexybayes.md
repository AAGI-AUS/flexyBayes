# Number of observations a flexyBayes fit was fitted to

A fit whose missing responses were augmented rather than dropped was
handed more rows than it has observations: the design cell of a lost
plot is still a row, carried as a latent quantity so the index set a
structured covariance is built over survives. The two counts are
therefore different numbers and the argument says which one is wanted.

## Usage

``` r
# S3 method for class 'flexybayes'
nobs(object, type = c("design", "observed"), ...)
```

## Arguments

- object:

  A `flexybayes` fit of any engine.

- type:

  Which count to return: `"design"` (default) or `"observed"`.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A single integer.

## Details

`type = "design"` (the default, and the historical behaviour) is the
number of rows the engine saw. `type = "observed"` is how many of those
carried an observed response, read from the record the missing-response
layer left on the fit.

## See also

[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md),
whose `n_design` and `n_observed` slots are the same two numbers.
