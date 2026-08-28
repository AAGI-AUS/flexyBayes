# Print a compact description of a flexyBayes fit

Reports the call the fit came from, the model the representation
describes, the engine that ran it, the design and observation counts,
and the headline posterior quantities. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for the full
coefficient and variance-component tables.

## Usage

``` r
# S3 method for class 'flexybayes'
print(x, ...)
```

## Arguments

- x:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the description it prints.

## Details

The header is shared with
[`print.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes_inla.md)
and the brms print, so the three cannot disagree about what the fit is.
Sampler lines print only on an engine that sampled: a nested Laplace
approximation has no chains and no warmup, and reporting them over one
tells a reader the fit is stochastic when it is not.
