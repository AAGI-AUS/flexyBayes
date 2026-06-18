# Extract canonical-named posterior means from an fb_greta() fit

Returns the posterior means of the target greta_arrays under their
canonical names. Differs from
[`coef.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes.md)
only in the source of the names (canonical map rather than fixed-effect
formula terms).

## Usage

``` r
# S3 method for class 'flexybayes_direct_greta'
coef(object, ...)
```

## Arguments

- object:

  A `flexybayes_direct_greta` object.

- ...:

  Additional arguments (ignored).

## Value

Named numeric vector of posterior mean target parameters.
