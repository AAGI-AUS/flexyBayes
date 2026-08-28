# Extract variance-covariance matrix of fixed effects

Extract variance-covariance matrix of fixed effects

## Usage

``` r
# S3 method for class 'flexybayes'
vcov(object, ...)
```

## Arguments

- object:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

The posterior covariance matrix of the fixed effects, square with one
row and column per coefficient and dimnames taken from
[`coef()`](https://rdrr.io/r/stats/coef.html). This is a posterior
covariance, not a sampling-theory variance estimate, though the
downstream packages that consume it treat it as one.
