# Extract model formula

Extract model formula

## Usage

``` r
# S3 method for class 'flexybayes'
formula(x, ...)
```

## Arguments

- x:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

The fixed-effect (population-level) formula as a `formula` object.
Random-effect and residual-structure terms are not part of it: read them
from the fit's `fb_terms` intermediate representation, or from
[`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md)
before fitting.
