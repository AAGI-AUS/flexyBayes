# Fixed-effect model matrix of a flexyBayes fit

Rebuilds the population-level design matrix from the formula and data
the fit carries. Random-effect and residual-structure terms are not part
of it: this is the basis the fixed-effect coefficients are expressed in.

## Usage

``` r
# S3 method for class 'flexybayes'
model.matrix(object, ...)
```

## Arguments

- object:

  A `flexybayes` fit carrying a fixed-effect formula and the data it was
  fitted to.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A numeric matrix with one row per observation and one column per
fixed-effect basis column, carrying the usual `assign` attribute.

## Details

The formula and data are resolved from whichever slots the object holds.
The brms-shaped emit keeps both under `$glm`; an INLA fit keeps its data
at `$data` and recovers its fixed-effect formula through
[`formula.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/formula.flexybayes_inla.md).
An object supplying neither is refused by name, since
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) on a
`NULL` formula silently returns the intercept-only matrix of the calling
frame's data.
