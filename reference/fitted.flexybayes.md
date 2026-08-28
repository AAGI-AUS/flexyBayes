# Extract fitted values

Extract fitted values

## Usage

``` r
# S3 method for class 'flexybayes'
fitted(object, ...)
```

## Arguments

- object:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A numeric vector of in-sample fitted values on the response scale, one
per row of the fitted data, each the posterior mean of that
observation's conditional expectation. Rows carried as latent under
`na_action = "augment"` receive a fitted value like any other.
