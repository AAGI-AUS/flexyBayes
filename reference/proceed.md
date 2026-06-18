# Advance a deferred-execution object into its fit

Generic for the inspect-then-fit pattern. The `<flexybayes_review>`
method restores the RNG snapshot captured at review-object construction,
runs the deferred fit via the backend driver, caches the result
in-place, and returns the fit. A second call returns the cached fit.

## Usage

``` r
proceed(x, ...)

# S3 method for class 'flexybayes_review'
proceed(x, ...)
```

## Arguments

- x:

  a `<flexybayes_review>` object.

- ...:

  reserved for future deferred-execution classes (e.g., deferred
  triangulation).

## Value

the fit object the originating call would have returned (class
`flexybayes`).
