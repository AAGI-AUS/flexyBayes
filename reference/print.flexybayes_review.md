# Print method for a deferred-execution review object

Two-line summary of a `<flexybayes_review>` object returned when
`review_code = TRUE` is passed to
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).
The first line reports the backend and the IR (intermediate
representation) dimensions; the second is a prompt directing the user to
[`cat_code()`](https://aagi-aus.github.io/flexyBayes/reference/cat_code.md)
for inspection and
[`proceed()`](https://aagi-aus.github.io/flexyBayes/reference/proceed.md)
to advance the fit.

## Usage

``` r
# S3 method for class 'flexybayes_review'
print(x, ...)
```

## Arguments

- x:

  A `<flexybayes_review>` object as returned by
  `flexybayes(review_code = TRUE)`.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the two summary lines it prints.
