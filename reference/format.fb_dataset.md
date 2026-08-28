# Format method for an internal `<fb_dataset>` wrapper

One-line summary of the internal dataset wrapper used by the preflight
layer.

## Usage

``` r
# S3 method for class 'fb_dataset'
format(x, ...)
```

## Arguments

- x:

  An `<fb_dataset>` object, the wrapper the preflight layer carries in
  place of the raw data.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

A length-one character string summarising the wrapped dataset.
