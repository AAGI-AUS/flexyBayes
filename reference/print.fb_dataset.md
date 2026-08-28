# Print method for an internal `<fb_dataset>` wrapper

Diagnostic print: one-line header from
[`format.fb_dataset()`](https://aagi-aus.github.io/flexyBayes/reference/format.fb_dataset.md),
then a per-column types row and a per-dictionary level-count row.

## Usage

``` r
# S3 method for class 'fb_dataset'
print(x, ...)
```

## Arguments

- x:

  An `<fb_dataset>` object, the wrapper the preflight layer carries in
  place of the raw data.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the header and the per-column rows
it prints.
