# Format method for fb_terms – one-line summary

Internal S3 method, registered for dispatch only. Used by R's default
print/format machinery for compact display in lists and `data.frame`
columns.

## Usage

``` r
# S3 method for class 'fb_terms'
format(x, ...)
```

## Arguments

- x:

  An `fb_terms` object, the parsed intermediate representation of a
  model.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

A length-one character string naming the ingest grammar, the response,
and the fixed and random term counts.
