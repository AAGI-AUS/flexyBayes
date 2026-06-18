# Print method for an internal `<fb_preflight>` summary

Diagnostic print of the design-memory preflight result: per-term
`design_memory_bytes` (formatted with thousand-separator " "),
`representation_class`, and the aggregate ceiling check. On refusal the
binding term + numeric ceiling appear below the per-term table.

## Usage

``` r
# S3 method for class 'fb_preflight'
print(x, ...)
```

## Arguments

- x:

  an `<fb_preflight>` object.

- ...:

  unused.

## Value

invisibly returns `x`.
