# Print method for an internal `<fb_preflight_refusal>` object

Three-line diagnostic: reason code, binding term + its byte estimate,
and the active ceiling with the suggested override.

## Usage

``` r
# S3 method for class 'fb_preflight_refusal'
print(x, ...)
```

## Arguments

- x:

  An `<fb_preflight_refusal>` object, carrying the reason code and the
  binding term.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the three diagnostic lines it
prints.
