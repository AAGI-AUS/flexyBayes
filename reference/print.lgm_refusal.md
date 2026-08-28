# Print method for lgm_refusal – LGM feasibility refusal

Internal S3 method, registered for dispatch only. Emits the structured
refusal template: rule id + one-line gloss + diagnostic + re-route
hint + override hint + docs pointer.

## Usage

``` r
# S3 method for class 'lgm_refusal'
print(x, ...)
```

## Arguments

- x:

  An `lgm_refusal` object, as returned by the feasibility gate when a
  model falls outside the latent Gaussian class.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the refusal template it prints.
