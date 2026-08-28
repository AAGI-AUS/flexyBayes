# Print method for flexybayes_inla

Internal S3 method. Brief one-screen description of an INLA fit produced
via `fb(... backend = "inla")` or `emit_inla()`.

## Usage

``` r
# S3 method for class 'flexybayes_inla'
print(x, ...)
```

## Arguments

- x:

  A `flexybayes_inla` object, the fit an INLA run returns.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the one-screen summary it prints.

## Details

Opens with the header every engine's print shares, so the three prints
cannot disagree about what the fit is or how many rows it saw, then adds
what belongs to this engine alone: the formula as it reached
`INLA::inla()`, the post-fit numerical-confirm verdict, and – on a fit
carrying a latent autoregressive field – the field's own parameters on
the correlation and standard-deviation scales.
