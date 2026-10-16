# Print method for lgm_refusal: LGM feasibility refusal

Internal S3 method, registered for dispatch only. Prints why flexyBayes
cannot fit the model on INLA and what to do about it, and nothing else.

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

Invisibly, `x` unchanged. Called for the refusal it prints.

## Details

The refusal used to print sixteen lines: a per-rule diagnostic string,
four lines explaining what `backend = "auto"` already does, an
unconditional sentence about a roadmap item for mixtures, printed
whether or not the model contained one, and a two-argument override
incantation. The R Journal review of 0.10.0 read one of these and said
so. A refusal is read by someone who wanted a fit, so it now carries the
reasons and one route onward, `rule_id` and `diagnostic` stay on the
object for anything that wants to branch on them, and `?lgm_gate`
documents the override.

## See also

[`vignette("flexyBayes-03-syntax")`](https://aagi-aus.github.io/flexyBayes/articles/flexyBayes-03-syntax.md)
lists every spelling the package covers and what each refusal means.
