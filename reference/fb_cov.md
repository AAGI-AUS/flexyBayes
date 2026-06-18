# Construct a structured-covariance carrier

Wraps a known covariance, precision, Cholesky factor, block-diagonal
list, or low-rank factor in a classed carrier object that the structured
random-effect terms consume, most directly
`vm(geno, cov = fb_cov(G, type = "chol"))`. It is the v0.4.0
constructor-noun replacement for the legacy `vm()` keyword forms
(`chol = `, `precision = `, `blocks = `, `low_rank_factor = `), which
deprecate at this release and are removed at v0.5.0.

## Usage

``` r
fb_cov(M, type = "dense", levels = NULL, scheme = NULL, ...)
```

## Arguments

- M:

  The covariance carrier. For `type = "blocks"` a base-R list of K
  square covariance matrices; otherwise a square numeric matrix (base-R
  or Matrix). Inside a formula this argument names a matrix resolved
  through `known_matrices`.

- type:

  Character(1): the carrier type, one of `"dense"`, `"chol"`,
  `"precision"`, `"blocks"`, `"low_rank"`. Defaults to `"dense"`.

- levels:

  Optional character vector of grouping-factor level labels the
  carrier's rows / columns align to. Carried as metadata; the fit-time
  validator checks alignment against the fitted factor.

- scheme:

  Character(1), required when `type = "low_rank"`: the name of a
  registered approximation scheme (see
  [`validate_approximation()`](https://aagi-aus.github.io/flexyBayes/reference/validate_approximation.md)).
  Ignored for the exact types.

- ...:

  Reserved for forward carrier-type options; currently carried verbatim
  on the object.

## Value

An `fb_cov` object: a classed list with the matrix `M`, the `type`,
optional `levels`, and (for low-rank) `scheme` as elements, plus
`representation_class` and `validation_summary` attributes.

## Details

The carrier's `type` selects how the downstream emit path derives a
covariance square root. `"dense"` takes the lower Cholesky of the
supplied covariance; `"chol"` uses the supplied factor directly;
`"precision"` inverts via the Cholesky of the precision matrix;
`"blocks"` assembles a block-diagonal covariance from a list of
per-block matrices; `"low_rank"` pairs a rank-K factor with a registered
approximation scheme (reserved at v0.4.0 – the carrier vocabulary is
active and validated, while the approximate-covariance fit route
activates in a subsequent release).

Inside a model formula the carrier is written inline – the matrix
argument names a matrix passed through `known_matrices`, exactly as the
legacy keyword forms did. The construction-time check is a light
structural probe (shape, lower-triangular pattern, symmetry, block
count); the full level-aware validation runs at fit time against the
grouping factor.

## See also

[`fb_approx()`](https://aagi-aus.github.io/flexyBayes/reference/fb_approx.md),
[`fb_engine()`](https://aagi-aus.github.io/flexyBayes/reference/fb_engine.md),
[`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md),
[`validate_approximation()`](https://aagi-aus.github.io/flexyBayes/reference/validate_approximation.md)

## Examples

``` r
G <- crossprod(matrix(rnorm(9L), 3L, 3L))
fb_cov(G, type = "dense")
#> <fb_cov> type = "dense"
#>   representation: dense_cov
#>   carrier: dense: 3x3
L <- t(chol(G))
fb_cov(L, type = "chol")
#> <fb_cov> type = "chol"
#>   representation: chol_cov
#>   carrier: chol: 3x3 lower-triangular
```
