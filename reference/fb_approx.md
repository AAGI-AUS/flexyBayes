# Construct an approximation-scheme specification

Names an approximation scheme and its tuning, validated against the
package's locked approximation registry. The result is consumed where an
approximate representation is accepted – most directly in a smooth term,
`s(x, representation = fb_approx("low_rank_smooth", rank = 5L))`.

## Usage

``` r
fb_approx(scheme, ...)
```

## Arguments

- scheme:

  Character(1): the approximation-scheme name. Must be a scheme
  registered in the approximation registry; an unregistered name is
  refused with the supported vocabulary.

- ...:

  Scheme-specific tuning carried on the object, for example `rank` for
  `"low_rank_smooth"`.

## Value

An `fb_approx` object: a classed list with the `scheme` and the tuning
kwargs as elements, and a `bias_bound_promise` attribute describing the
scheme's bias bound.

## Details

The only scheme with a smooth-basis fitting path in this release is
`"low_rank_smooth"`, a rank-`K` principal-component truncation of the
smooth basis. Its bias is the relative squared Frobenius error of the
truncation (Wood, 2017, chapter 5);
[`validate_approximation()`](https://aagi-aus.github.io/flexyBayes/reference/validate_approximation.md)
reports the realised capture against the threshold for a fitted model.

## See also

[`validate_approximation()`](https://aagi-aus.github.io/flexyBayes/reference/validate_approximation.md),
[`fb_engine()`](https://aagi-aus.github.io/flexyBayes/reference/fb_engine.md)

## Examples

``` r
a <- fb_approx("low_rank_smooth", rank = 5L)
a$scheme
#> [1] "low_rank_smooth"
a$rank
#> [1] 5
```
