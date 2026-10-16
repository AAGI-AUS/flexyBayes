# Test whether an object is an `fb_engine` specification

Test whether an object is an `fb_engine` specification

## Usage

``` r
is_fb_engine(x)
```

## Arguments

- x:

  Any R object. The test is a class check, not a structural one.

## Value

`TRUE` if `x` is an `fb_engine` object, `FALSE` otherwise.

## Deprecated

`is_fb_engine()` is deprecated in flexyBayes 0.13.0 and will be removed
in a later release. Use `inherits(x, "fb_engine")`. A call warns with a
condition of class `flexybayes_deprecated`.

## Examples

``` r
inherits(fb_engine("inla"), "fb_engine")
#> Warning: `fb_engine()` is deprecated as of flexyBayes 0.13.0 and will be removed in a later release.
#> Pass `backend` and engine settings directly to `flexybayes()`, for example `flexybayes(..., backend = "brms", chains = 4L)`.
#> [1] TRUE
inherits("inla", "fb_engine")
#> [1] FALSE
```
