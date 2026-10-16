# Test whether an object is an `fb_approx` specification

Test whether an object is an `fb_approx` specification

## Usage

``` r
is_fb_approx(x)
```

## Arguments

- x:

  Any R object. The test is a class check, not a structural one.

## Value

`TRUE` if `x` is an `fb_approx` object, `FALSE` otherwise.

## Deprecated

`is_fb_approx()` is deprecated in flexyBayes 0.13.0 and will be removed
in a later release. Use `inherits(x, "fb_approx")`. A call warns with a
condition of class `flexybayes_deprecated`.

## Examples

``` r
a <- fb_approx("low_rank_smooth")
inherits(a, "fb_approx")
#> [1] TRUE
inherits("low_rank_smooth", "fb_approx")
#> [1] FALSE
```
