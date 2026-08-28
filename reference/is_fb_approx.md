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

## Examples

``` r
a <- fb_approx("low_rank_smooth")
is_fb_approx(a)
#> [1] TRUE
is_fb_approx("low_rank_smooth")
#> [1] FALSE
```
