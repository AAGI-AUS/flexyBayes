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

## Examples

``` r
is_fb_engine(fb_engine("inla"))
#> [1] TRUE
is_fb_engine("inla")
#> [1] FALSE
```
