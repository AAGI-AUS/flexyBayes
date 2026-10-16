# Test whether an object is an `fb_cov` carrier

Test whether an object is an `fb_cov` carrier

## Usage

``` r
is_fb_cov(x)
```

## Arguments

- x:

  Any R object. The test is a class check, not a structural one.

## Value

`TRUE` if `x` is an `fb_cov` object, `FALSE` otherwise.

## Deprecated

`is_fb_cov()` is deprecated in flexyBayes 0.13.0 and will be removed in
a later release. Use `inherits(x, "fb_cov")`. A call warns with a
condition of class `flexybayes_deprecated`.

## Examples

``` r
K <- diag(3)
dimnames(K) <- list(letters[1:3], letters[1:3])
inherits(fb_cov(K, type = "dense"), "fb_cov")
#> [1] TRUE
inherits(K, "fb_cov")
#> [1] FALSE
```
