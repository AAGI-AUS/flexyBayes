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

## Examples

``` r
K <- diag(3)
dimnames(K) <- list(letters[1:3], letters[1:3])
is_fb_cov(fb_cov(K, type = "dense"))
#> [1] TRUE
is_fb_cov(K)
#> [1] FALSE
```
