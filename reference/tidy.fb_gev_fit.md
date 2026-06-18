# Tidy a GEV fit

Returns the GEV parameter summary as a `broom`-style `data.frame`, one
row per parameter (`location`, `scale`, `shape`), with the canonical
columns.

## Usage

``` r
# S3 method for class 'fb_gev_fit'
tidy(x, ...)
```

## Arguments

- x:

  An `fb_gev_fit` object from
  [`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md).

- ...:

  Currently unused; present for generic compatibility.

## Value

A `data.frame` with `term`, `estimate`, `std.error`, `conf.low`, and
`conf.high`.

## See also

[`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md)

## Examples

``` r
set.seed(1)
fit <- fb_gev(rgev(200L, 10, 2, 0.1))
tidy(fit)
#>       term    estimate std.error    conf.low  conf.high
#> 1 location 10.19181337 0.1453225  9.90698659 10.4766402
#> 2    scale  1.82586952 0.1085912  1.61303475  2.0387043
#> 3    shape  0.08191751 0.0526013 -0.02117914  0.1850142
```
