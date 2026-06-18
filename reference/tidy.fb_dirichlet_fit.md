# Tidy a Dirichlet fit

Returns the concentration-parameter summary as a `broom`-style
`data.frame`, one row per simplex component, with the canonical columns.

## Usage

``` r
# S3 method for class 'fb_dirichlet_fit'
tidy(x, ...)
```

## Arguments

- x:

  An `fb_dirichlet_fit` object from
  [`fb_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_dirichlet.md).

- ...:

  Currently unused; present for generic compatibility.

## Value

A `data.frame` with `term`, `estimate`, `std.error`, `conf.low`, and
`conf.high`.

## See also

[`fb_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_dirichlet.md)

## Examples

``` r
set.seed(1)
fit <- fb_dirichlet(rdirichlet(300L, c(2, 5, 3)))
tidy(fit)
#>    term estimate std.error conf.low conf.high
#> c1   c1 1.886471 0.1117041 1.667535  2.105407
#> c2   c2 4.673664 0.2811266 4.122666  5.224662
#> c3   c3 2.792916 0.1667429 2.466106  3.119727
```
