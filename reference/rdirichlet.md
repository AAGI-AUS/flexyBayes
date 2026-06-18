# Simulate from a Dirichlet distribution

Draws compositional rows from a Dirichlet with the given concentration
vector, via the gamma-normalisation construction.

## Usage

``` r
rdirichlet(n, alpha)
```

## Arguments

- n:

  Integer. The number of compositions (rows) to draw.

- alpha:

  Numeric vector of positive concentration parameters; its length sets
  the number of simplex components \\K\\.

## Value

A numeric `n` by `K` matrix whose rows sum to one.

## See also

[`fb_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_dirichlet.md)

## Examples

``` r
set.seed(1)
X <- rdirichlet(5L, alpha = c(2, 5, 3))
rowSums(X)
#> [1] 1 1 1 1 1
```
