# Simulate from a generalised extreme value distribution

Draws block-maxima observations from a GEV with the given location,
scale, and shape, by inverting the GEV cumulative distribution function.
The Gumbel limit (\\\xi \to 0\\) is handled exactly.

## Usage

``` r
rgev(n, location, scale, shape = 0)
```

## Arguments

- n:

  Integer. The number of observations to draw.

- location:

  Numeric. The location parameter \\\mu\\.

- scale:

  Numeric. The positive scale parameter \\\sigma\\.

- shape:

  Numeric. The shape parameter \\\xi\\. Defaults to `0` (Gumbel).

## Value

A numeric vector of length `n`.

## See also

[`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md)

## Examples

``` r
set.seed(1)
y <- rgev(100L, location = 10, scale = 2, shape = 0.15)
summary(y)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   7.375   9.758  10.680  11.417  12.939  24.111 
```
