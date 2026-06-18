# Generalised extreme value (GEV) family object

Constructs the family descriptor for the generalised extreme value
distribution, the limiting law of block maxima. The object mirrors the
shape of a base [`stats::family()`](https://rdrr.io/r/stats/family.html)
descriptor – a named list carrying the family name, the parameter names,
and the natural support – so it sits alongside the other flexyBayes
families, while signalling that GEV is fitted through the dedicated
[`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md)
entry point rather than the GLM-link emit path.

## Usage

``` r
fb_family_gev()
```

## Value

An object of class `c("fb_family_gev", "fb_family")`: a list with
`family` (the canonical string `"gen_extreme_value"`), `parameters` (the
character vector `c("location", "scale", "shape")`), `n_par` (the
integer `3`), `link` (`"identity"` on the location), and `fitter`
(`"fb_gev"`).

## Details

The GEV is parameterised by a location \\\mu\\, a positive scale
\\\sigma\\, and a shape \\\xi\\ (the extreme value index). The shape
governs the tail: \\\xi \> 0\\ gives the heavy-tailed Frechet type,
\\\xi \< 0\\ the bounded Weibull type, and \\\xi \to 0\\ the Gumbel
limit.

## See also

[`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md)

## Examples

``` r
fam <- fb_family_gev()
fam$parameters
#> [1] "location" "scale"    "shape"   
```
