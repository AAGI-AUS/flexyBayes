# Dirichlet family object

Constructs the family descriptor for the Dirichlet distribution, the
natural model for compositional (simplex) data. The object mirrors the
shape of the other flexyBayes family descriptors – a named list carrying
the family name, the parameter description, and the support – while
signalling that the Dirichlet is fitted through the dedicated
[`fb_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_dirichlet.md)
entry point rather than the GLM-link emit path.

## Usage

``` r
fb_family_dirichlet()
```

## Value

An object of class `c("fb_family_dirichlet", "fb_family")`: a list with
`family` (the canonical string `"dirichlet"`), `parameters`
(`"concentration (alpha), one per simplex component"`), `support`
(`"K-simplex"`), and `fitter` (`"fb_dirichlet"`).

## Details

The Dirichlet on the \\K\\-simplex is parameterised by a vector of \\K\\
positive concentration parameters \\\alpha_1, \ldots, \alpha_K\\. Their
sum controls the concentration (large sum gives compositions tightly
clustered around the mean), and the normalised vector \\\alpha / \sum_k
\alpha_k\\ is the mean composition.

## See also

[`fb_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_dirichlet.md)

## Examples

``` r
fam <- fb_family_dirichlet()
fam$family
#> [1] "dirichlet"
```
