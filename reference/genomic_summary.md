# Genomic summary of a fitted relationship model

Extract the breeder-facing genomic quantities: narrow-sense heritability
\\h^2\\, genomic estimated breeding values (GEBVs) with posterior
reliability, and the genetic / residual variances, from a fitted `vm()`
(genomic / GBLUP) or `ped()` (pedigree) model. The quantities are read
from the posterior draws engine-agnostically: an INLA or brms GBLUP fit
returns the same summary object, so the two engines' genomic summaries
can be read side by side.

## Usage

``` r
genomic_summary(fit, term = NULL)
```

## Arguments

- fit:

  A fitted `flexybayes` object carrying at least one `vm()` or `ped()`
  relationship term.

- term:

  Optional grouping-factor name selecting which relationship term to
  summarise when the model has more than one. Defaults to the first.

## Value

An `fb_genomic_summary` object: `heritability`, `genetic_variance`,
`residual_variance` (each a posterior summary), `gebv` (data frame of
breeding values with reliability), and metadata.

## Details

The heritability is computed per draw as \\h^2 = \sigma_g^2 /
(\sigma_g^2 + \sigma_e^2)\\ on the genotype-mean basis. The kinship
scaling convention is the analyst's (state it when reporting).
Reliability is \\1 - \mathrm{PEV}\_i / \sigma_g^2\\ from the posterior
variance of each breeding value. GEBVs are available on both active
backends (brms and INLA) natively.

## See also

[`vpredict()`](https://aagi-aus.github.io/flexyBayes/reference/vpredict.md)
for posterior functionals of the variance components, and
[`ranef()`](https://aagi-aus.github.io/flexyBayes/reference/ranef.md)
for the group-level effects themselves.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- flexybayes(
  yield ~ 1, random = ~ vm(geno, Gmat), data = met,
  known_matrices = list(Gmat = G), backend = "brms"
)
gs <- genomic_summary(fit)
gs$heritability
head(gs$gebv)
} # }
```
