# Genomic summary of a fitted relationship model

Extract the breeder-facing genomic quantities – narrow-sense
heritability \\h^2\\, genomic estimated breeding values (GEBVs) with
posterior reliability, and the genetic / residual variances – from a
fitted `vm()` (genomic / GBLUP) or `ped()` (pedigree) model. The
quantities are read from the posterior draws engine-agnostically: a
greta, INLA, or brms GBLUP fit returns the same summary object, so a
multi-backend genomic analysis is directly triangulatable.

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
(\sigma_g^2 + \sigma_e^2)\\ on the genotype-mean basis; the kinship
scaling convention is the analyst's (state it when reporting).
Reliability is \\1 - \mathrm{PEV}\_i / \sigma_g^2\\ from the posterior
variance of each breeding value. GEBVs are available on the brms and
INLA backends natively and on the greta backend (the breeding-value
vector is monitored).

## See also

[`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md)
for factor-analytic MET covariance,
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
for cross-engine agreement.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- flexybayes(
  yield ~ 1, random = ~ vm(geno, Gmat), data = met,
  known_matrices = list(Gmat = G), backend = "greta"
)
gs <- genomic_summary(fit)
gs$heritability
head(gs$gebv)
} # }
```
