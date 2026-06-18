# Triangulate genomic model outputs across engines or against a field lens

Compare two genomic analyses – heritability, variance components, and
genomic estimated breeding values – on a common footing. Each argument
is either a flexyBayes GBLUP / pedigree fit (greta, INLA, or brms) or a
generic *genomic lens* (`list(h2, var_g, var_e, gebv, label)`), the form
a field-standard oracle such as sommer's REML supplies. Cross-engine use
checks that the Bayesian backends agree; the lens form lets the koine
fourth opinion (REML / established tools) cross-check the Bayesian
answer – the orchestra's signature value in a field with decades of
established methods.

## Usage

``` r
triangulate_genomic(a, b, term = NULL, data_independence = NA)
```

## Arguments

- a, b:

  A flexyBayes GBLUP / pedigree fit, or a genomic-lens list. A lens
  entry for `h2` / `var_g` / `var_e` may be a numeric vector of draws or
  `list(estimate, se)` (a REML point); `gebv` is a named numeric vector
  keyed by genotype.

- term:

  Optional grouping-factor name selecting the relationship term when a
  fit has more than one.

- data_independence:

  `TRUE` if the two lenses used independently sourced data (suppresses
  the caveat); `NA` (default) or `FALSE` attaches it.

## Value

A `triangulate_genomic_result`: `components` (a data frame comparing
heritability and the variance components – value and interval per lens,
the difference, and an interval-overlap flag), `gebv` (the Pearson /
Spearman correlation of the matched breeding values and the number in
common), the lens labels, and the caveat.

## Details

Like
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md),
this measures inter-lens *agreement*, not correspondence: two lenses fit
to the same data share any fabricated upstream data fact, so the result
carries a shared-upstream caveat unless `data_independence = TRUE`.

## See also

[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md),
[`genomic_summary()`](https://aagi-aus.github.io/flexyBayes/reference/genomic_summary.md),
[`triangulate_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_gwas.md).

## Examples

``` r
# A Bayesian posterior (draws) against a field-standard REML point lens.
set.seed(1)
bv <- c(2, 1, 0, -1, -2)
bayes <- list(
  h2 = rnorm(500, 0.5, 0.05),
  gebv = stats::setNames(bv + rnorm(5, 0, 0.2), paste0("g", 1:5))
)
reml <- list(
  h2 = list(estimate = 0.48, se = 0.04),
  gebv = stats::setNames(bv + rnorm(5, 0, 0.2), paste0("g", 1:5))
)
triangulate_genomic(bayes, reml)
#> <triangulate_genomic_result>  lens vs lens
#>   heritability      : 0.501 vs 0.48  (diff 0.021, intervals overlap)
#>   genetic_variance  : (not supplied)
#>   residual_variance : (not supplied)
#>   GEBVs (5 common): r = 0.986 (Pearson), 1 (Spearman)
#>   ! triangulate_genomic measures inter-lens agreement, not correspondence: data independence was not declared. A fabricated upstream data fact (e.g. a mis-scaled relationship matrix) is shared by both and would not be caught by their agreement.
```
