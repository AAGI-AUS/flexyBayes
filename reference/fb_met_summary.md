# Breeder summary of a factor-analytic multi-environment-trial fit

For a `fa(env, k):gen` factor-analytic G x E fit, summarise the
quantities a plant breeder acts on: each genotype's overall performance
(the across-environment mean of its realised effects) and stability (the
across-environment spread), the genotype-by-environment BLUPs, and the
environment genetic-correlation matrix (the crossover structure). The
realised effects are identified – invariant to the rotation and sign
ambiguity of the raw loadings – so their posterior summaries are
interpretable; judge convergence on these and on the identified
covariance
([`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md))
rather than on the raw loadings.

## Usage

``` r
fb_met_summary(fit, genotype_levels = NULL, environment_levels = NULL)
```

## Arguments

- fit:

  A `flexybayes` greta fit with a `fa()` factor-analytic G x E term.

- genotype_levels, environment_levels:

  Optional character labels for the inner (genotype) and outer
  (environment) factors; default to positional labels.

## Value

An `fb_met_summary` object (one entry is built per `fa()` term, the
function returns the first / named): `op` (data frame of overall
performance per genotype with credible interval), `stability` (data
frame of across-environment spread per genotype), `gxe_blup` (the
posterior-mean genotype-by-environment effect matrix), `env_cor` (the
environment genetic-correlation matrix), `loadings` (posterior-mean
factor loadings), and metadata.

## See also

[`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md)
for the identified environment covariance and its convergence
diagnostic.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- flexybayes(
  yield ~ env, random = ~ fa(env, 2):gen, data = met, backend = "greta"
)
ms <- fb_met_summary(fit)
head(ms$op[order(-ms$op$mean), ]) # best genotypes on average
ms$env_cor # environment crossover structure
} # }
```
