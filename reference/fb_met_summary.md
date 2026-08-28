# Breeder summary of a factor-analytic multi-environment-trial fit

For a `fa(env, k):gen` factor-analytic G x E fit, summarise the
quantities a plant breeder acts on: each genotype's overall performance
(the across-environment mean of its realised effects) and stability (the
across-environment spread), the genotype-by-environment BLUPs, and the
environment genetic-correlation matrix (the crossover structure). The
realised effects are identified – invariant to the rotation and sign
ambiguity of the raw loadings – so their posterior summaries would be
interpretable.

## Usage

``` r
fb_met_summary(fit, genotype_levels = NULL, environment_levels = NULL)
```

## Arguments

- fit:

  A flexybayes fit.

- genotype_levels, environment_levels:

  Optional character labels for the inner (genotype) and outer
  (environment) factors; unused (see Lifecycle).

## Value

Does not return: raises the classed `met_summary_not_available` refusal.

## Lifecycle

No active engine emits an `fa(env, k):gen` term – both INLA and brms
refuse a factor-analytic structured-covariance term before a fit object
exists – so this function always abstains (`met_summary_not_available`).
What an active engine does report for a multi-environment trial is the
variance components, through
[`summary()`](https://rdrr.io/r/base/summary.html), and on brms the
genotype-by-environment covariance of a
[`diag()`](https://rdrr.io/r/base/diag.html) or `us()` term through
`brms::VarCorr()`.

## See also

[`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md)
for the identified environment covariance and its convergence
diagnostic.
