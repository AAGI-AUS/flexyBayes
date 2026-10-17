# Breeder summary of a factor-analytic multi-environment-trial fit

For a `fa(env, k):gen` factor-analytic G x E fit, summarise the
quantities a plant breeder acts on: each genotype's overall performance
(the across-environment mean of its realised effects) and stability (the
across-environment spread), the genotype-by-environment BLUPs, and the
environment genetic-correlation matrix (the crossover structure). The
realised effects are identified, invariant to the rotation and sign
ambiguity of the raw loadings, so their posterior summaries would be
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
  (environment) factors. Unused (see Lifecycle).

## Value

Does not return: raises the classed `met_summary_not_available` refusal.

## Lifecycle

No active engine fits an `fa(env, k):gen` term: both refuse it with
`factor_analytic_not_repeatable`. This function therefore always
abstains (`met_summary_not_available`). What an active engine does
report for a multi-environment trial is the variance components, through
[`summary()`](https://rdrr.io/r/base/summary.html), and on brms the
genotype-by-environment covariance of a
[`diag()`](https://rdrr.io/r/base/diag.html) or `us()` term through
[`brms::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html).

## See also

[`vpredict()`](https://aagi-aus.github.io/flexyBayes/reference/vpredict.md)
for posterior functionals of the variance components, including genetic
correlations.
