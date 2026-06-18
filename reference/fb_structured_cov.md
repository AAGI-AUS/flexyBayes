# Identified covariance for factor-analytic structured-covariance terms

For each `fa(outer, k)` term in a greta fit, reconstruct the implied
outer-factor covariance \\G = \Lambda\Lambda^\top +
\mathrm{diag}(\psi)\\ from the posterior draws and summarise it. Unlike
the raw loadings \\\Lambda\\ – which are identified only up to rotation
and sign, so their per-entry Rhat is meaningless – the covariance \\G\\
and the correlation derived from it are rotation- and sign-invariant.
They are therefore the identified quantities whose posterior is
interpretable and whose Rhat is a genuine convergence diagnostic.
Consult this in preference to the raw-loading Rhat when judging whether
a factor-analytic fit has converged.

## Usage

``` r
fb_structured_cov(fit)
```

## Arguments

- fit:

  A `flexybayes` fit produced on the greta backend with at least one
  `fa()` structured-covariance term.

## Value

A named list with one entry per factor-analytic term (named by the
term's outer factor). Each entry is a list with: `levels` (the
outer-factor levels labelling the rows/columns), `cov_mean`,
`cov_lower`, `cov_upper` (posterior mean and 95% interval of \\G\\),
`cor_mean` (posterior-mean correlation), `rhat` (entrywise Rhat of
\\G\\, `NA` if fewer than two chains), `max_rhat`, and `k` (the number
of factors). Returns an empty list (with a message) when the fit carries
no factor-analytic term. Non-factor-analytic structured terms (`us`,
`ar1`) are reported as not-yet-reconstructed.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- flexybayes(
  y ~ 1, random = ~ fa(env, 2):id(geno), data = met_data,
  family = "gaussian", backend = "greta"
)
sc <- fb_structured_cov(fit)
sc$env$cov_mean   # identified genetic covariance across environments
sc$env$max_rhat   # convergence of the identified quantity
} # }
```
