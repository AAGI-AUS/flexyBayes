# Identified covariance for factor-analytic structured-covariance terms

For each `fa(outer, k)` term in a fit, reconstruct the implied
outer-factor covariance \\G = \Lambda\Lambda^\top +
\mathrm{diag}(\psi)\\ from the posterior draws and summarise it. Unlike
the raw loadings \\\Lambda\\ – which are identified only up to rotation
and sign, so their per-entry Rhat is meaningless – the covariance \\G\\
and the correlation derived from it are rotation- and sign-invariant.

## Usage

``` r
fb_structured_cov(fit)
```

## Arguments

- fit:

  A flexybayes fit.

## Value

An empty list (with a message): no active backend fits an `fa()` term,
so no fit ever carries one to reconstruct. Non-factor- analytic
structured terms (`us`, `ar1`) are reported as not-yet-reconstructed.

## Examples

``` r
# Reports on factor-analytic structured covariance. A fit carrying no
# such term is told so rather than returned an empty object.
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
                    verbose = FALSE)
  flexyBayes:::fb_structured_cov(fit)
}
# }
```
