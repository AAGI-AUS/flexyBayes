# Genome-wide association scan (EMMAX / P3D)

Test each marker for association with a phenotype while correcting for
polygenic background and population / family structure via a genomic
relationship matrix. `fb_gwas()` fits the polygenic null mixed model
once to estimate the variance components by REML, then — holding those
components fixed (the P3D approximation of Kang et al. 2010) — tests
every marker by exact generalised least squares in the eigenbasis of the
relationship matrix. The single eigendecomposition is the shared
spectral primitive; each marker is then an `O(n)` weighted-least-squares
score test, so a whole-genome scan is feasible without a per-marker
model fit.

## Usage

``` r
fb_gwas(formula, data, markers, K = NULL, marker_map = NULL, tol = 1e-08)
```

## Arguments

- formula:

  A two-sided formula `y ~ covariates` for the response and fixed-effect
  covariates (an intercept is included unless removed with `- 1`). Use
  `y ~ 1` for an intercept-only background.

- data:

  A data frame with the response and covariates; one row per individual
  (the association unit).

- markers:

  A numeric `n x m` matrix of marker genotypes (allele dosages), one row
  per individual aligned to `data`, one column per marker. Column names
  are used as marker identifiers.

- K:

  An `n x n` genomic (or pedigree) relationship matrix. When `NULL`
  (default) it is built from the centred markers as \\G = Z_c Z_c^\top /
  m\\, the average-allele-frequency genomic relationship (Astle &
  Balding 2009). This normalises by the marker count `m`, unlike
  VanRaden's (2008) first method, which divides by \\2\sum_j p_j(1 -
  p_j)\\; the two coincide only when every marker has minor-allele
  frequency 0.5 and otherwise differ by a constant scale that the REML
  variance ratio absorbs, so the scan statistics are unaffected. Note
  that building `K` from the same markers being tested causes *proximal
  contamination*: a marker's own signal enters the polygenic term and
  deflates its test statistic (conservative, but a power loss). For a
  well-calibrated scan supply `K` from an independent background panel
  or a leave-one-chromosome-out construction.

- marker_map:

  Optional data frame with columns `marker`, `chr`, `pos` to annotate
  the results and order the Manhattan plot.

- tol:

  PSD / symmetry tolerance passed to the spectral primitive.

## Value

An `fb_gwas` object: a list with `results` (a data frame of `marker`,
`effect`, `se`, `statistic` (the chi-square), `p_value`, `p_bonferroni`,
`q_value` (Benjamini-Hochberg FDR), plus any map columns), `lambda_gc`
(the genomic-control inflation factor), `var_g` / `var_e` / `h2` (the
null REML components), and metadata.

## Details

The backends (greta / INLA / brms) are not used by the scan itself — it
is a deterministic frequentist fast path. They enter only when a handful
of significant loci are re-fit as full Bayesian models for credible
effect sizes, which is affordable at that reduced scale.

## References

Kang, H. M., et al. (2010). Variance component model to account for
sample structure in genome-wide association studies. *Nature Genetics*,
42(4), 348-354. Astle, W., & Balding, D. J. (2009). Population structure
and cryptic relatedness in genetic association studies. *Statistical
Science*, 24(4), 451-471. VanRaden, P. M. (2008). Efficient methods to
compute genomic predictions. *Journal of Dairy Science*, 91(11),
4414-4423.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 100L
m <- 200L
M <- matrix(rbinom(n * m, 2L, 0.3), n, m)
colnames(M) <- paste0("snp", seq_len(m))
y <- 2 * scale(M[, 50L]) + rnorm(n)
dat <- data.frame(y = y)
scan <- fb_gwas(y ~ 1, data = dat, markers = M)
head(scan$results[order(scan$results$p_value), ])
#>     marker     effect        se statistic      p_value p_bonferroni
#> 50   snp50  2.7711252 0.3074583 81.234363 2.004756e-19 4.009512e-17
#> 98   snp98 -0.5708480 0.3474396  2.699495 1.003801e-01 1.000000e+00
#> 199 snp199 -0.5404380 0.3300467  2.681273 1.015345e-01 1.000000e+00
#> 159 snp159  0.5256615 0.3279511  2.569178 1.089643e-01 1.000000e+00
#> 129 snp129 -0.5322155 0.3374802  2.487015 1.147892e-01 1.000000e+00
#> 179 snp179  0.5263342 0.3437879  2.343916 1.257730e-01 1.000000e+00
#>          q_value
#> 50  4.009512e-17
#> 98  9.851188e-01
#> 199 9.851188e-01
#> 159 9.851188e-01
#> 129 9.851188e-01
#> 179 9.851188e-01
scan$lambda_gc
#> [1] 0.5870781
# }
```
