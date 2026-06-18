# Triangulate two genome-wide association scans

Compare two GWAS results – typically a flexyBayes
[`fb_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gwas.md)
scan against a field-standard scan (GEMMA / rrBLUP, supplied through the
koine oracle) – by the agreement that matters for a scan: do the same
loci come up? Reports the Jaccard overlap of the genome-wide-significant
marker sets, the overlap among the top markers, the correlation of the
marker effects on the markers in common, and each scan's genomic-control
inflation factor.

## Usage

``` r
triangulate_gwas(a, b, alpha = 0.05, top_k = 10L)
```

## Arguments

- a, b:

  An `fb_gwas` object, or a GWAS lens
  `list(results = <data frame with marker, p_value, effect>, lambda_gc)`.

- alpha:

  Genome-wide significance threshold on the Bonferroni- adjusted p-value
  for the hit-set comparison (default 0.05).

- top_k:

  Size of the top-marker overlap comparison (default 10).

## Value

A `triangulate_gwas_result`: `jaccard` (significant-set overlap),
`n_sig_a` / `n_sig_b` / `n_sig_common`, `top_k_overlap`,
`effect_correlation` (on common markers), and `lambda_gc_a` /
`lambda_gc_b`.

## See also

[`fb_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gwas.md),
[`triangulate_genomic()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_genomic.md).

## Examples

``` r
mk <- function(sig) data.frame(
  marker = paste0("snp", 1:6),
  p_value = ifelse(seq_len(6) %in% sig, 1e-9, 0.5),
  p_bonferroni = ifelse(seq_len(6) %in% sig, 6e-9, 1),
  effect = c(2, 0, 0, -1.5, 0, 0)
)
a <- list(results = mk(c(1, 4)), lambda_gc = 1.01)
b <- list(results = mk(c(1, 4)), lambda_gc = 0.99)
triangulate_gwas(a, b)
#> <triangulate_gwas_result>
#>   significant-set Jaccard: 1  (2 common of 2 / 2)
#>   top-10 overlap: 6/10
#>   effect correlation (common markers): 1
#>   lambda_GC: 1.01 vs 0.99
```
