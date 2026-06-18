# Genomic-prediction accuracy by cross-validation

Estimate how accurately a genomic BLUP predicts the performance of
genotypes whose phenotypes are held out. For each fold the variance
components are estimated by REML on the training genotypes, the held-out
breeding values are predicted from the genomic relationship matrix
(\\\hat{u}\_V = K\_{VT}(K\_{TT} + \delta I)^{-1}(y_T - X_T\hat\beta)\\),
and the predictions are compared with the observed phenotypes. The
headline metric is *prediction accuracy*, the correlation between
predicted and observed across the held-out genotypes; *bias* is the
slope of observed on predicted (1 is unbiased).

## Usage

``` r
fb_gblup_cv(y, K, X = NULL, folds = 5L, repeats = 1L, seed = NULL, tol = 1e-08)
```

## Arguments

- y:

  Numeric vector of phenotypes, one per genotype, aligned to the rows of
  `K`.

- K:

  An `n x n` genomic (or pedigree) relationship matrix.

- X:

  Optional `n x p` fixed-effect design matrix (an intercept is used when
  `NULL`).

- folds:

  Number of cross-validation folds (default 5).

- repeats:

  Number of independent fold partitions to average over (default 1);
  repeated CV stabilises the estimate.

- seed:

  Optional integer seed for reproducible fold assignment.

- tol:

  PSD / symmetry tolerance for the spectral decomposition.

## Value

An `fb_gblup_cv` object: `accuracy` (pooled predicted-observed
correlation), `accuracy_per_fold`, `bias` (pooled observed-on- predicted
slope), `predicted` (out-of-fold predictions in input order, averaged
over repeats), `h2` (REML heritability on the full data), and metadata.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 80L; m <- 400L
Z <- matrix(rbinom(n * m, 2L, 0.3), n, m)
Zc <- scale(Z, scale = FALSE)
G <- tcrossprod(Zc) / m + diag(n) * 1e-3
u <- as.vector(t(chol(G)) %*% rnorm(n))
y <- u + rnorm(n)
cv <- fb_gblup_cv(y, G, folds = 5L, seed = 1L)
cv$accuracy
#> [1] -0.3181149
# }
```
