# Example Multi-Environment Trial (MET) Dataset

A synthetic multi-environment trial dataset suitable for demonstrating
mixed model formulations. Contains yield data for 10 genotypes evaluated
across 4 environments with 3 replicates per environment, plus spatial
layout, genomic markers, and multiple response types.

## Usage

``` r
met_example
```

## Format

A list with the following components:

- `dat`:

  A data frame with 120 rows and the following columns:

  `geno`

  :   Factor: genotype identifier (G1–G10)

  `env`

  :   Factor: environment identifier (E1–E4)

  `rep`

  :   Factor: replicate (R1–R3)

  `block`

  :   Factor: incomplete block within rep (B1–B6)

  `row`

  :   Factor: row position in field layout

  `col`

  :   Factor: column position in field layout

  `yield`

  :   Numeric: grain yield (continuous, primary response)

  `x_cov`

  :   Numeric: a continuous covariate

  `bin_y`

  :   Integer: binary response (0/1)

  `count_y`

  :   Integer: count response (Poisson-like)

- `G_mat`:

  10x10 genomic relationship matrix (G = MM'/2pq)

- `A_mat`:

  10x10 pedigree-based relationship matrix

- `M_geno`:

  10x50 marker matrix (50 SNP markers for 10 genotypes)

- `n_geno`:

  Integer: number of genotypes (10)

- `n_env`:

  Integer: number of environments (4)

## Source

Generated synthetically for package demonstration purposes.

## Details

The dataset is generated with known variance components:

- Genotype variance: sigma_g^2 = 4

- Environment variance: sigma_e^2 = 9 (varies by environment)

- GxE variance: sigma_ge^2 = 2

- Residual variance: sigma^2 = 1

The genomic relationship matrix is computed from 50 simulated SNP
markers using the VanRaden method 1. The pedigree matrix uses a simple
half-sib structure.

## Examples

``` r
data(met_example)
str(met_example$dat)
#> 'data.frame':    120 obs. of  10 variables:
#>  $ geno   : Factor w/ 10 levels "G1","G10","G2",..: 1 3 4 5 6 7 8 9 10 2 ...
#>  $ env    : Factor w/ 4 levels "E1","E2","E3",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ rep    : Factor w/ 3 levels "R1","R2","R3": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ block  : Factor w/ 6 levels "B1","B2","B3",..: 1 2 3 4 5 6 1 2 3 4 ...
#>  $ row    : Factor w/ 10 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
#>  $ col    : Factor w/ 12 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
#>  $ yield  : num  51.9 45 49.4 48.3 47.8 ...
#>  $ x_cov  : num  -0.23 0.837 -1.745 1.689 0.865 ...
#>  $ bin_y  : int  0 1 0 0 0 0 1 1 0 0 ...
#>  $ count_y: int  10 7 5 9 9 9 4 6 5 8 ...
#>  - attr(*, "out.attrs")=List of 2
#>   ..$ dim     : Named int [1:3] 10 4 3
#>   .. ..- attr(*, "names")= chr [1:3] "geno" "env" "rep"
#>   ..$ dimnames:List of 3
#>   .. ..$ geno: chr [1:10] "geno=G1" "geno=G2" "geno=G3" "geno=G4" ...
#>   .. ..$ env : chr [1:4] "env=E1" "env=E2" "env=E3" "env=E4"
#>   .. ..$ rep : chr [1:3] "rep=R1" "rep=R2" "rep=R3"
head(met_example$dat)
#>   geno env rep block row col    yield      x_cov bin_y count_y
#> 1   G1  E1  R1    B1   1   1 51.90921 -0.2297781     0      10
#> 2   G2  E1  R1    B2   2   2 44.97827  0.8366191     1       7
#> 3   G3  E1  R1    B3   3   3 49.44719 -1.7450559     0       5
#> 4   G4  E1  R1    B4   4   4 48.26289  1.6894589     0       9
#> 5   G5  E1  R1    B5   5   5 47.77688  0.8647780     0       9
#> 6   G6  E1  R1    B6   6   6 49.89760 -0.1507760     0       9

# Simple model (small budget for example purposes)
if (FALSE) { # \dontrun{
# live fit -- needs a backend (greta Python/TF, INLA, or brms/Stan)
fit <- flexybayes(
  fixed  = yield ~ env,
  random = ~ geno,
  data   = met_example$dat,
  n_samples = 100, warmup = 100, chains = 1, verbose = FALSE
)
summary(fit)
} # }
```
