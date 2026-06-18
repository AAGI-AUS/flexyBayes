# Specify priors via the PC-canonical hybrid DSL

The flexyBayes prior DSL (domain-specific language) lives on the
standard-deviation scale (never precision / variance) and accepts two
canonical idioms: distributional (`half_normal(scale = 1)`) and
tail-quantile / PC (penalised complexity) (`pc(upper = 1, prob = 0.01)`
meaning `Pr(sigma > 1) = 0.01`). The PC idiom is the cross-engine
interlingua – it survives every backend because it is a probability
statement, not a distributional name.

## Usage

``` r
fb_prior(...)
```

## Arguments

- ...:

  one or more two-sided formulas of the form
  `target ~ distribution(args)`. Examples:

  - `sigma ~ pc(upper = 2, prob = 0.05)`

  - `sd(group = "subject") ~ half_normal(scale = 1)`

  - `b("treatment") ~ student_t(df = 4, scale = 2.5)`

  - `cor(group = "subject") ~ lkj(eta = 2)`

  - `sd(group = "subject") ~ uniform(lower = 0, upper = 5)`

  Supported distribution families: `pc`, `half_normal`, `half_cauchy`,
  `student_t`, `normal`, `exponential`, `lkj`, `cauchy`, `gamma`,
  `uniform`. Note that `uniform()` on a variance component sits outside
  the PC-canonical interlingua, but both backends represent it
  faithfully on the SD scale: the INLA backend via an expression-prior
  on the log-precision, and the greta backend as a bounded
  `greta::uniform()` on each simple random-effect (and residual) SD.
  Structured-covariance terms (`us`, `fa`, `ar1`, `vm`, `ped`) on greta
  fall back to the legacy scale prior.

## Value

an `fb_prior` object (S3, inherits from list) with `$specs` carrying the
parsed `target` / `spec` pairs.

## Details

v0.1 supports the targets and distributions listed above in the file
header. Calls outside the supported set raise a structured error with
the supported list.

## Examples

``` r
p <- fb_prior(
  sigma                 ~ pc(upper = 2, prob = 0.05),
  sd(group = "subject") ~ half_normal(scale = 1),
  b("treatment")        ~ student_t(df = 4, scale = 2.5)
)
p
#> <fb_prior> 3 specifications
#>   sigma ~ pc(upper = 2, prob = 0.05)
#>   sd(group = "subject") ~ half_normal(scale = 1)
#>   b("treatment") ~ student_t(df = 4, scale = 2.5)
```
