# The marginal covariance of the response of a fitted model

Returns the covariance of the response over the rows of a fit's data,
given its fitted variance components, \$\$V = R + \sum_k \sigma_k^2 Z_k
G_k Z_k^\top + K,\$\$ as a structured object with one block per variance
component. [`as.matrix()`](https://rdrr.io/r/base/matrix.html) assembles
it as a dense matrix.

## Usage

``` r
fb_marginal_cov(fit, rows = NULL, variances = TRUE)

# S3 method for class 'fb_marginal_cov'
as.matrix(x, rows = NULL, cols = rows, ...)
```

## Arguments

- fit:

  A fitted `flexybayes` model, or an `fb_marginal_cov` object to take a
  subset of or to evaluate at other variances.

- rows:

  Integer positions of the rows to cover: rows of the fit's data, or of
  `fit` when it is an `fb_marginal_cov` object. `NULL` (the default)
  covers every row.

- variances:

  `TRUE` (the default) to use the fit's own variance components, `FALSE`
  for the structure alone, with every variance `NA`. or a list with
  `residual` (the residual variance, or one per section named by level),
  `terms` (a variance per block, named by block) and, for an AR1 x AR1
  field, `rho` (a list of `c(row = , col = )` named by block), all on
  the variance scale.

- x:

  An `fb_marginal_cov` object.

- cols:

  Integer positions of the columns, for a cross-covariance between two
  sets of rows. Defaults to `rows`.

- ...:

  Unused.

## Value

An object of class `fb_marginal_cov`, a list with

- `n`, `rows`:

  The number of rows covered, and their positions in the fit's data.

- `residual`:

  A list: `variance` (one value, or one per section named by level),
  `by` (the sectioning column, or `NULL`), `level` (each row's section,
  or `NULL`), `weights` (each row's observation weight, or `NULL`) and
  `known_variance` (each row's known variance, or `NULL`).

- `blocks`:

  A named list, one entry per variance component: `name`, `label`,
  `component` (the names the variance may carry in
  `summary(fit)$varcomp`), `variance`, `index`, `n_levels`, `value`,
  `key`, `columns` and `structure` as described above, and for an AR1 x
  AR1 field `rho` and `kronecker`, the factors `n_at` (trials), `row`
  and `col` (AR1 correlation matrices) with \\G = I\_{n\_{at}} \otimes
  R\_{row} \otimes R\_{col}\\.

- `known_covariance`:

  The known covariance matrix over the rows, whose diagonal is
  `residual$known_variance`, or `NULL`.

- `groups`:

  Whether V is block-diagonal over a grouping column: `applicable`,
  `column`, `blocks` (row positions per group, the grouping with the
  most groups that every block nests within) and `reason`.

- `varcomp`:

  The variance-component table the values were read from, or `NULL`.

## Details

Each variance enters at its posterior mean, \\E\[\sigma^2\] =
\mathrm{estimate}^2 + \mathrm{std.error}^2\\, read from the fit's
variance-component table, `summary(fit)$varcomp`, whose `estimate` and
`std.error` are the posterior mean and standard deviation of each
standard deviation. For a structure without a correlation parameter this
is the posterior mean of V itself. The two correlations of an AR1 x AR1
field enter at their posterior means (`estimate`). On brms, the rows of
a sectioned residual report the posterior median of each standard
deviation, so there the residual variances are median^2 + std.error^2.

The structures built are random intercepts (`(1 | g)`, `ide(g)`),
two-way interactions (`a:b`), uncorrelated random slopes (`(x || g)`),
`diag(f):g`, an `ar1(row):ar1(col)` field on one trial or
`id(trial):ar1(row):ar1(col)` on each trial, and `vm()` or `ped()` known
covariances, over an independent residual or one variance per level of a
factor (`dsum(~ units | f)`). Observation weights divide the residual
variance (\\\sigma^2 / w_i\\), a known variance per row is added to it,
and a known covariance matrix is added whole.

A model with any other structure is refused with the reason code
`marginal_cov_not_buildable`, naming the terms.

## See also

[`fb_variance_draws()`](https://aagi-aus.github.io/flexyBayes/reference/fb_variance_draws.md)
for posterior draws of the variance components.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(site = factor(rep(1:6, each = 5)))
  d$y <- rnorm(6)[d$site] + rnorm(30)
  fit <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
                    verbose = FALSE)
  v <- fb_marginal_cov(fit, rows = 1:10)
  v$blocks$site$variance
  as.matrix(v)
}
#>           [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#>  [1,] 2.091745 1.178568 1.178568 1.178568 1.178568 0.000000 0.000000 0.000000
#>  [2,] 1.178568 2.091745 1.178568 1.178568 1.178568 0.000000 0.000000 0.000000
#>  [3,] 1.178568 1.178568 2.091745 1.178568 1.178568 0.000000 0.000000 0.000000
#>  [4,] 1.178568 1.178568 1.178568 2.091745 1.178568 0.000000 0.000000 0.000000
#>  [5,] 1.178568 1.178568 1.178568 1.178568 2.091745 0.000000 0.000000 0.000000
#>  [6,] 0.000000 0.000000 0.000000 0.000000 0.000000 2.091745 1.178568 1.178568
#>  [7,] 0.000000 0.000000 0.000000 0.000000 0.000000 1.178568 2.091745 1.178568
#>  [8,] 0.000000 0.000000 0.000000 0.000000 0.000000 1.178568 1.178568 2.091745
#>  [9,] 0.000000 0.000000 0.000000 0.000000 0.000000 1.178568 1.178568 1.178568
#> [10,] 0.000000 0.000000 0.000000 0.000000 0.000000 1.178568 1.178568 1.178568
#>           [,9]    [,10]
#>  [1,] 0.000000 0.000000
#>  [2,] 0.000000 0.000000
#>  [3,] 0.000000 0.000000
#>  [4,] 0.000000 0.000000
#>  [5,] 0.000000 0.000000
#>  [6,] 1.178568 1.178568
#>  [7,] 1.178568 1.178568
#>  [8,] 1.178568 1.178568
#>  [9,] 2.091745 1.178568
#> [10,] 1.178568 2.091745
# }
```
