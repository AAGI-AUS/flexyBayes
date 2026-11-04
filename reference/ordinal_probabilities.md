# Category probabilities of an ordinal fit

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns, for
each row the model was fitted to, the posterior mean probability of each
category, at the fitted values of the random effects.
[`predict()`](https://rdrr.io/r/stats/predict.html) returns the same
probabilities for the rows of `newdata`, with the random effects set to
zero, or with `marginal = TRUE` integrated over their fitted variances,
which is the probability for a new unit in a new block.

## Usage

``` r
# S3 method for class 'flexybayes_ordinal'
fitted(object, n_draws = 1000L, seed = NULL, summary = TRUE, ...)

# S3 method for class 'flexybayes_ordinal'
predict(
  object,
  newdata = NULL,
  marginal = FALSE,
  n_draws = 1000L,
  seed = NULL,
  summary = TRUE,
  ...
)
```

## Arguments

- object:

  An ordinal fit, from
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  with `family = "cumulative"` or `"cratio"`.

- n_draws:

  The number of posterior draws. Default `1000`. On brms it is at most
  the number the sampler stored.

- seed:

  An integer seed for the draws, or `NULL` (the default) to use the
  session's random number stream.

- summary:

  `TRUE` (the default) returns the posterior mean. `FALSE` returns every
  draw.

- ...:

  Ignored by [`fitted()`](https://rdrr.io/r/stats/fitted.values.html).
  In [`predict()`](https://rdrr.io/r/stats/predict.html), an argument of
  the other [`predict()`](https://rdrr.io/r/stats/predict.html) methods
  that has no meaning here is refused by name.

- newdata:

  A data frame holding the fixed-effect variables, or `NULL` (the
  default) for the rows the model was fitted to.

- marginal:

  `FALSE` (the default) sets every random effect to zero. `TRUE`
  integrates over them.

## Value

A matrix with one row per row of the data and one column per category,
named by the category labels. With `summary = FALSE`, an array, draws by
rows by categories.

## Details

Each probability is averaged over `n_draws` joint posterior draws of the
thresholds, the coefficients and the variances, so it carries the
uncertainty of all three. With `marginal = TRUE` the random terms of a
new unit enter as one normal deviation whose variance is the sum of the
fitted variances (random intercepts, their interactions and an AR1
field, whose marginal variance is its fitted variance), and the
probability is integrated over it by 40-node Gauss-Hermite quadrature in
each draw. On brms,
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) reads
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html).

The probability of a set of categories is the sum of their columns, for
example `rowSums(p[, c("Z5", "Z6")])` for a stage of Z5 or later.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(
    trt = factor(rep(c("a", "b"), each = 60)),
    block = factor(rep(1:6, times = 20))
  )
  eta <- 0.8 * (d$trt == "b") + rnorm(6, 0, 0.3)[d$block]
  d$score <- cut(eta + rlogis(120), c(-Inf, -0.5, 0.5, 1.5, Inf),
    labels = c("low", "mid", "high", "top"), ordered_result = TRUE)
  fit <- flexybayes(score ~ trt, random = ~block, data = d,
    family = "cumulative", backend = "inla", verbose = FALSE)
  head(fitted(fit, seed = 1L))
  predict(fit, newdata = data.frame(trt = c("a", "b")), marginal = TRUE,
    seed = 1L)
}
#> Error in inla.require("sn", stop.on.error = TRUE): Package 'sn' is required to proceed, but is not installed. Please install.
# }
```
