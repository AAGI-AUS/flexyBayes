# Joint posterior draws of the variance components of an INLA fit

Draws the hyperparameters of an INLA fit jointly and returns them on the
scale and under the names of the fit's variance-component table
(`summary(fit)$varcomp`): a standard deviation for the residual
(`sigma`) and for each random term (`sd_<term>`), and a correlation as
it is.

## Usage

``` r
fb_variance_draws(object, n, seed = NULL)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- n:

  The number of draws.

- seed:

  An integer seed, or `NULL` (the default) to use the session's random
  number stream.

## Value

A data frame with `n` rows and one column per variance component.

## Details

The draws come from `INLA::inla.hyperpar.sample()` and are mapped to the
response's original scale. A hyperparameter held fixed, such as the
precision of a known variance, is not drawn. With `seed`, the draws from
a given fit are repeatable on one INLA thread
(`options(flexyBayes.inla_threads = "1:1")`). A refit of the same model
in another R session can reach a different mode (see
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md),
section *INLA threads*).

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(site = factor(rep(1:6, each = 5)))
  d$y <- rnorm(6)[d$site] + rnorm(30)
  fit <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
                    verbose = FALSE)
  head(fb_variance_draws(fit, 100L, seed = 1L))
}
# }
```
