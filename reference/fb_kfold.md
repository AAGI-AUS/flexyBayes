# Exact K-fold cross-validation of an ordinal fit

Refits the model once per fold with that fold's responses held out, and
scores each held-out group by its log predictive density, the group's
rows taken jointly. With `group` naming the plot, a plot is scored by
the probability of all of its plants' categories together, which is how
ordinal models of sub-sampled field trials are compared when
leave-one-out approximations are unreliable.

## Usage

``` r
fb_kfold(fit, folds, group = NULL, n_draws = 1000L, seed = NULL, ...)
```

## Arguments

- fit:

  An ordinal fit, from
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  with `family = "cumulative"` or `"cratio"`.

- folds:

  The fold of each row of the data the fit was given: a vector with one
  label per row and at least two distinct labels.

- group:

  The name of a column of the fit's data whose levels are scored
  jointly, for example the plot, or `NULL` (the default) to score each
  row on its own. The rows of a group must share a fold.

- n_draws:

  The number of posterior draws per fold. Default `1000`.

- seed:

  An integer seed, or `NULL` (the default) to use the session's random
  number stream. Fold `f` uses `seed + f`.

- ...:

  Further arguments of
  [`fb_refit()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refit.md)
  for each fold's refit, for example `control` on brms.

## Value

An object of class `fb_kfold`, a list with `estimates`, a data frame
holding `elpd` and its standard error `se`, `pointwise`, one row per
scored group with `group`, `fold`, `n_rows` and `elpd`, and `folds`,
`engine` and `group`.

## Details

On INLA a fold's responses are set to missing and its rows kept, so
every latent node, an AR1 field's included, stays in the model, and the
held-out linear predictors are drawn from the refit's joint posterior.
On brms the model is refitted to the other folds' rows, and a held-out
group's random effect at a level the refit did not see is drawn from its
fitted variance (brms's `allow_new_levels = TRUE` with
`sample_new_levels = "gaussian"`).

The score of group `g` is
`log(mean over draws s of prod over its rows i of P(y_i | draw s))`, the
`elpd` is the sum over groups, and its standard error is
`sqrt(G) * sd(score)` over the `G` groups. Rows with a missing response
are not scored.

Every category must keep a training observation in every fold, or the
refit would have fewer thresholds than the model it is scored against.

## References

Cao Z, Megirian R, Nguyen M, Sparks AH (2026). Bayesian ordinal analysis
of cereal growth stage data from oat trials in Western Australia.
arXiv:2610.37384.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(
    trt = factor(rep(c("a", "b"), each = 60)),
    plot = factor(rep(1:24, each = 5))
  )
  eta <- 0.8 * (d$trt == "b") + rnorm(24, 0, 0.5)[d$plot]
  d$score <- cut(eta + rlogis(120), c(-Inf, -0.5, 0.5, 1.5, Inf),
    labels = c("low", "mid", "high", "top"), ordered_result = TRUE)
  fit <- flexybayes(score ~ trt, random = ~plot, data = d,
    family = "cumulative", backend = "inla", verbose = FALSE)
  folds <- rep(rep(1:4, 6), each = 5)
  fb_kfold(fit, folds, group = "plot", seed = 1L)
}
#> Error in inla.require("sn", stop.on.error = TRUE): Package 'sn' is required to proceed, but is not installed. Please install.
# }
```
