# Refit a model on new data, with added fixed columns

Refits the model of `fit` on `data`, with the same random structure,
family, priors, engine and engine settings, known variance or
covariance, observation weights and offsets. A package built on
flexyBayes uses it to refit a user's model on a fold of the data, with
candidate columns added to the fixed part, or with brms prior rows and
Stan code of its own.

## Usage

``` r
fb_refit(
  fit,
  data = NULL,
  add_fixed = NULL,
  extra_prior = NULL,
  stanvars = NULL,
  formula_hook = NULL,
  readout = c("full", "varcomp", "fitted", "draws"),
  fixed = NULL,
  ...
)
```

## Arguments

- fit:

  A flexyBayes fit.

- data:

  A data frame with the columns the model uses and any added ones.
  `NULL`, the default, refits on the data the fit carries.

- add_fixed:

  Fixed-effect terms to append to the fixed formula: a one-sided formula
  such as `~ z1 + z2`, or a character vector of term labels. `NULL` adds
  none.

- extra_prior:

  A `brmsprior` from `brms::set_prior()` with rows to add to the prior,
  or `NULL`.

- stanvars:

  A `stanvars` object from `brms::stanvar()`, or `NULL`.

- formula_hook:

  A function of `formula` and `prior` returning
  `list(formula = , prior = )`, or `NULL`. See "Extending the brms
  model".

- readout:

  `"full"`, the default, or `"varcomp"`, `"fitted"` or `"draws"`. See
  "Readout".

- fixed:

  A two-sided formula that replaces the fit's fixed formula, for a refit
  whose fixed part differs by more than added terms. `NULL`, the
  default, keeps the fit's.

- ...:

  Named arguments of
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  that change how this refit runs, for example `seed`, `backend`,
  `n_samples`, `control`, `na_action`, `verbose`, `return_code`,
  `weights` or `known_variance`. On the brms route of "Extending the
  brms model", arguments that
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  does not have go to `brms::brm()`. Arguments that state the model
  (`random`, `residual`, `family`, `link`, `known_matrices`, `prior`,
  `prior_fixed_sd`, `prior_vc_sd`) are refused,
  [`update()`](https://rdrr.io/r/stats/update.html) refits a changed
  model.

## Value

A fitted flexyBayes object of the engine the refit ran on, or, with
`return_code = TRUE`, the code that engine would run.

## What the refit keeps

- The fixed formula is the one `fit` was called with, or `fixed`, with
  the terms of `add_fixed` appended.

- A prior written with
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md),
  `prior_fixed_sd` or `prior_vc_sd` is passed on as written. The package
  default prior is built again over the refitted model, as
  [`update()`](https://rdrr.io/r/stats/update.html) does, so on the same
  data every term keeps the prior it had and an added term gets the
  default.

- The engine is the one `fit` was fitted with, `"inla"` or `"brms"`,
  even when the original call asked for `"auto"`. Its recorded sampler
  settings, `seed`, `control`, `na_action`, `fixed_alias` and
  `aggregate` are passed on.

- Weights and a known variance belong to rows. A refit on the fit's own
  data keeps them. A known variance is read from the column of that name
  in the fit's data, so rows of `fit$data` carry it. For any other
  `data`, pass `weights = ` or `known_variance = ` with one value per
  row. without them a weighted fit is refused, as is a fit with a known
  variance whose `data` lacks the column.

## Extending the brms model

`extra_prior`, `stanvars` and `formula_hook` are brms objects, and the
refit then runs on brms: `fit` must be a brms fit, or `backend = "brms"`
must be given. The refit still starts from the brms formula and prior
rows the package builds for the model.

`formula_hook` is called as `formula_hook(formula, prior)` with that
formula and the `brmsprior` of those rows, and returns
`list(formula = , prior = )`, the formula and prior rows to fit. A
non-linear formula usually moves the package's rows onto one of its
parts, which is why the hook receives them. `extra_prior` rows are then
merged in, replacing any row with the same `class`, `coef`, `group`,
`resp`, `dpar` and `nlpar`, and `stanvars` is passed to `brms::brm()`.
The fitted object reads under the linear model's names.

## Readout

`readout` names what the caller reads from an INLA refit. `"varcomp"`
(variance components), `"fitted"` (fitted values) and `"draws"`
(posterior draws, for example on held-out rows) ask INLA for that
readout only: no DIC, WAIC or marginal densities, posterior
configurations only for `"draws"`, and no second fit to test an AR1
field. The posterior mode, and the readout taken from it, are those of
the full refit. A brms refit is always returned in full.

## See also

[`update.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/update.flexybayes.md)
to refit with any argument changed.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(site = factor(rep(1:12, each = 10)), x = rnorm(120),
                  z = rnorm(120))
  d$y <- rnorm(12)[d$site] + 0.5 * d$x + rnorm(120)
  fit <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                    verbose = FALSE)
  # The same model with z added, on the first nine sites.
  fold <- fb_refit(fit, data = d[d$site %in% 1:9, ], add_fixed = ~z,
                   verbose = FALSE)
  coef(fold)
}
# }
```
