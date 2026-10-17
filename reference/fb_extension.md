# Reading a fit, for packages that extend flexyBayes

A package built on flexyBayes reads a fitted model through these
functions,
[`fb_refit()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refit.md),
[`fb_marginal_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_marginal_cov.md),
[`fb_isolated_lapply()`](https://aagi-aus.github.io/flexyBayes/reference/fb_isolated_lapply.md)
and the refusal registry
([`fb_register_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_register_refusals.md),
[`fb_refuse()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refuse.md)),
and never through the package's internal functions. `fb_extension_api()`
gives the version of that set. A package checks it in its `.onLoad()`
and refuses to load against a version it was not written for.

Builds the fixed-effect design matrix of `formula` over `data`, as
flexyBayes builds it for its aliasing check, with the term of each
column and the order in which ASReml 4.2 enters the columns.

## Usage

``` r
fb_extension_api()

fb_model_terms(fit)

fb_response_values(response, data)

fb_offset_values(offsets, data, newdata = FALSE)

fb_model_weights(terms)

fb_known_variance(terms)

fb_engine_settings(engine)

fb_convergence(fit)

fb_inla_posterior_sample(fit, n, ...)

fb_check_refit(fit, caller = "")

fb_fixed_design(formula, data, train_rows = NULL)
```

## Arguments

- fit:

  A flexyBayes fit.

- response:

  The response as written in the model: a column name or an R expression
  in the columns of `data`.

- data:

  A data frame.

- offsets:

  The offset terms of a model as written, for example
  `"offset(log(area))"`, `fb_model_terms(fit)$offsets`.

- newdata:

  `TRUE` when `data` holds new rows for a prediction, so that a missing
  offset column is refused with a message that says so.

- terms:

  A model description from `fb_model_terms()`, or a fit.

- engine:

  `NULL`, an `fb_engine` object, or a named list of sampler settings.

- n:

  The number of posterior samples.

- ...:

  Further arguments to
  [`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html),
  for example `seed`, `num.threads` or `selection`.

- caller:

  Text placed in front of a refusal's message, naming the function that
  asked.

- formula:

  A one-sided formula of fixed terms.

- train_rows:

  Row numbers of the training rows, or `NULL` (default) to learn
  transforms from every row.

## Value

- `fb_extension_api()`: an integer, the version of the set. It is 1.

- `fb_model_terms()`: a list, or `NULL` for an object that carries none.

- `fb_response_values()`: a vector with one value per row of `data`, or
  `NULL` when the response cannot be evaluated there.

- `fb_offset_values()`: a numeric vector with one value per row of
  `data`, or `NULL` when the model has no offset.

- `fb_model_weights()`: the numeric weights, or `NULL` for an unweighted
  model.

- `fb_known_variance()`: `NULL` for a model without a known variance.
  Otherwise a list with `values` (one variance per row of the fitted
  data), `matrix` (the known covariance, or `NULL` for a vector),
  `column` (the data column that carries each row's variance) and
  `inla_index` (the index of the row effect that carries it in an INLA
  formula).

- `fb_engine_settings()`: a named list, empty for `NULL`.

- `fb_convergence()`: a list with `engine`, `ok` (`TRUE`, `FALSE` or
  `NA`), `reasons` and `evidence`.

- `fb_inla_posterior_sample()`: the list
  [`INLA::inla.posterior.sample()`](https://rdrr.io/pkg/INLA/man/posterior.sample.html)
  returns, with the linear predictor and the hyperparameters on the
  scale of the data.

- `fb_check_refit()`: `fit`, invisibly.

A list: `X` (the matrix), `terms`, `frame` (the model frame), `xlev`
(the factor levels), `assign` (the term index of each column),
`column_terms` (the term label of each column) and `asreml_order` (the
column positions in ASReml's order: term by term, and within a term by
descending factor levels with the first factor varying slowest).

## Details

- `fb_model_terms()` returns the model description of a fit: the parsed
  fixed, random, residual and addition terms, the response, the family
  and the link, as the engines received them.

- `fb_response_values()` evaluates a response, a column name or an
  expression such as `"log(y)"`, over `data`.

- `fb_offset_values()` sums the offsets of a model over the rows of
  `data`, refusing an offset that is not a finite number on every row.

- `fb_model_weights()` returns the observation weights a model was
  fitted with.

- `fb_known_variance()` describes the known variance or covariance of a
  model and where it is carried in the data.

- `fb_engine_settings()` checks and returns sampler settings.

- `fb_convergence()` returns the convergence verdict that `summary(fit)`
  reports, with its evidence.

- `fb_inla_posterior_sample()` draws joint posterior samples of an INLA
  fit on the scale of the data.

- `fb_check_refit()` refuses a fit that cannot be refitted as the model
  it was fitted as.

- `fb_fixed_design()` builds the design matrix of a fixed formula.

Factor levels and contrasts are those of all of `data`. With
`train_rows`, a transform that learns from the rows it is evaluated on
([`scale()`](https://rdrr.io/r/base/scale.html),
[`poly()`](https://rdrr.io/r/stats/poly.html),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html),
[`splines::bs()`](https://rdrr.io/r/splines/bs.html)) learns from the
training rows alone and is applied with those parameters to every row,
the rule of [`stats::predict()`](https://rdrr.io/r/stats/predict.html).
Missing values are kept as rows of `NA`.

## Examples

``` r
fb_extension_api()
#> [1] 1
d <- data.frame(y = c(1, 4, 9), area = c(1, 2, 3))
fb_response_values("sqrt(y)", d)
#> [1] 1 2 3
fb_offset_values("offset(log(area))", d)
#> [1] 0.0000000 0.6931472 1.0986123
fb_engine_settings(list(chains = 2L))
#> $chains
#> [1] 2
#> 
fb_fixed_design(~ area, d)$X
#>   (Intercept) area
#> 1           1    1
#> 2           1    2
#> 3           1    3
#> attr(,"assign")
#> [1] 0 1
```
