# Extract coefficients from a flexyBayes fit

Returns the fixed effects by default, and the random-effect predictions
or the unobserved design cells on request. An ASReml user reaches for
`coef(fit)$random`; the equivalent here is `coef(fit, what = "random")`,
which is what
[`ranef()`](https://aagi-aus.github.io/flexyBayes/reference/ranef.md)
calls.

## Usage

``` r
# S3 method for class 'flexybayes'
coef(object, what = c("fixed", "random", "missing", "all"), ...)
```

## Arguments

- object:

  A fitted `flexybayes` object of any backend.

- what:

  Which part of the fit to return. `"fixed"` (the default) is the
  fixed-effect posterior means. `"random"` is one data frame per
  grouping factor, with the columns `group`, `level`, `estimate`,
  `std.error`, `conf.low` and `conf.high`; a grouping factor carrying
  more than one effect names the effect in `level`. `"missing"` is the
  table of unobserved design cells, in the columns `row`, `estimate`,
  `std.error`, `conf.low`, `conf.high` plus any design index variables
  the fit recorded. `"all"` is the named list of all three.

- ...:

  Ignored, present for compatibility with the generic.

## Value

For `what = "fixed"`, a named numeric vector of the fixed effects'
posterior means on the treatment-contrast basis, named for the
design-matrix columns, empty when the model carries no fixed effects.
For `"random"`, a named list of data frames, empty when the model
carries no random terms. For `"missing"`, a data frame with zero rows
when every response was observed. For `"all"`, a list with the elements
`fixed`, `random` and `missing`.

## Details

The default is the historical return – a named numeric vector of the
fixed effects – so every caller that treats `coef(fit)` as a numeric
vector keeps working, including the emmeans and marginaleffects seams,
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and the tidiers.

## See also

[`ranef()`](https://aagi-aus.github.io/flexyBayes/reference/ranef.md)
for the random-effect table on its own,
[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md)
whose `random` and `missing` slots are the same two objects, and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) for the design and
observed counts.
