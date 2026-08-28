# Random-effect predictions from a flexyBayes fit

The random-effect table, one data frame per grouping factor, in the same
columns on every engine: `group`, `level`, `estimate`, `std.error`,
`conf.low`, `conf.high`. Identical to `coef(object, what = "random")`,
and present under this name because it is the name the mixed-model
ecosystem uses.

## Usage

``` r
ranef(object, ...)

# S3 method for class 'flexybayes'
ranef(object, ...)
```

## Arguments

- object:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

A named list of data frames, one per grouping factor. Empty when the
model carries no random terms.

## Details

The method is registered on both this package's own `ranef()` generic
and on nlme's, which is the generic lme4 and brms re-export. Attaching
any of those masks the local generic, and the registration on nlme's is
what keeps a bare `ranef(fit)` dispatching there. In a session where the
masking is ambiguous, call `flexyBayes::ranef(fit)`.

## See also

[`coef.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes.md),
which this delegates to.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
                    verbose = FALSE)
  re <- ranef(fit)
  str(re, max.level = 1)
}
# }
```
