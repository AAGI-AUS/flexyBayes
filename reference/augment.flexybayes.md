# Augment a flexyBayes fit with fitted values and residuals

Returns the model frame with two observation-level columns added: the
posterior-mean fitted value and the response residual.

## Usage

``` r
# S3 method for class 'flexybayes'
augment(x, data = NULL, ...)

# S3 method for class 'flexybayes_inla'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  A flexyBayes fit (`flexybayes` or `flexybayes_brms`).

- data:

  Optional `data.frame` to augment. Defaults to the data the model was
  fitted to.

- ...:

  Currently unused; present for generic compatibility.

## Value

The supplied (or original) `data.frame` with `.fitted` and `.resid`
columns appended.

## See also

[`tidy.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
augment(fit)
} # }
```
