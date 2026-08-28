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
# Not run: augment() appends fitted values and residuals to the data
# a sampled fit was built from, so the example needs a Stan compile.
# The fragment is complete and runs as written with a C++ toolchain.
if (FALSE) { # \dontrun{
set.seed(1)
dat <- data.frame(yield = rnorm(40), env = factor(rep(1:4, each = 10)))
fit <- flexybayes(yield ~ env, data = dat, backend = "brms",
                  chains = 1L, n_samples = 200L, warmup = 100L)
augment(fit)
} # }
```
