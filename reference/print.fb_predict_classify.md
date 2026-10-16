# Print a classify means table

Prints the table with a banner naming what the numbers are. The estimate
is a posterior mean. The banner also states that the table has no
standard errors of difference, which ASReml's
[`predict()`](https://rdrr.io/r/stats/predict.html) prints beside its
means.

## Usage

``` r
# S3 method for class 'fb_predict_classify'
print(x, digits = 4L, ...)
```

## Arguments

- x:

  An `fb_predict_classify` table, as returned by
  `predict(fit, classify = )`.

- digits:

  Number of significant digits for the printed table.

- ...:

  Ignored, present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the table it prints.

## Details

When the table includes random effects, the banner names the random
terms included, the variables averaged over and the variables left out,
as ASReml's own notes do. Otherwise, on an INLA fit, the banner says the
interval comes from the Gaussian approximation of the joint fixed-effect
posterior, not from INLA's own marginal densities.
