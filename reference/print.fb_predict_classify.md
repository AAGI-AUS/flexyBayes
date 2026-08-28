# Print a classify means table

Prints the table with a banner naming what the numbers are. The estimate
is a posterior mean of a marginal mean, and the table is not a
standard-error-of-difference table – ASReml's
[`predict()`](https://rdrr.io/r/stats/predict.html) prints one beside
the means and this does not, so the banner says so rather than leaving
the absence to be noticed.

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

On an INLA fit the banner carries a second line. The interval there
comes from the Gaussian approximation of the joint fixed-effect
posterior, which is what the estimation basis supplies to emmeans, and
not from INLA's own marginal densities.
