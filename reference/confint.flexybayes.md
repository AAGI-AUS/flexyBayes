# Credible intervals for fixed effects

Returns posterior quantile-based credible intervals, not frequentist
confidence intervals.

## Usage

``` r
# S3 method for class 'flexybayes'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A flexybayes object

- parm:

  Parameter names (NULL for all fixed effects)

- level:

  Credible level (default 0.95)

- ...:

  Additional arguments (ignored)

## Value

Matrix with lower and upper credible bounds
