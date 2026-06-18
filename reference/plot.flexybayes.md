# Plot diagnostics for a flexybayes model

Plot diagnostics for a flexybayes model

## Usage

``` r
# S3 method for class 'flexybayes'
plot(
  x,
  type = c("diagnostics", "residuals", "effects", "variance", "blups", "pp_check"),
  ...
)

# S3 method for class 'flexybayes_inla'
plot(x, ...)

# S3 method for class 'flexybayes_brms'
plot(x, ...)

# S3 method for class 'flexybayes_aggregated'
plot(x, ...)

# S3 method for class 'flexybayes_direct_greta'
plot(x, ...)

# S3 method for class 'flexybayes_glm'
plot(x, ...)
```

## Arguments

- x:

  A flexybayes object

- type:

  Character: type of plot to produce.

  - `"diagnostics"`: trace plots + density (requires bayesplot)

  - `"residuals"`: residuals vs fitted + QQ plot

  - `"effects"`: forest plot of fixed effects with credible intervals

  - `"variance"`: bar chart of variance components with credible
    intervals

  - `"blups"`: caterpillar plot of BLUPs

  - `"pp_check"`: posterior predictive check (observed vs replicated)

- ...:

  Additional arguments passed to plotting functions
