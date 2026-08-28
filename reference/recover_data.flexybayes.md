# emmeans support: recover model data (default method)

The bare `flexybayes` method: reached by any fit whose engine has no
more specific override, which today is the brms backend.

## Usage

``` r
# S3 method for class 'flexybayes'
recover_data(object, ...)

# S3 method for class 'flexybayes_inla'
recover_data(object, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- ...:

  Passed to `emmeans::recover_data()`.

## Value

A data frame of predictors for the reference grid.
