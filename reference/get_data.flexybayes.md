# Model data accessor (greta backend)

Registered for `insight::get_data()` so marginaleffects (which discovers
a model's data through insight) can build reference grids and average
predictions without an explicit `newdata`.

## Usage

``` r
# S3 method for class 'flexybayes'
get_data(x, ...)

# S3 method for class 'flexybayes_inla'
get_data(x, ...)
```

## Arguments

- x:

  A `flexybayes_inla` fit.

- ...:

  Ignored.

## Value

The model data frame.
