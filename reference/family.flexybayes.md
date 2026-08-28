# Extract model family

Extract model family

## Usage

``` r
# S3 method for class 'flexybayes'
family(object, ...)
```

## Arguments

- object:

  A fitted `flexybayes` object of any backend.

- ...:

  Ignored, present for compatibility with the generic.

## Value

The response family the model was fitted under, as the `family` object
the fit recorded at emit time – so it names the family and link the
engine actually used, not the string the caller passed.
