# Print method for fb_terms (intermediate representation)

Internal S3 method, registered for dispatch only. Used during
development and inside the
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
/
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
flow to inspect the parsed model object before backend dispatch.

## Usage

``` r
# S3 method for class 'fb_terms'
print(x, ...)
```

## Arguments

- x:

  An `fb_terms` object, the parsed intermediate representation of a
  model.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

Invisibly, `x` unchanged. Called for the term listing it prints.
