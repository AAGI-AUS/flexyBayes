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

  an `fb_terms` object.

- ...:

  unused.

## Value

invisibly returns `x`.
