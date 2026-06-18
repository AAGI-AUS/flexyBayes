# Emit the generated backend code for a deferred review object

Writes the backend code carried by a `<flexybayes_review>` object (greta
R code for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
/
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md);
Stan code via `brms::stancode()` for
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md))
to a connection.

## Usage

``` r
cat_code(x, ...)

# S3 method for class 'flexybayes_review'
cat_code(x, file = stdout(), ...)
```

## Arguments

- x:

  a `<flexybayes_review>` object.

- ...:

  method-specific arguments. The `flexybayes_review` method accepts
  `file` (connection; default
  [`stdout()`](https://rdrr.io/r/base/showConnections.html)).

- file:

  connection to write to (default
  [`stdout()`](https://rdrr.io/r/base/showConnections.html)).

## Value

invisibly returns the code string.
