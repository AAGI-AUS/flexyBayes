# Emit the generated backend code for a deferred review object

Writes the backend code carried by a `<flexybayes_review>` object (Stan
code via `brms::stancode()` for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
/
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

  A `<flexybayes_review>` object carrying the generated backend code.

- ...:

  Method-specific arguments. The `flexybayes_review` method accepts
  `file`, a connection defaulting to
  [`stdout()`](https://rdrr.io/r/base/showConnections.html).

- file:

  A connection to write the code to, defaulting to
  [`stdout()`](https://rdrr.io/r/base/showConnections.html).

## Value

Invisibly, the code as a character vector. Called for the code it writes
to `file`.

## Examples

``` r
# review_code = TRUE defers the fit, but generating the code still
# goes through brms::make_stancode(), so brms must be present.
if (requireNamespace("brms", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  rv <- flexybayes(y ~ x + (1 | g), data = d, backend = "brms",
                   review_code = TRUE)
  cat_code(rv)
}
```
