# Predict from a flexybayes_direct_greta fit

Direct-greta fits do not encode a Wilkinson-Rogers formula on the IR, so
[`predict()`](https://rdrr.io/r/stats/predict.html) cannot mechanically
apply the model matrix to `newdata`. The user must supply an explicit
predictor function `f(theta, newdata)` mapping a length-`p` named
parameter vector and a data frame to a length-`nrow(newdata)` vector of
predicted values. Posterior-mean prediction is computed by applying
`f()` to each posterior draw and averaging.

## Usage

``` r
# S3 method for class 'flexybayes_direct_greta'
predict(
  object,
  newdata,
  predictor,
  type = c("response", "link"),
  n_draws = 500L,
  ...
)
```

## Arguments

- object:

  A `flexybayes_direct_greta` object.

- newdata:

  A data frame.

- predictor:

  A function `function(theta, newdata) -> numeric` where `theta` is a
  named vector of canonical-named target parameters and `newdata` is the
  data frame passed in.

- type:

  Character; `"response"` (default; on the response scale) or `"link"`
  (on the linear-predictor scale). Currently identity link only;
  non-Gaussian families queued for v0.3.

- n_draws:

  Integer; number of posterior draws to use (default 500). The predictor
  function is applied to each draw and the mean is returned.

- ...:

  Additional arguments (ignored).

## Value

Numeric vector of length `nrow(newdata)`.
