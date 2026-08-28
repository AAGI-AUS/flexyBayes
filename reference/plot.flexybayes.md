# Plot diagnostics for a flexyBayes model

Draws one of seven standard displays for a fitted model, chosen by
`type`. Every backend reaches the same method, and a display that needs
a slot a given engine does not populate – posterior draws on the
deterministic INLA path, for instance – prints a message naming what is
missing and returns invisibly rather than erroring or drawing an empty
panel. The one exception is `"pp_check"`, which raises the refusal
[`pp_check.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/pp_check.flexybayes.md)
raises, so a caller can catch it by class.

## Usage

``` r
# S3 method for class 'flexybayes'
plot(x, type = NULL, ...)

# S3 method for class 'flexybayes_inla'
plot(x, ...)

# S3 method for class 'flexybayes_brms'
plot(x, ...)

# S3 method for class 'flexybayes_aggregated'
plot(x, ...)

# S3 method for class 'flexybayes_glm'
plot(x, ...)
```

## Arguments

- x:

  A fitted `flexybayes` object of any backend. Aggregated and
  generalised-linear fits reach the same method through their own
  registrations.

- type:

  A single string naming the display to draw. One of: `"diagnostics"`,
  trace plots and marginal densities per parameter (needs a sampled
  posterior); `"residuals"`, residuals against fitted values beside a
  normal quantile-quantile plot; `"effects"`, a forest plot of the fixed
  effects with their credible intervals, available on every backend that
  supplies [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html); `"variance"`, a
  bar chart of the variance components with credible intervals;
  `"blups"`, a caterpillar plot of the random-effect predictions ordered
  by magnitude; `"pp_check"`, a posterior predictive check overlaying
  replicated datasets on the observed response, which needs predictive
  draws and so refuses by name on a fit that carries none (see
  [`pp_check.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/pp_check.flexybayes.md));
  and `"variogram"`, the empirical semivariance of the residuals over
  the design index. `NULL` (the default) resolves to `"diagnostics"` on
  a sampled fit and `"residuals"` otherwise.

- ...:

  Further arguments passed to the underlying plotting call, for example
  `variable` to restrict which parameters a diagnostic display covers,
  or `type` and `ndraws` for `"pp_check"`.

## Value

Invisibly, the object the underlying plotting call returns – a ggplot2
object for the displays built with it, the semivariance table for
`"variogram"`, and `NULL` for those drawn on the base graphics device.
Called for the plot it draws.

## Details

The default display depends on what the fit carries. A sampled posterior
gets `"diagnostics"`, because the first question about a sampler is
whether it converged. A deterministic approximation has no chains to
trace, so an INLA fit gets `"residuals"` instead of a message declining
to draw the display it was asked for.

## See also

[`pp_check.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/pp_check.flexybayes.md)
for the posterior predictive check `type = "pp_check"` draws, and the
refusal it raises where a fit has no predictive draws.
