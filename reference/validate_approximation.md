# Validate an approximate model fit against its bias bound

`validate_approximation()` reports how much of a fitted model's
structure was lost to its approximation scheme, measured against the
scheme's declared pass threshold. It is the user-facing entry to the
per-scheme validation procedure registered for every approximate route;
the contract surfaces the realised error number while the user keeps the
accept / re-fit judgement.

## Usage

``` r
validate_approximation(fit, ...)
```

## Arguments

- fit:

  A fitted `flexybayes` object.

- ...:

  Passed to the per-scheme validation procedure (e.g. `threshold` for
  `low_rank_smooth`).

## Value

An `<fb_approximation_validation>` object: the scheme, the overall pass
/ fail verdict, the pass threshold, one result row per approximated
smooth (realised capture, bias bound, per-smooth pass flag), and the
registry's fallback hint.

## Details

Dispatch is on the fit's registered approximation scheme. At present the
only registered scheme is `low_rank_smooth` (the rank-K
principal-component truncation of an `s()` smooth basis): for such a
fit, the procedure reports the realised Frobenius capture \\\sum\_{i \le
K} d_i^2 / \sum_i d_i^2\\ of each truncated smooth against the default
pass threshold of `0.99`, where \\d_i\\ are the singular values of the
full smooth basis.

A fit carrying no recognised approximation (an exact fit) is refused
rather than returned as a vacuous pass.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for fitting; the approximation registry records each scheme's bias bound
and fallback.

## Examples

``` r
# An exact fit (no registered approximation) refuses rather than
# returning a vacuous pass. \donttest{} + an INLA guard: this needs
# a live fit and does not fire on a machine without INLA installed.
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
                     verbose = FALSE)
  tryCatch(validate_approximation(fit), error = function(e) {
    message(conditionMessage(e)) # "carries no recognised approximation"
  })
}
# }

# The low_rank_smooth verdict shape, on a plan-only object carrying
# the scheme (no active engine reaches low_rank_smooth on the
# current backends -- it was built for a since-withdrawn engine and
# now refuses at dispatch; see R/emit_smooth_low_rank.R). This part
# needs no engine and always runs.
demo_fit <- structure(
  list(
    exactness = "approximate_low_rank_smooth",
    extras = list(parse_info = list(approx = list(
      x = list(
        scheme = "low_rank_smooth", rank = 4L, k = 9L,
        frobenius_capture = 0.995,
        V_K = matrix(0, 9L, 4L), singular_values = rep(1, 9L)
      )
    )))
  ),
  class = c("flexybayes", "list")
)
v <- validate_approximation(demo_fit)
print(v)
#> <fb_approximation_validation>
#>   scheme:    low_rank_smooth
#>   threshold: Frobenius capture >= 0.99
#>   verdict:   PASS
#>     ok s(x): rank 4/9  capture 0.995  (bias bound 0.005)
v$pass
#> [1] TRUE
```
