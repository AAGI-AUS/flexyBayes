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
