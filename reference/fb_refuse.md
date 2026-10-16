# Raise a registered refusal

`fb_refuse()` stops with a refusal condition for a reason registered by
flexyBayes or by a package through
[`fb_register_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_register_refusals.md).
The condition has the classes `flexybayes_refusal_<reason_code>`, any
classes given in `class`, and `flexybayes_refusal`, so users catch a
plug-in's refusals exactly as they catch the core's.

## Usage

``` r
fb_refuse(reason_code, message = NULL, ..., class = character(0), call = NULL)
```

## Arguments

- reason_code:

  A registered reason code (see
  [`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md)).

- message:

  The message shown to the user. Defaults to the reason's registered
  description. A message that names the user's own variables or values
  is usually more helpful.

- ...:

  Further named fields to carry on the condition object, for handlers
  that need the offending value. The names `message`, `reason_code`,
  `package` and `call` are set by `fb_refuse()` and cannot be used.

- class:

  Additional condition classes, placed between the per-code class and
  `flexybayes_refusal`.

- call:

  The call to record on the condition. Defaults to `NULL`, as for the
  core's refusals.

## Value

Does not return. It signals an error. The condition carries `message`,
`reason_code`, `package` (the package that registered the code), `call`
and the fields given in `...`.

## See also

[`fb_register_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_register_refusals.md),
[`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md).

## Examples

``` r
res <- tryCatch(
  fb_refuse("unsupported_family"),
  flexybayes_refusal = function(e) e
)
res$reason_code
#> [1] "unsupported_family"
res$package
#> [1] "flexyBayes"
```
