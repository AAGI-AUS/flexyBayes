# Ingest an ASReml-format model specification into the flexyBayes IR

Parses an ASReml-style `fixed` / `random` / `rcov` specification into a
`fb_terms` object – flexyBayes's backend-agnostic intermediate
representation (IR) of a model. The IR is what every engine emits from,
so building it explicitly lets a power user inspect the parsed model,
cache it, or hand it to a fitting verb. Most users never call this
directly:
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
build the IR internally from the same arguments. Argument names and
defaults mirror
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
one-for-one, so this is a drop-in for that function's parsing step.

## Usage

``` r
fb_from_asreml(
  fixed,
  random = NULL,
  rcov = NULL,
  data,
  family = "gaussian",
  link = NULL,
  weights = NULL,
  known_matrices = list(),
  prior = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1
)
```

## Arguments

- fixed:

  Two-sided formula: `response ~ fixed_effects`.

- random:

  One-sided formula `~ random_terms` (ASReml syntax), or `NULL`.

- rcov:

  One-sided formula `~ residual_structure`, or `NULL`. `NULL` defaults
  to iid residuals (`list(list(type = "units"))`), matching
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md).

- data:

  A data.frame containing every variable referenced.

- family:

  Character family name (`gaussian`, `binomial`, `binary`, `poisson`,
  `negative_binomial`, `negbinom`, `gamma`, `beta`).

- link:

  Character link override, or `NULL` for the family default.

- weights:

  Optional numeric vector of length `nrow(data)`. When non-`NULL`,
  mapped to a single addition term with `type = "weights"`.

- known_matrices:

  Named list of known matrices (e.g. `list(Gmat = G)`); names are
  recorded in `data_summary$known_matrices`.

- prior:

  An
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  object, or `NULL`. `NULL` falls back to the scalar priors
  `prior_fixed_sd` / `prior_vc_sd`.

- prior_fixed_sd:

  Numeric SD for fixed-effect normal priors when `prior` is `NULL`.

- prior_vc_sd:

  Numeric hyperparameter for the variance-component priors when `prior`
  is `NULL`.

## Value

An `fb_terms` object with `source = "asreml"`.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal fitting entry;
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
and
[`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md)
for the other ingest dialects.

Other flexyBayes ingest adapters:
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md),
[`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md)

## Examples

``` r
df <- data.frame(
  yield = rnorm(20),
  geno  = factor(rep(letters[1:4], 5)),
  env   = factor(rep(c("a", "b"), 10))
)
ir <- fb_from_asreml(yield ~ env, random = ~ geno, data = df)
class(ir)
#> [1] "fb_terms" "list"    
```
