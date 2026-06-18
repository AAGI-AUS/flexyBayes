# Ingest a brms-format formula into the flexyBayes IR

Parses a brms / lme4-style two-sided formula
(`response ~ fixed + (1 | g)`) into a `fb_terms` object – flexyBayes's
backend-agnostic intermediate representation (IR). The IR is what every
engine emits from, so building it explicitly lets a power user inspect
the parsed model, cache it, or hand it to a fitting verb. Most users
never call this directly:
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
detect a brms-style `formula` argument and build the IR internally.

## Usage

``` r
fb_from_brms(
  formula,
  data,
  family = "gaussian",
  link = NULL,
  prior = NULL,
  weights = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  carry_n_rows = NULL,
  ...
)
```

## Arguments

- formula:

  A base-R formula (recommended) or a `brmsformula` object.
  `brmsformula` support requires brms to be installed.

- data:

  A data.frame containing every referenced variable. May be `NULL` only
  on the advanced metadata-only path (see `carry_n_rows`).

- family:

  Character family name or a base-R
  [`family()`](https://rdrr.io/r/stats/family.html) object.

- link:

  Character link override, or `NULL` for the family default.

- prior:

  An
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  object, or `NULL`.

- weights:

  Optional numeric weights vector of length `nrow(data)`, mapped to a
  single addition term with `type = "weights"`.

- prior_fixed_sd:

  Numeric SD for fixed-effect normal priors when `prior` is `NULL`.

- prior_vc_sd:

  Numeric hyperparameter for the variance-component priors when `prior`
  is `NULL`.

- carry_n_rows:

  Advanced: a positive integer enabling the metadata-only IR path
  (`data = NULL`). The formula's variables are realised as a one-row
  placeholder and the row count is recorded as `carry_n_rows`, for
  stress-testing the preflight layer without materialising the full
  data. Leave `NULL` for ordinary use.

- ...:

  Reserved for future brms-ingest options (specials, autocor handling);
  currently unused.

## Value

An `fb_terms` object with `source = "brms"`.

## Details

Linear fixed effects, random intercepts (`(1 | g)`, `(1 || g)`) and
uncorrelated random intercept+slope (`(x || g)`) are supported;
correlated random slopes (`(x | g)`), smoothers, Gaussian processes and
autocorrelation terms refuse at ingest with a structured message.

## See also

[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal fitting entry;
[`fb_from_asreml()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_asreml.md)
and
[`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md)
for the other dialects.

Other flexyBayes ingest adapters:
[`fb_from_asreml()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_asreml.md),
[`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md)

## Examples

``` r
df <- data.frame(
  y = rnorm(30),
  x = rnorm(30),
  g = factor(rep(letters[1:5], 6))
)
ir <- fb_from_brms(y ~ x + (1 | g), data = df)
class(ir)
#> [1] "fb_terms" "list"    
```
