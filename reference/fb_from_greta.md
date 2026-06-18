# Ingest a user-built greta model into the flexyBayes IR

Wraps a `greta_model` (built with `greta::model(...)`) in a `fb_terms`
object – flexyBayes's backend-agnostic intermediate representation (IR)
– so a natively-specified greta graph can flow through the same
downstream machinery (summaries,
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md),
canonical-name mapping) as a formula-ingested model. Unlike the formula
adapters this is a post-hoc wrapper around an already-built graph, so
the resulting IR is greta-only by construction: it carries no fixed /
random / rcov term lists, only the populated `greta_meta` slot.

## Usage

``` r
fb_from_greta(
  model,
  data = NULL,
  prior = NULL,
  canonical_names = NULL,
  known_matrices = list()
)
```

## Arguments

- model:

  A `greta_model` returned by `greta::model(...)`.

- data:

  Optional data.frame used to build the model, recorded on the IR for
  downstream methods. Not required for fitting – greta has already
  captured the data into its TensorFlow graph.

- prior:

  Optional
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  object. When supplied, every name must resolve to a target
  greta_array; semantic agreement with the graph's encoded priors is the
  caller's responsibility.

- canonical_names:

  Optional named character vector mapping greta-side parameter names to
  canonical names. `NULL` falls back to the verbatim greta names (with a
  one-time silenceable note, since
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  against a non-greta backend needs canonical names).

- known_matrices:

  Named list, mirroring the ASReml entry; recorded in
  `data_summary$known_matrices`.

## Value

An `fb_terms` object with `source = "greta"`, `intercept = NA`, empty
fixed / random / rcov term lists, and the populated `greta_meta` slot
(carrying the model graph, target arrays, and the canonical-name map).
Pass the returned IR straight to
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
to fit the graph via greta while keeping a canonical-name map for
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md):
`fb(fb_from_greta(model, canonical_names = ...))`.

## See also

[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
and
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
for the universal fitting entry (which accept the returned IR directly);
[`fb_from_asreml()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_asreml.md)
and
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
for the formula dialects.

Other flexyBayes ingest adapters:
[`fb_from_asreml()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_asreml.md),
[`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(greta)
mu    <- normal(0, 10)
sigma <- normal(0, 5, truncation = c(0, Inf))
y     <- as_data(rnorm(20))
distribution(y) <- normal(mu, sigma)
m  <- model(mu, sigma)
ir <- fb_from_greta(m, canonical_names = c(mu = "(Intercept)"))
fit <- fb(ir)   # fit the native graph via greta, keeping the map
} # }
```
