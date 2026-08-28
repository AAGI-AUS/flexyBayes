# Canonical parameter-name view for a flexyBayes fit

Returns the backend-native -\> canonical parameter-name map for a fit,
plus per-parameter value transforms where applicable (e.g., the INLA
precision-to-SD `sqrt(1/prec)` transform applied to hyperparameters
before triangulation). The canonical convention follows brms
(`(Intercept)`, `<term>`, `sd_<group>`, `sigma`, `r_<group>[<level>]`).

## Usage

``` r
canonical_names(fit, drop = FALSE, ...)

# S3 method for class 'flexybayes_inla'
canonical_names(fit, drop = FALSE, ...)

# S3 method for class 'flexybayes_brms'
canonical_names(fit, drop = FALSE, ...)
```

## Arguments

- fit:

  A `flexybayes_inla` or `flexybayes_brms` object.

- drop:

  Logical: if `TRUE` (default `FALSE`), drop backend-native names that
  are not in the registered map (e.g., INLA's `Predictor.<i>`
  latent-predictor draws). When `FALSE`, un-mapped names appear in the
  returned `$unmapped` element.

- ...:

  Additional arguments (ignored by current methods).

## Value

A list holding the canonical-name map and its provenance.

- `map`:

  Named character vector keyed by backend-native parameter name with
  canonical name as the value.

- `transform`:

  Named list of `function(x) -> x'` transforms keyed by canonical name.
  Empty list when no transforms apply.

- `source`:

  Character: `"registry"`, `"user"`, `"registry_fallback_verbatim"`, or
  `"legacy_inferred"`.

- `unmapped`:

  Character vector of backend-native names not in the map (when
  `drop = FALSE`).

- `prior_parametrization`:

  Character, present only on aggregated fits, naming which prior the fit
  ran under. `"per_row_equivalent"` is the legacy scalar bridge, whose
  precision prior makes the aggregated posterior match the per-row
  posterior to numerical precision. `"package_default"` is the automatic
  bounded-uniform-on-SD prior, which claims no such equivalence on this
  route. `"custom"` is an explicit prior from the caller (see the
  "Matched priors" note on
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)).

## Details

On `flexybayes_inla` and `flexybayes_brms` fits, the per-backend mapper
registered at package load (`inla` or `brms`) drives the resolution; the
returned map is cached on `fit$extras$canonical_map` for fast repeated
access.

## Examples

``` r
# The map is engine-specific, so it needs a fit to read.
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
                    verbose = FALSE)
  cn <- canonical_names(fit)
  head(cn$map)
}
# }
```
