# Canonical parameter-name view for a flexyBayes fit

Returns the backend-native -\> canonical parameter-name map for a fit,
plus per-parameter value transforms where applicable (e.g., the INLA
precision-to-SD `sqrt(1/prec)` transform applied to hyperparameters
before triangulation). The canonical convention follows brms
(`(Intercept)`, `<term>`, `sd_<group>`, `sigma`, `r_<group>[<level>]`).

## Usage

``` r
canonical_names(fit, drop = FALSE, ...)

# S3 method for class 'flexybayes'
canonical_names(fit, drop = FALSE, ...)

# S3 method for class 'flexybayes_inla'
canonical_names(fit, drop = FALSE, ...)

# S3 method for class 'flexybayes_brms'
canonical_names(fit, drop = FALSE, ...)

# S3 method for class 'flexybayes_direct_greta'
canonical_names(fit, drop = FALSE, ...)
```

## Arguments

- fit:

  A `flexybayes` (or `flexybayes_inla`) object.

- drop:

  Logical: if `TRUE` (default `FALSE`), drop backend-native names that
  are not in the registered map (e.g., INLA's `Predictor.<i>`
  latent-predictor draws). When `FALSE`, un-mapped names appear in the
  returned `$unmapped` element.

- ...:

  Additional arguments (ignored by current methods).

## Value

A list with components:

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

  Character, present only on aggregated-gaussian fits:
  `"per_row_equivalent"` when the default precision prior is in force
  (the aggregated posterior then matches the per-row posterior to
  numerical precision) or `"custom"` when an explicit prior was supplied
  (see the "Matched priors" note on
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)).

## Details

On `flexybayes` and `flexybayes_inla` fits, the per-backend mapper
registered at package load (`greta` or `inla`) drives the resolution;
the returned map is cached on `fit$extras$canonical_map` for fast
repeated access. On `flexybayes_direct_greta` fits (built via
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md))
the map comes from the user-supplied `canonical_names` argument, with a
verbatim-greta-name fallback when the argument is omitted.
