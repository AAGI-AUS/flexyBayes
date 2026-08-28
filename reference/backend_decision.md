# Backend dispatch trace for a flexyBayes fit

Returns the dispatch trace recorded at fit time: which backend was
selected, which gate checks ran, and why. Engine-agnostic – works on an
`flexybayes_inla` or `flexybayes_brms` fit alike, reading
`fit$extras$backend_decision` as populated by
`.build_routing_decision()` during dispatch (see the top of this file).

## Usage

``` r
backend_decision(fit)
```

## Arguments

- fit:

  A `flexybayes` object.

## Value

A list holding the recorded dispatch trace. The first four components
are present on every fit. The four routing-trace fields are present on
fits from v0.3.6 onwards and `NULL` on earlier ones, for backward
compatibility.

- `backend`:

  Character; one of `"inla"`, `"brms"`.

- `path`:

  Character; the dispatch path token.

- `gate_checks`:

  List or NULL; the `lgm_gate()` check trail (failures on refusal;
  capabilities on accept).

- `reason`:

  Character; the dispatch-decision rationale.

- `preflight_summary`:

  An `<fb_preflight>` object or NULL. Populated when `.fb_preflight()`
  ran (\>1e5-row path); NULL on the small-data fast path.

- `representation_plan`:

  Named list of slim per-term entries
  `(term_id, representation_class, justification)` derived from
  `preflight_summary`; NULL when no preflight.

- `rejected_routes`:

  List of `(backend, reason)` pairs for the backends considered but not
  chosen. Empty for explicit user requests (the routing policy is
  bypassed when the backend is named directly).

- `routing_policy_version`:

  Character; e.g. `"stage5a_v1"`. The audit-anchor for reproducibility –
  a policy change bumps this string.

## Examples

``` r
# backend_decision() reads a slot recorded on an already-fitted
# object, so (unlike fb_plan()) an engine has to run first. Wrapped
# in \donttest{} and guarded on INLA so this does not fire on a
# machine without it.
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
                     verbose = FALSE)
  bd <- backend_decision(fit)
  bd$backend   # "inla"
  bd$path      # the dispatch path token
  bd$reason    # why this backend was chosen
}
# }
```
