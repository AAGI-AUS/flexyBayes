# Backend dispatch trace for a flexyBayes fit

Returns the dispatch trace recorded at fit time: which backend was
selected, which gate checks ran, and why. On fb_greta() fits the trace
is trivial (the user bypassed the gate by entering on the greta-direct
path).

## Usage

``` r
backend_decision(fit)
```

## Arguments

- fit:

  A `flexybayes` object.

## Value

A list with the following components. The first four are present on
every fit; the four routing-trace fields are present on v0.3.6+ fits and
NULL on earlier fits for backward compatibility.

- `backend`:

  Character; one of `"greta"`, `"inla"`, `"brms"`, `"gretaR"`.

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
