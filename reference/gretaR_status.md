# Introspect the gretaR backend slot

Returns the activation state of the gretaR backend slot plus the
procedure to activate it. The gretaR slot is provisioned at v0.2 but not
yet a live backend; activation lands in v0.3 when the gretaR package is
publicly available and has cleared the audit gate. Until then, this
helper is the canonical discoverability surface — users who see "gretaR"
referenced in the backend matrix or release notes can call
`gretaR_status()` to inspect why the slot is dormant and how to wake it.

## Usage

``` r
gretaR_status()
```

## Value

A list with components:

- `activated`:

  Logical: `TRUE` if `options(flexyBayes.gretaR_activated)` is set AND a
  future audit mechanism has cleared. v0.2: always `FALSE`.

- `gretaR_installed`:

  Logical: result of `nzchar(system.file(package = "gretaR"))`.

- `audit_clean`:

  Logical or `NA`: audit-status indicator. v0.2: always `NA` (no
  audit-status mechanism shipped).

- `dormancy_reason`:

  Character: one of `"slot_provisioned_not_activated"`,
  `"gretaR_not_installed"`, `"gretaR_not_audit_clean"`, or
  `"gretaR_dispatch_eligible"` (the last meaning the slot is fully
  active).

- `activation_procedure`:

  Character vector: numbered steps to activate the slot, in order.

## Examples

``` r
gs <- gretaR_status()
gs$activated         # v0.2: FALSE
#> [1] FALSE
gs$dormancy_reason   # v0.2: "slot_provisioned_not_activated" or
#> [1] "gretaR_not_installed"
                     #       "gretaR_not_installed"
cat(gs$activation_procedure, sep = "\n")
#> Install the gretaR package when it goes public on CRAN.
#> Verify the local gretaR build against the audit checklist.
#> Run options(flexyBayes.gretaR_activated = TRUE) in the session.
```
