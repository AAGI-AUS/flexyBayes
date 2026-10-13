# O02: which structures return the same answer twice

**Run** 2026-09-11, `data-raw/asreml_oracle/screen_repeatability.R`, three
identical fits per structure, `agridat::besag.met` for the multi-environment
rows and a fixed-seed simulation for the slope rows.

A single paired comparison against ASReml is not evidence for a structure
whose repeated fits disagree with each other: the number it produces is one
draw from a distribution over modes. So the screen runs before any further
oracle scenario, and it decides which rows can be compared at all.

## Criterion, fixed before the screen ran

- **collapsed**: a component's SD below 1% of the residual SD, the
  package's own `.FB_COLLAPSE_FRACTION`.
- **UNSTABLE**: any component collapsed in one run and not another, or any
  non-collapsed component whose maximum across runs exceeds twice its
  minimum. A factor of two is generous on purpose: INLA's integration is not
  deterministic and small movement is expected, but a component that doubles
  between identical calls is not small movement.
- **ALL_COLLAPSED**: every run pinned the same component at zero.
- **STABLE**: none of the above.

## Amendment, made during the run

The first pass returned **STABLE** for `diag_gxe`, the structure already
measured as bistable across 24 fits. All three runs had collapsed alike, so
there was no flip and no live component to compare, and a consistently wrong
answer scored as a repeatable one.

Repeatability is not correctness. `ALL_COLLAPSED` was added as a third
verdict, and the restart probe folded into the screen as the sharper signal:
one fit plus one restart settles in seconds what three repeats cannot settle
at all.

This is worth recording as a property of the instrument. A screen that can
return STABLE for the worst structure in the set would have sent the oracle
suite on to compare it against ASReml as though the number meant something.

## Result, before adoption was implemented

| Structure | Verdict | Evidence |
|---|---|---|
| `dsum_residual` | STABLE | max ratio 1.36 |
| `us_dim2` | STABLE | max ratio 1.01 |
| `us_dim3` | STABLE | max ratio 1.02 |
| `cor_slope` | STABLE | max ratio 1.00 |
| `uncor_slope` | STABLE | max ratio 1.00 |
| `fa1` | UNSTABLE | max ratio **15.36** |
| `combined_met` | UNSTABLE | max ratio **8.36** |
| `diag_gxe` | UNSTABLE | collapse flip |
| `iid_intercept` | UNSTABLE | collapse flip |
| `nested` | ALL_COLLAPSED | restart probe flagged, delta **+132.1** |

`iid_intercept` is the simplest random effect the package fits, and it flips
on this dataset. That places the cause in the data and not in exotic
structures: `besag.met` is low-signal, with a residual SD near 19.7 against
genotype SDs near 5. Weak identification is what makes a second mode
reachable, and the more hyperparameters a structure carries, the more often
the optimiser finds it, which is why the factor-analytic row swings by
fifteen.

## Result, after mode adoption

Same screen, same data, same criterion, re-run against the rebuilt tarball
carrying the fix recorded in O01.

| Structure | Before | After | Probe |
|---|---|---|---|
| `iid_intercept` | UNSTABLE (flip) | **STABLE** (1.00) | not needed |
| `nested` | ALL_COLLAPSED | **STABLE** (1.10) | flagged, +132.3 |
| `dsum_residual` | STABLE (1.36) | STABLE (1.71) | not needed |
| `us_dim2` | STABLE (1.01) | STABLE (1.01) | not needed |
| `us_dim3` | STABLE (1.02) | STABLE (1.04) | not needed |
| `cor_slope` | STABLE (1.00) | STABLE (1.00) | not needed |
| `uncor_slope` | STABLE (1.00) | STABLE (1.00) | not needed |
| `diag_gxe` | UNSTABLE (flip) | UNSTABLE (2.29) | flagged, +194.2 |
| `combined_met` | UNSTABLE (8.36) | UNSTABLE (37.51) | flagged, +26.5 |
| `fa1` | UNSTABLE (15.36) | UNSTABLE (15.46) | flagged, +6.2 |

Seven of ten stable, against five before.

**What the fix does and does not do.** It closes the collapse class: a
component pinned at zero by mode selection is now found and replaced, and the
two structures whose instability was entirely collapse (`iid_intercept` and
`nested`) are stable. `diag_gxe` no longer flips between zero and real, and
its remaining spread of 2.29 sits just above a threshold set at 2.

It does not close instability in the high-dimensional structures.
`combined_met` carries twelve components and `fa1` thirteen, and both remain
badly unstable; `fa1` is unchanged at a factor of fifteen. Their ratios are
now computed between real values and not between zero and a real value,
which is a different and more ordinary kind of disagreement, but three
identical calls still return answers that differ by more than an order of
magnitude.

Adoption picks the better of two modes. Where the surface has many modes and
none is reached reliably, one restart from one perturbed start is not enough.
Those three rows stay out of the paired comparison.
## Consequence

Five of ten rows could be taken to the paired ASReml comparison as they
stood. The other five could not, and the answer was not a wider tolerance but
the mode-adoption fix recorded in O01.
