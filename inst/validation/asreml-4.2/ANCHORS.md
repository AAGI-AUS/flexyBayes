# Declared anchors for the banked ASReml evidence

Scenarios 01 to 15 were produced by `data-raw/asreml_oracle/run_tranche.sh`,
which refuses a dirty tree and builds the tarball from `git archive <sha>`.
Scenarios 20 to 22 are ASReml-only reference fits produced by
`data-raw/asreml_fa_reference.R`, which also refuses a dirty tree. The
commits below are the ones that evidence names.

`tests/testthat/test-oracle-fixture-provenance.R` fails when a fixture
names a commit that is not declared here, when a declared commit is not an
ancestor of HEAD, or when a fixture was produced on a dirty tree. It also
fails when its own checking function stops rejecting those things.

| Commit | Scenarios | ASReml | INLA | R | Run |
|---|---|---|---|---|---|
| `8638306fd8166dda4a73bfdef8987f20dfe7991a` | 01_iid_intercept, 02_nested_stratum, 03_diag_gxe, 04_dsum_residual, 05_combined_met, 06_us_dim2, 07_us_dim3, 15_weights | 4.2.0.392 | 25.10.19 | 4.5.2 | 2026-09-11 |
| `76de0074382930634779153097f6d4f107fa4b26` | 20_fa_besag_met, 21_fa_besag_met_30, 22_fa_dasilva_maize | 4.2.0.482 | not used | 4.6.1 | 2026-10-09 |

A superseded anchor stays listed once its fixtures are replaced only if
some fixture still names it. Removing a row while a fixture still points
at it is what the test exists to catch.
