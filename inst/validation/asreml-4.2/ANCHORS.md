# Declared anchors for the banked ASReml evidence

Every `.rds` beside this file was produced by
`data-raw/asreml_oracle/run_tranche.sh`, which refuses a dirty tree and
builds the tarball from `git archive <sha>` rather than from the working
directory. The commits below are the ones that evidence names.

`tests/testthat/test-oracle-fixture-provenance.R` fails when a fixture
names a commit that is not declared here, when a declared commit is not an
ancestor of HEAD, or when a fixture was produced on a dirty tree. It also
fails when its own checking function stops rejecting those things.

| Commit | Scenarios | ASReml | INLA | R | Run |
|---|---|---|---|---|---|
| `8638306fd8166dda4a73bfdef8987f20dfe7991a` | 01_iid_intercept, 02_nested_stratum, 03_diag_gxe, 04_dsum_residual, 05_combined_met, 06_us_dim2, 07_us_dim3, 15_weights | 4.2.0.392 | 25.10.19 | 4.5.2 | 2026-09-11 |

A superseded anchor stays listed once its fixtures are replaced only if
some fixture still names it. Removing a row while a fixture still points
at it is what the test exists to catch.
