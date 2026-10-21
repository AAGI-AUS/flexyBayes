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

| Commit | In history as | Scenarios | ASReml | INLA | R | Run |
|---|---|---|---|---|---|---|
| `8638306fd8166dda4a73bfdef8987f20dfe7991a` | `3c34a086f393375d61845f286cb79e28d07bc256` | 01_iid_intercept, 02_nested_stratum, 03_diag_gxe, 04_dsum_residual, 05_combined_met, 06_us_dim2, 07_us_dim3, 15_weights | 4.2.0.392 | 25.10.19 | 4.5.2 | 2026-09-11 |
| `76de0074382930634779153097f6d4f107fa4b26` | `88b1b9afc4ec1cbc621b3f313c3d1aaa79174088` | 20_fa_besag_met, 21_fa_besag_met_30, 22_fa_dasilva_maize | 4.2.0.482 | not used | 4.6.1 | 2026-10-09 |

The repository history was rewritten on 2026-10-15, before its first public
push of these commits. The fixtures keep the commit each run was built from,
in the first column; the second column is the same commit as it appears in
the rewritten history, taken from the rewrite's commit map. Between each pair
the package source (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`) and the scripts
that produced the evidence are identical. The test checks that the second
column is an ancestor of HEAD, and, where the original commit is still in the
local object store, that the package source of the two commits is identical.

A superseded anchor stays listed once its fixtures are replaced only if
some fixture still names it. Removing a row while a fixture still points
at it is what the test exists to catch.
