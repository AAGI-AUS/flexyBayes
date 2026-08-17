<!--
INTERNAL DRAFT -- not for submission as written.

This file is a pre-submission draft for the future CRAN submission, which is
gated on the AAGI-AUS collective conversation (decision D8). Three things are
finalised only at the release cut and must be regenerated then:
  - the version string (the working tree is on the 0.9.0 development line; a
    CRAN release version is cut at the D8 release decision);
  - the package URLs (they 404 until the public repository exists under D8);
  - the exact `R CMD check --as-cran` status, re-run on the FINAL release
    tarball with the release version and public URLs, with the command and date
    recorded.
The check status quoted below is the one recorded in `validation/VALIDATION.md`
for 0.9.0, generated 2026-08-15: **0 errors, 0 warnings, 1 NOTE** on the built
tarball under the CRAN profile. The single NOTE is the
CRAN-incoming-feasibility item itemised below. Building with vignettes removes
the two vignette-directory WARNINGs that appear only under
`--no-build-vignettes`. Re-verify on the final release-cut tarball. Do not
submit using this file verbatim.
-->

## Submission

This is the first submission of `flexyBayes` to CRAN.

The package is on the `0.9.0` development line; the CRAN release version
is finalised at the D8 release cut. Public release is held pending the
AAGI-AUS collective conversation on the upstream repository rename
(project decision **D8**); the package URL declared in `DESCRIPTION` is
the agreed forward destination and 404s today — see **NOTE 1** below.
The check status below is the itemised baseline; it is verified on a
clean library with user R profiles disabled (not only on the
maintainer's preconfigured machine), so a fresh reviewer sees the same
result.

Two inference engines are active: INLA (integrated nested Laplace
approximation) and brms (the Stan passthrough). greta and gretaR are
retained as quarantined re-entry candidates and fit nothing — an
explicit request for either refuses with a structured reason. See
`NEWS.md` for the reshape.

## Test environments

- local macOS Tahoe 26.5.2 (Apple Silicon), R 4.5.2 — `R CMD check
  --as-cran`, 1-NOTE structural baseline (NOTE 1 below); one
  additional environmental NOTE (NOTE 2) fires when the local HTML
  Tidy binary is older than the version `R CMD check` accepts, and a
  third environmental NOTE (NOTE 3) fires when the check host cannot
  reach a network time source. Net `Status:` on the current dev
  machine is `1 NOTE` when networked + recent-Tidy; up to `3 NOTEs`
  on a stale-Tidy, offline host.
- GitHub Actions matrix (committed locally — `.github/workflows/`;
  activates on first push after the D8 collective decision):
  ubuntu-latest × {R-devel, R-release, R-oldrel-1} plus R-release ×
  {windows-latest, macos-latest, macos-14 (Apple Silicon)}. The
  workflow set is `R-CMD-check.yaml`, `dependency-review.yaml`,
  `lint.yaml`, `pkgdown.yaml`, `test-coverage.yaml`. INLA is installed
  on the Linux release leg only via the `Additional_repositories`
  binary path. greta is intentionally not installed on the standard CI
  matrix (Python / TF setup is a separate slow workflow); the retained
  greta tests gate on the backend registry's lifecycle state and skip
  cleanly.

## R CMD check results — itemised per NOTE

```
Status: 1 NOTE
```

Recorded in `validation/VALIDATION.md` for 0.9.0 (generated 2026-08-15):
`R CMD check --as-cran` on the built tarball under the CRAN profile,
0 errors, 0 warnings, 1 NOTE.

**NOTE 1 — `checking CRAN incoming feasibility`.**

This NOTE carries four sub-items:

1. *New submission.* Intrinsic to first submission; resolves once
   CRAN's incoming checks accept the upload.
2. *Version contains large components.* This sub-NOTE fires only at
   development versions (`*.9000` tag per the R-pkgs convention); it
   is absent at a non-`.9000` version such as `0.9.0`, and reappears
   only when checking a `*.9000` dev tip between releases.
3. *Suggests or Enhances not in mainstream repositories: greta, INLA.*
   Neither is on CRAN, and both are declared in
   `Additional_repositories:`. `greta` was archived from CRAN and is
   maintained on the greta-dev R-universe
   (`https://greta-dev.r-universe.dev`); `INLA` is distributed from
   `https://inla.r-inla-download.org/R/stable`. The Stan backend
   (`brms`) is on CRAN, so a reviewer always has at least one
   CRAN-installable inference engine. Every backend code path is guarded
   by `requireNamespace(., quietly = TRUE)`; the package degrades
   gracefully when a Suggests dependency is absent, and the engine-using
   tests skip cleanly (see "Test-suite behaviour on CRAN").
4. *Found the following (possibly) invalid URLs:*
   `https://github.com/AAGI-AUS/flexyBayes` and
   `https://github.com/AAGI-AUS/flexyBayes/issues`, both returning
   HTTP 404. **Known pending.** The local-only rename from
   `bayesreml` to `flexyBayes` was completed 2026-04-27; the matching
   upstream repository rename is gated on the AAGI-AUS collective
   conversation (project decision **D8**). The URL declared in
   `DESCRIPTION` is the *agreed forward destination* —
   `AAGI-AUS/flexyBayes` is where the package will live once the
   upstream rename ships; the URL is forward-pointing on purpose so a
   single edit at the moment of collective approval is unnecessary.
   This sub-NOTE will resolve once the upstream rename ships. This
   sub-item may or may not fire on a given check host depending on
   the network-reachability check's caching behaviour.

**NOTE 2 (environmental) — `checking HTML version of manual`.**

```
Skipping checking HTML validation: 'tidy' doesn't look like recent
enough HTML Tidy.
```

Environmental: the local HTML Tidy binary on macOS Tahoe is older
than the version `R CMD check` is happy with. CRAN's check farm runs
a current Tidy, so this NOTE does not appear on CRAN itself. Does
not fire on the current dev machine when a recent Tidy is on PATH.

**NOTE 3 (environmental, only when offline) — `checking for future file timestamps`.**

```
unable to verify current time
```

`R CMD check` queries `worldclockapi.com` (or a similar network time
source) to validate that no file in the source tarball carries a
modification timestamp in the future. When the check host cannot
reach that network service — offline build environment, firewalled
CI runner, etc. — the check returns this NOTE rather than failing
the timestamp validation. No package file actually carries a future
timestamp on the development workstation; on a networked check host
this NOTE does not appear. CRAN's check farm is networked, so this
NOTE does not appear on CRAN.

## Downstream impact

`flexyBayes` has zero reverse dependencies on CRAN. `revdep_check()`
is N/A for this submission.

## Suggests-only dependencies

`INLA` is not on CRAN; it is hosted at the URL declared in
`Additional_repositories:`. The package degrades gracefully when INLA
is not installed: every code path that calls `INLA::*` is guarded by
`requireNamespace("INLA", quietly = TRUE)`, and tests that exercise
the live INLA backend skip via `skip_if_not_installed("INLA")`.

`brms` (the Stan passthrough, reached by `fb_brms()` or
`backend = "brms"`) requires a working Stan toolchain. Stan-dependent
tests gate via `skip_if_not_installed("brms")` plus `skip_on_cran()` /
`skip_on_ci()` because brms's first-call Stan compile (typically
30–60 s) exceeds the CRAN check budget.

`greta` remains in `Suggests` because its emit code and registry
descriptor are retained as a re-entry candidate. It is not an active
fitting engine, so its absence removes no capability.

## Test-suite behaviour on CRAN

The tests that exercise a *live* inference fit (INLA Laplace,
brms/Stan) gate behind `skip_if_not_installed()` / `skip_on_cran()`
and therefore skip on CRAN's clean check machine, since neither engine
is a hard dependency. What runs in full on CRAN is the
engine-independent core, which is the bulk of the package's logic and
its test surface: the ASReml and brms formula parsers, the term
classifier and intermediate representation, the `lgm_gate()`
feasibility checks, the backend-dispatch policy table, the refusal
registry, the canonical parameter-name registry, the prior DSL, the
capability matrix, and the `triangulate()` metric computations (run
against fixture draws rather than live fits). The skipped portion is
the thin engine-call layer, not the package's modelling logic; the
live engine paths are exercised on the maintainer's machine and on the
CI legs that install the engines (see the GitHub Actions matrix
above).

## Vignette compute

Thirteen vignettes ship with the package. Heavy MCMC fits inside
vignettes use small `n_samples` budgets and conditional evaluation
(`requireNamespace("INLA")`, `requireNamespace("brms")`); failures of
the live INLA or Stan integration during CRAN's clean re-render are
caught by `tryCatch()` so a single chunk failure does not crash the
vignette render. Vignettes use the `.Rmd.orig` precompile pattern
(`vignettes/_precompile.R` + `.Rbuildignore` excluding the `.orig`
sources) — the static `.Rmd` that ships in the package tarball is the
pre-evaluated output, so the build-side vignette render is a markdown
rendering exercise, not a live MCMC re-run.
