## Submission

This is a new submission of flexyBayes 0.18.1.

flexyBayes fits mixed models written in ASReml-style or brms-style formula
syntax with one of two inference engines, INLA or brms (Stan). The fitted
object has the same accessors whichever engine ran. A model structure that
neither engine can represent is refused by name, and the refusal names the
nearest structure the package does fit.

## Test environments

- Local: macOS Tahoe 26.6.2, aarch64-apple-darwin23 (Apple Silicon),
  R 4.6.1 (2026-06-24). `R CMD check --as-cran` on the built tarball, twice
  with INLA, brms and cmdstanr installed and once with all three hidden.

Windows R-devel (`devtools::check_win_devel()`) and the R-hub
multi-platform check have not yet been run on this version. They will be run
immediately before submission.

## R CMD check results

0 errors | 0 warnings | 1 note

All three local checks gave the same NOTE:

```
* checking CRAN incoming feasibility ... [3s/22s] NOTE
Maintainer: ‘Max Moldovan <max.moldovan@adelaide.edu.au>’

New submission

Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/stable

Found the following (possibly) invalid URLs:
  URL: https://aagi-aus.github.io/flexyBayes/
    From: DESCRIPTION
          man/flexyBayes-package.Rd
          README.md
          NEWS.md
    Status: 404
    Message: Not Found
```

It has three parts.

**New submission.**

**INLA is not in a mainstream repository.** INLA is distributed from its own
repository, `https://inla.r-inla-download.org/R/stable`, which is declared in
`Additional_repositories:` and which the check resolves. INLA is in
`Suggests`, not `Imports`, and the package does not require it. Examples
and tests that need INLA are guarded, so without it the examples run, those
tests skip and the vignettes build. Requesting the INLA engine without INLA
installed stops with a message naming the repository to install it from.
brms is on CRAN, so one engine can always be installed.

**A 404 URL.** The documentation site named in `DESCRIPTION` is not yet
published. It will resolve before the package is submitted.

A local check can also report `checking for future file timestamps ...
unable to verify current time` when the time server cannot be reached from
the maintainer's network. It did not occur in the three checks above.

## Cores

`Config/testthat/parallel` is `false`, so test files run in one process.
INLA fits run on one thread (`num.threads = "1:1"`), and the accessors that
draw from an INLA posterior also use one thread. By default each INLA fit is
run three times, once in the session and twice in separate R processes
started through callr, and these run at the same time. While `R CMD check`
limits cores (`_R_CHECK_LIMIT_CORES_` set), each fit runs once, so the
examples use one process. The test suite sets
`options(flexyBayes.inla_best_of = 1)`, so each test fits once. flexyBayes
passes no `cores` argument to `brms::brm()` unless the user supplies one,
and the tests that set `mc.cores` set it to 2.

## Reverse dependencies

None. This is a new submission.

## Vignettes

The six vignettes ship pre-computed. Their model fits are run by the
maintainer with `vignettes/_precompile.R`, and the tarball contains the
resulting markdown with the output already in it. Building the vignettes at
check time fits no model, so it needs neither INLA nor a Stan toolchain.
