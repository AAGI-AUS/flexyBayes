## Submission

This is a new submission of `flexyBayes`, version 0.9.1.

The package specifies mixed models in ASReml-style or brms-style formula
syntax and estimates them through one of two inference engines, INLA or
brms (a Stan passthrough), returning one object whose accessors behave
the same whichever engine ran. Model structures neither engine can
represent are refused by name, with the nearest implemented alternative
given in the refusal, rather than silently translated into something
else.

## Test environments

- local: macOS Tahoe 26.5.2, aarch64-apple-darwin20 (Apple Silicon),
  R 4.5.2 (2025-10-31). `R CMD check --as-cran` on the built tarball.

Windows R-devel (`devtools::check_win_devel()`) and the R-hub
multi-platform check are to be run by the maintainer immediately before
submission. They are not reported here because they had not been run
when this file was written, and a check result is not worth quoting
unless it happened.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the CRAN incoming-feasibility note, in full:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Max Moldovan <max.moldovan@adelaide.edu.au>'

New submission

Suggests or Enhances not in mainstream repositories:
  INLA
Availability using Additional_repositories specification:
  INLA   yes   https://inla.r-inla-download.org/R/stable

Found the following (possibly) invalid URLs:
  URL: https://aagi-aus.github.io/flexyBayes/
    From: DESCRIPTION
          man/flexyBayes-package.Rd
    Status: 404
    Message: Not Found
```

It carries three items.

**New submission.** Expected, and it resolves on acceptance.

**INLA is not in a mainstream repository.** INLA is distributed from its
own repository at `https://inla.r-inla-download.org/R/stable`, which is
declared in `Additional_repositories:` and which the check resolves
successfully on the line above. This is the established arrangement for
packages that offer INLA as an inference engine, and INLA is in
`Suggests`, never in `Imports`.

The package does not require it. Every call into INLA sits behind
`requireNamespace("INLA", quietly = TRUE)`, and each of the three emit
entry points opens with an installation check that stops with a typed
message naming the repository to install from. Requesting the INLA
backend without INLA installed therefore gives a refusal that says what
to do, not a "there is no package called 'INLA'" error from somewhere
inside the call stack. Every example, test and vignette path that would
reach INLA is guarded, so on a machine without it the examples run, the
INLA-dependent tests skip, and the vignettes build. brms is on CRAN, so
a reviewer always has one installable engine.

`greta` is in `Suggests` as a dormant engine and installs from CRAN:
it fits nothing in this version, `backend = "auto"` never selects it,
and an explicit request for it is refused by name.

**A 404 URL.** The documentation site declared in `DESCRIPTION` does not
resolve yet, because the repository it is built from is not yet public.
The maintainer will not submit until that URL resolves. It is quoted
here so the check output above is reproduced faithfully rather than
edited.

## Cores

`Config/testthat/parallel` is `true`. testthat's own worker default is
a hard-coded two (verified against the installed testthat 3.3.2
sources), which is within the CRAN two-core limit. Nothing in the
package calls `parallel::detectCores()`, sets `mc.cores`, or opens a
cluster, and both calls into INLA's posterior sampler pin
`num.threads = "1:1"`. No sampler is given a `cores` argument anywhere
in the tests, examples or vignettes, so brms runs on its default of one.

## Reverse dependencies

None. This is a new submission.

## Vignettes

The twelve vignettes ship pre-evaluated. Their sources carry live model
fits, and those are run by the maintainer through
`vignettes/_precompile.R`; what the tarball contains is the resulting
static markdown, with the fitted output already in it. Building the
vignettes at check time renders markdown and fits nothing, so the
vignette build needs neither INLA nor a Stan toolchain.
