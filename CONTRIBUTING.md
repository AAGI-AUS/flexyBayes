# Contributing to flexyBayes

Thanks for your interest in contributing. The package is an `AAGI-AUS`
collective effort; we welcome bug reports, fixes, documentation
improvements, and discussion of new features.

## Bug reports and feature requests

Please use the GitHub issue tracker:
<https://github.com/AAGI-AUS/flexyBayes/issues>.

Useful information to include:

- a minimal reproducible example;
- the output of
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html);
- the output of `packageVersion("flexyBayes")` and, if relevant,
  `packageVersion("greta")` and `packageVersion("INLA")`;
- whether the issue reproduces on a fresh R session.

## Pull requests

Small fixes (typos, doc clarifications, minor bug fixes) are welcome as
direct PRs against `main`. For non-trivial changes, please open an issue
first to discuss the approach.

### Development setup

``` r

# clone the repo, then in R:
install.packages(c("devtools", "testthat", "roxygen2"))
devtools::load_all()
devtools::test()
devtools::check()
```

For release-grade checks, **build the tarball first and check the
tarball** rather than running `R CMD check` against the source
directory:

``` sh
R CMD build flexyBayes
R CMD check --as-cran --no-manual flexyBayes_*.tar.gz
```

`R CMD check` against the raw source directory fails immediately with
`Required fields missing or empty: 'Author' 'Maintainer'`. The `Author`
/ `Maintainer` fields are materialised from `Authors@R` by
`R CMD build`; this is the canonical check path for the project.
`devtools::check()` calls `R CMD build` internally, so the R-side
workflow is unaffected.

The package depends on `greta` (default MCMC backend) and `INLA`
(Laplace backend). `greta` requires a working Python + TensorFlow
install — see `greta::install_greta_deps()`. `INLA` is hosted at
<https://inla.r-inla-download.org/R/stable>.

### Code style

- snake_case for function and argument names.
- Australian / British English in user-facing prose (`optimise`,
  `behaviour`, `parameterise`, `centred`).
- ISO-8601 dates (`YYYY-MM-DD`).
- Two-space indentation; line width ≤ 80 chars where practicable.

### Tests

Add tests under `tests/testthat/`. Tests that depend on `greta`, `INLA`,
or `brms` should `skip_if_not_installed(.)`. The test suite uses
testthat edition 3.

#### Three-tier test discipline

The test suite is organised into three implicit tiers via
`testthat::skip_*` guards. Pick the tier that matches what you can warm
in your environment.

**Tier 1 — CRAN-fast (no heavy engines).**

``` r

devtools::test()           # ~30 s; engine-free path only
```

Guards: `skip_on_cran()`, `skip_on_ci()`, `skip_if_not_installed()`.
Covers IR parsing, prior DSL, refusal templates, registry lookups,
canonical-name transforms, dispatch trace shape, review-code workflow,
gretaR dormant scaffold. Expected: roughly **PASS 700+ / FAIL 0** with
many SKIPs.

**Tier 2 — local integration (greta + INLA warm; brms gated).**

``` r

Sys.setenv(NOT_CRAN = "true")
devtools::test()           # ~3-5 min on a warm TF backend
```

Same suite, but `skip_on_cran()` is now `FALSE`, so greta + INLA
round-trip tests fire. `skip_if_not_installed("brms")` still gates the
Stan passthrough tests. Expected: **PASS 850+ / FAIL 0** with a small
handful of SKIPs (brms round-trip, vdiffr).

**Tier 3 — full triangulation (greta + INLA + brms / Stan).**

``` r

Sys.setenv(NOT_CRAN = "true")
Sys.setenv(FLEXYBAYES_RUN_BRMS = "true")  # opt-in flag
devtools::test()           # ~10-20 min; brms first-call Stan
                           # compile is 30-60 s per backend test
```

Adds the Stan passthrough round-trips and three-engine triangulation
tests on top of Tier 2. The gate is **FAIL 0**; PASS counts grow across
the 0.8.x line (the 90+-file suite is well above the early-release
floors quoted above). The residual SKIPs are vdiffr snapshots and
asreml-dependent reference comparisons.

#### Full `R CMD check --as-cran` rehearsal

``` bash
cd /path/to/flexyBayes_dev
R CMD build flexyBayes
_R_CHECK_CRAN_INCOMING_REMOTE_=true R CMD check --as-cran \
  flexyBayes_<version>.tar.gz
```

Expected: **Status: 2 NOTEs** baseline preserved (CRAN-incoming
feasibility with URL 404 sub-item gated on the upstream rename
conversation; HTML Tidy environmental). A third NOTE
(`unable to verify current time`) fires on offline / firewalled check
hosts where `R CMD check` cannot reach a network time source; it does
not appear on networked machines or on CRAN. See `cran-comments.md` for
the per-NOTE itemisation.

#### Re-precompile vignettes

When a `.Rmd.orig` source changes — or when the `DESCRIPTION` version
bumps and you want the
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) chunks to
refresh — re-precompile via:

``` bash
R CMD build flexyBayes
R CMD INSTALL flexyBayes_<version>.tar.gz
cd flexyBayes && Rscript vignettes/_precompile.R
```

Pre-requisites: greta + TF warm, INLA installed, brms installed. The
driver knits each `.Rmd.orig` into its sibling `.Rmd`; failures are
reported per-vignette and the script exits non-zero on any. Expected
wall-time on a warm M1 / M2 Mac: 10-15 minutes for the full 16-vignette
deck.

### Documentation

Roxygen2 (`@param`, `@return`, `@examples`, `@family`) is the source of
truth for `man/`; do not hand-edit `.Rd` files. Long-running examples
should use `\donttest{}`, not `\dontrun{}`.

### Vignettes

Vignettes follow the brms / tidymodels register: motivation → theory →
implementation → interpretation → pitfalls → extensions. See
`vignettes-stage1-plan.md` (workspace-level) for the deck-wide
conventions.

## Community

The project follows the [Contributor Covenant Code of
Conduct](https://aagi-aus.github.io/flexyBayes/CODE_OF_CONDUCT.md). By
participating, you agree to abide by its terms.

For security-relevant disclosures, please follow
[SECURITY.md](https://aagi-aus.github.io/flexyBayes/SECURITY.md) rather
than the public issue tracker.
