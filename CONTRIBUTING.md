# Contributing to flexyBayes

flexyBayes is an AAGI-AUS collective effort. Bug reports, fixes,
documentation improvements and discussion of new features are welcome.

## Bug reports and feature requests

Please use the GitHub issue tracker:
<https://github.com/AAGI-AUS/flexyBayes/issues>.

A report is easiest to act on when it includes:

- a minimal reproducible example
- the output of `sessionInfo()`
- the output of `packageVersion("flexyBayes")` and, where relevant,
  `packageVersion("INLA")` or `packageVersion("brms")`
- whether the issue reproduces in a fresh R session

## Pull requests

Small fixes (typos, documentation, minor bugs) are welcome as direct pull
requests against `main`. For anything larger, please open an issue first
so the approach can be agreed before the work is done. A pull request that
closes an issue says so in its body with `Fixes #<issue-number>`.

### Development setup

```r
# clone the repository, then in R:
install.packages(c("devtools", "testthat", "roxygen2"))
devtools::load_all()
devtools::test()
devtools::check()
```

For a release-grade check, build the tarball first and check the tarball:

```sh
R CMD build flexyBayes
R CMD check --as-cran --no-manual flexyBayes_*.tar.gz
```

`R CMD check` on the raw source directory stops with
`Required fields missing or empty: 'Author' 'Maintainer'`, because those
fields are written from `Authors@R` by `R CMD build`. `devtools::check()`
builds first, so the R-side workflow is unaffected.

The package fits through two engines, INLA (Laplace approximation) and
brms (Stan MCMC, on CRAN). INLA is installed from
<https://inla.r-inla-download.org/R/stable>, which `DESCRIPTION` names in
`Additional_repositories`. Both engines are optional and the package
installs without either.

### Code style

- snake_case for function and argument names.
- Australian/British English in user-facing prose (`optimise`,
  `behaviour`, `parameterise`, `centred`).
- ISO 8601 dates (`YYYY-MM-DD`).
- Formatting follows `air.toml` (two-space indentation, 80-character
  lines). Run `air format .` before committing, and leave code you did not
  change as it is.
- `lintr::lint_package()` reads `.lintr` and should return nothing new.

### Tests

Tests live under `tests/testthat/` and use testthat edition 3. A test that
needs INLA or brms starts with `skip_if_not_installed()`. Files run
sequentially (`Config/testthat/parallel: false`), because INLA and brms
fits do not share a process safely.

The suite has two tiers. `devtools::test()` and `devtools::check()` run
the short tier, which takes about 15 minutes on a ten-core Apple Silicon
machine with INLA held to one thread. The blocks that dominate the wall
clock run only when `FLEXYBAYES_RUN_LONG` is set:

```sh
# the short tier, the default
Rscript -e 'devtools::test()'

# every tier, as a release runs it
FLEXYBAYES_RUN_LONG=true Rscript -e 'devtools::test()'
```

A block joins the long tier by calling `skip_if_not_long_tier()` after
its `skip_on_cran()` and `skip_if_not_installed()` guards. Choose a block
for the long tier on its measured time in the short tier as it runs,
because several files build an expensive fixture in one block and reuse
it in later ones. Stress tests on very large inputs run only when
`FLEXYBAYES_RUN_STRESS=true`.

Every release runs every tier. `tools/run_test_suite.R` records one row
per block with the reason it was skipped, and criterion D11 of
`tools/definition_of_done.R` fails when a block was held back, so a
release cannot be evidenced by a short-tier run.

### `R CMD check --as-cran`

```sh
R CMD build flexyBayes
_R_CHECK_CRAN_INCOMING_REMOTE_=true R CMD check --as-cran flexyBayes_<version>.tar.gz
```

The expected result is one NOTE, on CRAN incoming feasibility: a new
submission, INLA from an additional repository, and the documentation
URL until the site is published. On a machine that cannot reach a network
time service a second NOTE, `unable to verify current time`, can appear.
It is a property of that host. `cran-comments.md` lists the NOTE for the
current release.

### Vignettes

The six vignettes ship pre-evaluated. Each `.Rmd.orig` source holds live
model fits and is knitted into its sibling `.Rmd` by the maintainer, so
building the package renders static markdown and fits nothing. When a
`.Rmd.orig` changes, or the version bumps and the `sessionInfo()` chunks
should refresh:

```sh
R CMD build flexyBayes
R CMD INSTALL flexyBayes_<version>.tar.gz
cd flexyBayes && Rscript vignettes/_precompile.R
```

This needs INLA and brms installed with a working Stan toolchain. The
driver reports each failure per vignette and exits non-zero on any.
`Rscript vignettes/_precompile.R --only 01` refreshes one vignette.

### Documentation

roxygen2 (`@param`, `@return`, `@examples`, `@family`) is the source of
truth for `man/`. Do not hand-edit `.Rd` files. Long-running examples use
`\donttest{}`, not `\dontrun{}`.

## Community

The project follows the [Contributor Covenant Code of
Conduct](CODE_OF_CONDUCT.md). By participating, you agree to abide by its
terms.

For security disclosures, please follow [SECURITY.md](SECURITY.md) instead
of the public issue tracker.
