# Getting help with flexyBayes

## Where to ask

| Question | Where |
|---|---|
| How do I do X with flexyBayes? | [GitHub Issues](https://github.com/AAGI-AUS/flexyBayes/issues) |
| I think this is a bug. | [Bug report](https://github.com/AAGI-AUS/flexyBayes/issues/new?template=bug_report.yml), with a [reprex](https://reprex.tidyverse.org) |
| I would like a new feature, engine or model structure. | [Feature request](https://github.com/AAGI-AUS/flexyBayes/issues/new?template=feature_request.yml) |
| A security problem | See [`SECURITY.md`](SECURITY.md). Please do **not** open a public issue. |
| Bayesian mixed models in general | [Stan Forums](https://discourse.mc-stan.org/) or the [INLA discussion group](https://www.r-inla.org/contact-us), which are larger and more active than this project. |

## Maintenance

flexyBayes is maintained by a small team within Analytics for the Australian
Grains Industry (AAGI). This file is in the GitHub repository and on the
documentation site, not in the package tarball.

| Item | Current value |
|---|---|
| Primary maintainer | Max Moldovan, Adelaide University (email in `DESCRIPTION`) |
| Number of maintainers | Not yet set |
| Weekly maintenance hours | Not yet set |
| Issue response time | None committed for 0.18.1 |
| Security response time | Acknowledgement within 7 days (see [`SECURITY.md`](SECURITY.md)) |
| If the package is no longer maintained | Not yet set. The default is handover to AAGI. |

## A good reproducible example

Please include:

1. A small `data.frame` that can be copied and pasted.
2. The exact `flexybayes()` call.
3. The full output of `sessionInfo()`, including the INLA, brms and Stan
   versions.
4. The error or unexpected output, copied as printed.

For a problem during fitting, run with `verbose = TRUE` and include the
generated model code it prints: the Stan program for `backend = "brms"`, or
the INLA formula for `backend = "inla"`.

## What we cannot help with

- General R questions. See [Stack Overflow `[r]`](https://stackoverflow.com/questions/tagged/r).
- ASReml-R syntax outside the subset flexyBayes implements. See VSNi
  support.
- Tuning INLA beyond the settings flexyBayes exposes. See
  [r-inla.org](https://www.r-inla.org/).
- Installing a Stan toolchain for the brms engine. See
  [mc-stan.org](https://mc-stan.org/) and the `rstan` getting-started guide.
