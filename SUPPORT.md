# Getting help with flexyBayes

Thanks for using `flexyBayes`. The fastest path to a useful answer
depends on what kind of question you have.

## Where to ask

| Question type | Where |
|----|----|
| “How do I do X with flexyBayes?” | [GitHub Discussions — Q&A](https://github.com/AAGI-AUS/flexyBayes/discussions/categories/q-a) |
| “I think this is a bug.” | [GitHub Issues — bug report](https://github.com/AAGI-AUS/flexyBayes/issues/new?template=bug_report.yml) — please include a [reprex](https://reprex.tidyverse.org). |
| “I’d like a new feature / backend / capability.” | [GitHub Issues — feature request](https://github.com/AAGI-AUS/flexyBayes/issues/new?template=feature_request.yml) |
| Security-sensitive disclosure | See [`SECURITY.md`](https://aagi-aus.github.io/flexyBayes/SECURITY.md) — please do **not** open a public issue. |
| General Bayesian-mixed-model methodology | [Stan Forums](https://discourse.mc-stan.org/) or [INLA mailing list](https://www.r-inla.org/contact-us) — broader and more active than this project. |

## Maintenance capacity

`flexyBayes` is maintained by a small team within the Australian
Agricultural Genomics Institute (AAGI-AUS) collective. The capacity
block below is updated alongside each release.

> **Note.** The bus-factor / SLA / abandonment-protocol entries below
> are pending the AAGI-AUS collective conversation that also governs
> project decisions D5 (paper coauthors), D8 (upstream rename), and D10
> (paper venue). Until that conversation lands, do not infer commitments
> from this file beyond the items already labelled with concrete values.
> This file is `.Rbuildignore`’d and ships only on the GitHub repository
> / pkgdown site, not in the CRAN tarball.

| Field | Current value | Notes |
|----|----|----|
| Primary maintainer | Max Moldovan | Adelaide University; contact via package author email. |
| Bus factor | Pending AAGI-AUS collective decision | Will be set jointly with D5 / D8 / D10. |
| Weekly maintenance hours | Pending AAGI-AUS collective decision | Honest median across recent quarters will be reported once measured. |
| Issue response SLA | Pending AAGI-AUS collective decision | No SLA committed under the current development release (`0.8.3`). |
| Critical security SLA | 7 days | See [`SECURITY.md`](https://aagi-aus.github.io/flexyBayes/SECURITY.md). |
| Abandonment protocol | Pending AAGI-AUS collective decision; default fallback = AAGI-AUS institutional handover | Concrete handover / archive plan ratified with D5 / D8. |

## Reprex etiquette

A reproducible example (reprex) lets us help you in one round-trip
instead of three. The minimum:

1.  A small `data.frame` we can copy-paste.
2.  The exact
    [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
    /
    [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
    call.
3.  The full output of
    [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) —
    especially the greta / TensorFlow / INLA versions.
4.  The error / unexpected output, copied verbatim.

For runtime issues, please run with `verbose = TRUE` and include the
emitted greta code (or, for `backend = "inla"`, the printed INLA
formula).

## What we cannot help with

- General R help — see [Stack Overflow
  `[r]`](https://stackoverflow.com/questions/tagged/r).
- ASReml-R proprietary syntax outside the subset that `flexyBayes`
  implements — see VSNi support.
- INLA model-tuning beyond the `lgm_gate()` / `priors_to_inla()` surface
  — see [r-inla.org](https://www.r-inla.org/).
- greta’s own Python / TensorFlow installation — see
  [greta-dev/greta](https://github.com/greta-dev/greta).
