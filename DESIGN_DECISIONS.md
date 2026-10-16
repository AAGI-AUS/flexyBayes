# flexyBayes design decisions

flexyBayes records its main design decisions as numbered architecture
decision records (ADRs). Source comments and tests cite them by number,
for example `ADR 0019`, and this file lists the decision behind each
number. Entries in the refusal and representation registries record the
decision that introduced them in a `registered_in_adr` field, as an ADR
number or, for some entries, as the release and topic. The full records
are kept by the maintainers and are not distributed with the package.

An ADR number is never reused. A decision that replaces an earlier one
refers to it, and the earlier record is not edited.

| ADR | Decision |
|----|----|
| 0003 | Imports budget and the PC-prior default |
| 0004 | Uniform-on-SD default prior supersedes the PC-prior default |
| 0005 | Canonical parameter-name registry across the fitting engines (three at the time, two active today, see `NEWS.md`, 0.9.3) |
| 0006 | `backend = "auto"` argument and the backend decision trace |
| 0011 | `review` / `code = FALSE` planning options on the entry points |
| 0012 | Direct native-model entry semantics for the engine withdrawn in 0.9.3 (see `NEWS.md`) |
| 0013 | Backend slot for the engine withdrawn in 0.9.3 (scaffold only, see `NEWS.md`) |
| 0014 | [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) and brms-formula entry semantics |
| 0015 | Stan passthrough emit on the [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) backend |
| 0017 | Gate truth equals emit truth: INLA emit refusals are folded into the feasibility gate |
| 0018 | Smooth-basis retention on the internal model representation for [`predict()`](https://rdrr.io/r/stats/predict.html) |
| 0019 | Factor-continuous indexed-interaction emit |
| 0020 | Uncorrelated random slopes `(x || g)` |
| 0021 | Dataset and preflight contract |
| 0022 | Gaussian exact aggregation by sufficient statistics |
| 0023 | Chunked indexed prediction |
| 0024 | Backend routing trace and preflight metadata |
| 0025 | Known-covariance input formats |
| 0027 | Approximation registry contract and exactness conditions |
| 0029 | Triangulation independence-axis vocabulary |
| 0030 | Architectural contract: the seven truth conditions |
| 0031 | Backend registry: engine pins and the Stan / brms naming |

## The backend contract

The backend contract is the interface an inference engine must meet to
be a flexyBayes backend.
`tests/testthat/test-backend-conformance-open-oracle.R` checks both
engines against it, comparing each with base R
[`lm()`](https://rdrr.io/r/stats/lm.html) and with lme4 REML fits of the
same models.
