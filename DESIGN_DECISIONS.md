# flexyBayes design decisions

flexyBayes records its architecture decisions as **ADRs** (Architecture
Decision Records). Each structural refusal and each representation in
the package carries a `registered_in_adr` tag naming the ADR that
introduced it, so a reader can trace any behaviour back to the decision
(and the alternatives weighed) that produced it. Source comments and
tests sometimes cite an ADR by number (for example `ADR 0019`); this
file is the index that resolves those numbers. The full decision records
are maintained in the project’s decision ledger.

The numbering is stable: an ADR number is never reused, and a
superseding decision references the one it replaces rather than editing
it.

| ADR | Decision |
|----|----|
| 0003 | Imports budget and the PC-prior default |
| 0004 | Uniform-on-SD default prior supersedes the PC-prior default |
| 0005 | Canonical parameter-name registry across greta / INLA / brms |
| 0006 | `backend = "auto"` argument and the backend decision trace |
| 0011 | `review` / `code = FALSE` planning options on the entry points |
| 0012 | [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md) and direct greta-model entry semantics |
| 0013 | gretaR backend slot (scaffold / provision-only) |
| 0014 | [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) and brms-formula entry verb semantics |
| 0015 | Stan passthrough emit on the [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md) backend |
| 0017 | Gate-truth = emit-truth: fold INLA emit refusals into the feasibility gate |
| 0018 | Smooth-basis retention on the IR for [`predict()`](https://rdrr.io/r/stats/predict.html) |
| 0019 | Factor-continuous indexed-interaction emit |
| 0020 | Uncorrelated random slopes `(x || g)` |
| 0021 | `fb_dataset()` / `fb_preflight()` contract |
| 0022 | Gaussian exact aggregation by sufficient statistics |
| 0023 | Chunked indexed prediction |
| 0024 | Backend routing-trace and preflight metadata |
| 0025 | Known-covariance input formats |
| 0027 | Approximation registry contract and exactness conditions |
| 0029 | Triangulation independence-axis vocabulary |
| 0030 | Architectural contract: the seven truth-conditions |
| 0031 | Backend registry: engine pins and the Stan / brms naming |

A small number of registry entries carry the placeholder `9999` for a
behaviour whose decision record has not yet been ratified.

## Contracts

A few cross-cutting interfaces are referred to by a short contract code
in the code and tests:

- **the backend contract** – the interface any inference engine must
  honour to be a first-class flexyBayes backend (the conformance battery
  in `tests/testthat/test-backend-conformance.R` is its executable
  form);
- **the surrogate / ensemble / data-source contracts** – the seams
  through which flexyBayes consumes externally-produced inputs without
  depending on a particular producer.

The guiding principle throughout is *depend on a versioned contract,
never on a producer’s internals*.
