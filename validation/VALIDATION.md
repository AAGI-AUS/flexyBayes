# flexyBayes validation record

- package version: **0.9.0**
- commit: release tree tagged `v0.9.0` (squash of the 0.9.x development line; the granular audit trail is in the project ledger)
- generated: 2026-08-17
- archetype: A1 (statistical-method) · audit schema `rpkg-audit-v1`

Assembled by `validation/assemble_report.py` from the artefacts listed at the foot. Re-run any gate and regenerate; nothing here is hand-entered.

## Gates

| Gate | What it checks | Result | Detail |
| --- | --- | --- | --- |
| R CMD check --as-cran | built tarball, CRAN profile | PASS | `pass` |
| F19 live tests | devtools::test() actually run | PASS | 3238 pass / 0 fail / 199 skip |
| F16 vibe hygiene | E13-E20 detector suite | FAIL | score 0.8571 -- 6 pass / 1 fail (all should_fix, no blocker) |
| F17 conformance | pkg-validation conformance evidence | na | absent |
| Audit blockers | I1-I14 + coherence C1-C3 | PASS | 0 blocker / 5 should-fix / 15 nice-to-have |

### Blockers

None.

### Elegance and vibe-hygiene detectors

None of these is a release blocker; all are graded `should_fix` or lower, and the failures are recorded rather than waived. E10's three hits are a detector false positive -- the flagged comments explain why a code path exists, which is the opposite of decorative.

| Detector | Result | Evidence |
| --- | --- | --- |
| E13 Comment-paraphrase ratio | PASS | 22/907 paraphrasing; e.g. emit_count_aggregated.R:509, plot.R:257, plot.R:290 |
| E14 Validation-scattering index | PASS | delegated 0.80 (131/163 delegated+duplicated) ≥ floor 0.75 for A1; 229 single-use self-validation guard(s) exempt |
| E15 Roxygen sentence completeness | PASS | 497/497 full sentences |
| E16 LLM-filler density | PASS | 0 fillers / 85338 words = 0.0/1k |
| E17 Banner coverage | PASS | 52/73 long units bannered |
| E18 Error-message actionability | PASS | 54/54 actionable |
| E19 Human-pass evidence | FAIL | no De-Vibe / human-pass record (critical_review_*.md / de_vibe*.md / report/internal/ / NEWS / recipe 40) |
| E20 Conformance coverage of generated surfaces | N/A | no evidence.json (no pkg-validation conformance study; looked under manifest, package-root, and sibling report/) |

## Quality ladder

- current rung: **pilot**  (target releasable)
- next driver: raise Vibe Hygiene to >= 90% (recipe 40 De-Vibe pass)
- why: hygiene floor gates 'shaped'

The rung is capped by the vibe-hygiene floor, not by correctness: audit blockers are 0 and every numerical gate below passes. Raising it needs a De-Vibe pass (recipe 40), and E19 in particular asks for a human review record, which is not something this report can supply for itself.

**E19 candidates, pending ratification.** Two external review passes over the 0.9.0 development line exist in the project's private review archive: a structured review briefing and an independent adversarial release review whose findings were verified and closed before this release. Their conclusions and the decisions they forced are recorded in the project decision ledger (cairn `2026-08-17-release-nine-review-record`). They are recorded as candidates only: the detector stays FAIL until a maintainer ratifies a review as the human pass, and nothing in this file should be read as a completed human pass.

## Numerical validation, against independent oracles

Each row is a claim checked against a reference the package did not author. `n` is the number of paired comparisons.

| Claim | Oracle | n | Statistic | Result |
| --- | --- | --- | --- | --- |
| Missing-plot device equals the observed-data REML (p <= 20% missing) | independently written dense-V REML | 4000 | 0 exceed 1e-2, largest 0.0097. | PASS |
| ...same, at 30% missing (estimand degrades, not the device) | as above | 7125 | 0.6% exceed 1e-2 | reported |
| ...same, at 50% missing (estimand degrades, not the device) | as above | 1200 | 1.8% exceed 1e-2 | reported |
| ...same, at 70% missing (estimand degrades, not the device) | as above | 1185 | 14.3% exceed 1e-2 | reported |
| diag(site):gen == brms (0 + site || gen) | ASReml + independent lme4 expansion | 3 arms | REML arms agree to 3.51e-05 | PASS |
| ...and all three fit k variances, 0 covariances | parameter count | 3 arms | diagonal structure confirmed | PASS |
| Head-to-head on 5 models (see head-to-head.pdf) | ASReml 4.2 | 17 components | parameter counts equal on every model | PASS |

## Computational cost

| Item | Setting | Wall time |
| --- | --- | --- |
| Simulation sweep | 170 distinct cells x 25 seeds x 5 methods (185 runs incl. slice overlap) | 6.0 h |
| ASReml diagnostics pass | 4250 refits | 1.4 min |

Per-model head-to-head timings are in `head-to-head/head-to-head.pdf` Table 2. Summary: ASReml returns a REML point estimate in 0.02-0.82 s; flexyBayes returns a posterior in 2.9-3.6 s through INLA and 25-37 s through brms.

## Artefacts this reads

- `validation/audit2.json` -- rpkg audit, machine-readable (yes)
- `validation/audit2.md` -- rpkg audit, scorecard (yes)
- `design-preserving-missingness/results/report_sweep.txt` -- simulation sweep analysis (yes)
- `design-preserving-missingness/results/report_diagnostic.txt` -- variance-component reliability diagnostic (yes)
- `design-preserving-missingness/results/oracle_het.log` -- heterogeneity three-arm oracle (yes)
- `head-to-head/results/comparison.rds` -- head-to-head estimates (yes)
- `head-to-head/head-to-head.pdf` -- head-to-head report (PDF) (yes)
- `head-to-head/head-to-head.html` -- head-to-head report (HTML) (yes)

## Known non-blocking notes

- `R CMD check --as-cran` reports one NOTE: new submission, and `Suggests`/`Enhances` not in a mainstream repository. Both are expected for a package not yet on CRAN whose optional backends live on r-universe.
- 199 tests skip. They are the greta / gretaR quarantine guards and backend-availability guards, not silent omissions.
