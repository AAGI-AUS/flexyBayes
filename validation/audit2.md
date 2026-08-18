# R Package Audit: flexyBayes

## Package
- Package: flexyBayes
- Version: 0.9.1
- License: MIT + file LICENSE

## R CMD check --as-cran
- Status: **Verified-Pass**
- Log: `/var/folders/l4/dcx9vkrx2bj0w0kxfvlyrc7r0000gp/T/rpkg_audit_2eozaup8/flexyBayes.Rcheck/00check.log`
- ERRORs: 0 · WARNINGs: 0 · NOTEs: 1

## Blocker (0)
- None

## Should-fix (5)
- **Tests** `tests/testthat/test-tooling-skip-ledger.R`: Trivial assertion detected.
- **Security** `tests/testthat/test-fb-preflight-stress.R`: Audit shell calls for user-input injection risk.
- **Security** `R/emit_gretaR.R`: Audit shell calls for user-input injection risk.
- **Security** `R/fb_preflight.R`: Audit shell calls for user-input injection risk.
- **Security** `R/emit_greta.R`: Avoid eval(parse()); it is hard to secure and test.

## Nice-to-have (15)
- **Tests** `tests/testthat/test-tidy-hub.R`: 1 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-gretaR-slot.R`: 1 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-inla-verification-simple-slope-uncor.R`: 1 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-validation-lmer.R`: 4 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-emit-smooth-low-rank.R`: 1 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-fb-plan.R`: 1 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-known-covariance-inputs.R`: 2 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-emit-gaussian-aggregated.R`: 2 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-factor-continuous-inla-verification.R`: 2 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-tooling-skip-ledger.R`: 2 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-backend-conformance-open-oracle.R`: 6 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-emit-inla.R`: 1 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-backend-conformance.R`: 6 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Tests** `tests/testthat/test-fb-preflight-stress.R`: 1 bare skip() call(s); a manual unconditional skip can hide core behaviour (the conditional skip_*() guards are exempt).
- **Build-check** `/var/folders/l4/dcx9vkrx2bj0w0kxfvlyrc7r0000gp/T/rpkg_audit_2eozaup8/flexyBayes.Rcheck/00check.log`: R CMD check NOTE (env): * checking CRAN incoming feasibility ... [3s/33s] NOTE Maintainer: ‘Max Moldovan <max.moldovan@adelaide.edu.au>’ New submission Suggests or Enhances not in main

## Human Review Questions
- What existing package or workflow does this improve on?
- What behavior would be dangerous if wrong?
- Which tests would fail if the core assumption were false?
- Who will maintain the package and answer issues?
- Are license, ownership, and citation decisions approved?

## Live test suite — Verified-Pass

- live devtools::test() green: 4298 PASS / 199 SKIP / 0 FAIL. 243 skip_on_cran() site(s) — a cached --as-cran log would hide these.

## Vibe Hygiene Model (E13–E20, v1.0.0) — flexyBayes

| Row | Status | Severity | Evidence | Note |
|---|---|---|---|---|
| E13 Comment-paraphrase ratio | **PASS** | should_fix | R source | 23/1038 paraphrasing; e.g. emit_count_aggregated.R:532, plot.R:335, plot.R:379 |
| E14 Validation-scattering index | **PASS** | should_fix | R source | delegated 0.80 (138/172 delegated+duplicated) ≥ floor 0.75 for A1; 237 single-use self-validation guard(s) exempt |
| E15 Roxygen sentence completeness | **PASS** | should_fix | R source | 539/539 full sentences |
| E16 LLM-filler density | **PASS** | should_fix | R source | 0 fillers / 101975 words = 0.0/1k |
| E17 Banner coverage | **PASS** | should_fix | R source | 57/79 long units bannered |
| E18 Error-message actionability | **PASS** | should_fix | R source | 54/54 actionable |
| E19 Human-pass evidence | **FAIL** | should_fix | flexyBayes | no De-Vibe / human-pass record (critical_review_*.md / de_vibe*.md / report/internal/ / NEWS / recipe 40) |
| E20 Conformance coverage of generated surfaces | **N/A** | should_fix | flexyBayes | no evidence.json (no pkg-validation conformance study; looked under manifest, package-root, and sibling report/) |

**Tally**: FAIL=1, N/A=1, PASS=6

**Vibe Hygiene Score** (Pass / (Pass + Fail), excluding N/A): 6/7 = **86%** (target ≥ 90%)

rpkg-audit-v1 JSON -> flexyBayes/validation/audit2.json (blockers=0, F17=na, F19=pass, collect=full)
