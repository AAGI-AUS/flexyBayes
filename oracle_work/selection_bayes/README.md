# Feasibility probes for the Bayesian selection module

These are the scripts that decided the architecture of
`SPEC_bayesian_selection_v5_2026-09-12.md`. They were run before the
specification was written, and the section-13 table of that document is their
result. Nothing in the contract rests on a capability that was not executed
here.

**Run 2026-09-12** on an Apple-silicon macOS machine, under R 4.6.1,
brms 2.23.0, CmdStan 2.38.0 via cmdstanr 0.9.0, projpred 2.10.0, glmnet 5.0,
posterior 1.7.0.

The scripts are copied verbatim from the session in which they were run. They
are evidence, not package code: they are not sourced by anything, they are not
part of `R CMD check`, and they are kept so that a reader can re-run the
measurement rather than take the table's word for it.

## What each script establishes

### `f0_probe.R` -- what brms 2.23.0 permits (probes A1-A8)

| Probe | Result |
|---|---|
| A1 class-wide horseshoe with `(1 \| g)` | compiles; non-centred `zb .* sdb` with `hs_local`, `hs_global`, `hs_slab` |
| A2 coefficient-specific `double_exponential` on candidates plus `normal` on a protected coefficient | compiles; exact per-coefficient targeting |
| A3 class-wide horseshoe plus a protected coefficient row | refused by brms |
| A4 non-linear-parameter partition, horseshoe on `cand`, normal on `prot`, `(1 \| g)` | compiles and samples |
| A5 class-wide R2D2 | compiles; `R2D2_phi` simplex, `R2D2_R2` beta |
| A6 elastic net as `target +=` over `b[cand_idx]` via `stanvars` and `check = FALSE` | compiles; the protected coefficient keeps its own row; samples cleanly (0 divergences, R-hat 1.00, bulk ESS above 1800) |
| A7 horseshoe with `(1 \| gr(g, cov = K))` | compiles |
| A8 compile-and-sample smoke test, 60 rows | 10 s wall |

A3 is the reason `regularized_horseshoe()` and `r2d2()` refuse a model with
protected fixed effects beyond the intercept, and A2 and A6 are the reason
`ridge()`, `laplace()` and `elastic_net()` do not.

### `bridge_probe.R`, `bridge_probe2.R`, `map_probe.R` -- the optimiser conventions (probes B1-B4)

| Probe | Result |
|---|---|
| B1 ridge: closed form against `glmnet(alpha = 0, lambda = 1 / (n s^2), standardize = FALSE, intercept = FALSE)` | 6e-2 off; recovered to 3e-5 only after multiplying `lambda` by `sd(y)` |
| B2 lasso: `glmnet(alpha = 1, lambda = 1 / (n b))`, KKT residual on the canonical objective | 3e-10 |
| B3 ridge: Stan L-BFGS mode against the closed form, `sigma` held at 1 | 1.5e-7 |
| B4 lasso and elastic net: Stan L-BFGS mode against `glmnet` | 0.16 and 0.15 off while reporting convergence |

B1 is why the ridge exact-limit oracle is the closed form
`(A'A + I / s^2)^{-1} A'z` rather than `glmnet`, and why the installed
`glmnet`'s conventions are pinned by a test that measures them. B4 is why
L-BFGS on the Stan program is not used as a mode oracle for a non-smooth
prior.

### `projpred_probe.R` -- the projection structure gate (probes C1-C7)

| Probe | Result |
|---|---|
| C1 `get_refmodel()` on a horseshoe `brmsfit` with `(1 \| g)` | accepted |
| C2 `varsel(forward, nterms_max = 6)`, n = 120, 9 terms | 41 s; ranking `x1, x2, (1 \| g), p1, ...`; elpd plateau at size 3, matching the planted truth |
| C3 `cv_varsel(kfold, custom grouped folds, validate_search = TRUE)` | 89 s including four brms refits |
| C4 `search_terms` forcing a protected term | works; the forced term is first in every path |
| C5 `get_refmodel()` with `gr(g, cov = K)` | refused |
| C6 `get_refmodel()` on the non-linear partition | refused |
| C7 `posterior_epred(newdata, allow_new_levels = TRUE, sample_new_levels = "gaussian")` | works |

C5 and C6 are the red rows of the structure gate. C6 together with A4 is why
the non-linear-parameter route is not shipped: exact targeting that cannot be
projected does not serve the primary decision.

### `timing_gates.R` and `timing_gates.log` -- cost at working shape (F1, F2)

F1, at barley shape (n = 48, p = 176, 48 groups, two chains of 2000,
`adapt_delta = 0.95`): 20 s wall, 8 divergences, max R-hat 1.013, minimum bulk
ESS 226, minimum tail ESS 86, top-10 recovery 0.6 at SNR 1. That is the
measurement behind the reference fit's defaults being a floor rather than a
recommendation, and behind the convergence contract rather than the settings
being what decides adequacy. F2, at deposit shape (n = 3,427, p = 758, 400
groups), was still running when the contract was written; the log carries F1
only.

### `f0b_stanvar_horseshoe.R` and its log -- the block-targeted horseshoe spike

Written for this unit rather than copied. A regularised horseshoe density
written by hand over `b[cand_idx]` in `stanvars` would give both exact block
targeting and a linear formula that projpred accepts, which is the one gap the
five shipped families leave. The script fits it and the native class-wide
horseshoe on identical simulated data (n = 100, p = 20, three true non-zero
coefficients, one random intercept over ten groups) with identical sampler
settings, and states its exit rule before the run: the route ships only if
divergences, R-hat, ESS and posterior means match the native arm within
Monte-Carlo error.

**Outcome: the exit rule is not met, and the route is not wired.** The two
arms agree on the posterior -- the largest difference between their
coefficient means is 0.0056, which is 1.45 Monte-Carlo standard errors, and
both arms rank the three true non-zero coefficients first -- but the sampler
does not agree with itself. The native arm ran 0 divergences, maximum R-hat
1.0022, minimum bulk ESS 2,790 and minimum tail ESS 3,025 in 13 s; the
hand-written arm ran 434 divergences out of 4,000 transitions, maximum R-hat
1.0122, minimum bulk ESS 665 and minimum tail ESS 138 in 16 s. That is the
funnel a centred global-local prior has and a non-centred one does not, and
it is the reason brms writes `b = zb .* sdb` rather than putting a normal
with a random scale on `b` directly. The same targeting cannot be had
non-centred without giving up the coefficient rows that were the point.

The five shipped families therefore stand: `ridge()`, `laplace()` and
`elastic_net()` target blocks exactly, and `regularized_horseshoe()` and
`r2d2()` are class-wide with `scope = "all_fixed"` as the named opt-in.
