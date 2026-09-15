# S01 -- selection validation against ASReml and glmnet

**Run** 2026-09-13. **Programme** V1 to V4 below, each item's criterion fixed
before the item ran and recorded with it. **Systems** ASReml 4.2.0.392 and
lme4 2.0.1 under R 4.5.2 with INLA 25.10.19 and glmnet 4.1.10; brms 2.23.0 via
CmdStan, projpred 2.10.0, INLA 26.8.7 and glmnet 5.0 under R 4.6.1.

## Provenance, and a correction to it

The brief named source SHA `8777be4` and tarball SHA-256 `6d23fcde...`. Both
moved during the session, and the record has to say so rather than repeat the
figure it was given.

| Artefact | State when the work began | State now | What ran on it |
|---|---|---|---|
| package HEAD | `8777be4` | `9b9fd44` (mine) on top of `3641b5b` and `3670c4c` (a concurrent process) | -- |
| R 4.6 library | build of `8777be4` | rebuilt 00:51:01 from a tree carrying `3670c4c` | V2 (first run), V3, V4 |
| `oracle_work/oraclelib` (R 4.5) | empty of 0.11.0 | installed 00:50 from the `8777be4` tarball, replaced 00:59 by `git archive 9b9fd44` | V2 (R 4.5 run) on the first; the tranche and every later R 4.5 script on the second |
| root tarball | `6d23fcde...` | `d13aca46...`, rebuilt 00:57 by the concurrent process | the R 4.5 install at 00:50 used the first |

`3670c4c` is the only code commit between `8777be4` and the build every result
here rests on. Its whole diff is in `R/fb_select.R`: `scoring_unit` becomes
`NA` rather than the string `not scored out of fold` when a route has no
out-of-fold table, the calibration approximation gains the words `units of`,
and two print methods word the absent case. No estimand, no arithmetic and no
decision changes, so no number below moves under it. That was checked rather
than assumed: **V2 was re-run on the current R 4.6 build and every one of its
33 measurements came back bit-identical** to the run on the earlier build.

**Tranche invocation**, recorded as SPEC.md's manifest requires:

```
FLEXYBAYES_ASREML_LIB=$HOME/Library/R/arm64/4.5/library \
  data-raw/asreml_oracle/run_tranche.sh <pkg>/oracle_work \
  17_v1_alpha_nested 18_v1_sim_met_candidates 19_v1_random_slope
```

source SHA `9b9fd44896d5cd6c4661c0ed357f491c40e8c50f`, clean tree, tarball
`b77bf009cfbaea0c9e08f927c712e9f439b4b780c7d57bbb322cd1721bec5ddf`, three
fixtures in `oracle_work/oracle_out/`, five seconds each.

---

## V1. The mixed-model side, against ASReml

Three designs, each fitted twice in each system -- protected-only and
protected plus candidates. Criteria are SPEC.md Amendment 5, written before
the tranche ran. lme4's REML fit is carried alongside ASReml on every
component, because two of the three ASReml formulations were written for this
suite and a disagreement they produced could be a defect in the comparator.

### V1 layer totals

| Row | Substrate | Layer 1 structural | Layer 2 (existing tolerances) | Layer 3 (Amendment 5) |
|---|---|---|---|---|
| 17 `alpha_nested` | `agridat::john.alpha`, `~ rep + rep:block` | 4 pass / 0 fail | 2 pass / 0 fail / 1 not asserted | 11 pass / 6 fail / 1 not asserted |
| 18 `sim_met_candidates` | 20 sites by 10 plots, seed 20260913 | 4 pass / 0 fail | 10 pass / 0 fail | 23 pass / 2 fail / 4 not asserted |
| 19 `random_slope` | `lme4::sleepstudy`, `(Days | Subject)` | 4 pass / 0 fail | 3 pass / 0 fail | 16 pass / 2 fail / 1 not asserted |

Every fail is a variance component and every one is characterised below. No
fixed effect and no name-matched BLUP missed.

### V1-F. Fixed effects against the REML BLUE

| Criterion | Statistic | Worst measured | Verdict |
|---|---|---|---|
| V1-F1 `abs(b_fb - b_asr) <= 2 * sd_post` | ratio of the difference to the posterior sd | 0.0816 (row 19, protected-only intercept) | PASS, 22 of 22 |
| V1-F2 relative, where `abs(z_asr) >= 2` | relative difference | 3.18e-3 (row 19, `Days`) | PASS, 11 of 11 asserted |

Four candidate coefficients on row 18 are not asserted at V1-F2 because
ASReml's own estimate is not resolved from zero, which is the gate the
criterion declares rather than a miss. On row 18 the two systems agree on
every one of the eight fixed effects to five significant figures: ASReml
`c1 = 0.55781`, flexyBayes `0.55773`; ASReml `c4 = -0.57515`, flexyBayes
`-0.57524`.

### V1-B. BLUP agreement, on levels matched by name

| Row | Group | Levels matched | Correlation (criterion 0.99) | Largest absolute difference | Verdict |
|---|---|---|---|---|---|
| 17 | `rep` | 3 of 3 | 1.000000 | 0.00214 | PASS |
| 17 | `rep:block` | 0 of 21 | -- | -- | not asserted |
| 18 | `site` | 20 of 20 | 1.000000 | 0.00718 | PASS |
| 19 | `Subject` intercept and slope | 0 of 18 | -- | -- | not asserted |

The two blocked rows are the reporting defect in section "Defects found"
below, not a numerical disagreement.

### V1-V. Variance components, posterior mean and median against REML

Every value on the variance scale. The posterior mean, median **and** mode
come from the INLA marginal transformed to the variance scale and summarised
afterwards, which is the only order that gives the moments of the variance
itself; squaring a reported standard deviation would compare a quantity
neither system reports. `cv` is ASReml's own `std.error / component`, the
identification gate Amendment 3 declared.

| Row | Component | ASReml | lme4 | cv | Identified | Posterior mean | median | mode | rel mean | rel median | rel mode |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 17 | `rep` | 0.113961 | 0.113941 | 1.118 | no | 1.93562 | 0.456283 | 0.062413 | 15.99 | 3.004 | 0.452 |
| 17 | `rep:block` | 0.061952 | 0.061944 | 0.590 | no | 0.078044 | 0.067238 | 0.047885 | 0.260 | 0.085 | 0.227 |
| 17 | `units!R` | 0.085236 | 0.085226 | 0.256 | yes | 0.091342 | 0.087904 | 0.081267 | 0.072 | 0.031 | 0.047 |
| 18 | `site` | 0.466270 | 0.466218 | 0.397 | yes | 0.556329 | 0.510870 | 0.433550 | 0.193 | 0.096 | 0.070 |
| 18 | `units!R` | 1.016091 | 1.015980 | 0.108 | yes | 1.030950 | 1.023950 | 1.010220 | 0.015 | 0.008 | 0.006 |
| 19 | `us(2)_1:1` intercept | 612.170 | 612.100 | 0.471 | yes | 720.755 | 701.522 | 655.052 | 0.177 | 0.146 | 0.070 |
| 19 | `us(2)_2:2` slope | 35.0763 | 35.0717 | 0.421 | yes | 52.0975 | 50.1869 | 46.3709 | 0.485 | 0.431 | 0.322 |
| 19 | `units!R` | 655.028 | 654.940 | 0.118 | yes | 661.134 | 657.965 | 652.353 | 0.009 | 0.004 | 0.004 |

**lme4 confirms ASReml on all eight components** to five significant figures,
including the correlated-slope formulation this suite wrote. The comparator
is sound, so every disagreement below is between REML and the posterior.

The six V1-V fails sort into three kinds.

**Kind 1, a component REML barely identifies (row 17 `rep`, `rep:block`).**
`rep` has three levels and ASReml's own `std.error` exceeds its estimate by
twelve per cent, so the identification gate Amendment 3 declared excludes it
and the Layer 2 variance-ratio check on the same row correctly recorded
`not asserted -- weakly identified`. Amendment 5's boundary test is looser
than that gate and asserted anyway. The posterior mean of a three-level
variance under a uniform prior on the standard deviation is dominated by the
prior; the median is 4.0 times the REML estimate and the mode 0.55 times it.
This is a property of the data, and the fix belongs to Amendment 5's gate
rather than to the package: **Layer 3 should carry Amendment 3's
identification gate as well as the boundary test.**

**Kind 2, the mean-against-mode estimand difference (row 18 `site`, row 17
`rep:block` and `units!R`).** On every one of these the posterior median is
inside the 15 per cent criterion and only the mean is outside it: `site`
0.193 on the mean against 0.096 on the median, `rep:block` 0.260 against
0.085. A posterior over a non-negative parameter is right-skewed and its mean
sits above its mode by an amount that grows with the skew, so a mean that
misses while the median holds is the difference Amendment 3 recorded, not a
translation defect. Amendment 5 asserted both precisely so this could be
told apart, and it was.

**Kind 3, a well-identified component where the mean, the median AND the mode
all miss (row 19 slope variance).** ASReml 35.08 and lme4 35.07 against a
posterior mean of 52.10, a median of 50.19 and a **mode of 46.37**, which is
32 per cent high. The intercept variance and the residual on the same fit
agree. This is a defect, it is finding O04's defect on a second route, and it
is described under "Defects found" below.

### V1-I. The out-of-fold increment, against ASReml refits on the same folds

`fb_select(base, ~ c1..c6, penalty = NULL, workflow = "bayes")` on design (b),
five site folds, `new_cluster`. The oracle refits both models by REML on each
training fold and scores the held-out sites at the marginal predictive
`N(x'beta_hat, sigma_e^2 + sigma_site^2)`, aggregating by the package's own
rule re-implemented rather than called. The partition is read off the
result's own fold objects, never rebuilt from the seed.

| Quantity | ASReml plug-in | flexyBayes Bayesian | Difference |
|---|---|---|---|
| protected-only score (loss) | 18.21685 | 18.22399 | 0.00713 |
| full score (loss) | 16.41723 | 16.42386 | 0.00663 |
| **increment** | **1.79963** | **1.80013** | **0.00050** |
| unit-level standard error, full arm | 0.81863 | 0.71626 | -0.10237 |
| independent units | 20 | 20 | 0 |

| Criterion | Measured | Verdict |
|---|---|---|
| V1-I1 the increments agree in sign | both positive | PASS |
| V1-I2 within one unit-level standard error | 0.00050 against 0.71626 | PASS |

The increment agrees to five parts in ten thousand of a standard error. The
package's fold loop cost 204.7 s of INLA fitting and posterior sampling; the
ten ASReml refits cost 0.2 s in total.

---

## V2. Exact limits, against a closed form, a KKT residual and glmnet

Whitened synthetic designs at three sizes, one of them wider than it is long.
Every expected value is computed in the validation script: the ridge solution
by an LU solve of the normal equations where the package uses a Cholesky, the
optimality residual by the subgradient condition written out, the unpenalised
limit by QR. The scaled max norm is `max|a - b| / max(max|b|, 1)`, because a
lasso solution has exact zeros and elementwise relative difference is
undefined on them.

| Check | Criterion | p=6, n=50 | p=30, n=100 | p=120, n=60 | Verdict |
|---|---|---|---|---|---|
| package ridge against the closed form | 1e-6 | 2.64e-16 | 3.93e-16 | 1.02e-14 | PASS |
| glmnet ridge against the closed form, at the measured convention | 1e-4 | 7.36e-11 | 5.78e-09 | 7.74e-07 | PASS |
| glmnet ridge at the naive `lambda / n` map | reported | 1.75e-2 | 1.68e-2 | 3.02e-2 | -- |
| lasso, KKT residual of glmnet's point | 1e-6 | 6.38e-09 | 5.32e-07 | 9.60e-07 | PASS |
| elastic net at alpha 0.25 / 0.5 / 0.75, KKT residual | 1e-4 | 5.66 / 6.87 / 5.92 e-09 | 6.29 / 4.22 / 2.42 e-07 | 1.92 / 1.24 / 1.34 e-06 | PASS |
| proximal mode against glmnet | 1e-6 | 6.36e-11 | 8.64e-09 | **3.40e-06** | **FAIL at p > n** |
| unpenalised limit against QR least squares | 1e-6 | 5.15e-16 | 1.72e-15 | not applicable | PASS |
| package KKT function against the one written here | 1e-12 | 0 | 0 | 0 | PASS |

The naive-map row is the hypothesis the convention measurement rejects, kept
so the size of the correction is visible: without the response-scale term the
ridge answer is wrong in the second decimal place.

**The one fail is the comparator, and it can be shown to be.** At `p = 120`,
`n = 60` the two points differ by 3.40e-6 while the package's point carries a
KKT residual of 2.6e-14 and glmnet's carries 9.6e-7 -- forty-five thousand
times larger. The canonical objective at the package's point is lower by
9.6e-12, which is the numerical floor of the objective itself, so the two
points sit in a flat direction of a near-degenerate wide problem and the
package's is the one nearer the optimum. Recorded as a comparator-precision
result, not a package defect. The criterion was declared at 1e-6 before the
run and is left where it is; what changed is the reading of it.

At `p = 120`, `n = 60` the unpenalised limit is not applicable: with more
columns than rows the minimiser is a set, not a point, so there is no referent
to compare against. Recorded rather than skipped.

### The convention, at two independent glmnet versions

This is the check the package's own test cannot make, because it measures
whichever glmnet is installed.

| Rule | glmnet 5.0, L2 max abs | glmnet 4.1.10, L2 max abs | glmnet 5.0, L1 KKT | glmnet 4.1.10, L1 KKT |
|---|---|---|---|---|
| `none` | 1.888850e-02 | 1.888850e-02 | 9.758056e-09 | 9.758058e-09 |
| `sd_population` | 7.363719e-03 | 7.363719e-03 | 2.949693 | 2.949693 |
| `sd_sample` | 7.185723e-03 | 7.185723e-03 | 2.995683 | 2.995683 |
| **`rms`** | **1.284008e-10** | **1.284006e-10** | 4.863440 | 4.863440 |

Both versions select `response_scale_rule = rms` and `l1_scale_rule = none`,
and the winning residual is seven orders of magnitude below the runner-up in
both. The two versions agree to 1e-14 absolute on every rule's residual, and
every one of the 33 checks returns the same verdict under both. The measured
L2 convention is therefore a property of glmnet across a major version
boundary rather than of the installed copy.

---

## V3. Real data, side by side (Tier 3, reported), and the one gated identity

`agridat::john.alpha` with eight simulated candidates, four of them in a
compound-symmetric block at correlation 0.7 and two of those four carrying
signal (`z1 = 0.35`, `z3 = -0.30`). Base `yield ~ 1 + (1 | rep) +
(1 | plot_block)` on brms, folds on `plot_block` (18 units, 6 folds), four
chains, seed 20260913.

### Certificate stanzas

Both arms: `mode_oracle$available = FALSE`, with the reason that the identity
is defined only with the covariance plugged in and the residual scale held
fixed, and the Bayesian route does neither. That is the gated criterion V3-2
and both arms **PASS**.

Both arms: `predictive_side_by_side$available = TRUE`, metric `nlpd`, at the
cluster unit, on 18 units. Equivalence levels
`gls = no_identity_expected [reml_plugin_covariance]` and
`glmnet = posterior_predictive_comparable [reml_plugin_covariance]`. The
audit's finding 3 asked for exactly this split and for
`posterior_predictive_comparable` not to be printed until both arms emit a
shared unit-level table; under `benchmark = TRUE` they do, and the stanza
carries the table.

### Shared unit-level nlpd, 18 clusters

| Arm | laplace run | elastic-net run |
|---|---|---|
| bayes (the primary) | 0.9368 (se 0.0737) | 0.9315 (se 0.0712) |
| gls, REML plug-in covariance | 0.8809 (se 0.0497) | 0.8809 (se 0.0497) |
| glmnet | 0.9143 (se 0.0573) | 0.9195 (se 0.0570) |

The three losses sit inside one standard error of each other. The gls arm is
identical across the two runs, as it must be: it carries no penalty.

### Three supports, three procedures

| Procedure | Support over the eight candidates | True kept | False kept |
|---|---|---|---|
| projection, laplace arm | `z1 + z3` | 2 | 0 |
| projection, elastic-net arm | `z1 + z3` | 2 | 0 |
| posterior threshold, effect size 0.15, probability 0.5 | `z1 + z2 + z3 + z4` | 2 | 2 |
| glmnet cross-validated, `lambda.1se` = 0.16486 | empty | 0 | 0 |
| glmnet cross-validated, `lambda.min` = 0.028147 | `z1 z2 z3 z4 z6 z7 z8` | 2 | 5 |

The two false positives the threshold rule keeps are `z2` and `z4`, both
inside the correlated block, which is where a coordinate-wise rule is least
able to separate a signal from its neighbour. glmnet's own one-standard-error
rule empties the model on this trial. Reported, not gated: the four
procedures optimise different quantities.

Wall time: base fit 8.1 s; laplace arm with benchmark 318.5 s; elastic-net arm
with benchmark 276.1 s; posterior-threshold decision 9.1 s.

### The MAP identity, and the gap the audit asked for

Substrate: the whitened, protected-residualised design, 72 rows by 8 columns.
The plug-in covariance is built from an lme4 REML fit (`rep` 0.1132,
`plot_block` 0.04183, residual 0.3711), so the substrate is not one the
package built for itself.

| Quantity | Value | Criterion | Verdict |
|---|---|---|---|
| KKT residual, package proximal MAP | 3.997e-14 | -- | -- |
| KKT residual, glmnet at the measured convention | 7.577e-07 | -- | -- |
| **MAP against glmnet, scaled max norm** | **8.909e-09** | 1e-4 | **PASS** |
| sampled posterior mean against the MAP, same design, same penalty, `sigma` held at 1 | 1.073e-2 | reported | -- |
| KKT residual of the sampled posterior mean | 1.895 | reported | -- |
| exact zeros in the MAP / in the posterior mean | 0 / 0 | -- | -- |

The identity holds on the substrate where it is defined, to nine decimal
places, and the sampled posterior mean on that same substrate sits a million
times further from the mode than glmnet does -- its KKT residual is 1.895
against 4e-14. The package claims the identity nowhere: `mode_oracle` is
empty on both sampled arms, and the glmnet arm is labelled
`posterior_predictive_comparable`, never `map_exact_fixed_covariance`.

One caveat on this run: at a Laplace scale of 0.5 the canonical lasso solution
on this design keeps all eight columns, so the contrast between a sparse mode
and a dense mean is not exercised here. V2 exercises it at three sizes, where
the lasso leaves 5 of 6, 28 of 30 and 56 of 120 coefficients non-zero and the
posterior mean would leave none.

**A measurement error of this suite's own, recorded because it changes how the
first attempt reads.** The first run of this comparison called `glmnet()` at
its default `thresh = 1e-7`, and the returned point sat 0.0202 from the
canonical optimum in the KKT norm; the identity then measured 1.96e-4 and
failed its 1e-4 criterion. Tightening the solver to the settings V2 uses
(`thresh = 1e-16`, `maxit = 1e6`) moved the measurement to 8.9e-9. The failing
number measured where glmnet stopped.

---

## V4. Regularised horseshoe

Design (b), 20 sites by 10 plots, six candidates of which `c1` and `c4` carry
signal, seed 20260913.

### V4(i) intercept-only base, so the candidate block is the whole class

| Quantity | Value | Criterion | Verdict |
|---|---|---|---|
| `calibration_n` | 20 | equals the 20 sites, not the 200 rows | **PASS** |
| `calibration_basis` / `calibration_unit` | `resampling_units` / `site` | -- | -- |
| decision resolved | `projection` from `auto`, structure gate green on all three features | -- | -- |
| chosen terms | `(1 | site) + c4 + c1` | -- | -- |
| chosen candidates against the planted truth | `c4 + c1`; 2 of 2 true, 0 false | -- | -- |
| convergence verdict | pass | -- | -- |
| divergent transitions | 0 | -- | -- |
| max Rhat / min bulk ESS / min tail ESS | 1.0033 / 1125 / 1066 | -- | -- |
| min EBFMI / transitions at max treedepth | 0.626 / 0 | -- | -- |
| posterior predictive check, p(mean) and p(sd) | 0.460, 0.502 | -- | -- |
| shrinkage factors, labelled approximate | c1 0.141, c2 0.440, c3 0.526, c4 0.034, c5 0.570, c6 0.564 | -- | -- |
| effective non-zero at kappa 0.5 | 3 against `p0 = 2` | -- | -- |
| wall time | base fit 7.0 s, `fb_select` 136.6 s | -- | -- |

`(1 | site)` appears in the chosen list because the search forces the base's
random terms into every submodel, as the documentation says. The candidate
half of the support is exactly the planted truth.

**The first attempt refused, and the refusal was right.** At the module's
default `adapt_delta = 0.95` the reference fit ended two transitions of six
thousand in a divergence and the run refused with
`penalty_sampler_diagnostics_failed`, saying that a selection read off that
posterior would be a statement about the sampler. The numbers above are from
a re-run at `adapt_delta = 0.999`, passed through
`fb_engine("brms", adapt_delta = 0.999)`. Two observations follow, both in
"Smaller findings" below.

### V4(ii) the same request with a protected covariate present

| Criterion | Measured | Verdict |
|---|---|---|
| refuses with `penalty_block_targeting_unsupported` | that code, from a `flexybayes_refusal` condition | **PASS** |
| the refusal names a block-targeting family | names `ridge()`, `laplace()` and `elastic_net()` | **PASS** |
| the refusal names `scope = "all_fixed"` | yes, with what it costs -- protected coefficients under the global-local prior, forced into every submodel, and an `approximations` entry | **PASS** |

### V4(iii) the alternative the refusal names, run

`laplace(scale = 0.5)` on the same base with the protected covariate present
ran to completion in 109.9 s and chose `x_prot + (1 | site) + c4 + c1`:
the protected covariate retained, both true candidates recovered, no false
candidate. `calibration_n` 20.

---

## Defects found

### D1. The correlated-slope route on INLA is not scale invariant, and its slope variance is high by a third

**Locus.** The `iid2d` field the `(x | g)` route emits, reached through
`.fb_inla_cor_slope_index()` / `.fb_inla_cor_slope_slope_index()` in
`R/emit_inla.R` (around lines 3218 to 3223) and its hyperparameter block. This
is finding **O04's defect on a second route**: O04 measured it on `us()`
through `iidKd` and wrote that whether the joint Wishart prior is the cause
was not established. It now has two loci, so a fix has to cover both.

**Evidence 1, against REML.** Row 19, `lme4::sleepstudy`, 18 subjects. ASReml
35.076 and lme4 35.072 for the slope variance against a posterior mean of
52.10, median 50.19 and mode 46.37. The mode is the summary a REML estimate is
comparable to and it is 32 per cent high. The intercept variance on the same
fit is 7 per cent high at the mode and the residual 0.4 per cent, so this is
not a whole-block scaling error.

**Evidence 2, the metamorphic property, which needs no oracle.** Multiplying
the response by a constant is exact, and a variance ratio and a correlation
are invariant to it. REML returns the same numbers at every scale.

| Response scaled by | `sd_intercept / sigma` | `sd_slope / sigma` | intercept-slope correlation |
|---|---|---|---|
| 0.1 | 1.0363 | 0.2827 | -0.0215 |
| 1 | 0.9395 | 0.2807 | 0.0233 |
| 10 | 1.9301 | 0.3339 | **0.9689** |
| ASReml, any scale | 0.9667 | 0.2314 | 0.0656 |

At ten times the response the reported intercept-slope correlation moves from
0.02 to 0.969 -- across essentially the whole of its range -- on a
transformation that changes no model. The ratio of the intercept standard
deviation to the residual doubles. The tenfold fit also raises the package's
own numerical-confirm gate (`INLA mode.status = 8`), which is the one
protection a user has here and it fires only at the extreme.

**Why the row-19 comparison passed its correlation check anyway.** The true
correlation on `sleepstudy` is near zero, so a collapse toward zero lands
inside the 0.10 tolerance. A substrate with a real correlation would not.

### D2. `ranef()` on an INLA fit loses the level labels of a structured group, and returns two list elements under one name

**Locus.** `.fb_random_inla()` in `R/summary_asreml.R`.

- Line 1002, `level_vals <- tab$ID %||% seq_len(nrow(tab))`, followed by an
  `is.character(level_vals)` guard on the relabelling. For an interacted or a
  correlated-slope group INLA's `ID` column is **integer**, so the guard never
  fires and the levels come back as `1 ... n`. Measured: `rep:block` on
  `john.alpha` returns levels `"1" ... "18"` while ASReml returns
  `rep_R1:block_B1 ...`, and `(Days | Subject)` returns `"1" ... "36"` while
  the subject identifiers sit unused in `fit$level_labels$Subject`.
- Line 1019, `names(out) <- vapply(names(sr), relabel, ...)`. For a
  correlated slope, `R/canonical_names.R` lines 302 to 304 deliberately map
  both INLA index names to one label, so the returned list has **two elements
  named `"Days | Subject"`** and `ranef(fit)[["Days | Subject"]]` reaches only
  the first. The two blocks are the 18 intercept deviations and the 18 slope
  deviations; a caller cannot tell them apart by name.

**Consequence measured here.** Two of the four V1-B1 BLUP comparisons could
not be made at all. Aligning by row order was refused: `interaction()` orders
`rep` fastest and ASReml orders it slowest, so a positional match would have
been silently wrong on exactly the row where it looked easiest.

### D3. `?fb_engine` does not document the option that unblocks a refused selection

**Locus.** `R/fb_engine.R` line 52 (and the generated `man/fb_engine.Rd`
line 12) say the recognised option names are `n_samples`, `warmup` and
`chains`. Line 36 of the same file has
`.fb_engine_opt_names <- c("n_samples", "warmup", "chains", "adapt_delta")`,
and its own comment says `adapt_delta` exists because the convergence contract
refuses on a single divergent transition and "without a lever a caller who
hits one has nothing to reach for but a different seed".

Both selection runs in this programme hit that refusal -- V4(i) on two
divergences of six thousand, V3's laplace arm on **one** -- and the lever that
clears it is the one option the help page omits. A reader following the
documentation would conclude the refusal cannot be cleared.

---

## Smaller findings, reported rather than filed as defects

1. **The convergence contract refuses on one divergent transition in six
   thousand.** V3's laplace arm refused at 0.017 per cent divergences. The
   contract is strict by design and the refusal message is exact about why;
   the observation is that a heavy-tailed prior on a 72-row trial reaches it
   routinely, and that the remedy is D3's undocumented option.
2. **Amendment 5's Layer 3 gate should carry Amendment 3's identification
   gate.** Three of the six V1-V fails are on components ASReml itself reports
   with `std.error / component` above 0.5. Layer 2 declines to assert on them
   and Layer 3 asserted. This is a defect in the criterion written for this
   suite, not in the package, and it is recorded here rather than fixed
   quietly by moving a tolerance.
3. **The elastic-net arm's `gls` benchmark is byte-identical across the two
   penalties**, as it must be, which is a small confirmation that the
   benchmark arm is not reading the primary arm's penalty.

---

## What could not be run, and why

| Item | Reason |
|---|---|
| BLUP agreement on `rep:block` (row 17) and on `Subject` (row 19) | D2: one system returns integer level indices, so no level matches by name, and positional matching is a guess this suite refuses to make |
| The unpenalised limit at `p = 120, n = 60` (V2) | with more columns than rows the unpenalised minimiser is a set, not a point |
| Stan's optimiser as a MAP oracle for the Laplace family (V3) | SPEC section 4.3 disqualifies L-BFGS for non-smooth priors on probe B4's evidence; the proximal solver with its KKT certificate was used instead, and glmnet is the independent check on it |
| A sparse-mode against dense-mean contrast on the V3 substrate | at Laplace scale 0.5 the canonical solution on that 72 by 8 design keeps all eight columns; the contrast is exercised in V2 instead |
| `r2d2`, `scope = "all_fixed"`, `elastic_net` under a class-wide block, and the `stanvar_density` horseshoe spike | outside the four V-items as briefed |
| A repeatability screen on rows 17 to 19 | Amendment 2 requires one before a comparison, and the tranche ran the comparison without it; the run-to-run spread is visible instead in the two independent fits of row 19 taken here, whose slope-variance posterior means were 49.24 and 52.10 |
| `R CMD check --as-cran`, the full test suite, and any vignette rebuild | outside this programme, and a concurrent process held the 4.6 library |

## Files

- Scenarios, comparators and the runner: `data-raw/asreml_oracle/scenarios.R`,
  `compare.R`, `run_one.R`; criteria in `SPEC.md` Amendment 5.
- Fixtures: `oracle_work/oracle_out/17_v1_alpha_nested.rds`,
  `18_v1_sim_met_candidates.rds`, `19_v1_random_slope.rds`.
- Validation scripts, logs and one csv per item:
  `oracle_work/selection_validation/`.
