# O04: `us()` on INLA disagrees with ASReml, and is not scale-invariant

**Found** 2026-09-11, oracle scenarios 06 and 07, on `agridat::besag.met`
restricted to two and three counties. ASReml 4.2.0.392, INLA 25.10.19.

## The disagreement

Structural layer exact on both scenarios: same parameter count, same
estimable fixed effects, same observations. The estimates are another matter.

**Two environments.**

| Quantity | ASReml | flexyBayes |
|---|---|---|
| genotype variance, C1 | 23.39 | 124.56 |
| genotype variance, C2 | 33.43 | 133.06 |
| genetic correlation | **0.930** | **0.147** |
| residual variance | 296.99 | 281.46 |

**Three environments**, correlations only:

| Pair | ASReml | flexyBayes |
|---|---|---|
| C1, C2 | 0.898 | 0.171 |
| C1, C3 | 0.915 | 0.220 |
| C2, C3 | 0.955 | 0.309 |

The residual agrees closely in both. The genotype variances are four to five
times too large and the correlations are collapsed toward zero. For a
multi-environment trial the genetic correlation between environments is not a
nuisance parameter, it is the question the model is asked, and it is the
number that decides whether a breeding programme treats two environments as
one target or two.

This is not instability. The repeatability screen rated both rows STABLE, at
ratios of 1.01 and 1.04 across three identical fits. The answer is
reproducible and different.

## The harder finding: the fit is not scale-invariant

Multiplying the response by a constant is an exact transformation. Variance
ratios and correlations are invariant to it, and ASReml returns 0.930
whatever the scale.

| Response scaled by | correlation | sd_gen / sigma |
|---|---|---|
| 1 | 0.150 | 0.669, 0.688 |
| 0.1 | 0.144 | 0.664, 0.684 |
| 10 | **-0.622** | **3.465**, 0.899 |

At ten times the scale the correlation changes sign and a variance ratio
moves by a factor of five.

This needs no oracle. It is a metamorphic property the package could have
tested against itself at any point, and it is on the audit's own work-package
4 list ("response and covariate unit rescaling") which has not been built.
An independent oracle found it first only because the oracle was built first.

## What is not yet established

Whether the cause is the joint Wishart prior the `iidKd` field carries. The
summary reports `engine default` on the `us()` components while `sigma`
reports the package's own prior, and the release notes record the Wishart
scale being taken from the response variance for `iid2d` and `iidKd`. On this
data the response variance is about 1,483 against genotype variances of 23 to
70, so a scale set from the response would sit far above the components it
governs. That is a hypothesis with a plausible mechanism, and it is written
here as a hypothesis: the emitted `iidKd` hyper block was not read, because
this route does not populate the formula slot the other routes do.

The first version of finding O01 was exactly this shape, a plausible
mechanism, written up before the emitted model was read, and wrong. The
mechanism here is not to be reported as the cause until the hyper block is
read directly.

## Bearing on the registry

The rows *Unstructured covariance over a factor* on INLA and on brms are
classified `covariance-family`, meaning the fitted covariance set is the same
after a documented parameter map. On INLA that claim is not supported by this
measurement and the row must not keep it unqualified. The brms route was not
exercised here and is untouched by this finding.

## Disposition, 2026-09-11

Committed, and the claim surfaces brought into line with the measurement.

- The registry row *Unstructured covariance over a factor, on INLA* is
  reclassified from `covariance-family` to **`analogue`**, with the
  measurement in the row and not in a note beneath it. The change
  propagates to `inst/TRANSLATION.md`, the shipped CSV, and (through the
  capability matrix) to `README.md` and `inst/KNOWN_ISSUES.md`, where it is
  open issue 1. The capability-matrix note previously said the Wishart
  "shrinks the correlations toward each other", which is the wrong direction
  as well as the wrong claim.
- Scenarios `06_us_dim2` and `07_us_dim3` are committed, and the tranche has
  been re-run from a committed SHA so the banked evidence names the code it
  came from.
- The scale-invariance check is now a test and not an exercise, in
  `tests/testthat/test-metamorphic-scale-invariance.R`. It is two-sided: an
  invariance **gate** on the validated structures, which must pass, and a
  **tripwire** on this defect, which asserts the gap and names what to close
  when the gap goes away. A defect left as a silent skip becomes an unknown
  defect.
- Decision recorded 2026-09-11: `us()` on INLA demoted to an analogue.

## One thing measured while writing that test

The rescaled fit does not converge, and says so. At a tenfold response scale
the fit returns `Mode status: 6` and `Numerical confirm: FAIL`, and the
package warns. **The unscaled fit reports `Mode status: 0` and passes its
numerical confirm**, and it is the one a user actually runs, and it is wrong.

That is the same shape as O01b, one level along. The numerical gate is not
blind here, it fires on the harder problem and stays quiet on the easier one
that happens to be misspecified. A gate that passes the fit you will use and
fails the fit you will not is worse than no gate for this purpose, because it
reads as assurance. Whatever the eventual cause turns out to be, a repair has
to be judged on the unscaled fit, not on whether the warning appears.

## The banked run, and what three refits actually show

The numbers at the top of this finding come from the first run, which was made
from an uncommitted tree. The evidence that now ships was produced at commit
`8638306` on a clean tree, from `git archive`, with the scenarios asserting
the correlation and not only the fixed effects.

| Quantity | ASReml | flexyBayes (banked) | Three refits |
|---|---|---|---|
| correlation, 2 environments | 0.930 | **0.142** | 0.141 to 0.147 |
| correlation, C1 vs C2 (3 env) | 0.898 | **0.198** | 0.171 to 0.198 |
| correlation, C1 vs C3 | 0.915 | **0.231** | 0.220 to 0.238 |
| correlation, C2 vs C3 | 0.955 | **0.320** | 0.309 to 0.320 |
| genotype variances, 2 env | 23.4, 33.4 | 124.3, 132.6 | |
| residual, 2 env | 297.0 | 281.6 | |

**The gap is the finding; the third digit is one run's.** ASReml returns the
same correlation on every refit. Ours moves in the second decimal (0.171 to
0.198 on one pair is a 16 per cent spread) so the earlier sentence in this
finding, that three identical fits agree "to 1 per cent", was true of the
*variance components* the repeatability screen measures and is not true of the
correlations derived from them. The claim surfaces now quote the banked value
with the observed range beside it instead of a single figure that the next
refit would contradict.

## The instrument was also at fault, and that is the third time

Scenarios 06 and 07 were added to test `us()`. On their first two runs they
banked `numeric 1 pass / 0 fail` (the intercept) because the scenario
schema had slots for variance ratios and BLUPs and none for a correlation, so
`cmp_correlation()` and its declared 0.10 tolerance, both present in
`compare.R` since SPEC.md was written, were never called. An oracle reporting
a pass on the row its own finding says is wrong is worse than no oracle for
that row.

With the correlation asserted, scenario 06 reports `1 pass / 1 fail` and
scenario 07 reports `0 pass / 4 fail`, which is what the measurement always
said. This joins O03, where the comparator read sectioned residuals on the
wrong scale, as the second defect found in the instrument and not in the package, and the third instance overall of the rule that **the instrument is
an artefact too**.
