# O04 -- `us()` on INLA disagrees with ASReml, and is not scale-invariant

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
nuisance parameter -- it is the question the model is asked, and it is the
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
4 list -- "response and covariate unit rescaling" -- which has not been built.
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

The first version of finding O01 was exactly this shape -- a plausible
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
  measurement in the row rather than in a note beneath it. The change
  propagates to `inst/TRANSLATION.md`, the shipped CSV, and -- through the
  capability matrix -- to `README.md` and `inst/KNOWN_ISSUES.md`, where it is
  open issue 1. The capability-matrix note previously said the Wishart
  "shrinks the correlations toward each other", which is the wrong direction
  as well as the wrong claim.
- Scenarios `06_us_dim2` and `07_us_dim3` are committed, and the tranche has
  been re-run from a committed SHA so the banked evidence names the code it
  came from.
- The scale-invariance check is now a test rather than an exercise, in
  `tests/testthat/test-metamorphic-scale-invariance.R`. It is two-sided: an
  invariance **gate** on the validated structures, which must pass, and a
  **tripwire** on this defect, which asserts the gap and names what to close
  when the gap goes away. A defect left as a silent skip becomes an unknown
  defect.
- Decision recorded at
  `ledger/r-package/2026-09-11-us-inla-demoted-to-analogue.cairn.md`.

## One thing measured while writing that test

The rescaled fit does not converge, and says so. At a tenfold response scale
the fit returns `Mode status: 6` and `Numerical confirm: FAIL`, and the
package warns. **The unscaled fit reports `Mode status: 0` and passes its
numerical confirm**, and it is the one a user actually runs, and it is wrong.

That is the same shape as O01b, one level along. The numerical gate is not
blind here -- it fires on the harder problem and stays quiet on the easier one
that happens to be misspecified. A gate that passes the fit you will use and
fails the fit you will not is worse than no gate for this purpose, because it
reads as assurance. Whatever the eventual cause turns out to be, a repair has
to be judged on the unscaled fit, not on whether the warning appears.
