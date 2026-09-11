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
