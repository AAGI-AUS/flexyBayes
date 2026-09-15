# O01 -- `diag(f):g` is bistable, and its summary misreports the prior

**Found** 2026-09-11 by the ASReml oracle suite, scenario `03_diag_gxe`, at
source SHA `5e4e0ef`, tarball SHA-256 `6b8ae8aa…`. Substrate
`agridat::besag.met` (6 counties, 64 genotypes, 3 replicates, 1,152 rows).
ASReml 4.2.0.392, INLA 25.10.19, R 4.5.2.

## Correction to the first version of this finding

The first draft of this document claimed that a valid `sd(gen)` prior on a
`diag(f):g` term was accepted and silently discarded, and that the route had
no prior branch at all. **Both claims were false**, and they were reached the
wrong way: from the `prior` column of `summary(fit)$varcomp`, which is a
report about the fit rather than the fit itself.

Reading the emitted INLA formula settles it. Without a user prior the six
`f()` terms carry the package's own uniform-SD default:

~~~
hyper = list(prec = list(prior = "expression: U=192.5870372718752; ..."))
~~~

With `fb_prior(sd(gen) ~ uniform(lower = 0, upper = 40))` they carry:

~~~
hyper = list(prec = list(prior = "expression: U=40; ..."))
~~~

The prior layer reaches this term, the user's prior is applied, and the
default is applied when none is given. The lesson is the project's own:
a display column is a proxy, and a proxy is not the fact.

## What is actually defective

**O01a -- the summary misreports the prior on `at_simple` rows.**
`summary(fit)$varcomp$prior` reads `engine default` for every
`sd_gen__env<level>` row while the emitted formula carries either the package
default or the user's prior. The `sigma` row reports its prior correctly. This
is a reporting defect on a claim surface, and it is how the false version of
this finding was reached, so it has cost something already.

**O01b -- the fit is bistable, and an informative prior does not settle it.**
Eight identical calls under the default prior:

| Run | genotype SD range across six environments | collapsed to ~0 |
|---:|---|---:|
| 1 | 2.76 to 7.30 | 0 of 6 |
| 2 | 9.8e-10 to 0.056 | 2 of 6 |
| 3 | 5.2e-15 to 2.7e-08 | 6 of 6 |
| 4 | 6.0e-14 to 1.5e-10 | 6 of 6 |
| 5 | 1.2e-16 to 3.0e-10 | 6 of 6 |
| 6 | 2.8e-10 to 8.6e-07 | 6 of 6 |
| 7 | 6.4e-13 to 8.7e-04 | 5 of 6 |
| 8 | 5.1e-07 to 2.27 | 1 of 6 |

Every run reported `Mode status: 0 (converged)` and `Numerical confirm: PASS`.
Tightening the prior to `uniform(0, 40)` did not prevent collapse: the six
upper credible bounds came back between 2.9e-10 and 3.9e-06.

ASReml on the same rows returns genotype variances of 0.00, 3.28, 34.11, 0.00,
3.77, 29.69 -- two at the boundary and four substantial, and stable across
refits. The disagreement with the incumbent is therefore real and is not
explained by the prior.

## The prior is not the lever -- measured, not assumed

Twenty-four fits, four prior settings, six identical runs each, same data.

| Prior on `sd(gen)` | runs collapsing all six | max SD on a run that did not |
|---|---:|---|
| default `uniform(0, 192.6)` | 5 of 6 | 6.42 |
| `uniform(0, 40)` | 5 of 6 | 7.71 |
| `uniform(0, 10)` | 5 of 6 | 1.28 |
| `pc(upper = 10, prob = 0.1)` | 3 of 6 | 3.44, 3.13 |

ASReml, stable across refits: 0.00, 1.81, 5.84, 0.00, 1.94, 5.45 on the SD
scale.

Tightening the uniform bound by a factor of twenty changes nothing. The
penalised-complexity prior roughly halves the collapse rate and does not
remove it. Nineteen of twenty-four fits report no genetic variance in any
environment. The five that do not land near ASReml's largest components,
which says the non-collapsed mode is the right one and the optimiser reaches
it about one time in five.

This separates the finding from its two predecessors. The factor-analytic
specific variance at 0.007 and the `iid2d` Wishart at unit scale were both
our own defaults, and setting the prior correctly fixed both. Here the prior
was verified to reach the model, varied across two orders of magnitude and
two families, and the behaviour is unchanged. On the evidence this is mode
selection in the engine's hyperparameter optimisation, not a prior this
package chose.

Two qualifications, both material. The structure is not universally unstable:
a small balanced simulation (4 environments, 12 genotypes, 2 replicates)
returns well-identified components on every run. And `besag.met` is a
low-signal case -- residual SD near 19.7 against genotype SDs near 5.8 -- so
the components are weakly identified even where they are real. What is not
acceptable is that the fit reports `Mode status: 0` and
`Numerical confirm: PASS` on both modes, so nothing in the returned object
distinguishes the run that found the variance from the four that did not.

## What the package already does well here

The boundary-collapse warning fires on every collapsed run, names each
component with its upper credible bound against the residual scale, states
that the optimiser converged to a solution with the term at zero, says in
terms that **this is not always a property of the data**, and offers three
routes including an informative prior. A user is neither left in silence nor
told a falsehood about the cause. The earlier draft of this document said the
warning misattributed the cause; that was unfair and is withdrawn.

## Resolution -- the better mode is adopted, not merely flagged

Ratified by Max 2026-09-11. `R/inla_restart_probe.R` detects components that
came back pinned at zero variance, restarts the fit from a mode where they
carry real variance, and compares marginal likelihoods. Where the restart is
better by more than two log units, **the restarted fit is the one returned**,
the numerical gate is re-run against it, and the adoption is stated in the
convergence block and as a message.

Flagging alone was considered and rejected. The measured gap on `besag.met`
ran from 65 to 154 marginal-log-likelihood units, and the inferior mode
reported no genetic variance in any of six environments where ASReml found
substantial variance in two. Returning that fit with a caveat attached leaves
the user holding the wrong answer.

Measured after adoption, five identical calls:

| Environment | flexyBayes after adoption (SD) | ASReml (SD) |
|---|---|---|
| C1 | 4.07 -- 5.83 | 0.00 |
| C2 | 4.55 -- 5.29 | 1.81 |
| C3 | 5.76 -- 7.77 | 5.84 |
| C4 | 3.42 -- 4.82 | 0.00 |
| C5 | 3.78 -- 5.57 | 1.94 |
| C6 | 6.03 -- 7.96 | 5.45 |

Every run now returns substantial variance in every environment, and the
run-to-run spread falls inside the screen's factor-of-two stability
criterion. The two systems agree closely where the component is
well-identified (C3, C6) and differ where ASReml has pinned the component at
the boundary (C1, C4) or near it (C2, C5).

**That remaining difference is expected and is not a translation defect.**
REML can and does place a variance estimate exactly at zero; the mean of a
posterior over a non-negative parameter cannot be zero, and a uniform prior on
the standard deviation puts mass away from it. This is the documented
estimation difference between the two systems, which is why the suite's
identification gate declines to assert on a component whose ASReml standard
error is not small relative to its estimate -- and C1 and C4, estimated at
exactly zero, are precisely such components.

## The cost

One additional fit, and only when a component comes back pinned at zero. A
fit with nothing at the boundary never enters the probe and pays nothing.

## Bearing on the translation registry

The row *Diagonal heterogeneous variance over a factor* is classified
`covariance-family`. Its structural layer passes exactly against ASReml here:
7 variance parameters both sides, 6 estimable fixed effects, 1,152
observations. The numerical layer cannot be asserted while repeated identical
fits do not agree with each other, which is a stronger statement than
disagreeing with ASReml.
