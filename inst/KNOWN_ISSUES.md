# flexyBayes -- known issues and open computational problems

This file ships with the package (`system.file("KNOWN_ISSUES.md", package =
"flexyBayes")`). It is an honest, current statement of what flexyBayes can and
cannot do, and an open invitation to the team to help close the gaps. It is
written for the internal release: the package validates cleanly on the model
classes it supports, and it is candid about the ones it does not yet support.

The single most important thing to know: **the full multi-environment-trial
(MET) model -- genotype-by-environment interaction random effects plus a
heteroscedastic per-environment residual -- is not fittable by any trustworthy
flexyBayes backend in this release.** ASReml and `lme4` fit it; flexyBayes does
not yet. The reasons are specific and, in most cases, addressable. They are
below, with the intended route to a fix, so that they can be solved
collectively rather than rediscovered.

## Backend support by model class

Status reflects the 2026-06 empirical validation (`barrero.maize`,
`crowder.seeds`, `yates.oats`, and a simulation with known truth).

| Model class | greta | INLA | brms (Stan) |
|---|---|---|---|
| Gaussian LMM, simple random intercepts | validated | validated | validated |
| GLMM (binomial / Poisson / NB), simple RE | validated | validated | validated |
| Random slopes, structured covariance (`us`/`fa`/`ar1`), GBLUP, pedigree, separable spatial, univariate splines | supported (syntax emits; convergence is model-specific) | supported | supported |
| Streaming exact aggregation (factor fixed + simple RE) | supported | supported | n/a |
| **Interaction / nested random effects (GxE)** -- `gen:loc`, `gen:loc:yearf` | expressible, does **not** converge | **refused** | **refused** |
| **Heteroscedastic / per-stratum residual** -- `dsum(~ units \| env)` | expressible, does **not** converge | **refused** | **refused** |
| Interaction *fixed* effects, binomial path -- `y ~ a*b` | supported | **known bug** | supported |

"Supported" means the syntax emits and the model is tested; it does not promise
that every fit converges at vignette-scale budgets. Convergence is always
model- and backend-specific, and the package reports it -- treat a printed fit
with a high R-hat badge as a diagnostic, not a result.

## The full-MET boundary, stated plainly

The realistic Barrero Model 1 needs two things this release lacks on the
trustworthy backends:

1. **Interaction / nested random effects** (`gen:loc`, `gen:loc:yearf`, ...).
2. **A heteroscedastic residual** (one error variance per environment).

What each backend does today:

- **INLA** refuses both. Interaction random effects are addressable (an `iid`
  effect over the combined factor). The 107-stratum residual is **not** cleanly
  addressable: INLA integrates over hyperparameters numerically, which is
  tractable only for ~15-20 hyperparameters; 5 + 107 = 112 is far beyond that.
  This is a property of INLA, not of our code.
- **brms** refuses both, but Stan supports both natively (`(1 | a:b)` and a
  distributional `sigma ~ env`). This is the **most promising path**: an emit
  extension, not new methodology, and Stan's NUTS converges where greta does
  not. The effective random-effect dimension on observed cells is a few
  thousand (not the 170k of the full factorial), well within Stan.
- **greta** expresses both (the validation produced per-environment residual
  terms) but does **not converge** (R-hat was infinite at a 300-iteration
  budget and 302 at 2000). Fixing this is open sampler research.

The recommended direction (see the roadmap in the project's release docs):
**make brms the full-MET workhorse**, give **INLA** the interaction random
effects plus a *hierarchically-shrunk* residual (a prior
`log sigma^2_env ~ Normal(mu, tau)` that replaces 107 free hyperparameters with
two -- a model ASReml cannot fit and INLA can), and **scope greta out** of this
class rather than letting it be a silent fallback.

## Open issues (contributions welcome)

Each is reproducible on this release. Priority is for the MET use case.

1. **INLA: interaction / nested random effects.** Map `a:b` to
   `f(interaction(a, b), model = "iid")`. Emit-layer change; INLA represents
   this natively. *Unblocks GxE on the fast engine.*
2. **brms: interaction / nested random effects.** Emit `(1 | a:b)` /
   `(1 | a:b:c)`. Likely the quickest high-value backend unlock. *Unblocks GxE
   on a convergent HMC engine.*
3. **Heteroscedastic residual (`rcov` / `dsum`).** brms: a distributional
   `sigma ~ env`. INLA: a hierarchical / few-stratum reformulation (the full
   many-independent-variance form is outside INLA's envelope -- do not promise
   it). *Completes full-MET on brms.*
4. **INLA bug: interaction fixed effect on the binomial path.**
   `y ~ a*b, family = "binomial", backend = "inla"` fails inside INLA with
   `object 'a_b' not found`. The additive model works; the interaction column
   is not constructed safely. *A real bug, not a refusal.*
5. **Aggregated-binomial input.** No clean `cbind(success, failure)` or
   `trials =` on the main fit entry; the streaming path has `trials` but the
   modelling entry does not. Today the only working form is Bernoulli long
   expansion. *Usability.*
6. **greta convergence / auto-routing safety.** greta does not converge on
   high-dimensional crossed / interaction random effects, and it is a
   `backend = "auto"` fallback. Add a preflight risk classifier, route
   high-risk models to INLA/brms, and escalate catastrophic R-hat from a badge
   to a hard, unmissable diagnostic. *Safety; do not gate a release on solving
   greta's sampler.*

## Minor / environment-specific notes

- **greta readiness probe noise.** `fb_backend_status()` now captures the
  Python / TensorFlow discovery output, and `fb_backend_status(deep = FALSE)`
  skips that probe entirely (a fast, non-invasive check). On a misconfigured
  Python stack a *subprocess* launcher may still write to the OS console below
  the level R can capture; `deep = FALSE` avoids triggering it.
- **INLA hardware-probe chatter.** On some operating systems / INLA builds, the
  INLA *binary* prints harmless hardware-probe lines (for example `/bin/kstat:
  No such file or directory`) to the console during a fit. This is the INLA
  subprocess, not flexyBayes, and is below the level R can capture; it is
  cosmetic and does not affect results. Not reproducible on every platform.

## How to help

- Pick an issue above; the emit-layer ones (1, 2, 3-brms, 4) are the
  highest-value and lowest-risk. The IR / gate / emit / refusal architecture is
  already in place -- these are new term types and emit branches, not new
  subsystems.
- A fix is "done" only when `summary()`, `prior_summary()`, `canonical_names()`,
  `triangulate()`, and prediction all understand the new term. A half-plumbed
  capability is worse than a refusal.
- Validate against ASReml / `lme4` on a `barrero.maize` subset before claiming
  a class is supported, and add the result to the validation study.

The detailed issue write-ups (full repros, acceptance criteria, implementation
paths) and the strategy roadmap live in the project's release documentation.
