# flexyBayes -- known issues and open computational problems

This file ships with the package (`system.file("KNOWN_ISSUES.md", package =
"flexyBayes")`). It is the current statement of what flexyBayes fits, what it
refuses, and where the open computational problems sit, together with an
invitation to the team to help close them. It is written for the internal
release: the package validates on the model classes it supports, and states
plainly the ones it does not.

The thing to know before relying on flexyBayes for a multi-environment trial:
**the MET model -- interaction random effects together with a per-environment
residual -- fits on brms, and the evidence behind that row is small-scale.**
The combination generates both blocks in the emitted Stan program and samples
cleanly on simulated data (the live test) with the pieces corroborated on
`agridat::besag.met`; INLA refuses both halves structurally, and a national
trial series has not been run. The detail is below, so the scale question can
be answered collectively rather than rediscovered.

## Backend support by model class

The table below is generated from a single R-level source
(`.fb_capability_matrix()`) by `tools/generate_capability_matrix.R`, and every
verdict in it is re-derived from the gate and emit code by
`tests/testthat/test-capability-matrix.R`. A hand edit to the rows fails that
test. This replaces the hand-maintained matrix that shipped through 0.8.x,
which had drifted from the emit layer in four places.

<!-- capability-matrix:begin -->
| Model class | Spelling | INLA | brms | Notes |
|---|---|:-:|:-:|---|
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | fits | fits | The certified overlap class; both engines emit it and `triangulate()` compares them. |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | fits | fits | INLA's likelihood allowlist is read from `INLA::inla.models()` when INLA is installed. |
| Uncorrelated random slope | `(x \|\| g)` | refuses | fits | The INLA mapping named greta as one of its three verification arbitrators, so it stays deferred until the criterion is rebuilt around the active engines. The deferral is host-independent -- no local artefact lifts it. `auto` routes to brms. |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | refuses | fits | The indexed-slope INLA mapping shares the deferred three-arbitrator verification with the uncorrelated random slope, and refuses on every host. `auto` routes to brms. |
| Correlated random slope | `(x \| g)` | refuses | refuses | Refused at ingest, before any engine is chosen. Fit `(x \|\| g)` when the correlation is not of inferential interest. |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | refuses | fits | INLA collapses the finest strata, so it refuses rather than reporting a zero. brms emits `(1 \| a:b)`. |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | refuses | fits | One variance per level of `f`, no covariance between levels. All three spellings emit identical code. |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | refuses | fits | The correlated sibling of the diagonal structure -- `k(k+1)/2` parameters against `diag()`'s `k`. At one observation per cell the residual variance is confounded with the diagonal of the covariance: the covariance block converges and `sigma` does not, and a longer chain does not help. Replicate within cell, or put an informative prior on the residual. |
| Heterogeneous variances with one shared correlation | `~ corh(f):g` | refuses | refuses | No active engine has an equicorrelation group-level structure. Use `diag(f):g` or `us(f):g`. |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | refuses | fits | Lowered to distributional regression on log sigma, `sigma ~ 0 + f`. Refused for families with no residual scale. |
| Combined interaction random effects and heterogeneous residual (full MET) | `random = ~ gen + gen:env` with the `dsum` residual | refuses | fits | The emit carries both the group-level term and the `sigma` predictor, and a live fit samples cleanly on simulated multi-environment data. `auto` reaches brms for this class. |
| Factor-analytic genotype-by-environment covariance | `~ fa(env, k):gen` | refuses | refuses | Parsed for the formula catalogue and refused at dispatch -- no active engine emits a factor-analytic covariance. |
| Multi-trait covariance | `~ us(trait):vm(gen)` | refuses | refuses | No active engine represents a trait-by-genotype unstructured covariance. |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits | INLA takes the sparse-precision, pedigree-precision and block carriers; brms additionally takes dense and Cholesky. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses | A latent AR1 field plus the Gaussian observation nugget -- four hyperparameters, one observation per grid node. This is not ASReml's three-parameter nugget-free residual, so the residual spelling refuses and names this one. |
| Univariate P-spline | `~ spl(x)` | fits | refuses | Mapped to INLA's second-order random walk. brms has no lowering for the smooth basis. |
| Observation weights | `weights = w` | refuses | refuses | Parsed and recorded, consumed by no active emitter. A non-constant vector refuses rather than returning the unweighted posterior. |
| Exact sufficient-statistic aggregation | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a | Exact cell-likelihood aggregation for iid exponential-family models with small cell count. The brms path has no aggregated emit. |

`fits` -- the engine emits the structure and a test exercises it. `emits` -- the engine generates the structure and no live fit has yet confirmed it samples acceptably. `refuses` -- the request raises rather than fitting something else. `n/a` -- the class does not apply to that engine's interface.

This block is generated from `.fb_capability_matrix()` by `tools/generate_capability_matrix.R`. Edit the R table, re-run the generator, and let `tests/testthat/test-capability-matrix.R` check that every verdict still matches the gate and emit code. Do not edit the rows here by hand.
<!-- capability-matrix:end -->

`fits` does not promise that every fit converges at vignette-scale budgets.
Convergence is model- and engine-specific, and the package reports it -- treat
a printed fit carrying a high R-hat badge as a diagnostic, not a result.

## The MET boundary, stated plainly

A realistic multi-environment-trial model (the Barrero Model 1 shape) needs
two things beyond a simple random intercept:

1. **Interaction / nested random effects** (`gen:loc`, `gen:loc:yearf`).
2. **A heteroscedastic residual** (one error variance per environment).

Where each stands today:

- **brms** fits both, separately and together, and each is checked against
  an ASReml oracle. The interaction random effects emit as `(1 | a:b)` and
  recover every variance component against the ASReml / lme4 REML reference
  on `agridat::besag.met`. The heteroscedastic residual lowers to
  distributional regression on log sigma, `sigma ~ 0 + f`; fitted through
  `flexybayes()` on 360 simulated plots its per-site posterior-mean
  variances are 0.1146 / 1.1516 / 4.6981 against ASReml's
  0.1093 / 1.1248 / 4.6079, within 4.8% and largest on the smallest of the
  three. The `diag()` / `idh()` / `at()` and `us()` genotype variances emit
  and were validated on the parameter count before the values, which is
  what caught `lme4`'s `||` quietly fitting a correlated block for a factor
  slope.
- **The combination fits.** A model carrying interaction random effects and
  a sectioned residual at once generates both the group-level term and the
  `sigma` predictor in the emitted Stan program, and samples with
  acceptable diagnostics: on simulated multi-environment data
  (`tests/testthat/test-met-combined.R`) the fit returns every R-hat below
  1.05, no divergent transitions, and per-environment residual variances in
  the order they were simulated in. The capability row reads `fits` on that
  evidence, and `auto` reaches brms for this class. What is still worth
  saying is that `us(env):gen` is a different matter -- at one observation
  per genotype-environment cell its covariance diagonal is confounded with
  the residual, and a longer chain does not separate them.
- **INLA** refuses both by name. Interaction random effects are addressable
  in principle (an `iid` effect over the combined factor), and INLA
  collapses the finest strata on this model class, which is why the gate
  refuses rather than reporting a zero. The 107-stratum residual is **not**
  cleanly addressable: INLA integrates over hyperparameters numerically,
  which is tractable for roughly 15 to 20 of them, and 5 + 107 = 112 is far
  beyond that. This is a property of INLA, not of this code.

The direction of travel: **make brms the full-MET workhorse** -- the emit and
the live fit are both in place, so the remaining work is scale (the verified
fit is 120 rows, not a national trial series) -- and give **INLA** the
interaction random effects
plus a *hierarchically shrunk* residual (a prior
`log sigma^2_env ~ Normal(mu, tau)` replacing 107 free hyperparameters with
two, a model ASReml cannot fit and INLA can).

## Open issues (contributions welcome)

Each is reproducible on this release. Priority is for the MET use case.

1. **Take the combined MET fit to realistic scale.** The emit and one live
   fit are in place (120 rows, 4 chains, clean diagnostics). What is missing
   is the same check on a trial series of realistic size, where
   `adapt_delta` and the rest of Stan's `control` list start to matter.
   Those are reachable -- `flexybayes()` and `fb_brms()` forward `seed` and
   `control` to `brms::brm()` -- so the gap is the run, not the route.
   *Closes the package's headline gap.*
2. **INLA: interaction / nested random effects.** Map `a:b` to
   `f(interaction(a, b), model = "iid")` and characterise the finest-stratum
   collapse rather than refusing categorically. *Unblocks GxE on the fast
   engine.*
3. **Heteroscedastic residual on INLA.** A hierarchical / few-stratum
   reformulation. The full many-independent-variance form is outside INLA's
   envelope -- do not promise it.
4. **Exact aggregation of a factor interaction.** The per-row INLA path
   fits `y ~ a*b` on any family (see the resolved list below), but the
   aggregated path cannot: it expands the fixed effects to a model matrix
   and names its columns in the INLA formula, and INLA does not resolve a
   column called `a2:b2` even backticked. `aggregate = TRUE` on such a
   model refuses by name and points at the per-row route. The fix is a
   column-naming scheme carried through `summary.fixed`,
   `marginals.fixed` and the latent field together, so the posterior is
   never labelled with a synthesised token. *Restores compression on the
   replicated factorial where it is most useful.*
5. **Aggregated-binomial input.** No clean `cbind(success, failure)` or
   `trials =` on the main fit entry; the streaming path has `trials` but the
   modelling entry does not. Today the only working form is Bernoulli long
   expansion. *Usability.*
6. **Observation weights.** Parsed, recorded, and consumed by no emitter. The
   refusal is correct; the fix is an inverse-variance Gaussian emit checked
   against an analytic oracle, with the other weight senses (frequency,
   likelihood-power, trials, exposure) kept distinct rather than folded into
   one argument.

## Closed since the previous revision

Listed rather than deleted, because a reader who met one of these needs to
know it went away. Each was re-checked live against this tree.

- **INLA interaction fixed effect on the binomial path.** `y ~ a*b` with
  `family = "binomial"` and `backend = "inla"` was recorded as failing
  inside INLA with `object 'a_b' not found`. The per-row path fits it: on
  a replicated 4-by-4 factorial the sixteen posterior means agree with a
  `glm()` maximum-likelihood oracle to within 0.11 on the log-odds scale,
  the largest gap sitting on the sparsest cell.
- **Unrecognised ASReml functions.** The parser's default branch used to
  read an unknown call as a variable of that name, and `ar2`, `cor` and
  `str` reached an untyped `stop()`. The grammar is a closed set now:
  `foo(g)` refuses as `asreml_function_not_recognised`, and `ar2`, `str`
  and the `cor` family each refuse under their own code.

## Minor / environment-specific notes

- **INLA hardware-probe chatter.** On some operating systems and INLA builds,
  the INLA *binary* prints harmless hardware-probe lines (for example
  `/bin/kstat: No such file or directory`) to the console during a fit. This
  is the INLA subprocess, not flexyBayes, and is below the level R can
  capture; it is cosmetic and does not affect results. Not reproducible on
  every platform.
- **greta readiness probe noise.** `fb_backend_status()` captures the Python /
  TensorFlow discovery output, and `fb_backend_status(deep = FALSE)` skips
  that probe entirely (a fast, non-invasive check). On a misconfigured Python
  stack a *subprocess* launcher may still write to the OS console below the
  level R can capture; `deep = FALSE` avoids triggering it. greta is
  quarantined as a fitting engine, so this affects the status report only.

## How to help

- Pick an issue above; 1 and 5 are the highest-value and lowest-risk. The IR /
  gate / emit / refusal architecture is already in place -- these are new
  emit branches, reason codes and verification runs, not new subsystems.
- A fix is "done" only when `summary()`, `prior_summary()`,
  `canonical_names()`, `triangulate()`, and prediction all understand the new
  term, and when the capability row that claims it is backed by an assertion
  in `tests/testthat/test-capability-matrix.R`. A half-plumbed capability is
  worse than a refusal, and an unbacked table row is how this file went stale
  the first time.
- Validate against ASReml or `lme4` on a `barrero.maize` subset before
  claiming a class is supported, and add the result to the validation record.

The detailed issue write-ups (full repros, acceptance criteria, implementation
paths) and the strategy roadmap live in the project's release documentation.
