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
`agridat::besag.met`. INLA refuses both halves structurally, and a national
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
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | fits | fits | The certified overlap class, which both engines emit and `triangulate()` compares. |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | fits | fits | INLA's likelihood allowlist is read from `INLA::inla.models()` when INLA is installed. |
| Hurdle gamma (zero mass plus a positive gamma part) | `family = "hurdle_gamma"` | refuses | fits | brms-native (`dpars` mu, shape, hu); the zero-mass probability `hu` keeps brms's own prior. INLA's likelihood roster carries no counterpart, so the family gate refuses it there and `auto` routes to brms. |
| Uncorrelated random slope | `(x \|\| g)` | refuses | fits | The three-arbitrator verification named a since-withdrawn engine as one arbitrator, so the INLA mapping stays deferred until the criterion is rebuilt around the active engines. The deferral is host-independent -- no local artefact lifts it. `auto` routes to brms. |
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
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits | INLA takes the sparse-precision, pedigree-precision and block carriers, and brms additionally takes dense and Cholesky. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses | A latent AR1 field plus the Gaussian observation nugget -- four hyperparameters, one observation per grid node. This is not ASReml's three-parameter nugget-free residual, so the residual spelling refuses and names this one. |
| Univariate P-spline | `~ spl(x)` | fits | refuses | Mapped to INLA's second-order random walk. brms has no lowering for the smooth basis. |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits | Precision weighting, Var(y_i) = sigma^2 / w_i (the ASReml / lme4 / glm(weights=) sense): INLA's `scale = w`; on brms a known offset on the sigma distributional parameter, NOT brms's own `weights()` addition term (a different, likelihood-power quantity per brms's own documentation). Both engines match lme4::lmer(weights=) closely on a shared simulated fixture. Any other family, or a non-identity link on Gaussian, refuses by name (`weights_requires_gaussian`); `aggregate = TRUE` alongside weights also refuses by name (`weights_not_aggregatable`). |
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
  distributional regression on log sigma, `sigma ~ 0 + f`, and fitted through
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
   `trials =` on the main fit entry. The streaming path has `trials` but the
   modelling entry does not. Today the only working form is Bernoulli long
   expansion. *Usability.*
6. **Observation weights.** Parsed, recorded, and consumed by no emitter. The
   refusal is correct, and the fix is an inverse-variance Gaussian emit checked
   against an analytic oracle, with the other weight senses (frequency,
   likelihood-power, trials, exposure) kept distinct rather than folded into
   one argument.
7. **A bounded uniform prior whose upper bound cuts into the posterior can
   crash the INLA binary.** *Symptom*: the `inla` subprocess exits with a
   segmentation fault, INLA retries once with improved initial values and
   the retry crashes too, and `backend = "inla"` then raises `INLA fit
   failed: The inla-program exited with an error`. The R session itself is
   unaffected -- the crash is inside the subprocess, not in R -- so the
   failure arrives as an ordinary error and nothing is lost. *Trigger*: a
   `fb_prior()` bounded uniform on the standard-deviation scale whose upper
   bound sits inside the region the data support. Reproduced on a 72-row
   randomised block model whose `sd_Block` posterior mean is 2.28 with a
   97.5% bound of 6.35, priored `sd(group = "Block") ~ uniform(0, 3)`.
   *Workaround*: widen the bound -- `uniform(0, 10)` fits the same model
   cleanly -- or use a prior with no upper bound at all, such as
   `sd(group = "Block") ~ half_normal(scale = 3)`, which expresses the same
   scale belief without a hard edge. The package's own auto-default is a
   bounded uniform, so a user narrowing that bound is an ordinary path into
   this rather than an exotic one. On `backend = "auto"` the crash is not
   surfaced as a refusal: dispatch falls back to brms and the user reads a
   Stan fit where an INLA failure happened. *The crash is upstream. What is
   ours to fix is the silence around it.*
8. **The autoregressive field can lose its signal on an incomplete grid.**
   *Symptom*: the field standard deviation runs to a floor near 0.01, both
   correlations sit at approximately zero with credible intervals spanning
   almost the whole of `[-1, 1]`, and the variance the field should have
   carried is absorbed into the nugget, which rises correspondingly. The
   fit is then an independent-errors model wearing a spatial model's
   name. *Trigger*: unobserved plots. On one 12 by 10 grid with 12 of 120
   responses missing under the default `augment` policy, ten different
   hole patterns were fitted and five lost the field while five returned
   ordinary correlations -- the same data, the same model, the same
   settings, differing only in which plots were missing. Repeating a
   single fit three times on one unchanged data frame reproduced the same
   answer on some patterns and moved between a collapsed and a recovered
   field on others, so the outcome is not always a function of the input
   alone. *What the fit says about it*: nothing, before this release. The
   convergence block reported a converged mode and a passing numerical
   confirm throughout, and the augmentation record reported that it had
   completed the design; the credible intervals were the only tell. A fit
   whose field does not identify now raises a warning naming the
   parameters that lost their identification, and
   `summary(fit)$varcomp` carries `collapsed` in the `note` cell of the
   `sd_spatial` row. *Workaround*: complete the grid, either by supplying
   the field-book rows for every sown plot with the response set to `NA`
   or by padding the array with `fb_complete_grid()`; fit the same model
   on the brms engine for a second reading; or give the field an
   informative prior with `fb_prior()` rather than leaving it to run to a
   boundary. *What is still open* is the underlying identification, not
   the reporting of it: a hyperparameterisation or a default field prior
   under which an incomplete grid does not put the solution at a boundary
   in the first place.

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

## Response families: where the boundary is, and whose it is

The entry allowlist (`.resolve_family()`) is narrower than either engine's
own family roster. That is deliberate -- a family is admitted when this
package has an emit for it and a test exercises it -- but the refusal used
to say only "unsupported", which left a reader unable to tell a flexyBayes
boundary from an engine one. The four families a field engagement reached
for, each checked against the installed engines rather than recalled:

| Family | brms 2.23.0 | INLA 25.10.19 | flexyBayes | Where the boundary is |
|---|:-:|:-:|:-:|---|
| `hurdle_gamma` | native (`dpars` mu, shape, hu) | absent | **fits on brms since 0.9.2** | -- |
| `zero_inflated_gamma` | absent (`brmsfamily()` refuses) | absent | refuses | Neither engine has it. Use `hurdle_gamma`: for a gamma positive part, which has no mass at zero, the hurdle and the zero-inflated model are the same model. |
| `tweedie` | absent | **present** in the likelihood roster | refuses | **flexyBayes**, not the engine. No INLA emit (no link / power parameter, no validated fit). A feature request, not a defect -- see below. |
| `compound_poisson` | absent | present as `tweedie` | refuses | Same as `tweedie`; the two name the same distribution. |

Each of these refusals now carries the boundary note in its message, so the
distinction reaches the user at the point of refusal rather than only here.

**Tweedie on INLA is the open feature ask.** INLA carries the likelihood;
flexyBayes has no emit for it. The work is a family row, the `p` power
parameter threaded through `control.family`, and a validated fit against an
independent implementation (`statmod`/`cplm` or `mgcv::Tweedie`) before the
row can read `fits`. Until then, for compound-Poisson gamma data: fit the
zero and positive parts separately and recombine on the posterior, or use
`family = "hurdle_gamma"` with `backend = "brms"`, which is the same
two-part decomposition with the parts fitted jointly.

## Minor / environment-specific notes

- **INLA hardware-probe chatter.** On some operating systems and INLA builds,
  the INLA *binary* prints harmless hardware-probe lines (for example
  `/bin/kstat: No such file or directory`) to the console during a fit. This
  is the INLA subprocess, not flexyBayes, and is below the level R can
  capture, and it is cosmetic and does not affect results. Not reproducible on
  every platform.

## How to help

- Pick an issue above. Items 1 and 5 are the highest-value and lowest-risk. The IR /
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
