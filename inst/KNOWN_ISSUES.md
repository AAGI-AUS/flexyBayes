# flexyBayes -- known issues and open computational problems

This file ships with the package (`system.file("KNOWN_ISSUES.md", package =
"flexyBayes")`). It is the current statement of what flexyBayes fits, what it
refuses, and where the open computational problems sit, together with an
invitation to the team to help close them. It is written for the internal
release: the package validates on the model classes it supports, and states
plainly the ones it does not.

The thing to know before relying on flexyBayes for a multi-environment trial:
**the baseline heteroscedastic MET model -- iid genotype, iid
genotype-by-environment and a per-environment residual -- fits on both
engines, and the evidence behind that row is small-scale.** On brms the
combination generates both blocks in the emitted Stan program and samples
cleanly on simulated data (the live test), with the pieces corroborated on
`agridat::besag.met`; on INLA it recovers the simulated components (see the
table). Unstructured and factor-analytic genotype-by-environment covariance
combined with a sectioned residual do not plan on either engine, and a
national trial series has not been run. The detail is below, so the scale
question can be answered collectively.

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
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | fits | fits | The class both engines emit, and the one their agreement is checked on in the test suite. |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | fits | fits | INLA's likelihood allowlist is read from `INLA::inla.models()` when INLA is installed. |
| Hurdle gamma (zero mass plus a positive gamma part) | `family = "hurdle_gamma"` | refuses | fits | brms-native (`dpars` mu, shape, hu); the zero-mass probability `hu` keeps brms's own prior. INLA's likelihood roster carries no counterpart, so the family gate refuses it there and `auto` routes to brms. |
| Uncorrelated random slope | `(x \|\| g)` | fits | fits | Two independent iid fields on INLA, one for the intercept and one for the slope, checked against lme4's REML fit of the same model. Reported as `sd_<g>` and `sd_<x>_<g>` on INLA and as `sd_<g>` and `sd_<g>__<x>` on brms; triangulate() maps both to the INLA spelling. |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | fits | fits | An ordinary treatment-coded design matrix on both engines. The INLA posterior means sit on the least-squares estimates of the same model under a weak prior. |
| Correlated random slope | `(x \| g)` | refuses | fits | The uncorrelated slope with the correlation freed: one intercept SD, one slope SD and the correlation between them. brms writes `(x \| g)` verbatim and reports `sd_<g>`, `sd_<g>__<x>` and `cor_<g>__Intercept__<x>`. Estimated correlation matrices are fitted on brms only: INLA refuses by name (`correlated_covariance_requires_brms`) and `backend = "auto"` routes to brms. The INLA route was one iid2d field under a joint Wishart prior, and paired against ASReml 4.2.0.392 and lme4 on `lme4::sleepstudy` it returned a slope variance of 46.4 at the posterior mode against 35.08 and 35.07, while rescaling the response by ten moved the intercept-slope correlation from 0.02 to 0.969 where REML holds at 0.066 (oracle row 19, defect D1 of finding S01). |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | fits | fits | An iid field over the interaction index on INLA, `(1 \| a:b)` on brms -- the same model, and the two agree with lme4's REML fit of it on `agridat::john.alpha`. Only observed cells get a level, as in lme4 and brms. |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | fits | fits | One variance per level of `f`, no covariance between levels. All three spellings emit identical code. INLA writes k independent iid fields, one per level; brms writes `(0 + f \|\| g)`. Both report the component as `sd_<g>__<f><level>` and agree with lme4's REML fit of the same model. `at(f, level)` conditions on selected levels -- a different model -- and still refuses. |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | refuses | fits | `(0 + f \| g)` on brms. Estimated correlation matrices are fitted on brms only: INLA refuses by name (`correlated_covariance_requires_brms`) and `backend = "auto"` routes to brms. The INLA route was one iidKd field under a joint Wishart prior, and paired against ASReml 4.2.0.392 on `agridat::besag.met` it returned a genetic correlation of 0.142 where ASReml returns 0.930, while eight identical calls on the response rescaled by ten none converged and returned correlations from -0.894 to +0.686 (finding O04). The correlated sibling of the diagonal structure -- `k(k+1)/2` parameters against `diag()`'s `k`. At one observation per cell the residual variance is confounded with the diagonal of the covariance: the covariance block converges and `sigma` does not, and a longer chain does not help. Replicate within cell, or put an informative prior on the residual. |
| Heterogeneous variances with one shared correlation | `~ corh(f):g` | refuses | refuses | No active engine has an equicorrelation group-level structure. Use `diag(f):g` or `us(f):g`. |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | fits | fits | One residual variance per section. INLA carries it through the multiple-likelihood interface, a matrix response with one Gaussian family per level; brms lowers it to distributional regression on log sigma, `sigma ~ 0 + f`. Both report `sigma_<level>` on the standard-deviation scale. Refused for families with no residual scale; `at(f, level):units` conditions on selected levels, which is a different model, and refuses. |
| Baseline heteroscedastic MET: iid genotype, iid genotype-by-environment and a sectioned residual | `random = ~ gen + gen:env` with the `dsum` residual | fits | fits | The baseline multi-environment-trial model, shown on small simulated data; unstructured and factor-analytic genotype-by-environment covariance do not combine with the sectioned residual on either engine. Both halves refused on INLA before 0.10.0.9000, so the combination was unreachable there. It now recovers the genotype, genotype-by-environment and all three per-section residual standard deviations on simulated data. |
| Factor-analytic genotype-by-environment covariance | `~ fa(env, k):gen` | fits | refuses | Each common factor is one latent field entered per level of `f` as a SCALED COPY, and INLA's copy scaling IS the loading, so the loadings carry a posterior rather than a point estimate. Factor j is anchored at level j, the zero upper triangle ASReml uses against rotation, so an anchored loading is absent from the summary rather than reported as 1. Verified at k = 1 and k = 2. brms group-level effects are either uncorrelated or fully unstructured, with no factor-analytic form between. |
| Multi-trait covariance | `~ us(trait):vm(gen)` | refuses | refuses | No active engine represents a trait-by-genotype unstructured covariance. |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits | INLA takes the sparse-precision and pedigree-precision carriers and refuses a dense, Cholesky or block carrier by name; brms takes dense, Cholesky and precision and refuses a block carrier. Under `fb_select()`: shrinkage yes, projection no. A joint shrinkage prior over candidate fixed effects fits on a base carrying one of these terms, but the term lowers to `(1 \| gr(g, cov = K))` and projpred refuses a function call on the right-hand side of a group-term `\|` (probe C5), so the projection decision is red. `decision = "auto"`, the default, resolves to `full` there; naming `"projection"` refuses before the reference model is sampled. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses | A latent AR1 field plus the Gaussian observation nugget -- four hyperparameters, one observation per grid node. This is not ASReml's three-parameter nugget-free residual, so the residual spelling refuses and names this one. |
| Per-trial separable AR1 field | `random = ~ at(trial):ar1(row):ar1(col)` | fits | refuses | One field realisation per level of `trial`, via INLA's `replicate =` mechanism, but the row correlation, column correlation and field SD are shared across every level -- not estimated per trial. `at(trial, level):ar1(row):ar1(col)` (a level argument, asking for a single conditioned trial or for per-trial hyperparameters) refuses by name (`at_field_per_level_hyper_not_representable`). brms has no lowering for either spelling. |
| Grouped temporal AR1 field | `random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)` | fits | refuses | AR1 along `t` within each level of `g` and independent across levels: one correlation, one field SD and the Gaussian nugget, three hyperparameters shared by every level. The units-by-occasions repeated-measures covariance of Verbyla and Cullis (1992), written on the random side; an occasion a unit missed stays as a design row with the response NA under `na_action = "augment"`. Either factor order is accepted. brms has no lowering for the grouped autoregressive precision and refuses. |
| Univariate P-spline | `~ spl(x)` | fits | refuses | Mapped to INLA's second-order random walk. brms has no lowering for the smooth basis. An analogue of ASReml's spline rather than a translation of it -- basis, scaling and prediction rule differ, and no equivalence has been established (see TRANSLATION.md). |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits | Precision weighting, Var(y_i) = sigma^2 / w_i (the ASReml / lme4 / glm(weights=) sense): INLA's `scale = w`; on brms a known offset on the sigma distributional parameter, NOT brms's own `weights()` addition term (a different, likelihood-power quantity per brms's own documentation). Both engines match lme4::lmer(weights=) closely on a shared simulated fixture. Any other family, or a non-identity link on Gaussian, refuses by name (`weights_requires_gaussian`); `aggregate = TRUE` alongside weights also refuses by name (`weights_not_aggregatable`). |
| Exact sufficient-statistic aggregation (binomial, Poisson) | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a | Exact cell-likelihood aggregation for binomial (logit) and Poisson (log) models with a small cell count; the posterior is the per-row posterior. The marginal likelihood is not reported (`aggregated_evidence_unavailable`): the engine's value omits a constant that depends on the data. The brms path has no aggregated emit. |
| Exact sufficient-statistic aggregation (Gaussian) | `aggregate = TRUE`, `flexybayes_stream(family = "gaussian")` | refuses | n/a | Refused by name (`aggregation_residual_prior_not_composed`). The aggregated Gaussian route composes the within-cell sum of squares with a fixed prior on the residual precision, so a stated or default residual prior would not be the prior fitted. Gaussian models are fitted per row, which is what `aggregate = "auto"` does; `flexybayes_stream(fit = FALSE)` still returns the sufficient statistics. |
| Ridge shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ ridge(scale))` | n/a | fits | A joint shrinkage prior over a declared block of candidate fixed effects, compiled to one coefficient-specific `normal` row per candidate so the block is targeted exactly and no protected coefficient is touched. Only `fb_select()` takes a shrinkage prior, and it fits its reference model on brms whatever engine fitted the base, carrying the base's structure and priors through and recording the re-emission on the certificate. No call asks INLA to fit one, and INLA has no implementation of the block. |
| Laplace shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ laplace(scale))` | n/a | fits | One coefficient-specific `double_exponential` row per candidate. The posterior mean is not a lasso estimate and the draws carry no exact zeros; the alias `lasso()` says so once per session. Fitted on brms, as for ridge. A Laplace prior is Gaussian given one local scale per coefficient, so the model can be written as a latent Gaussian model, but it then carries a hyperparameter for every coefficient, and INLA's numerical integration is built for models with a few hyperparameters. |
| Elastic-net shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ elastic_net(lambda, alpha))` | n/a | fits | The canonical elastic-net density written as a `target +=` contribution over the candidate indices through `stanvars` and `check = FALSE`, so the block is targeted exactly while a protected coefficient keeps its own prior row. Fixed hyperparameters only in this release, because the normalising constant is then irrelevant. Fitted on brms, as for ridge. |
| Regularised horseshoe shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ regularized_horseshoe(p0))` | n/a | fits | brms's native `horseshoe()`, with `par_ratio` set from `p0` and calibrated against the number of resampling units rather than the number of rows. Exact block targeting ONLY when the base's fixed part is intercept-only: it is a class-wide prior, brms refuses a coefficient-specific row beside a class-wide special prior, so a base carrying protected fixed effects beyond the intercept refuses with `penalty_block_targeting_unsupported`. The two ways forward are `laplace()`, `ridge()` or `elastic_net()`, which target the block exactly, or `scope = "all_fixed"`, an opt-in under which the protected coefficients also sit under the global-local prior, the projection search forces them into every submodel, and the certificate records the named approximation `protected_under_global_shrinkage`. A candidate block containing a design column the re-emitted reference cannot standardise -- a factor contrast, an interaction product -- refuses with `penalty_candidate_columns_unscaled`. Fitted on brms, as for ridge; like the Laplace prior, it is Gaussian only given a local scale per coefficient. |
| R2D2 shrinkage prior on candidate fixed effects | `fb_penalty(term(...) ~ r2d2(mean_R2, prec_R2, cons_D2))` | n/a | fits | brms's native `R2D2()`, whose scales are stated on the explained-variance grammar rather than on a coefficient scale, so `p0` is refused for it. Class-wide, with the same targeting rule as the regularised horseshoe: exact only when the base's fixed part is intercept-only, otherwise `laplace()`, `ridge()` or `elastic_net()`, or `scope = "all_fixed"` with the named approximation. Unstandardised candidate columns refuse as for the horseshoe. Fitted on brms, as for ridge. |

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
  evidence. What is still worth
  saying is that `us(env):gen` is a different matter -- at one observation
  per genotype-environment cell its covariance diagonal is confounded with
  the residual, and a longer chain does not separate them.
- **INLA** fits both, separately and together, on the small simulated
  multi-environment data of `tests/testthat/test-sectioned-residual-inla.R`,
  where it recovers the genotype, genotype-by-environment and per-section
  residual standard deviations. Scale is the limit. INLA integrates over its
  hyperparameters numerically, and `lgm_gate()` refuses a model with more
  than 15 of them by name (`hyperparam_budget`), so a residual with one
  variance for each of 107 environments is refused on INLA.

The direction of travel: take the baseline MET model to realistic scale on
**brms**, where the emit and the live fit are both in place (the verified fit
is 120 rows, not a national trial series), and give **INLA** a
*hierarchically shrunk* residual, a prior `log sigma^2_env ~ Normal(mu, tau)`
that replaces 107 free hyperparameters with two.

## Open issues (contributions welcome)

Each is reproducible on this release. Priority is for the MET use case.

1. **Estimated correlation matrices are fitted on brms only.** `us(f):g`
   and the correlated random slope `(x | g)` refuse on INLA by name
   (`correlated_covariance_requires_brms`), and `backend = "auto"` routes
   both to brms. On INLA each was one `iidKd` or `iid2d` field under a joint
   Wishart prior, and neither reproduced REML's covariance. `us(env):gen` on
   `agridat::besag.met` returned a genetic correlation of 0.142 where ASReml
   returns 0.930 (finding O04), and `(Days | Subject)` on `sleepstudy` a
   slope variance of 46.4 at the posterior mode against 35.08 in ASReml and
   35.07 in lme4 (defect D1 of finding S01). Rescaling the response by ten
   moved both correlations a long way. Under an identical Wishart prior INLA
   reproduces Stan to within 1.5 to 6.5 per cent, and no Wishart scale gave
   REML-sized variances with a correlation left to the data, so the joint
   prior couples the variances and the correlations whatever the response
   scale. Readmitting either structure on INLA would take a prior that
   separates the two, which INLA's built-in models do not provide.
2. **Take the combined MET fit to realistic scale.** The emit and one live
   fit are in place (120 rows, 4 chains, clean diagnostics). What is missing
   is the same check on a trial series of realistic size, where
   `adapt_delta` and the rest of Stan's `control` list start to matter.
   Those are reachable -- `flexybayes()` and `fb_brms()` forward `seed` and
   `control` to `brms::brm()` -- so the gap is the run, not the route.
   *Closes the package's headline gap.*
3. **INLA: interaction / nested random effects.** Map `a:b` to
   `f(interaction(a, b), model = "iid")` and characterise the finest-stratum
   collapse rather than refusing categorically. *Unblocks GxE on the fast
   engine.*
4. **Heteroscedastic residual on INLA.** A hierarchical / few-stratum
   reformulation. The full many-independent-variance form is outside INLA's
   envelope -- do not promise it.
5. **Exact aggregation of a factor interaction.** The per-row INLA path
   fits `y ~ a*b` on any family (see the resolved list below), but the
   aggregated path cannot: it expands the fixed effects to a model matrix
   and names its columns in the INLA formula, and INLA does not resolve a
   column called `a2:b2` even backticked. `aggregate = TRUE` on such a
   model refuses by name and points at the per-row route. The fix is a
   column-naming scheme carried through `summary.fixed`,
   `marginals.fixed` and the latent field together, so the posterior is
   never labelled with a synthesised token. *Restores compression on the
   replicated factorial where it is most useful.*
6. **Aggregated-binomial input.** No clean `cbind(success, failure)` or
   `trials =` on the main fit entry. The streaming path has `trials` but the
   modelling entry does not. Today the only working form is Bernoulli long
   expansion. *Usability.*
7. **Observation weights outside the Gaussian/identity case.** `weights =`
   now fits (0.9.3) for the Gaussian family on an identity link, on both
   engines, in the ASReml / lme4 / `glm(weights=)` precision sense --
   checked against an analytic oracle and matching `lme4::lmer(weights=)`
   closely on a shared simulated fixture (see the capability matrix
   above). What remains open: any other family, and a non-identity link
   on Gaussian, still refuse by name (`weights_requires_gaussian`)
   rather than lowering to that family's own weight sense (frequency,
   likelihood-power, trials, exposure), and `aggregate = TRUE` alongside
   weights refuses (`weights_not_aggregatable`) because the closed-form
   sufficient statistics do not carry a per-observation weight through
   the compression.
8. **A bounded uniform prior whose upper bound cuts into the posterior can
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
   this rather than an exotic one. On `backend = "auto"` the crash is
   still not a typed refusal: dispatch re-routes to brms and returns a
   Stan fit where an INLA failure happened. It is no longer silent --
   the re-route emits a one-shot message naming the numerical fallback,
   and the returned object's class and `backend` field record which
   engine produced it, so a fit can always be traced to its engine.
   *The crash is upstream. What remains ours is that the re-route is a
   message rather than a refusal you have to acknowledge.* Two changes in
   0.12.0 bear on it. On `backend = "inla"` the refusal now says when the
   model states a uniform bound, and names the two remedies above. A bound
   tighter than the value INLA starts its search from is now enforced:
   before 0.12.0 the search began outside the bound, where the prior is
   flat, and never entered it, so `sd(group = "geno") ~ uniform(0, 0.05)`
   on a Poisson model returned a standard deviation of 0.61. The search
   now starts inside the bound, and the same fit returns 0.05.
9. **The autoregressive field can lose its signal on an incomplete grid.**
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
10. **Large Gaussian data has no aggregated route in this release.**
    `aggregate = TRUE` and `flexybayes_stream()` refuse a Gaussian model
    (`aggregation_residual_prior_not_composed`), because the aggregated
    Gaussian likelihood combined the within-cell sums of squares with a
    fixed prior on the residual precision and so did not fit the residual
    prior `prior_summary()` reported. A Gaussian model is fitted per row,
    which on a single random intercept completed at five million rows in
    19.5 GB (see `inst/validation/benchmark_scaling.md`). *What is open*:
    an aggregated Gaussian likelihood that composes the stated or default
    residual prior exactly and restores the marginal likelihood, planned
    for 0.13.0.
11. **Closed: held-out clusters under `prediction = "new_cluster"` are
    relabelled for every grouping term, or the term refuses by name.**
    `fb_select()` gives the assessment rows new labels for every column a
    random term's effect is indexed by: the group of a random intercept, of
    a random slope and of an `ide()` term, every factor of a two-way or
    multi-way interaction, and the inner factor of `diag(f):g`, `us(f):g`
    and `fa(f, k):g`. A `vm()`, `ped()` or autoregressive term cannot take
    a new level and refuses with `selection_new_cluster_not_representable`;
    `prediction = "known_cluster"` keeps its levels. A base written in the
    brms grammar, such as `y ~ x + (x || g)`, carries its random terms into
    every fold refit. On the benchmark route a two-way interaction is one
    covariance term, and every random term other than an intercept or a
    two-way interaction (a random slope, a multi-way interaction,
    `diag(f):g`, `vm()`) refuses with `selection_covariance_not_buildable`.
12. **Closed: on INLA, a fold refit's default priors read its training
    rows.** Each INLA fold refit in `fb_select()` receives the training rows
    and the assessment rows with the response set to missing. The default
    bound on an uncorrelated random slope's standard deviation and the
    default fixed-effect prior are scaled by covariate standard deviations
    taken over the rows with an observed response, which in a fold refit are
    its training rows, the rows the response-scaled defaults already read. A
    fit with `na_action = "augment"` reads the same rows, so the covariates
    of rows without a response leave its defaults unchanged.
13. **The selection certificate identifies columns by name.** The score
    contract that decides whether two arms' out-of-fold scores are
    comparable hashes column names: rescaling a covariate between arms
    (`x1` to `1000 * x1 + 17`) leaves its digest unchanged, and a
    transformed response is summarised without its values. *What is open*:
    hash the values the arms were scored on.
14. **The boundary warning on a random slope depends on the covariate's
    units.** The warning for a variance component pinned at zero compares
    the slope's standard deviation with a threshold that does not account
    for the scale of its covariate, so the same model can warn or not
    depending on the units the covariate is recorded in. *What is open*: a
    term-aware threshold.
15. **Blocks-format `vm()` and `ped()` carriers are refused on INLA.** A
    blocks carrier is one block-diagonal covariance with one standard
    deviation, and INLA fits each block as its own field with its own
    standard deviation, so with two or more blocks the fit was a different
    model, reported as one variance-component row per block. With one block
    the fit equals the precision carrier of the same matrix. Checked on a
    30-level pedigree against brms and against the exact posterior, its
    interval for the relationship standard deviation is 1 to 1.5 per cent
    narrower at each end and the residual standard deviation sits 0.5 to
    1.5 per cent low; that is INLA's approximation, and the precision and
    pedigree-precision carriers show the same. `backend = "inla"` refuses a
    blocks carrier with `inla_gate_refused` (rule `random_term_type_inla`).
    *Workaround*: the precision carrier of the block-diagonal inverse,
    `vm(geno, cov = fb_cov(Q, type = "precision"))` with
    `Q = Matrix::bdiag(lapply(blocks, solve))`, or the dense carrier of the
    block-diagonal matrix on brms. *What is open*: a blocks carrier fitted
    with one shared standard deviation, and an accuracy bound for INLA's
    variance components.
16. **A Poisson response that is almost all zeros, with a weakly identified
    group, can crash INLA's program, and `backend = "inla"` refuses it.** On
    a Poisson model with a random intercept over 28 clusters of two and 3
    non-zero counts in 56 rows, INLA's program stops with a segmentation
    fault after the mode is found, while it integrates over the group
    precision. It happens on every release from 0.10.0 to 0.12.0 and under
    INLA 26.08.07 and 25.10.19, and the package refuses with
    `inla_program_failed`. The same design with 7 non-zero counts fits. It
    is not a regression: the execution-grid cell that recorded 0.10.0 as
    fitting drew the 7-count response, and a fixture added to the grid
    later moved the random-number state so the cell drew the 3-count one.
    The cause inside INLA is not established. *Workaround*: state
    `sd(group = ...) ~ half_normal(scale = 1.53)`, which fitted the 3-count
    design in 3 of 3 runs, or fit on brms.
17. **On INLA, a hyperparameter marginal can lie far outside a stated
    uniform bound, and the fit is refused.** On a Poisson model with two
    random intercepts, one of which splits the rows exactly as a fixed
    factor does, under the default uniform bound of 3 on each standard
    deviation, the duplicated term's mode lies on the bound. INLA estimates
    curvature by stepping across the bound, where the prior's log-density
    is -1e10. Under INLA 26.08.07, in 45 of 100 identical runs the marginal
    for the other term's precision had its mode near 5e-07, a standard
    deviation near 1,400, while INLA's optimiser mode and brms (median 0.53,
    95% interval 0.17 to 2.11) were inside the bound; under INLA 25.10.19,
    100 of 100 runs stayed inside. 0.12.0 returned such a fit with a
    warning. From 0.12.1 every INLA fit is checked against INLA's own record
    of its priors, a fit whose mode or reported interval lies outside a
    bound is refused with `prior_support_violated`, and the check is kept in
    `fit$extras$support_audit`. *Workaround*: drop a random term that
    duplicates a fixed factor, state
    `sd(group = ...) ~ half_normal(scale = 1.53)` (no runaway fit in 24 runs
    of the same design; posterior median 0.505 against 0.489 on brms), or
    fit on brms.

    The same check refuses a fit whose mode is inside the bound when INLA's
    reported interval reaches past it, which happens on a weakly identified
    component with few levels. On `agridat::john.alpha` with
    `random = ~ rep + rep:block`, the default bound on `rep` (three levels)
    is 3.05 and INLA's interval runs from 0.17 to 3.65, so the default INLA
    fit is refused in 8 of 8 repeats. *Workaround*: state a half-normal prior
    on the standard deviations, for example with scale the standard
    deviation of the response, or fit on brms. *What is open*: the default
    prior (the release plan moves it to a half-normal prior on both engines
    in 0.13.0).
18. **Formula transforms that read other rows are refused on the Bayesian
    selection route.** `fb_select()` refits the model formula on every fold,
    and INLA evaluates it over the fold's training and assessment rows
    together, so `scale(z1)` in a candidate formula was centred and scaled
    with the assessment rows included: in one fold refit a single assessment
    value moved the coefficient of `scale(z1)` from 0.10 to 0.52. The
    Bayesian route (the penalty-free fold refits, the out-of-fold score of
    `benchmark = TRUE`, and the projection's k-fold search) now refuses any
    formula variable whose value on a row depends on other rows with
    `selection_formula_transform_not_fold_local`, before anything is fitted.
    The benchmark workflows learn `scale()`, `poly()`, `splines::ns()` and
    `splines::bs()` on each fold's training rows and apply them to its
    assessment rows, and refuse any other such transform, such as
    `I(x - mean(x))`. *Workaround*: name the untransformed column (a
    penalty's candidate standardisation is learned inside each fold), or
    write the transform as a rule of one row at a time, such as `I(x^2)`.
    *What is open*: carrying the training rows' transform parameters into
    the INLA and brms fold refits.
19. **A fold drop the route cannot write is refused.** `fb_select()` fixes
    each fold's candidate space from the covariates on its training rows
    and drops a candidate column those rows leave unidentified, with `NA`
    for its coefficient and the drop on the certificate. Two routes cannot
    drop a single column, and refuse with
    `selection_fold_candidate_drop_not_representable` before anything is
    fitted. The Bayesian fold refits write the model as formula text, so a
    factor or interaction term with only some of its columns unidentified
    on a fold is refused. The penalised reference's out-of-fold score
    (`benchmark = TRUE` with a penalty) compiles its prior over the whole
    candidate block, so any drop there is refused. *Workaround*: score the
    candidates with `workflow = "unpenalized_benchmark"` or the `glmnet`
    arm, which drop single columns; cut the folds on a unit whose training
    rows identify every candidate; or remove the unidentified columns from
    `candidates`. *What is open*: single-column drops in the formula
    refits, and a reference prior compiled over each fold's kept columns.

## Closed since the previous revision

Listed rather than deleted, because a reader who met one of these needs to
know it went away. Each was re-checked live against this tree.

- **Gaussian fits on INLA at large response scale.** Every Gaussian model on
  INLA was measured wrong once the response standard deviation reached the
  low hundreds: on `lme4::sleepstudy` the plain random intercept returned
  ratios of the subject variance to the residual between 0.03 and 6.99
  against REML's 1.198, several of those fits passing the numerical-confirm
  gate. Every Gaussian INLA fit now runs on a standardised response, and the
  result is mapped back to the original units exactly. The default priors
  the package chooses are scaled by the response and by each covariate. On
  `lme4::sleepstudy` the default-prior fit agrees across response
  multipliers from 0.01 to 1000. The refusal of a Gaussian INLA fit whose
  response standard deviation exceeds 150
  (`inla_response_scale_unvalidated`) is gone: an explicit
  `backend = "inla"`, `fb_inla()`, `backend = "auto"`, the aggregated route
  and the streamed route each fit a response on a thousandfold scale and
  agree with the same model fitted on the original one. The checks are in
  `tests/testthat/test-response-standardisation.R`,
  `tests/testthat/test-inla-scale-aware-defaults.R` and
  `tests/testthat/test-inla-containment.R`.
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
- **Per-trial separable AR1 field.** The single-lattice AR1 field
  (`random = ~ ar1(row):ar1(col)`, one trial, one grid) had no
  multi-environment spelling. `random = ~ at(trial):ar1(row):ar1(col)`,
  with no `level` argument on `at()`, now fits on INLA via the
  `replicate =` mechanism: one field realisation per level of `trial`,
  but the field's row correlation, column correlation and SD are
  **shared across every level** -- not estimated per trial. The
  boundary is where it was: `at(trial, level):ar1(row):ar1(col)` (a
  `level` argument, asking for a single conditioned trial or for
  per-trial hyperparameters) refuses
  (`at_field_per_level_hyper_not_representable`), and brms has no
  lowering for either spelling.
- **Units followed through time.** `random = ~ ar1(t):g`, in either
  factor order, fits one AR1 series per level of `g` on INLA,
  independent across levels, with the correlation and field SD shared
  by every level: the Verbyla-Cullis repeated-measures covariance on
  the random side. The model line, the field block and the note naming
  the ASReml residual the field is not used to describe this field as
  `ar1(t):ar1(g)` with four parameters; they now describe it as what it
  is (AR1 along `t`, independent across `g`, three hyperparameters).
  brms refuses.

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
  `canonical_names()` and prediction all understand the new term, and when
  the capability row that claims it is backed by an assertion
  in `tests/testthat/test-capability-matrix.R`. A half-plumbed capability is
  worse than a refusal, and an unbacked table row is how this file went stale
  the first time.
- Validate against ASReml or `lme4` on a `barrero.maize` subset before
  claiming a class is supported, and add the result to the validation record.

The detailed issue write-ups (full repros, acceptance criteria, implementation
paths) and the strategy roadmap live in the project's release documentation.
