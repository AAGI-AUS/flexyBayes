# flexyBayes: known issues and open computational problems

This file ships with the package (`system.file("KNOWN_ISSUES.md", package =
"flexyBayes")`). It states what flexyBayes fits, what it refuses, and where
the open computational problems are. Contributions to any of them are
welcome.

Before relying on flexyBayes for a multi-environment trial, note that **the
baseline heteroscedastic MET model (iid genotype, iid
genotype-by-environment and a per-environment residual) fits on both
engines, but the evidence for that row is small-scale.** On brms the emitted
Stan program contains both blocks and samples cleanly on simulated data (the
live test), and the separate parts are corroborated on `agridat::besag.met`.
On INLA it recovers the simulated components (see the table). `us(env):gen`
with a sectioned residual is written out on brms but has not been checked
against a reference fit, `fa()` is refused on both engines, and a national
trial series has not been run.

## Backend support by model class

The table below is generated from a single R-level source
(`.fb_capability_matrix()`) by `tools/generate_capability_matrix.R`, and
`tests/testthat/test-capability-matrix.R` re-derives every verdict in it from
the gate and emit code. A hand edit to the rows fails that test. A note that
says INLA or brms cannot do something is checked against the engine itself
in `tests/testthat/test-engine-claims.R`. Every other refusal is
flexyBayes's own, and its note says so. The README shows the rows that at
least one engine fits, without the notes.

<!-- capability-matrix:begin -->
<!-- Generated from .fb_capability_matrix() in R/capability_matrix.R by tools/generate_capability_matrix.R. Edit the R table and re-run the generator; tests/testthat/test-capability-matrix.R fails on a hand edit here. -->
| Model class | Spelling | INLA | brms | Notes |
|---|---|:-:|:-:|---|
| Gaussian LMM, simple random intercept | `random = ~ g` / `(1 \| g)` | fits | fits | The class both engines fit, and the one on which the tests check that the two agree. |
| GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect | `(1 \| g)` with `family =` | fits | fits | On INLA each family is checked against `INLA::inla.models()` when INLA is installed. |
| Hurdle gamma (zero mass plus a positive gamma part) | `family = "hurdle_gamma"` | refuses | fits | A brms family (parameters mu, shape and hu); `hu` keeps brms's own prior. INLA has no hurdle-gamma likelihood, so the INLA route refuses it and `backend = "auto"` sends it to brms. |
| Uncorrelated random slope | `(x \|\| g)` | fits | fits | Two independent iid fields on INLA, one for the intercept and one for the slope, checked against lme4's REML fit of the same model. Reported as `sd_<g>` and `sd_<x>_<g>` on INLA and as `sd_<g>` and `sd_<g>__<x>` on brms; `canonical_names()` maps both to the INLA names. |
| Factor-by-numeric fixed interaction | `y ~ f * x` with numeric `x` | fits | fits | An ordinary treatment-coded design matrix on both engines. Under a weak prior the INLA posterior means sit on the least-squares estimates of the same model. |
| Correlated random slope | `(x \| g)` | refuses | fits | One intercept SD, one slope SD and the correlation between them. brms takes `(x \| g)` as written and reports `sd_<g>`, `sd_<g>__<x>` and `cor_<g>__Intercept__<x>`. flexyBayes fits estimated correlation matrices on brms only, and `backend = "auto"` sends them there. On INLA the structure was one iid2d field under a joint Wishart prior. Paired against ASReml 4.2 and lme4 on `lme4::sleepstudy` it gave a slope variance of 46.4 at the posterior mode, against 35.08 and 35.07, and rescaling the response by ten moved the intercept-slope correlation from 0.02 to 0.969, where REML stays at 0.066. |
| Nested / interaction random effects, multi-stratum | `~ gen:env`, `~ env:rep:block` | fits | fits | An iid field over the interaction index on INLA and `(1 \| a:b)` on brms. They are the same model, and both agree with lme4's REML fit of it on `agridat::john.alpha`. Only observed cells get a level, as in lme4 and brms. |
| Heterogeneous variance by factor level | `~ diag(f):g`, `~ idh(f):g`, `~ at(f):g` | fits | fits | One variance per level of `f`, no covariance between levels. The three spellings give identical code: k independent iid fields on INLA, `(0 + f \|\| g)` on brms. Both report `sd_<g>__<f><level>` and agree with lme4's REML fit of the same model. `at(f, level):g` keeps the named levels of `f` only, each with its own variance, and agrees with ASReml on both engines. |
| Unstructured genotype-by-environment covariance | `~ us(f):g` | refuses | fits | Written `(0 + f \| g)` on brms. flexyBayes fits estimated correlation matrices on brms only, and `backend = "auto"` sends them there. On INLA the structure was one iidKd field under a joint Wishart prior. On `agridat::besag.met` it gave a genetic correlation of 0.142 where ASReml gives 0.930, and eight identical calls on the response rescaled by ten did not converge and returned correlations from -0.894 to +0.686. The structure has k(k+1)/2 parameters against the k of `diag()`. With one observation per cell the residual variance is confounded with the diagonal of the covariance, and a longer chain does not separate them: replicate within cells, or put an informative prior on the residual. |
| Heterogeneous variances with one shared correlation | `~ corh(f):g` | refuses | refuses | flexyBayes has no lowering for this structure on either engine. With one variance for every level and a correlation of zero or more, the same covariance is `random = ~ g + f:g`, which fits on both engines. INLA has an exchangeable group model with one variance and one shared correlation, which flexyBayes does not use. Use `diag(f):g` for separate variances without correlation, or `us(f):g` on brms for every correlation. |
| Heterogeneous residual by factor level | `residual = ~ dsum(~ units \| f)` / `~ at(f):units` | fits | fits | One residual variance per section. INLA fits one Gaussian likelihood per level; brms fits `sigma ~ 0 + f`. Both report `sigma_<level>` on the standard-deviation scale. With `diag(f):g` beside it past 15 variance parameters, INLA fits by empirical Bayes; on 12 simulated trials every per-trial variance ASReml identifies is within 0.17 of its ASReml standard error. Refused for a family with no residual scale, and for `at(f, level):units`, which conditions on selected levels. |
| Baseline heteroscedastic MET: iid genotype, iid genotype-by-environment and a sectioned residual | `random = ~ gen + gen:env` with the `dsum` residual | fits | fits | The baseline multi-environment-trial model. On small simulated data both engines recover the genotype, genotype-by-environment and per-section residual standard deviations. `us(env):gen` with the sectioned residual is written out on brms but has not been checked against a reference fit, and `fa()` is refused on both engines. |
| Factor-analytic genotype-by-environment covariance | `~ fa(env, k):gen` | refuses | refuses | flexyBayes does not fit this structure on either engine. INLA can represent it, one latent field per common factor entered at each level as a scaled copy, but on simulated trials with a known covariance the genetic-correlation intervals of INLA fits covered the true value 134 of 150 times for k = 1 and 125 of 150 for k = 2, against 136 required. Each factor's loading is fixed at one at an anchor environment, chosen from the data as the environment with the largest common genetic variance. flexyBayes has no factor-analytic lowering for brms. `diag(env):gen` fits on both engines, and `us(env):gen` on brms. |
| Multi-trait covariance | `~ us(trait):vm(gen)` | refuses | refuses | flexyBayes has no lowering for this structure on either engine. brms writes it directly, as `(0 + trait \| gr(gen, cov = K))` on data with one row per trait or with `mvbind()` on data with one column per trait. INLA can build it from a known-covariance field for each trait joined by a scaled `copy`. flexyBayes emits neither. Fit it with brms directly, or fit each trait separately. |
| Known-covariance genomic / pedigree random effect | `~ vm(g, K)`, `~ ped(a, A)` | fits | fits | On INLA flexyBayes takes the sparse-precision and pedigree-precision carriers and refuses a dense, Cholesky or blocks carrier by name; on brms it takes dense, Cholesky and precision carriers and refuses a blocks carrier. |
| Separable AR1 spatial field | `random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)` | fits | refuses | A latent AR1 field plus the Gaussian observation nugget: four hyperparameters, one observation per grid node. ASReml's three-parameter residual without a nugget is written `residual = ~ ar1(row):ar1(col)` and fits on INLA for a Gaussian response; the 1-D residual `~ ar1(t)` refuses. flexyBayes has no brms lowering for either. |
| Per-trial separable AR1 field | `random = ~ id(trial):ar1(row):ar1(col)` | fits | refuses | One field realisation per level of `trial`, through INLA's `replicate =` mechanism. The row correlation, column correlation and field SD are shared by every level, not estimated per trial. `at(trial):ar1(row):ar1(col)` fits a field with its own SD and correlations on each trial, as ASReml does, and `at(trial, level):ar1(row):ar1(col)` on the named trials only; both agree with ASReml. flexyBayes has no brms lowering for the field. |
| Grouped temporal AR1 field | `random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)` | fits | refuses | AR1 along `t` within each level of `g`, independent across levels: one correlation, one field SD and the Gaussian nugget, shared by every level. This is the units-by-occasions repeated-measures covariance of Verbyla and Cullis (1992), written on the random side; an occasion a unit missed stays as a row with the response `NA` under `na_action = "augment"`. Either factor order is accepted. flexyBayes has no brms lowering for it. |
| Univariate P-spline | `~ spl(x)` | fits | refuses | Fitted on INLA as a second-order random walk. It is an analogue of ASReml's spline, not a translation: the basis, scaling and prediction rule differ, and no equivalence has been established (see TRANSLATION.md). flexyBayes has no brms lowering for it. |
| Observation weights (Gaussian, identity link) | `weights = w` | fits | fits | Precision weights, Var(y_i) = sigma^2 / w_i, as in ASReml, lme4 and `glm(weights = )`: INLA's `scale = w`, and on brms a known offset on sigma instead of brms's own `weights()` term, which multiplies each observation's log-likelihood instead. Both engines match `lme4::lmer(weights = )` closely on a shared simulated fixture. Any other family, or a non-identity link on the Gaussian, refuses by name, and so does `aggregate = TRUE` with weights. |
| Exact sufficient-statistic aggregation (binomial, Poisson) | `aggregate = TRUE`, `flexybayes_stream()` | fits | n/a | Exact cell-likelihood aggregation for binomial (logit) and Poisson (log) models with few cells; the posterior equals the per-row posterior. The marginal likelihood is not reported, because the engine's value omits a constant that depends on the data. flexyBayes has no aggregated route on brms. |
| Exact sufficient-statistic aggregation (Gaussian) | `aggregate = TRUE`, `flexybayes_stream(family = "gaussian")` | refuses | n/a | Refused by name. The aggregated Gaussian route combined the within-cell sums of squares with a fixed prior on the residual precision, so a stated or default residual prior would not be the prior fitted. Gaussian models are fitted per row, which is what `aggregate = "auto"` does; `flexybayes_stream(fit = FALSE)` still returns the sufficient statistics. |

`fits`: the engine fits it, and a test checks the fit. `refuses`: flexyBayes stops with a message saying why and what to use instead. `n/a`: the request does not apply to that engine.
<!-- capability-matrix:end -->

`fits` does not promise that every fit converges at the sampling budgets used
in the vignettes. Convergence depends on the model and the engine, and the
package reports it. Treat a printed fit with a high R-hat badge as a
diagnostic, not a result.

## The multi-environment trial boundary

A realistic multi-environment-trial model (the shape of Barrero Model 1)
needs two things beyond a simple random intercept:

1. **Interaction or nested random effects** (`gen:loc`, `gen:loc:yearf`).
2. **A heteroscedastic residual** (one error variance per environment).

Where each stands:

- **brms** fits both, separately and together, and each is checked against
  ASReml. The interaction random effects are emitted as `(1 | a:b)` and
  recover every variance component of the ASReml and lme4 REML reference on
  `agridat::besag.met`. The heteroscedastic residual becomes a distributional
  regression on log sigma, `sigma ~ 0 + f`. Fitted through `flexybayes()` on
  360 simulated plots, its per-site posterior-mean variances are
  0.1146 / 1.1516 / 4.6981 against ASReml's 0.1093 / 1.1248 / 4.6079, within
  4.8% and furthest apart on the smallest of the three. The `diag()`,
  `idh()`, `at()` and `us()` genotype variances are emitted and were checked
  on the number of parameters before the values. That check found that
  `lme4`'s `||` silently fits a correlated block for a factor slope.
- **The combination fits.** A model with interaction random effects and a
  sectioned residual together produces both the group-level term and the
  `sigma` predictor in the emitted Stan program, and samples with acceptable
  diagnostics. On simulated multi-environment data
  (`tests/testthat/test-met-combined.R`) every R-hat is below 1.05, there are
  no divergent transitions, and the per-environment residual variances come
  back in the order they were simulated in. The capability row reads `fits`
  on that evidence. `us(env):gen` is different. At one observation per
  genotype-environment cell its covariance diagonal is confounded with the
  residual, and a longer chain does not separate them.
- **INLA** fits both, separately and together, on the small simulated
  multi-environment data of `tests/testthat/test-sectioned-residual-inla.R`,
  where it recovers the genotype, genotype-by-environment and per-section
  residual standard deviations. INLA integrates over its hyperparameters
  numerically, and flexyBayes refuses an INLA model with more than 15 of
  them by name (`hyperparam_budget`), counting one per section of a
  sectioned residual. Past 15, a model whose excess is one variance per
  level of `diag(f):g` or `dsum(~ units | f)` is fitted by empirical Bayes
  (`int.strategy = "eb"`). The fixed and random effects are then conditional
  on the posterior mode of the variance parameters, each variance keeps its
  own posterior, and the printed fit says so. On 12 simulated trials with
  both per-trial structures (26 variance parameters,
  `inst/extdata/model_c_reference.rds`), every per-trial variance that ASReml
  identifies is within 0.17 of its ASReml standard error.

The next step is to take the baseline MET model to realistic scale on
**brms**, where the emit and the live fit are both in place (the verified fit
has 120 rows, not a national trial series). A *hierarchically shrunk*
per-trial variance, with a prior `log sigma_trial ~ Normal(mu, tau)`, is not
offered. Fitted through brms on the same 12 trials, it moved well-identified
residual variances by up to 0.8 of an ASReml standard error. That is partial
pooling, and a different model from ASReml's.

## Open issues (contributions welcome)

Each is reproducible on this release. The multi-environment trial use case
has priority.

1. **Estimated correlation matrices are fitted on brms only.** `us(f):g`
   and the correlated random slope `(x | g)` are refused on INLA by name
   (`correlated_covariance_requires_brms`), and `backend = "auto"` sends
   both to brms. On INLA each was one `iidKd` or `iid2d` field under a joint
   Wishart prior, and neither reproduced REML's covariance. `us(env):gen` on
   `agridat::besag.met` gave a genetic correlation of 0.142 where ASReml
   gives 0.930, and `(Days | Subject)` on `sleepstudy` gave a slope variance
   of 46.4 at the posterior mode against 35.08 in ASReml and 35.07 in lme4.
   Rescaling the response by ten moved both correlations a long way. Under an
   identical Wishart prior INLA reproduces Stan to within 1.5 to 6.5 per
   cent, and no Wishart scale gave REML-sized variances with the correlation
   left to the data. The joint prior therefore couples the variances and the
   correlations whatever the response scale. INLA's correlated iid models
   (iid2d to iid5d, and iidkd) have a Wishart prior only. Fitting either
   structure on INLA again would need a prior that separates the variances
   from the correlations, for example one built from scaled copies of
   independent fields, which flexyBayes does not build.
2. **Take the combined MET fit to realistic scale.** The emit and one live
   fit are in place (120 rows, 4 chains, clean diagnostics). What is missing
   is the same check on a trial series of realistic size, where
   `adapt_delta` and the rest of Stan's `control` list start to matter.
   Both can be set, because `flexybayes()` and `fb_brms()` pass `seed` and
   `control` to `brms::brm()`. What is missing is the run.
   *Closes the package's main gap.*
3. **Closed: INLA fits interaction and nested random effects.** `a:b` and
   `a:b:c` are fitted as an iid field over the interaction index, and agree
   with lme4's REML fit of the same model on `agridat::john.alpha` (see the
   table).
4. **A residual with many sections on INLA.** A residual per section fits on
   INLA (see the table). Past 15 variance parameters it is fitted by
   empirical Bayes (see above). *What is open*: checking fits with many more
   sections than the 12 trials compared with ASReml, and a hierarchical
   form, a prior `log sigma^2_env ~ Normal(mu, tau)`, that replaces many free
   variances with two.
5. **Exact aggregation of a factor interaction.** The per-row INLA path
   fits `y ~ a*b` on any family (see the closed items below), but the
   aggregated path cannot. It expands the fixed effects to a model matrix
   and names its columns in the INLA formula, and INLA does not resolve a
   column called `a2:b2`, even in backticks. `aggregate = TRUE` on such a
   model refuses by name and points to the per-row route. The fix is a
   column-naming scheme carried through `summary.fixed`, `marginals.fixed`
   and the latent field together, so that the posterior is never labelled
   with a made-up name. *Restores compression on the replicated factorial,
   where it is most useful.*
6. **Aggregated-binomial input.** The main fitting entry point has no
   `cbind(success, failure)` or `trials =` argument. The streaming path has
   `trials`, but the main entry point does not. The only working form is
   expansion to one Bernoulli row per trial. *Usability.*
7. **Observation weights outside the Gaussian/identity case.** `weights =`
   fits (since 0.9.3) for the Gaussian family on an identity link, on both
   engines, in the ASReml, lme4 and `glm(weights=)` precision sense. It is
   checked against an analytic reference and closely matches
   `lme4::lmer(weights=)` on a shared simulated data set (see the capability
   table above). Still open: any other family, and a non-identity link on
   the Gaussian, are refused by name (`weights_requires_gaussian`) instead of
   being given that family's own meaning of a weight (frequency,
   likelihood power, trials, exposure). `aggregate = TRUE` with weights is
   refused (`weights_not_aggregatable`), because the closed-form sufficient
   statistics do not carry a per-observation weight through the compression.
8. **A bounded uniform prior whose upper bound cuts into the posterior can
   crash the INLA program.** *Symptom*: the `inla` subprocess exits with a
   segmentation fault, INLA retries once with improved initial values, the
   retry also crashes, and `backend = "inla"` then raises `INLA fit failed:
   The inla-program exited with an error`. The R session itself is not
   affected, because the crash is inside the subprocess. The failure arrives
   as an ordinary error and nothing is lost. *Trigger*: a `fb_prior()`
   bounded uniform on the standard-deviation scale whose upper bound sits
   inside the region the data support. Reproduced on a 72-row randomised
   block model whose `sd_Block` posterior mean is 2.28 with a 97.5% bound of
   6.35, with the prior `sd(group = "Block") ~ uniform(0, 3)`.
   *Workaround*: widen the bound (`uniform(0, 10)` fits the same model
   cleanly), or use a prior with no upper bound, such as
   `sd(group = "Block") ~ half_normal(scale = 3)`, which expresses the same
   belief about scale without a hard edge. Since 0.13.0 the default prior is
   a half-normal on each standard deviation, which has no upper bound, so
   this needs a uniform bound stated with `fb_prior()`. On
   `backend = "auto"` the crash is still not a typed refusal. Dispatch
   refits the model on brms and returns a Stan fit where an INLA fit failed.
   This is not silent. The refit prints a one-time message naming the
   numerical fallback, and the returned object's class and `backend` field
   record which engine produced it, so every fit can be traced to its
   engine. *The crash is in INLA. The open part in flexyBayes is that the
   refit is reported by a message, not by a refusal the user must
   acknowledge.* Two changes in 0.12.0 bear on it. On `backend = "inla"` the
   refusal now says when the model states a uniform bound, and names the two
   remedies above. A bound tighter than the value INLA starts its search
   from is now enforced. Before 0.12.0 the search began outside the bound,
   where the prior is flat, and never entered it, so
   `sd(group = "geno") ~ uniform(0, 0.05)` on a Poisson model returned a
   standard deviation of 0.61. The search now starts inside the bound, and
   the same fit returns 0.05.
9. **The autoregressive field can lose its signal on an incomplete grid.**
   *Symptom*: the field standard deviation runs to a floor near 0.01, both
   correlations sit at about zero with credible intervals covering almost
   the whole of `[-1, 1]`, and the variance the field should have carried is
   absorbed into the nugget, which rises to match. The fit is then an
   independent-errors model reported as a spatial one. *Trigger*: unobserved
   plots. On one 12 by 10 grid with 12 of 120 responses missing under the
   default `augment` policy, ten different patterns of missing plots were
   fitted. Five lost the field and five returned ordinary correlations, with
   the same data, model and settings, differing only in which plots were
   missing. Repeating a single fit three times on one unchanged data frame
   gave the same answer on some patterns and moved between a collapsed and a
   recovered field on others, so the outcome is not always a function of
   the input alone. *What the fit reported*: before this release, nothing.
   The convergence block reported a converged mode and a passing numerical
   confirmation throughout, and the augmentation record reported that it had
   completed the design. The credible intervals were the only sign. A fit
   whose field is not identified now raises a warning naming the parameters
   that lost their identification, and `summary(fit)$varcomp` shows
   `collapsed` in the `note` cell of the `sd_spatial` row. *Workaround*:
   complete the grid, either by supplying the field-book rows for every sown
   plot with the response set to `NA` or by padding the array with
   `fb_complete_grid()`. Alternatively, fit the same model on the brms
   engine for a second reading, or give the field an informative prior with
   `fb_prior()`. *Since 0.14.5* the field SD has the half-normal default and
   no longer runs to a floor. A field with no signal reads about 0.9 beside
   a nugget of 0.8. The warning and the `collapsed` cell now come from a
   comparison with the same model without the field (NEWS, 0.14.5). The
   measurements above were taken under INLA's own field prior and have not
   been repeated.
10. **Large Gaussian data has no aggregated route in this release.**
    `aggregate = TRUE` and `flexybayes_stream()` refuse a Gaussian model
    (`aggregation_residual_prior_not_composed`). The aggregated Gaussian
    likelihood combined the within-cell sums of squares with a fixed prior on
    the residual precision, so it did not fit the residual prior that
    `prior_summary()` reported. A Gaussian model is fitted per row, which on
    a single random intercept completed at five million rows in 19.5 GB (see
    `inst/validation/benchmark_scaling.md`). *What is open*: an aggregated
    Gaussian likelihood that includes the stated or default residual prior
    exactly and restores the marginal likelihood, planned for 0.16.0 at the
    earliest.
11. **Withdrawn in 0.16.0.** This item concerned penalised selection, which is no longer part of flexyBayes.
12. **Withdrawn in 0.16.0.** This item concerned penalised selection, which is no longer part of flexyBayes.
13. **Withdrawn in 0.16.0.** This item concerned penalised selection, which is no longer part of flexyBayes.
14. **Closed: the boundary warning on a random slope no longer depends on
    the covariate's units.** The slope's standard deviation is multiplied by
    its covariate's before the comparison, in the fit-time warning and in
    `summary(fit)$varcomp`. The check is in
    `tests/testthat/test-weak-variance-flag.R`.
15. **Blocks-format `vm()` and `ped()` carriers are refused on INLA.** A
    blocks carrier is one block-diagonal covariance with one standard
    deviation. INLA fitted each block as its own field with its own standard
    deviation, so with two or more blocks the fit was a different model,
    reported as one variance-component row per block. With one block the fit
    equals the precision carrier of the same matrix. Checked on a 30-level
    pedigree against brms and against the exact posterior, its interval for
    the relationship standard deviation is 1 to 1.5 per cent narrower at each
    end, and the residual standard deviation sits 0.5 to 1.5 per cent low.
    That is INLA's approximation, and the precision and pedigree-precision
    carriers show the same. `backend = "inla"` refuses a blocks carrier with
    `inla_gate_refused` (rule `random_term_type_inla`). *Workaround*: the
    precision carrier of the block-diagonal inverse,
    `vm(geno, cov = fb_cov(Q, type = "precision"))` with
    `Q = Matrix::bdiag(lapply(blocks, solve))`, or the dense carrier of the
    block-diagonal matrix on brms. *What is open*: a blocks carrier fitted
    with one shared standard deviation, and an accuracy bound for INLA's
    variance components.
16. **A Poisson response that is almost all zeros, with a weakly identified
    group under a uniform bound, can crash INLA's program, and
    `backend = "inla"` refuses it.** On a Poisson model with a random
    intercept over 28 clusters of two and 3 non-zero counts in 56 rows,
    INLA's program stops with a segmentation fault after the mode is found,
    while it integrates over the group precision. This happens on every
    release from 0.10.0 to 0.12.0 under INLA 26.08.07 and 25.10.19, and on
    0.12.1 under INLA 26.08.07 (in the execution grid,
    `inst/validation/execution_grid/`), all under the uniform default, and
    the package refuses with `inla_program_failed`. The same design with 7
    non-zero counts fits. It is not a regression. The execution-grid cell
    that recorded 0.10.0 as fitting drew the 7-count response, and a fixture
    added to the grid later moved the random-number state so that the cell
    drew the 3-count one. The cause inside INLA is not known. Since 0.13.0
    the default half-normal prior fits the 3-count design (3 of 3 runs, with
    a 95% interval for the group standard deviation of 0.10 to 3.89, against
    0.04 to 3.27 on brms under the same prior, where the sampler reported 138
    divergent transitions in 4,000 draws). A stated `uniform(0, 3)` on the
    group standard deviation still crashes in 3 of 3 runs. *Workaround*: keep
    a half-normal prior on the group standard deviation, or fit on brms.
17. **On INLA, a hyperparameter marginal can lie far outside a stated
    uniform bound, and the fit is refused.** Take a Poisson model with two
    random intercepts, one of which splits the rows exactly as a fixed
    factor does, under the default uniform bound of 3 on each standard
    deviation. The duplicated term's mode lies on the bound. INLA estimates
    curvature by stepping across the bound, where the prior's log-density
    is -1e10. Under INLA 26.08.07, in 45 of 100 identical runs the marginal
    for the other term's precision had its mode near 5e-07, a standard
    deviation near 1,400, while INLA's optimiser mode and brms (median 0.53,
    95% interval 0.17 to 2.11) were inside the bound. Under INLA 25.10.19,
    100 of 100 runs stayed inside. 0.12.0 returned such a fit with a
    warning. From 0.12.1 every INLA fit is checked against INLA's own record
    of its priors, a fit whose mode or reported interval lies outside a
    bound is refused with `prior_support_violated`, and the check is kept in
    `fit$extras$support_audit`. *Workaround*: drop a random term that
    duplicates a fixed factor, state
    `sd(group = ...) ~ half_normal(scale = 1.53)` (no runaway fit in 24 runs
    of the same design, posterior median 0.505 against 0.489 on brms), or
    fit on brms.

    The same check refuses a fit whose mode is inside the bound when INLA's
    reported interval reaches past it, which happens on a weakly identified
    component with few levels. On `agridat::john.alpha` with
    `random = ~ rep + rep:block`, the default bound on `rep` (three levels)
    is 3.05 and INLA's interval runs from 0.17 to 3.65, so the default INLA
    fit is refused in 8 of 8 repeats. *Workaround*: state a half-normal prior
    on the standard deviations, for example with scale the standard
    deviation of the response, or fit on brms. Since 0.13.0 the default
    prior is `half_normal(scale = s)` on each standard deviation, on both
    engines, with `s` equal to the former bound. It has no bound for an
    interval to cross, and the default `john.alpha` fit returns. *What is
    open*: a user-stated uniform bound is still refused when INLA's reported
    interval reaches past it.
18. **Withdrawn in 0.16.0.** This item concerned penalised selection, which is no longer part of flexyBayes.
19. **Withdrawn in 0.16.0.** This item concerned penalised selection, which is no longer part of flexyBayes.
20. **Closed on INLA: a weakly identified variance component is flagged.**
    Under the default half-normal prior on a standard deviation, a component
    the data barely inform is reported at a value driven mostly by the prior.
    A component whose 95% interval covers 40% or more of its prior's
    probability is marked `weakly identified` in `summary(fit)$varcomp`,
    and the fit warns. A random term with the same grouping as a fixed
    factor is flagged this way. From 0.18.0 a second check is added: a 95%
    interval on the standard deviation 4.24 or more times its posterior
    median also flags the component. On the 80 calibration fits it adds one
    flag, a five-level group with no true variance. Every other flagged fit
    was already flagged by the first check. It does not flag the
    `besag.met` counties of item 24. *What is open*: brms fits are not
    checked, because that needs the prior each component received to be
    read back from Stan.
21. **The factor-analytic covariance is refused: on simulated trials its
    genetic-correlation intervals covered the true values less often than
    required.**
    `fa(env, k):gen` was fitted on INLA from 0.10.0.9000. From 0.13.0 both
    engines refuse it (`factor_analytic_not_repeatable`), and flexyBayes has
    no factor-analytic lowering for brms. *Why*: in a test whose thresholds
    were set before it ran, INLA fitted 10 simulated trials per k (60
    genotypes, 6 environments, 2 replicates) with a known covariance. The 95%
    intervals of the environment genetic variances covered the true value 58
    of 60 times for k = 1 and 55 of 60 for k = 2, against 53 required. The
    intervals of the genetic correlations covered it 134 and 125 of 150
    times, against 136. The posterior medians also missed the point
    tolerances, 15% on a variance and 0.10 on a correlation, but no
    estimator can meet those at 60 genotypes: the genetic effects drawn
    already differ from the population covariance by about 13% on a
    variance. For k = 2 the genetic correlations were overstated, most where
    the true loadings are negative. The loading prior, normal(1, 0.316),
    gives a negative loading a probability of 0.00078. On
    `agridat::omer.sorghum` the common part of the genetic variance at
    INLA's posterior mode differed from ASReml 4.2's by a relative 1.53
    (153%, k = 1, environment E4) and 1.65 (165%, k = 2, environment E5),
    where ASReml estimated one or two specific variances as zero. The loadings
    are anchored at one for factor *j* at an anchor level, which represents
    the family only where that level carries a common-factor loading. From
    0.17.0 the anchors are the levels with the largest common genetic
    variance in the data. Up to 0.16.0 they were the first levels in factor
    order.
    *Workaround*: `diag(env):gen`, one variance per environment with no
    covariance between them, fits on both engines. `us(env):gen` fits on brms
    and estimates the correlations directly.
    *What is open*: a loading prior under which genetic correlations can be
    negative. The emit is unchanged and can still be reached by a direct
    `emit_inla()` call, as the `us()` branch can.
    A rebuilt fit that estimates the loadings directly, as ASReml does, was
    tested in October 2026 against thresholds fixed before it ran, and failed
    for every k. In simulation-based calibration its genetic-correlation
    intervals were too narrow: each 5% tail held 6 to 11% of the true values.
    The Stan fits used as the reference did not converge on 36 of 60
    simulated data sets at k = 2. The rebuilt fit is not in this release.
22. **Separate variances with one shared correlation, `corh(f):g`, are refused
    on both engines.** The structure sits between `diag(f):g`, with k
    variances, and `us(f):g`, with k(k+1)/2 parameters: k variances and one
    correlation. flexyBayes has no lowering for it. INLA has an exchangeable
    group model with one variance and one shared correlation, which
    flexyBayes does not use. *Workaround*: with one variance for every level
    and a correlation of zero or more, `random = ~ g + f:g` is the same
    covariance and fits on both engines. Otherwise use `diag(f):g`, or
    `us(f):g` on brms. *What is open*: a lowering, checked against a
    reference fit. `cor(f)`, one correlation with no variances of its own, is
    refused for the same reason.
23. **A multi-trait covariance, `us(trait):vm(gen)`, is refused on both
    engines.** flexyBayes has no lowering for a covariance between traits
    crossed with a known relationship matrix. Both engines can represent it.
    brms writes it as `(0 + trait | gr(gen, cov = K))` on data with one row
    per trait, or with `mvbind()` on data with one column per trait, and INLA
    can build it from a known-covariance field for each trait joined by a
    scaled `copy` (`tests/testthat/test-engine-claims.R`). *Workaround*: fit
    the model with brms directly, or fit each trait separately. *What is
    open*: a lowering on either engine, which is a new model structure and
    waits for a release that admits one.
24. **The weak-identification flag can stay silent on a variance REML cannot
    identify.** On `agridat::besag.met` with `random = ~ diag(env):gen`,
    ASReml estimates the county `C2` genotype variance at 3.28 with a
    standard error of 23.6, and `C5` at 3.77 with 23.7. Under the default
    best of three INLA fits the posterior median standard deviations are
    4.42 and 4.42, variances of 19.5 and 19.6, with 95% intervals of 1.73
    to 9.18 and 1.68 to 9.42 on the standard deviation. Each interval is
    1.7 times its median and covers 3% of the prior's probability. Neither
    the prior-share check nor the interval-to-median check of item 20 flags
    them, so C2 and C5 are not marked `weakly identified`. Read the REML
    standard error beside the posterior for such a component. *What is
    open*: a flag that catches this case. brms fits are not checked (item
    20). The regression test is in `tests/testthat/test-asreml-parity.R`
    (DEF-201), skipped while this is open.
25. **Closed: `ar1(t):g` names its correlation `rho_t`.** For
    `random = ~ ar1(t):g`, in either factor order, the correlation along `t`
    is `rho_<factor>` in `summary(fit)$varcomp`, in
    `summary(fit)$spatial_field` and in the draws. A fit saved before this
    release still prints `rho_row` in its stored table until it is refitted,
    and asking `vpredict()` for `rho_row` on a new fit stops and names the
    `rho_<factor>` to use. The regression test is in
    `tests/testthat/test-ar1-varcomp-names.R` (DEF-225).
26. **A refit in a fresh R process can reach a different INLA mode.** INLA
    reports a converged mode in every case, so `summary()` does not show it.
    Measured with INLA 26.8.7 on one thread, with the same call and the same
    data in each process:
    - On a 48-trial barley series (10,890 plots), four processes reached two
      modes, 0.10 apart in log marginal likelihood. The variance component
      standard deviations differed by at most 3.4% between them.
    - On six trials of the same series, eight processes reached four modes,
      up to 0.45 apart.
    - Over 32 processes on those six trials, some with INLA's optimiser
      settings changed (reordering, search directions, initial values, the
      start-up benchmark below), the main-plot standard deviation differed by
      up to 26%, the locality standard deviation (four localities) by 14%,
      the others by 7% or less and the residual by 0.4%.

    flexyBayes passes INLA the same arguments in every run. INLA's log was
    identical across runs up to its start-up benchmark of linear-algebra
    routines, and the differences remained with that benchmark switched
    off. Setting `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
    `VECLIB_MAXIMUM_THREADS` and `MKL_NUM_THREADS` to 1 did not remove them
    either: six fresh processes on the six trials reached three modes, as
    six at 8 threads did. INLA's `restart = 3` made all four processes on the
    full series agree, but on the lower mode, so it does not fix this.

    Since 0.17.0 each INLA fit is run three times, once in the session and
    twice in fresh R processes, and the fit with the highest log marginal
    likelihood is returned (`options(flexyBayes.inla_best_of)`, see
    `?flexybayes`, section *INLA threads*). The spread of the three fits is
    kept in `$extras$reproducibility$best_of`, and `print()` names any
    standard deviation that differs by more than 5% across them. On the six
    trials, four processes each running three fits returned three distinct
    modes, at log marginal likelihoods -1459.696 (twice), -1459.523 and
    -1459.520. Within a process the three marginal likelihoods differed by
    up to 0.29, the main-plot standard deviation by up to 13.6%, and the
    others by 3.9% or less. Best of three narrows the spread but does not
    make fresh processes agree. The three fits run at once when free memory
    holds them (`options(flexyBayes.inla_best_of_workers)`), and on the six
    trials best of three took a median 29 s against 64 s one after another.
    *Workaround*: where a conclusion rests on a weakly informed variance,
    raise `flexyBayes.inla_best_of` or refit in several fresh processes (for
    example with `fb_isolated_lapply()`) and compare. *What is open*: no
    fitting route is yet known that reliably reaches the highest mode. Best
    of 3 never reached the highest mode seen on these data (-1459.362).

## Closed since the previous revision

These were open in the previous revision of this file.

- **Gaussian fits on INLA at large response scale.** Every Gaussian model on
  INLA was wrong once the response standard deviation reached the low
  hundreds. On `lme4::sleepstudy` the plain random intercept returned ratios
  of the subject variance to the residual between 0.03 and 6.99 against
  REML's 1.198, and several of those fits passed the numerical-confirmation
  gate. Every Gaussian INLA fit now runs on a standardised response, and the
  result is mapped back to the original units exactly. The default priors
  the package chooses are scaled by the response and by each covariate. On
  `lme4::sleepstudy` the default-prior fit agrees across response
  multipliers from 0.01 to 1000. The refusal of a Gaussian INLA fit whose
  response standard deviation exceeds 150
  (`inla_response_scale_unvalidated`) is gone. An explicit
  `backend = "inla"`, `fb_inla()`, `backend = "auto"`, the aggregated route
  and the streamed route each fit a response on a thousandfold scale and
  agree with the same model fitted on the original one. The checks are in
  `tests/testthat/test-response-standardisation.R`,
  `tests/testthat/test-inla-scale-aware-defaults.R` and
  `tests/testthat/test-inla-containment.R`.
- **INLA interaction fixed effect on the binomial path.** `y ~ a*b` with
  `family = "binomial"` and `backend = "inla"` was recorded as failing
  inside INLA with `object 'a_b' not found`. The per-row path fits it. On a
  replicated 4-by-4 factorial the sixteen posterior means agree with a
  `glm()` maximum-likelihood fit to within 0.11 on the log-odds scale, the
  largest gap being on the sparsest cell.
- **Unrecognised ASReml functions.** The parser used to read an unknown call
  as a variable of that name, and `ar2`, `cor` and `str` reached an untyped
  `stop()`. The grammar is now a closed set. `foo(g)` refuses as
  `asreml_function_not_recognised`, and `ar2`, `str` and the `cor` family
  each refuse under their own code.
- **Per-trial separable AR1 field.** The single-lattice AR1 field
  (`random = ~ ar1(row):ar1(col)`, one trial, one grid) had no
  multi-environment form. `random = ~ id(trial):ar1(row):ar1(col)`
  now fits on INLA through the `replicate =` mechanism, with one field
  realisation per level of `trial`. The field's row correlation, column
  correlation and SD are **shared across every level**, not estimated per
  trial. ASReml reads `at(trial):ar1(row):ar1(col)` as per-trial parameters,
  and since 0.16.0 that form fits on INLA with its own SD and both
  correlations on each trial. `at(trial, level):ar1(row):ar1(col)` fits the
  field on the named trials only. brms has no lowering for the field.
- **Units followed through time.** `random = ~ ar1(t):g`, in either
  factor order, fits one AR1 series per level of `g` on INLA,
  independent across levels, with the correlation and field SD shared
  by every level: the Verbyla-Cullis repeated-measures covariance on
  the random side. The model line, the field block and the note naming
  the ASReml residual the field is not used to describe this field as
  `ar1(t):ar1(g)` with four parameters. They now describe it correctly
  (AR1 along `t`, independent across `g`, three hyperparameters).
  brms refuses.

## Response families: where the boundary is, and whose it is

The families flexyBayes accepts (`.resolve_family()`) are fewer than either
engine supports. A family is admitted when the package has an emit for it
and a test exercises it. Four families users have asked for, each checked
against the installed engines by `tests/testthat/test-engine-claims.R`:

| Family | brms 2.23.0 | INLA 26.08.07 | flexyBayes | Where the boundary is |
|---|:-:|:-:|:-:|---|
| `hurdle_gamma` | native (`dpars` mu, shape, hu) | absent | **fits on brms since 0.9.2** | None |
| `zero_inflated_gamma` | absent (`brmsfamily()` refuses) | absent | refuses | Neither engine has it. Use `hurdle_gamma`. For a gamma positive part, which has no mass at zero, the hurdle and the zero-inflated model are the same model. |
| `tweedie` | absent | **present** in the likelihood list | refuses | **flexyBayes**, not the engine. No INLA emit (no link or power parameter, no validated fit). A feature request, not a defect (see below). |
| `compound_poisson` | absent | present as `tweedie` | refuses | Same as `tweedie`. The two name the same distribution. |

Each of these refusals carries the boundary note in its message.

**Tweedie on INLA is an open feature request.** INLA has the likelihood, and
flexyBayes has no emit for it. The work is a family row, the `p` power
parameter passed through `control.family`, and a fit validated against an
independent implementation (`statmod`/`cplm` or `mgcv::Tweedie`) before the
row can read `fits`. Until then, for compound-Poisson gamma data, fit the
zero and positive parts separately and recombine on the posterior, or use
`family = "hurdle_gamma"` with `backend = "brms"`, which is the same
two-part model with the parts fitted jointly.

## Minor and environment-specific notes

- **INLA hardware-probe messages.** On some operating systems and INLA
  builds, the INLA program prints harmless hardware-probe lines (for example
  `/bin/kstat: No such file or directory`) to the console during a fit. They
  come from the INLA subprocess, not flexyBayes, below the level R can
  capture. They are cosmetic and do not affect results, and they do not
  appear on every platform.

- **Missing host core count.** Some sandboxed macOS sessions return `NA` from
  `parallel::detectCores()`, from which INLA constructs `NA:1` and then prints
  a parser fallback. From 0.13.0, flexyBayes passes `1:1` when the caller gave
  no thread setting and the host count is unavailable. A valid INLA default or
  an explicit value from the caller is not replaced. This corrects the
  resource record. It does not suppress the INLA hardware-probe messages
  above.

## How to help

- Pick an issue above. Items 1 and 5 are the most valuable and the least
  risky. The parser, gate, emitters and refusal registry are already in
  place, so these need new emit branches, reason codes and verification
  runs, not new subsystems.
- A fix is done only when `summary()`, `prior_summary()`,
  `canonical_names()` and prediction all handle the new term, and when an
  assertion in `tests/testthat/test-capability-matrix.R` backs the
  capability row that claims it.
- Validate against ASReml or `lme4` on a `barrero.maize` subset before
  claiming a class is supported, and add the result to `inst/validation/`.

Issues are tracked at <https://github.com/AAGI-AUS/flexyBayes/issues>.
