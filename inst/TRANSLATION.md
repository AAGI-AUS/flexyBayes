# What flexyBayes promises when it translates an ASReml model

flexyBayes accepts a documented subset of ASReml-shaped formula syntax and
lowers it to INLA or brms. The lowerings are not all of one kind, and this
document says which kind each one is.

The distinction matters because "ASReml syntax" as a blanket description
covers three quite different things: structures where the model fitted is the
model ASReml would fit, structures where the covariance family is the same but
the identification and the priors are not, and structures that are related to
ASReml's but built from a different basis. A user choosing between this
package and ASReml needs to know which of the three they are getting.

## What no row claims

ASReml returns AI-REML estimates. flexyBayes returns posterior summaries from
INLA or brms. Even where the likelihood and the covariance family are
identical, a posterior mean need not equal a REML estimate, and no row in the
table below claims that it does. What the `exact` rows claim is that the
design matrices and the covariance algebra are the same.

Cross-engine agreement is also not evidence for these rows. Both engines
consume the same parsed representation, so a parser error that misreads a
structure produces the same wrong model twice, and the two agree perfectly.
Agreement between INLA and brms tests the emitters against each other; it
cannot test either against ASReml.

## Translation classes

| Class | What it promises |
|---|---|
| **exact** | Same response likelihood, fixed-effect design, random-effect design, covariance family, residual covariance, and missing-row estimand. Estimation method and prior differ and are stated separately. |
| **covariance-family** | The fitted covariance set is the same after a documented parameter map, while identification, priors, or the inferential target differ. |
| **conditional** | Covariance-family equivalent only where the stated precondition holds. The precondition is named in the table and is not currently enforced by the code. |
| **analogue** | Scientifically related, with a different basis, penalty, nugget, rank, or prediction rule. Not the same model. |
| **refused** | Cannot be represented faithfully. The request terminates with a typed condition naming an alternative where one exists. |

## Syntax vintages

ASReml 4 removed two spellings that this package still accepts. Residual-side
`at()` was replaced by `dsum()`, and `ped()` was replaced by `vm()`. Both are
kept here so that models written for ASReml 3 continue to run, and both are
labelled as the ASReml 3 forms they are. Current ASReml may reject them.

The canonical spellings are `dsum(~ units | env)` for a sectioned residual and
`vm(g, A)` for a known relationship matrix.

## The registry

The table below is generated from `.fb_translation_registry()` in
`R/translation_registry.R`. A test fails when this file drifts from a fresh
render, and a machine-readable copy ships at
`inst/translation/translation_registry.csv`.

The `Independent oracle` column names what grounds each row and the version it
was taken at. Several rows read `none at this anchor`. That is accurate rather
than provisional: the 0.10 development line added or materially changed the
INLA factor-analytic, `us()` through `iidKd`, correlated-slope, sectioned
residual and combined multi-environment paths, and the banked ASReml
comparisons were taken at flexyBayes 0.9.0, before those changes. Until a
paired ASReml run at the current anchor exists for a row, that row's
translation class is a statement about the code rather than a measured result.

The `Prior cardinality` column says how many engine slots a prior on the
structure reaches: one, several (with the expansion named), none because the
engine keeps its own prior, or none because the structure is refused. Each
fit records the same count for its own priors, and `prior_summary()` lists
the slots.

<!-- translation-registry:begin -->
| Requested structure | flexyBayes spelling | Syntax vintage | Engine | Translation class | Prior cardinality | Condition | Independent oracle |
|---|---|---|---|---|---|---|---|
| Fixed factors, covariates, means model | `fixed = ~ trt + x` | ASReml 4 | inla, brms | **exact** | one-to-one: a b() prior reaches its own coefficient | Ordinary R contrast rules; rank-deficiency policy is R's. | lm() / lme4 unit tests at this anchor |
| Simple iid random intercept | `~ g`, `id(g)`, `ide(g)`, `(1 \| g)` | ASReml 4 | inla, brms | **exact** | one-to-one | -- | ASReml 4.2.0.392 at flexyBayes 0.9.0; lme4 REML at this anchor |
| Nested and interacted iid terms | `~ rep:block`, `~ gen:env` | ASReml 4 | inla, brms | **exact** | one-to-one: keyed by the colon-joined group | Each component term remains iid. | ASReml 4.2.0.392 at flexyBayes 0.9.0 (RCBD / IBD / MET) |
| Diagonal heterogeneous variance over a factor | `diag(f):g`, `idh(f):g`, `at(f):g` | ASReml 4 | inla, brms | **covariance-family** | one-to-many: one field per level of the outer factor | Repeated identical fits on a low-signal trial do not agree: three runs on agridat::besag.met span a factor of 2.3 in a component. Mode adoption (R/inla_restart_probe.R) removed the collapse-to-zero failure this row used to show, and did not remove the remaining spread. ASReml is stable on the same rows. | ASReml 4.2.0.392 at this anchor on agridat::besag.met. Structural layer exact: 7 variance parameters, 6 estimable fixed effects, 1,152 observations, identical on both sides. Numerical layer FAILS its declared tolerance on the intercept -- 148.60 against 149.14, a difference of 0.55 against a tolerance of 0.35. ASReml pins two of the six genotype variances at exactly zero and the posterior does not, so the two systems weight the fixed effects differently. The discrepancy is recorded rather than excused, and the tolerance is not moved |
| Sectioned residual, one variance per section | `residual = ~ dsum(~ units \| env)` | ASReml 4 | inla, brms | **exact** | one-to-many: one likelihood per section | Inner term is plain `units`; a structured inner is refused. | ASReml 4.2.0.392 at this anchor on agridat::yan.winterwheat. Structural layer exact: 10 variance parameters, 9 estimable fixed effects, 162 observations, identical on both sides. Genotype variance ratio 0.0843 against 0.0854, a difference of 0.0012 against a declared tolerance of 0.05; intercept identical to four decimal places; BLUP agreement r = 1.0000. Repeatability screen stable. This row is validated at the current anchor |
| Sectioned residual, ASReml 3 spelling | `residual = ~ at(env):units` | ASReml 3 alias | inla, brms | **exact** | one-to-many: one likelihood per section | ASReml 4 removed residual-side at() in favour of dsum(), so current ASReml may reject this spelling. Accepted here for compatibility; dsum() is the canonical form. | Shares the dsum() route and its evidence |
| Per-trial genotype and residual variances | `random = ~ diag(trial):g`, `residual = ~ dsum(~ units \| trial)` | ASReml 4 | inla, brms | **covariance-family** | one-to-many: one field and one likelihood per trial | Past 15 variance parameters INLA fits the model by empirical Bayes (int.strategy = "eb"): the fixed and random effects are conditional on the posterior mode of the variance parameters, and each variance keeps its own posterior. Other models past 15 refuse with hyperparam_budget. | ASReml 4.2.0.482 at this anchor on 12 simulated trials, 26 variance parameters (data-raw/model_c_reference.R): every per-trial variance ASReml identifies within 0.17 of its ASReml standard error on INLA and 0.28 on brms |
| Unstructured covariance over a factor, on brms | `us(f):g` | ASReml 4 | brms | **covariance-family** | one-to-one: one standard-deviation row for the group | -- | Live emission and parameter-count tests; no current ASReml comparison |
| Unstructured covariance over a factor, on INLA | `us(f):g` | ASReml 4 | inla | **refused** | refused | Estimated correlation matrices are fitted on brms only. The INLA lowering, an iidKd field under a joint Wishart prior, refuses with correlated_covariance_requires_brms and backend = "auto" routes to brms. The row was covariance-family until 2026-09-11, analogue until 2026-09-13, and refused from then. | Refusal is the claim; asserted by name. The refusal rests on ASReml 4.2.0.392 on agridat::besag.met, oracle scenarios 06 and 07 (finding O04): a genetic correlation of 0.142 against 0.930 over two environments and 0.198 / 0.231 / 0.320 against 0.898 / 0.915 / 0.955 over three, with eight identical calls on the response rescaled by ten none converged and correlations from -0.894 to +0.686. Under every identical Wishart prior tried INLA and Stan agreed, and no Wishart scale gave REML-sized variances (INLA scale diagnosis, 2026-09-13), so the joint prior rather than the integration drives the disagreement |
| Factor-analytic covariance | `fa(f, k):g` | ASReml 4 | inla | **refused** | one-to-many: one field per common factor and one per level | Refused from 0.13.0 with factor_analytic_not_repeatable, on both engines, and flexyBayes has no factor-analytic lowering for brms. The loadings are anchored to one for factor j at an anchor level, which represents the factor-analytic family only where that level carries a common-factor loading. From 0.17.0 the anchors are the levels with the largest common genetic variance in the data; until 0.16.0 they were the first levels in factor order. The row was covariance-family until 2026-09-11, conditional from then, and refused from 2026-09-16. | Refusal is the claim; asserted by name. The refusal rests on INLA fits of 20 simulated trials with a known covariance, and on agridat::omer.sorghum against ASReml 4.2; the numbers are under Reading a row |
| Known covariance or precision carrier | `vm(g, K)`, `vm(g, cov = fb_cov(...))` | ASReml 4 | inla, brms | **exact** | one-to-one: one field; a blocks carrier is refused on both engines | Representation metadata (dense / chol / precision / blocks) must be correct. | Dense and sparse carrier tests; no current ASReml fixture for the inverse form |
| Pedigree relationship carrier | `ped(g, A)` | ASReml 3 alias | inla, brms | **exact** | one-to-one | ASReml 4 removed ped() in favour of vm(). Accepted here for compatibility; vm(g, A) is the current spelling. | Shares the vm() carrier machinery and its evidence |
| Correlated random slope, on brms | `(x \| g)` | lme4 / brms | brms | **covariance-family** | one-to-one: the intercept and slope standard deviations each take a row | brms writes the term verbatim and fits it natively. | Live emission and recovery tests; no current ASReml comparison on this engine |
| Correlated random slope, on INLA | `(x \| g)` | lme4 / brms | inla | **refused** | refused | Estimated correlation matrices are fitted on brms only. The INLA lowering, one iid2d field over a two-block index under a joint Wishart prior, refuses with correlated_covariance_requires_brms and backend = "auto" routes to brms. The row was covariance-family until 2026-09-13 and is refused from then. | Refusal is the claim; asserted by name. The refusal rests on ASReml 4.2.0.392 and lme4 on lme4::sleepstudy, oracle row 19 (finding S01, defect D1): a slope variance of 46.37 at the posterior mode against 35.076 and 35.072, and a tenfold rescaling of the response that moved the intercept-slope correlation from 0.02 to 0.969 while REML held at 0.066. Under an identical Wishart prior INLA matched Stan within 1.5 to 6.5 per cent, and a simulated response on a kilogram scale returned a residual standard deviation of 0.13 against a truth of 250 (INLA scale diagnosis, 2026-09-13) |
| Uncorrelated random slope | `(x \|\| g)` | lme4 / brms | inla, brms | **covariance-family** | one-to-one: the intercept and slope standard deviations each take a field | -- | lme4 REML at this anchor |
| Temporal AR1 on the random side | `~ ar1(t)`, grouped AR1 | ASReml 4 | inla | **exact** | engine default: the field's precision and correlation keep INLA's priors | Exact for the latent-field model only: an AR1 field plus an observation residual. This is not a translation of ASReml's residual-side AR1, which carries no nugget. | Independent dense-V REML and simulation at this anchor |
| Separable residual AR1 x AR1 | `residual = ~ ar1(row):ar1(col)` | ASReml 4 | inla | **covariance-family** | engine default: the field's precision and correlations keep INLA's priors | Gaussian response with the identity link, one plot per grid cell. The observation error is fixed at a variance of exp(-12) on the standardised response; the posterior medians move by less than 1e-4 between exp(-10) and exp(-16). brms refuses. | ASReml 4.2.0.482 at this anchor on a simulated 16 x 10 trial (data-raw/at_ar1_reference.R): every posterior median within 0.3 posterior SDs of the REML estimate |
| Separable residual AR1 x AR1 per section | `residual = ~ dsum(~ ar1(row):ar1(col) \| trial)` | ASReml 4 | inla | **covariance-family** | one-to-many: one field per section, each with its own hyperparameters | As for the single residual field, on each section's own grid. brms refuses. | ASReml 4.2.0.482 at this anchor on three simulated trials (data-raw/at_ar1_reference.R): every posterior median within 0.3 posterior SDs of the REML estimate |
| Random term on selected levels of a factor | `at(f, level):g` | ASReml 4 | inla, brms | **exact** | one-to-many: one field per named level | g is a plain factor. Levels are named by value or by position in levels(f); a level the data lack refuses. | ASReml 4.2.0.482 at this anchor on six simulated trials, three with main plots (data-raw/at_ar1_reference.R), on both engines: every posterior median within 0.3 posterior SDs of the REML estimate |
| Separable AR1 field per level of a factor | `at(trial):ar1(row):ar1(col)`, `at(trial, level):ar1(row):ar1(col)` | ASReml 4 | inla | **exact** | one-to-many: one field per level, each with its own hyperparameters | An AR1 field with its own SD and correlations on each named level, plus the observation residual; each level's plots fill its grid. brms refuses. | ASReml 4.2.0.482 at this anchor on three simulated trials (data-raw/at_ar1_reference.R): every posterior median within 0.3 posterior SDs of the REML estimate |
| Penalised spline | `spl(x)`, `s(x)` | flexyBayes extension | inla | **analogue** | one-to-one: smooth() on the random walk's precision | Lowered to INLA's second-order random walk. ASReml's spline carries its own knot count, duplicate-point handling, prediction points and spline.scale, and no comparison of basis, null space, penalty, effective degrees of freedom or prediction has been made. Related, not the same model. | Emission and fit tests only |
| Gaussian relative weights | `weights = w` | ASReml 4 | inla, brms | **exact** | one-to-one: the residual prior is unchanged | Relative weights only, giving variance proportional to sigma^2 / w. ASReml's absolute-weight / fixed-dispersion mode is a different estimand and is not carried. | lme4 at this anchor |
| Missing response retained as a prediction target | `na_action = "augment"` | flexyBayes extension | inla, brms | **covariance-family** | one-to-one: the priors are unchanged | A representation-preserving policy rather than a universal ASReml equivalence: the design row is kept as a latent target. The estimator assumptions are the user's to check. | Large missing-plot simulation study and independent dense-V REML |
| Prediction over a random classify factor | `predict(fit, classify = "gen")` with `gen` random | ASReml 4 | inla, brms | **exact** | one-to-one: the priors are unchanged | ASReml's default prediction rule: the fixed part averaged over the other fixed factors, plus the random effects of the classify terms; other random terms are ignored. present, average, ignore, sed and vcov refuse with classify_argument_not_supported. | ASReml 4.2.0.482 at this anchor on agridat::john.alpha and 16 genotypes of agridat::besag.met (data-raw/asreml_classify_random.R): within 0.11 SE for genotypes and 0.4 SE for genotype-by-county cells, and within 0.006 SE with the variances held at ASReml's |
| Prediction at a level the fit did not see | `predict(fit, newdata, allow_new_levels = TRUE, mean_by = ~ trial)` | flexyBayes extension | inla, brms | **covariance-family** | one-to-one: the priors are unchanged | A posterior predictive, with new levels drawn from their fitted variances and the observation noise included. Random slopes, splines, known covariances, factor-analytic and unstructured terms, observation weights, a known variance and non-Gaussian responses refuse by name. | ASReml 4.2.0.482 at this anchor (data-raw/asreml_predict_new_level.R): at ASReml's REML variances the variance of a trial, cell or locality mean agrees with ASReml's plug-in value to 2e-15 |
| Heterogeneous variance with one common correlation | `corh(f):g` | ASReml 4 | inla, brms | **refused** | refused | flexyBayes has no lowering for it on either engine; it is not the diagonal or the unstructured form. INLA has an exchangeable group model, with one variance and one shared correlation, which flexyBayes does not use. | Refusal is the claim; asserted by name |
| Tensor-product and multivariate smooths | `te(x, z)`, `ti()`, `t2()`, `s(x, z)` | lme4 / brms | inla, brms | **refused** | refused | Univariate smooths only. | Refusal is the claim; asserted by name |
| Structured inner term inside a sectioned residual | `dsum(~ ar1(a) \| b)`, any inner term but `units` and `ar1(a):ar1(b)` | ASReml 4 | inla, brms | **refused** | refused | Reducing it to a per-section variance would drop the structure written. | Refusal is the claim; asserted by name |
<!-- translation-registry:end -->

## Repeatability comes before comparison

A row's oracle column can only mean something if the fit it describes is
reproducible. On a low-signal trial several structures here are not: repeated
identical calls return answers that differ, in the worst case by more than an
order of magnitude.

Two consequences are visible in the table. Where a row reads
`ASReml 4.2.0.392 at this anchor`, the fit was screened for repeatability
first and compared second. Where a row says it was *not* compared at this
anchor, that is because a single paired comparison would have measured one
draw from a multi-modal optimisation rather than the translation.

flexyBayes now detects the commonest form of this. When a variance component
comes back pinned at zero, the fit is restarted from a mode where that
component carries real variance, and if the restart is materially better it
is the fit returned -- stated in the convergence block, not applied silently.
That closed the failure entirely for simple and nested random effects. It did
not close it for the factor-analytic and combined multi-environment rows,
which carry twelve and thirteen variance parameters.

A user fitting the combined multi-environment structure to a trial with weak
genetic signal should refit and compare before acting on a variance
component.

## Reading a row

*Factor-analytic covariance* is the row worth reading closely. The lowering
builds iid latent factors with INLA copy coefficients for the loadings and
separate specific variances, and it fixes the loading for factor *j* at level
*j* to one. Fixing a loading is a scale and sign reparameterisation only when
the true loading at that level is not zero. A covariance in which the anchor
level has no common-factor loading, while other levels share the factor, has
no finite representation in this form -- and the fit still runs. Until 0.16.0
the anchor was taken from factor-level order, so which covariances were
representable depended on how the levels happened to be ordered. From 0.17.0
the anchors are chosen from the data, as the levels with the largest common
genetic variance.

That was why the row was `conditional` rather than `covariance-family`, and
why the summary stopped printing `Representation: exact` for a factor-analytic
term. From 0.13.0 the row is `refused`. The emit remains in the source,
reachable by a direct `emit_inla()` call, as the `us()` branch is; what is
withdrawn is the claim attached to its output.

The refusal rests on a simulation test whose thresholds were set before it
ran. For each of k = 1 and k = 2, INLA fitted 10 simulated trials of 60
genotypes in 6 environments with 2 replicates. The posterior median and 95%
interval of every environment's genetic variance (G_ee) and every genetic
correlation were compared with the known covariance. All 20 fits converged.

| Check | k = 1 | k = 2 | Required |
|---|---|---|---|
| G_ee intervals covering the true value | 58 of 60 | 55 of 60 | at least 53 |
| Genetic-correlation intervals covering the true value | 134 of 150 | 125 of 150 | at least 136 |
| G_ee medians more than 15% from the true value | 30 of 60 | 29 of 60 | none |
| Correlation medians more than 0.10 from the true value | 26 of 150 | 61 of 150 | none |

No estimator could meet the two checks on posterior medians at 60
genotypes: the genetic effects drawn for 60 genotypes already differ from
the population covariance by about 13% on a G_ee, before any residual noise.
The coverage rows are the informative ones.

For k = 2 the fits overstated the genetic correlations. The loading prior,
normal(1, 0.316) on every free loading, gives a negative loading a
probability of 0.00078, so genetic correlations are positive under it. The
simulated k = 2 covariance has second-factor loadings of -0.5 and -0.4. The
correlation between E4 and E6, whose true value is 0.025, was overestimated
by 0.26 on average.

On `agridat::omer.sorghum` (18 genotypes, 6 environments, 4 replicates) the
comparison is with ASReml 4.2. At INLA's posterior mode, the common part of
the genetic variance, the diagonal of LL', differed from ASReml's by a
relative 1.53 (153%) for k = 1 in environment E4 and 1.65 (165%) for k = 2
in environment E5. ASReml estimated the specific variance as zero in one
environment for k = 1 and in two for k = 2.
