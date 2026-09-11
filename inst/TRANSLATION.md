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

<!-- translation-registry:begin -->
| Requested structure | flexyBayes spelling | Syntax vintage | Engine | Translation class | Condition | Independent oracle |
|---|---|---|---|---|---|---|
| Fixed factors, covariates, means model | `fixed = ~ trt + x` | ASReml 4 | inla, brms | **exact** | Ordinary R contrast rules; rank-deficiency policy is R's. | lm() / lme4 unit tests at this anchor |
| Simple iid random intercept | `~ g`, `id(g)`, `ide(g)`, `(1 \| g)` | ASReml 4 | inla, brms | **exact** | -- | ASReml 4.2.0.392 at flexyBayes 0.9.0; lme4 REML at this anchor |
| Nested and interacted iid terms | `~ rep:block`, `~ gen:env` | ASReml 4 | inla, brms | **exact** | Each component term remains iid. | ASReml 4.2.0.392 at flexyBayes 0.9.0 (RCBD / IBD / MET) |
| Diagonal heterogeneous variance over a factor | `diag(f):g`, `idh(f):g`, `at(f):g` | ASReml 4 | inla, brms | **covariance-family** | Repeated identical fits on a low-signal trial do not agree: three runs on agridat::besag.met span a factor of 2.3 in a component. Mode adoption (R/inla_restart_probe.R) removed the collapse-to-zero failure this row used to show, and did not remove the remaining spread. ASReml is stable on the same rows. | ASReml 4.2.0.392 at this anchor on agridat::besag.met. Structural layer exact: 7 variance parameters, 6 estimable fixed effects, 1,152 observations, identical on both sides. Numerical layer FAILS its declared tolerance on the intercept -- 148.60 against 149.14, a difference of 0.55 against a tolerance of 0.35. ASReml pins two of the six genotype variances at exactly zero and the posterior does not, so the two systems weight the fixed effects differently. The discrepancy is recorded rather than excused, and the tolerance is not moved |
| Sectioned residual, one variance per section | `residual = ~ dsum(~ units \| env)` | ASReml 4 | inla, brms | **exact** | Inner term is plain `units`; a structured inner is refused. | ASReml 4.2.0.392 at this anchor on agridat::yan.winterwheat. Structural layer exact: 10 variance parameters, 9 estimable fixed effects, 162 observations, identical on both sides. Genotype variance ratio 0.0843 against 0.0854, a difference of 0.0012 against a declared tolerance of 0.05; intercept identical to four decimal places; BLUP agreement r = 1.0000. Repeatability screen stable. This row is validated at the current anchor |
| Sectioned residual, ASReml 3 spelling | `residual = ~ at(env):units` | ASReml 3 alias | inla, brms | **exact** | ASReml 4 removed residual-side at() in favour of dsum(), so current ASReml may reject this spelling. Accepted here for compatibility; dsum() is the canonical form. | Shares the dsum() route and its evidence |
| Unstructured covariance over a factor, on brms | `us(f):g` | ASReml 4 | brms | **covariance-family** | -- | Live emission and parameter-count tests; no current ASReml comparison |
| Unstructured covariance over a factor, on INLA | `us(f):g` | ASReml 4 | inla | **analogue** | Lowered to an iidKd field carrying a joint Wishart prior scaled to the data. The prior geometry and the parameter coordinates are not ASReml's; a supplied SD prior is refused rather than dropped. This row claimed covariance-family until 2026-09-11, when a paired ASReml comparison measured a different covariance rather than a differently-parameterised one, and a metamorphic check found the fit is not invariant to rescaling the response. Correlations from this route are not to be read as ASReml's. | ASReml 4.2.0.392 at this anchor on agridat::besag.met, oracle scenarios 06 and 07 (finding O04). Structural layer exact on both. Numerical layer FAILS: genetic correlation 0.147 against ASReml's 0.930 over two environments, and 0.171 / 0.220 / 0.309 against 0.898 / 0.915 / 0.955 over three, with genotype variances four to five times larger while the residual agrees closely. The repeatability screen rates both rows stable, so the answer is reproducible and different rather than noisy. Multiplying the response by ten -- an exact transformation under which a correlation is invariant -- moves the correlation to -0.622. The cause is not established; the emitted iidKd hyper block has not been read |
| Factor-analytic covariance | `fa(f, k):g` | ASReml 4 | inla | **conditional** | Loadings are anchored to one for factor j at level j, taken from factor-level order. This spans the factor-analytic family only where the leading k-by-k loading block supports the anchoring: a covariance in which the anchor level has zero common-factor loading cannot be represented with finite values. There is no pivot selection and the precondition is not enforced. | Simulated recovery with favourable anchors; one stored ASReml illustration at k = 1. Not compared at this anchor: three identical fits span a factor of 15 in a component, so a single paired comparison would measure the draw rather than the translation |
| Known covariance or precision carrier | `vm(g, K)`, `vm(g, cov = fb_cov(...))` | ASReml 4 | inla, brms | **exact** | Representation metadata (dense / chol / precision / blocks) must be correct. | Dense and sparse carrier tests; no current ASReml fixture for the inverse form |
| Pedigree relationship carrier | `ped(g, A)` | ASReml 3 alias | inla, brms | **exact** | ASReml 4 removed ped() in favour of vm(). Accepted here for compatibility; vm(g, A) is the current spelling. | Shares the vm() carrier machinery and its evidence |
| Correlated random slope | `(x \| g)` | lme4 / brms | inla, brms | **covariance-family** | brms fits this natively. On INLA it is an iid2d block whose Wishart scale is set from the data; repeated identical fits have returned materially different modes where the correlation is weakly identified. | Live INLA recovery tests against lme4; none against ASReml us(2):g |
| Uncorrelated random slope | `(x \|\| g)` | lme4 / brms | inla, brms | **covariance-family** | -- | lme4 REML at this anchor |
| Temporal AR1 on the random side | `~ ar1(t)`, grouped AR1 | ASReml 4 | inla | **exact** | Exact for the latent-field model only: an AR1 field plus an observation residual. This is not a translation of ASReml's residual-side AR1, which carries no nugget. | Independent dense-V REML and simulation at this anchor |
| Separable residual AR1 x AR1 | `residual = ~ ar1(row):ar1(col)` | ASReml 4 | inla, brms | **refused** | The package's latent field plus nugget is a four-parameter neighbour, not ASReml's nugget-free residual, so it refuses rather than substituting. | Refusal is the claim; asserted by name |
| Penalised spline | `spl(x)`, `s(x)` | flexyBayes extension | inla | **analogue** | Lowered to INLA's second-order random walk. ASReml's spline carries its own knot count, duplicate-point handling, prediction points and spline.scale, and no comparison of basis, null space, penalty, effective degrees of freedom or prediction has been made. Related, not the same model. | Emission and fit tests only |
| Gaussian relative weights | `weights = w` | ASReml 4 | inla, brms | **exact** | Relative weights only, giving variance proportional to sigma^2 / w. ASReml's absolute-weight / fixed-dispersion mode is a different estimand and is not carried. | lme4 at this anchor |
| Missing response retained as a prediction target | `na_action = "augment"` | flexyBayes extension | inla, brms | **covariance-family** | A representation-preserving policy rather than a universal ASReml equivalence: the design row is kept as a latent target. The estimator assumptions are the user's to check. | Large missing-plot simulation study and independent dense-V REML |
| Heterogeneous variance with one common correlation | `corh(f):g` | ASReml 4 | inla, brms | **refused** | Neither active engine represents it; it is not the diagonal or the unstructured form. | Refusal is the claim; asserted by name |
| Tensor-product and multivariate smooths | `te(x, z)`, `ti()`, `t2()`, `s(x, z)` | lme4 / brms | inla, brms | **refused** | Univariate smooths only. | Refusal is the claim; asserted by name |
| Structured inner term inside a sectioned residual | `dsum(~ ar1(a) \| b)` | ASReml 4 | inla, brms | **refused** | Reducing it to a per-section variance would drop the structure written. | Refusal is the claim; asserted by name |
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
which carry twelve and thirteen variance parameters and remain unstable.

A user fitting those structures to a trial with weak genetic signal should
refit and compare before acting on a variance component.

## Reading a row

*Factor-analytic covariance* is the row worth reading closely. The lowering
builds iid latent factors with INLA copy coefficients for the loadings and
separate specific variances, and it fixes the loading for factor *j* at level
*j* to one. Fixing a loading is a scale and sign reparameterisation only when
the true loading at that level is not zero. A covariance in which the anchor
level has no common-factor loading, while other levels share the factor, has
no finite representation in this form -- and the fit still runs. The anchor is
taken from factor-level order, so which covariances are representable depends
on how the levels happen to be ordered.

That is why the row is `conditional` rather than `covariance-family`, and why
the summary no longer prints `Representation: exact` for a factor-analytic
term.
