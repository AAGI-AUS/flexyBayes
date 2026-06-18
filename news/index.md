# Changelog

## flexyBayes 0.8.3

A documentation-honesty and ergonomics release on the 0.8.x line. There
are no modelling-behaviour changes; the additions are accessor coverage,
clearer refusals, and metadata / vignette consistency for the first
public development release.

### New features

- **[`glance()`](https://generics.r-lib.org/reference/glance.html) /
  [`augment()`](https://generics.r-lib.org/reference/augment.html) gain
  explicit INLA methods.** Calling
  [`glance()`](https://generics.r-lib.org/reference/glance.html) or
  [`augment()`](https://generics.r-lib.org/reference/augment.html) on an
  INLA fit (`flexybayes_inla`) previously raised a bare “no applicable
  method” error, because that class does not inherit `flexybayes`. They
  now dispatch to informative refusals pointing to
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)
  (coefficient-level summaries) and
  [`summary()`](https://rdrr.io/r/base/summary.html) /
  [`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md)
  (variance components). Net:
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) covers all
  three backends;
  [`glance()`](https://generics.r-lib.org/reference/glance.html) /
  [`augment()`](https://generics.r-lib.org/reference/augment.html) cover
  the greta and brms classes, with a clean message on INLA.

### Minor improvements and fixes

- **[`fb_met_summary()`](https://aagi-aus.github.io/flexyBayes/reference/fb_met_summary.md)
  distinguishes the wrong-backend case.** Passing an INLA or brms fit
  now returns a backend-specific message (breeder summaries need a greta
  factor-analytic fit; the INLA / brms MET path reports variance
  components via [`summary()`](https://rdrr.io/r/base/summary.html) /
  [`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md))
  rather than a generic “not a flexybayes object” error.
- **Development-release signalling and metadata reconciliation.**
  `.onAttach()` and `README` state that all exports are experimental and
  the package is not on CRAN. The package-level documentation,
  `CITATION` / `codemeta` / Zenodo metadata, `_pkgdown.yml`,
  `API_STABILITY.md`, and the security / support / contributor docs are
  reconciled to the 0.8.x line and the three-backend (greta / INLA /
  brms) surface.
- **Vignette convergence disclaimers** propagated to the remaining
  small-budget reference vignettes; every vignette that prints a high
  R-hat now carries an “illustration of output shape, not inference”
  callout.
- `NEWS` broom-coverage wording corrected; package-level roxygen updated
  to three engines and sixteen vignettes.

## flexyBayes 0.8.2

### New features

- **[`fb_log_posterior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_log_posterior.md)
  – a log-posterior producer for downstream tools.** A new exported
  generic turns a fitted flexyBayes object into a vectorised,
  domain-safe, unnormalised log-posterior callable that
  `proxymix::from_fb_posterior()` compresses into a closed-form
  Gaussian-mixture proxy. This is the single inference-result outflow
  from flexyBayes, with a real posterior rather than a mock. The
  returned callable takes a numeric matrix (rows = parameter draws,
  columns = parameters on the natural / constrained scale) and returns
  one unnormalised `log p(theta | data)` per row; it carries
  `parameter_names`, an `NA` `log_normalizer` (a posterior’s marginal
  likelihood is generally unknown – honest, not fabricated), the
  parameters’ `support_lower` / `support_upper` bounds, and the fit’s
  posterior `draws` to seed the consumer’s proposal.

  The **greta** backend is the canonical real producer: it evaluates the
  retained model graph’s unadjusted joint density at the free-state
  image of the supplied natural-scale parameters, which is the
  unnormalised natural-scale log-posterior exactly (validated against an
  analytic conjugate posterior to machine precision – correlation 1,
  constant offset). Out-of-support rows return `-Inf` rather than
  erroring. The **brms** and **INLA** backends abstain with an
  informative, classed `fb_c4_unavailable` condition: brms’s log-density
  lives on the Stan unconstrained scale with a version-fragile name
  mapping, and INLA’s posterior is a deterministic Laplace / grid
  approximation rather than a sampling log-density, so an honest abstain
  is preferred to a plausible-but-wrong producer (the Independent Oracle
  Principle). flexyBayes does not depend on proxymix – the cross-package
  demonstration lives in a separate integration harness, preserving the
  acyclic dependency graph.

## flexyBayes 0.8.1

A hub-ergonomics release: a single broom dialect across the three
backends, two standalone maximum-likelihood fitters for distributions
outside the GLM emit path, and an INLA-led multi-environment-trial
vignette.

### New features

- **[`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  [`glance()`](https://generics.r-lib.org/reference/glance.html) /
  [`augment()`](https://generics.r-lib.org/reference/augment.html) – one
  broom dialect across the hub.** The broom-style
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) generic
  (re-exported from `generics`, now a lightweight `Imports`) returns a
  flat one-row-per-term `data.frame` with the canonical `term` /
  `estimate` / `std.error` / `conf.low` / `conf.high` columns for the
  greta (`flexybayes`), brms (`flexybayes_brms`), and INLA
  (`flexybayes_inla`) fit classes alike – the INLA class gains its own
  [`tidy.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes_inla.md)
  method (it does not inherit from `flexybayes`).
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and
  [`augment()`](https://generics.r-lib.org/reference/augment.html) cover
  the greta and brms classes; on an INLA fit they error (turned into an
  informative refusal in 0.8.3). A cross-engine comparison table is now
  an [`rbind()`](https://rdrr.io/r/base/cbind.html) of two
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) outputs
  rather than a hand-built reconciliation of three different backend
  layouts. Mirrors the
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) method
  kernR adopted, so the orchestra speaks one tidy dialect.
- **[`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md)
  – generalised extreme value (block-maxima) fitter.** Fits the
  location, scale, and shape of a GEV distribution to block maxima
  (annual maximum rainfall, peak yields) by dependency-free maximum
  likelihood, and reports return levels for the requested return
  periods. The family descriptor
  [`fb_family_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_family_gev.md)
  and the simulator
  [`rgev()`](https://aagi-aus.github.io/flexyBayes/reference/rgev.md)
  ship alongside. `family = "gen_extreme_value"` in
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  now routes to
  [`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md)
  with an explicit pointer rather than a generic refusal – block maxima
  have no GLM mean-link, so they do not belong on the formula emit path.
  (A scalable Bayesian GEV on INLA’s native `gev` family is planned.)
- **[`fb_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_dirichlet.md)
  – compositional (simplex) Dirichlet fitter.** Fits the concentration
  vector of a Dirichlet distribution to compositional rows (soil texture
  fractions, species abundance, allele frequencies) by maximum
  likelihood (the default, dependency-free) or via greta’s native
  `dirichlet` distribution (`method = "greta"`), and reports the fitted
  mean composition. The family descriptor
  [`fb_family_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_family_dirichlet.md)
  and the simulator
  [`rdirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/rdirichlet.md)
  ship alongside; `family = "dirichlet"` in
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  routes here.

### Minor improvements and fixes

- The MET and genomics vignette now leads with the scalable INLA MET
  path (the diagonal genotype-by-environment model fits in seconds with
  trustworthy posteriors) as the recommended route, and presents the
  greta factor-analytic route as the slower, harder-mixing alternative
  for the stability decomposition.
- `generics` moves from `Suggests` to `Imports`, so
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  [`glance()`](https://generics.r-lib.org/reference/glance.html) /
  [`augment()`](https://generics.r-lib.org/reference/augment.html) are
  available without attaching `broom`; the runtime
  [`registerS3method()`](https://rdrr.io/r/base/ns-internal.html) shim
  in `.onLoad()` is replaced by static `S3method()` registration.

## flexyBayes 0.8.0

### Genomics and MET expansion

- **[`triangulate_genomic()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_genomic.md)
  /
  [`triangulate_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_gwas.md)
  – cross-engine and field-standard genomic triangulation.**
  [`triangulate_genomic()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_genomic.md)
  compares two GBLUP / pedigree analyses on heritability, the variance
  components, and the breeding values (matched by genotype) – either two
  flexyBayes fits (does greta agree with brms and INLA?) or a flexyBayes
  fit against a generic *genomic lens* (`list(h2, var_g, var_e, gebv)`),
  the form a field-standard REML answer from sommer supplies.
  [`triangulate_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_gwas.md)
  compares two genome scans by the agreement that matters – the Jaccard
  overlap of the significant marker sets, the top-marker overlap, and
  the effect correlation. Both carry the same shared-upstream caveat as
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  (agreement is not correspondence). flexyBayes core never depends on
  the field tools; the lens form lets the companion build the koine
  fourth opinion. Breeding-value labels are now the genotype factor
  levels on every backend (previously the greta and INLA paths used
  positional labels), so GEBVs match across engines.
- **[`fb_met_summary()`](https://aagi-aus.github.io/flexyBayes/reference/fb_met_summary.md)
  – breeder summary of a factor-analytic MET fit.** For a
  `fa(env, k):gen` factor-analytic genotype-by-environment fit it
  reports the quantities a plant breeder acts on: each genotype’s
  overall performance (the across-environment mean of its realised
  effects) and stability (the across-environment spread), the
  genotype-by-environment BLUPs, the factor loadings, and the
  environment genetic-correlation matrix – the crossover structure
  (negative correlations are rank reversals across environments). These
  come from the *realised* effects, which are identified (invariant to
  the loadings’ rotation and sign ambiguity), so their posteriors are
  interpretable. The greta factor-analytic fit now monitors those
  realised genotype-by-environment effects (it previously monitored only
  the loadings and specific variances).
- **[`fb_gblup_cv()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gblup_cv.md)
  – genomic-prediction accuracy by cross-validation.** The payoff layer
  of genomic selection: how well does a GBLUP trained on phenotyped
  genotypes predict the held-out performance of genotypes seen only
  through their markers? Each fold estimates the variance components by
  REML on the training set and predicts the held-out breeding values
  from the relationship matrix (the exact GBLUP prediction equation,
  evaluated through the spectral primitive), then reports prediction
  accuracy (predicted-observed correlation), bias (the observed-on-
  predicted slope), and the per-fold spread, with repeated-CV averaging.
  Validated against an exact full-matrix prediction; accuracy is near
  zero for an unrelated / non-heritable trait, rises with heritability,
  and is substantial when genotypes are related – the realistic
  selection setting.
- **[`fb_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gwas.md)
  – genome-wide association scan.** A whole-genome EMMAX / P3D scan
  (Kang et al. 2010): the polygenic null mixed model is fit once by REML
  to estimate the variance components, then every marker is tested by
  exact generalised least squares under those fixed components. The
  shared spectral primitive turns the per-marker test into an `O(n)`
  weighted least squares after a single eigendecomposition, so the scan
  is feasible without a per-marker model fit, and it needs no MCMC
  backend – it is a deterministic frequentist fast path (the backends
  enter only at optional top-hit refinement). Returns marker effects,
  standard errors, chi-square statistics, p-values, Bonferroni and
  Benjamini-Hochberg FDR adjustments, the genomic-control inflation
  factor, and the null REML heritability, with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (Manhattan /
  QQ) methods. The REML variance components are validated against an
  exact full-matrix reference and against ’s independent REML; the
  per-marker statistic is validated against exact per-marker GLS.
- **Genomic BLUP is now three-engine triangulatable.** The genomic /
  pedigree relationship random effect `vm(geno, G)` / `ped(animal, A)`
  reaches all three backends: greta and brms via the dense relationship
  matrix (brms’s native `gr(geno, cov = G)` group term, which Cholesky-
  factors the covariance internally), and INLA via the precision carrier
  `vm(geno, cov = fb_cov(solve(G), type = "precision"))` (the `generic0`
  sparse-precision path – a dense GBLUP precision is INLA-feasible but
  forgoes INLA’s sparsity advantage, so it is opt-in by carrier rather
  than silent). A simulated heritability is recovered on all three
  engines with mutually-close posterior means, so
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  can cross-check a GBLUP fit across paradigms.
- **[`genomic_summary()`](https://aagi-aus.github.io/flexyBayes/reference/genomic_summary.md)**
  extracts the breeder-facing quantities from a fitted relationship
  model on any backend: narrow-sense heritability , genomic estimated
  breeding values (GEBVs) with posterior reliability, and the genetic /
  residual variances. The greta GBLUP fit now also monitors the
  breeding-value vector, so GEBVs are available on greta, brms, and INLA
  alike (previously the greta path reported only the variance
  components).
- **Spectral efficiency primitive (foundations).** A shared internal
  eigendecomposition primitive now underpins the genomics / MET work: it
  decomposes a relationship matrix once and exposes the rotation that
  turns a structured genetic random effect into an independent one. This
  is the machinery a genome-wide scan reuses across every marker (the
  rotated model has a diagonal residual covariance, so each marker is an
  `O(n)` score test rather than a fresh mixed-model fit), and that
  scales genomic BLUP and variance-component estimation. It is
  positive-semidefinite aware – numerical-noise negative eigenvalues are
  clamped and reported, while genuine indefiniteness is refused – and
  validated against full-covariance generalised least squares to machine
  tolerance. Internal at this stage; the fit routes that consume it
  follow.
- **Genomic output contract (foundations).** A standardised genomic
  summary (narrow-sense heritability , genomic estimated breeding values
  with posterior reliability, and – for whole-genome marker models –
  marker effects with posterior retention probabilities) is now computed
  engine- agnostically from posterior draws, so the greta / INLA / brms
  paths all feed the same triangulatable result. The fit-level accessor
  lands with the multi-backend GBLUP route.

### Triangulation honesty

- [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  gained a `data_independence` argument and a `shared_upstream_caveat`
  result field (Independent Oracle Principle). It measures inter-fit
  *agreement*, and the backend-independence registry certifies code (not
  data) independence – so if both fits share the same upstream data, a
  fabricated data fact is common-mode and their agreement cannot detect
  it. Unless the caller declares `data_independence = TRUE`, the result
  carries a caveat (surfaced prominently by
  [`print()`](https://rdrr.io/r/base/print.html)) that agreement does
  not test a shared upstream data fact. Pure metadata; the metrics are
  unchanged.

### Engine reliability and diagnostics

- **Faithful default prior on the INLA backend.** The package’s default
  `uniform(0, scale)` prior on every variance-component standard
  deviation (and any user-supplied `uniform()` / `half_normal()` prior)
  is now represented *exactly* on the INLA backend via an
  expression-prior on the log-precision, replacing the former PC-prior
  approximation. The PC approximation concentrated prior mass at zero
  and so shrank a small-group variance component more than the greta
  backend’s flat uniform did, producing a cross-engine prior mismatch
  that surfaced as spurious
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  disagreement on the variance components for models with few groups.
  The two engines now carry genuinely the same default prior, so their
  variance-component posteriors agree far more closely. This changes the
  default INLA variance-component posterior for random-effects models
  fit without an explicit prior.
- **Convergence warning.** MCMC fits (`greta`, `brms`) now emit a
  warning when the sampler may not have converged (a parameter with Rhat
  at or above 1.1, or a low effective sample size), rather than
  surfacing it only as a print-method badge. Silence it for
  intentionally short fits with
  `options(flexyBayes.silence_convergence_warning = TRUE)`. The
  deterministic INLA path carries no such warning.
- **New
  [`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md)**
  reports the *identified* covariance for factor-analytic `fa()` terms,
  with an entrywise Rhat. The raw loadings are identified only up to
  rotation and sign, so their per-entry Rhat is meaningless; is
  rotation- and sign-invariant and is the quantity whose convergence is
  interpretable. The convergence warning points to it when a
  factor-analytic / unstructured term is present.
- **Fixed** an intercept-only model on the aggregated INLA path
  (`y ~ 1 + (1 | g)`, `backend = "inla"`, `aggregate = TRUE`), which
  crashed in the fixed-effect summariser because a length-one variance
  vector was misread by [`diag()`](https://rdrr.io/r/base/diag.html) as
  a dimension. The one-coefficient covariance is now built correctly.

### Usability

- **New
  [`fb_backend_status()`](https://aagi-aus.github.io/flexyBayes/reference/fb_backend_status.md)**
  reports which inference backends are installed and usable in the
  current session (greta additionally needs a reachable Python /
  TensorFlow stack), with an actionable install hint per backend. It is
  read-only and runs without any backend present.

### Lean-core split

- The orchestra-composition layer has been extracted to the companion
  package **flexyBayesOrchestra** (`Imports: flexyBayes`): the surrogate
  emulators (`fit_surrogate()`, `fb_surrogate_ies()`,
  `fb_surrogate_ppc()`), ensemble sources (`fb_ensemble()`,
  `as_fb_ensemble()`, `read_fb_ensemble()`, `verify_fb_ensemble()`,
  `register_ensemble_source()`), the PESTO ensemble-derived priors
  (`fb_pesto()`, `fb_prior_from_ensemble()`), the surrogate conformers
  (`register_pesto_surrogates()`, `register_kernr_surrogates()`,
  `register_surrogate()`), and the dormant koine fourth-opinion slot
  (`koine_status()`). Install `flexyBayesOrchestra` for these;
  `flexyBayes` alone now ships the mixed-model and cross-engine
  triangulation core.
- Dropped `kernR`, `PESTO`, `S7`, and `gretaR` from `Suggests` and
  removed the `Remotes:` field, so the package installs cleanly from
  CRAN-style repositories with no GitHub- or r-universe-only
  dependencies. The `greta` (CRAN) and `INLA` backends remain; `INLA` is
  served by its own `Additional_repositories`. The `gretaR` R-native
  engine is still wired as a dormant, opt-in backend: install `gretaR`
  yourself and it is detected at run time, but it is no longer a
  declared dependency of the public core.

### Methodology

- **Identified factor-analytic loadings.** `fa(x, k)` now emits the
  loadings matrix `Lambda` with the standard lower-triangular,
  positive-diagonal identification (Lopes & West 2004): the strict upper
  triangle is zeroed and the diagonal is constrained positive.
  Previously `Lambda` was a free `normal(0, 1)` matrix, which is
  unidentified up to rotation, sign, and column permutation for `k > 1`,
  so its posterior summaries were not interpretable. Predictions and the
  g-side covariance are unchanged – the product `F %*% t(Lambda)` is
  rotation-invariant – so only the loadings interpretation is affected.
- **Factor-analytic rank upper bound.** `fa(x, k)` is now refused at fit
  time when `k` is not strictly below the number of levels of `x` (new
  refusal `fa_rank_exceeds_dim`). A factor-analytic covariance is
  identifiable only for `k < n_outer`: at `k = n_outer` the loadings and
  specific variances are an over-parameterised reparameterisation of the
  unstructured form, and at `k > n_outer` the lower-triangular loadings
  carry empty columns. The check is data-aware (the number of levels is
  known only once the term is matched against the data), complementing
  the existing data-free `k >= 1` floor (`fa_rank_invalid`). The refusal
  message points to `us(x)` for a full unstructured covariance.
- **Removed the `rhat_means` triangulation metric.**
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  no longer reports `rhat_means`. It pooled two *different* engines’
  posteriors as if they were chains of a single sampler, which conflates
  genuine between-engine approximation bias with within-sampler
  non-convergence – so it was not a valid convergence diagnostic, and it
  was mislabelled “rank-normalised R-hat” although it applied no
  rank-normalisation. Cross-engine discrepancy is reported by the
  distributional metrics: `wasserstein_1`, `sd_ratio`, `mean_diff`, and
  the quantile differences. (The per-fit, within-engine rank-normalised
  R-hat used as a convergence diagnostic is unaffected.)
- **Honest default-prior provenance.** The default variance-component
  prior (bounded uniform on each SD; `U = 5 * sd(y)` for Gaussian) was
  attributed “following Gelman (2006)”. Gelman (2006) in fact recommends
  a half-t / half-Cauchy for variance components with few groups and
  cautions against a flat prior there; the uniform-on-SD default is a
  weakly-informative choice for *moderate* group counts, and the
  `5 * sd(y)` bound is a flexyBayes heuristic. The attribution is
  corrected throughout (the one-time announcement message,
  [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md),
  the docstrings, the README, and the regression + priors vignettes).
  The default behaviour is unchanged; `fb_prior(half_cauchy(...))`
  remains the documented choice for small `J` (see the priors vignette).

### Documentation

- **Vignettes migrated to the engine-pin API.** Eight vignettes still
  showed the pre-0.5.0 universal
  `fb_brms(..., backend = "greta" / "inla" / "auto")` call surface,
  which the 0.5.0 engine-pin refactor removed
  ([`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
  now pins to Stan and rejects a conflicting `backend`). Every example
  is updated to the current surface:
  [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
  /
  [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
  for a fixed engine,
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  /
  [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  for the universal entry that takes a `backend` argument. The
  brms-style *grammar* is unchanged and still accepted by every entry;
  only the engine selection moved. Stale prose is corrected accordingly
  (the `backend` default is `"auto"`, not `"greta"`; the Stan/brms emit
  backend is live, not “queued”; the native-greta canonical-name map is
  supplied via `fb_from_greta(model, canonical_names = ...)`). Vignettes
  were re-precompiled.

## flexyBayes 0.7.0

> **Note (current state).** The orchestra-composition features described
> in this section — the ensemble-source data contract, simulator-derived
> priors, the surrogate emulators and predictive checks, and the dormant
> fourth-opinion slot — were subsequently extracted to the companion
> package **flexyBayesOrchestra** (see the lean-core split under
> “(development version)” above). They are no longer part of the
> `flexyBayes` exported surface; the entries below record the 0.7.0
> history.

This release reconciles three parallel development streams onto the
engine, data, and surrogate axes: gretaR activated as a fourth inference
engine (the engine axis), the ensemble-source data contract and its
first consumers (the data axis), and the surrogate carried end-to-end
from a reference emulator to a distribution-preserving predictive check
(the surrogate axis).

### New features

- **gretaR activated as a fourth inference engine.**
  `flexybayes(..., backend = "gretaR")` fits supported hierarchical
  models through an out-of-process torch-NUTS worker, returning
  canonical parameter names and a `draws_array` that flows into
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  like any other engine. The integration is governed by a versioned
  backend contract with an **executable conformance battery** – a new
  engine onboards by adding one descriptor, gated on
  [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)-agreement
  within an SBC-calibrated threshold. A **dormant koine backend slot**
  is provisioned (scaffolded and correct, switched on when koine’s
  programmatic model builder is confirmed against the contract).
- **`fb_ensemble()` and the `ensemble_source` data contract.** A
  canonical ingest shape for calibrated parameter ensembles, with PESTO
  (via its manifest), raw `apsimx`, and file as conforming producers –
  so flexyBayes is standalone-functional for methods work without
  depending on any one producer’s internals. **`fb_pesto()`** turns a
  calibrated ensemble into an informative `fb_prior` for a downstream
  fit (moments method; joint and KDE methods refuse cleanly behind their
  gates).
- **[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), and
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) methods for INLA
  fits.** [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) now return
  values for a `flexybayes_inla` fit (previously `NULL`);
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) returns an honest
  `NA` with the right shape where INLA does not expose a likelihood,
  rather than erroring.
- **`fb_surrogate_ppc(whiten = TRUE)` whitens the outputs before the
  MMD.** The predictive draws and `observed` are Mahalanobis-whitened by
  their pooled covariance – a label-agnostic representation that puts
  the outputs on a common footing, removing the cross-output scale
  heterogeneity that makes a single median bandwidth too coarse to
  resolve a distortion in a low-variance direction. The map depends only
  on the pooled cloud (not the draws-vs-observed labels), so the
  permutation test stays valid, and the same map is applied to both
  clouds, so a relative scale, shape, or mean difference is preserved.
  It lifts power substantially on localised / heterogeneous-scale
  distortions a raw-output check misses (with no blind spot found across
  a mean / variance / shape battery), at a valid level. Opt-in (a
  graduated ridge guards a near-singular pooled covariance); composes
  with `project` and `aggregate`.
- **`fb_surrogate_ppc(aggregate = TRUE)` runs an aggregated multi-kernel
  MMD.** An MMDAgg test (Schrab et al. 2023) over an RBF bandwidth grid
  around the median heuristic, sharing the permutations across
  bandwidths and aggregating by the weighted-quantile statistic
  `min_g p_g / w_g` – it removes the bandwidth choice and is robust when
  the median heuristic is mis-scaled, at a valid level. The kernel
  weights are set by `agg_weights` (`"uniform"` default, or
  `"increasing"` / `"decreasing"` / `"centred"`, or a numeric vector):
  up-weighting the scale a distortion is expected to live at recovers
  power, and `"uniform"` recovers the min p-value aggregation. It
  orchestrates kernR’s single-kernel `mmd_test` over the grid (the MMD
  machinery stays kernR’s) and composes with `project`. Returns an
  `fb_mmd_agg` object.
- **New vignette: distribution-preserving surrogates (emulate, invert,
  check).** A self-contained walkthrough of the surrogate workflow –
  ingest an ensemble, fit a surrogate, read its predictive distribution,
  invert it with `fb_surrogate_ies()`, and check it with
  `fb_surrogate_ppc()` – runnable on the built-in reference surrogate
  with no other package.
- **`fb_surrogate_ppc(project = TRUE)` runs a factor-space MMD for large
  output counts.** For a surrogate that stores a low-rank covariance
  (`cov_rank`), the predictive-check draws and `observed` are projected
  onto the surrogate’s `q` factor directions before the MMD, so the test
  runs in `q` dimensions, not `k` – cheaper and higher-powered on the
  correlated cross-output structure (the MMD loses power as the
  dimension grows). It is, by design, blind to distortions orthogonal to
  the factor space (idiosyncratic per-output errors, which the marginal
  and `joint = FALSE` checks cover), so it complements the full-space
  check.
- **`fit_surrogate(cov_method = "auto")` picks the low-rank factor
  method.** In addition to `"eigen"` (default) and `"fa"`, `"auto"`
  selects between them from the heteroscedasticity of the eigen-residual
  diagonal – `"fa"` when the residual is uneven (where eigen-truncation
  is biased), `"eigen"` otherwise. The resolved method is recorded in
  the fitted object.
- **`fit_surrogate(ard_refit = "cv")` now uses k-fold
  cross-validation.** The held-out predictive-density objective is
  summed over a k-fold partition of the subsample rather than a single
  split – a less noisy objective, so the Nelder-Mead lands on sharper
  lengthscales.
- **`fb_surrogate_ppc()` samples a low-rank covariance factored.** When
  the surrogate stores a low-rank covariance (`cov_rank`), the joint
  predictive-check draws are sampled from the factored form directly – a
  `q`-dimensional factor noise plus a `k`-dimensional idiosyncratic
  noise – so the `k x k` covariance slice is never formed (the native
  low-rank consumer for the PPC, `O(k q)` per draw). Draws from the same
  per-row covariance as the dense path.
- **`fb_surrogate_ies()` consumes a low-rank covariance natively
  (Woodbury).** When the surrogate stores a low-rank covariance
  (`cov_rank`), the ES-MDA Kalman update folds the factored form
  straight in via the Woodbury identity (a rank-sized solve) and samples
  the observation perturbation factored, so the `d x d` data covariance
  is never formed – the inversion cost scales with the covariance rank,
  not the output count. Produces the same Kalman gain as the dense path
  (a dense-joint or marginal surrogate keeps the dense path).
- **`fit_surrogate(cov_method = "fa")` learns the low-rank factor by
  EM.** In addition to the default eigen-truncation (`"eigen"`), the
  rank-`cov_rank` factor can be fitted by factor-analysis EM, which
  recovers the loadings more faithfully when the residual diagonal is
  heteroscedastic (the eigenvectors of the full covariance are otherwise
  contaminated by the uneven diagonal). The covariance diagonal stays
  exact either way.
- **`fit_surrogate(ard_refit = "cv")` regularises the fuller
  evidence-ARD.** The `ard_refit = TRUE` profiling maximises the
  in-sample evidence and can over-shorten lengthscales; `"cv"` instead
  maximises the held-out predictive log-density on a train/validation
  split, keeping the relevance sharpening without the in-sample
  over-shortening.
- **`fb_surrogate_ppc()` gains a joint draw path.** With `joint = TRUE`
  (default) the predictive-check draws are sampled from the full
  cross-output predictive covariance (`predict(cov = TRUE)` + a per-row
  Cholesky), so the emulator’s cross-output correlation is exercised –
  the second consumer of the joint covariance after
  `fb_surrogate_ies()`. On a correlated-output simulator, marginal draws
  (`joint = FALSE`) discard that structure and false-reject a faithful
  surrogate, so the joint path is the correct distribution-preserving
  check. A single output, or a marginals-only surrogate, behaves exactly
  as before.
- **`fit_surrogate(method = "rff")` gains `cov_rank` for large output
  counts.** Stores the cross-output noise covariance in a rank-`q`
  factor-analysis form (`U U' + diag`), `O(k q)` instead of `O(k^2)`.
  The per-output marginal sd and the covariance diagonal are reproduced
  exactly; only the off-diagonal covariance is approximated.
  `predict(<fb_surrogate>, newdata, cov = "lowrank")` returns the
  compact factored predictive (`cov_lowrank = list(u, d, infl)`), and
  the new `fb_cov_slice()` reconstructs a `k x k` covariance slice from
  either the dense or the factored form.
- **`fit_surrogate(method = "rff")` gains `ard_refit` (fuller
  evidence-ARD).** When `TRUE`, the prior/noise precisions `(a, t)` are
  re-estimated by evidence at each candidate lengthscale (the evidence
  profiled over `(a, t)`), a fuller – and at equal lengthscales never
  lower – marginal likelihood. It sharpens relevance ranking but can
  over-shorten lengthscales (a higher evidence does not guarantee better
  held-out recovery for finite-feature models), so it is opt-in
  (`ard_refit = FALSE` default).
- **The reference `rff` surrogate gains ARD lengthscales and a joint
  multi-output fit, both on by default.**
  `fit_surrogate(method = "rff")` now learns per-input (automatic
  relevance determination) lengthscales by maximising the model evidence
  (`ard = TRUE`), so irrelevant parameters are stretched out of the
  kernel; and it fits a joint matrix-normal – inverse-Wishart model
  across outputs (`joint = TRUE`), preserving the cross-output
  predictive covariance rather than only the marginals. Both default
  behaviours can be turned off (`ard = FALSE`, `joint = FALSE`) to
  recover the isotropic, per-output fit. The per-output `<o>_sd`
  contract is unchanged.
- **`predict(<fb_surrogate>, newdata, cov = TRUE)` returns the joint
  predictive distribution**: an `fb_surrogate_prediction` with the
  predictive mean matrix and a per-row cross-output covariance array. A
  marginals-only surrogate (independent `rff`, or a conformer reporting
  only `<o>_sd`) returns a diagonal covariance, so the call is always
  honest.
- **`fb_surrogate_ies()` consumes the joint predictive covariance.**
  When the surrogate is joint, the full cross-output predictive
  covariance enters the ES-MDA data covariance
  (`use_surrogate_cov = TRUE`, default), so correlated emulator errors
  are assimilated coherently; a marginals-only surrogate falls back to
  the diagonal path, identical to before.
- **Sibling surrogate backends register as conformers.** PESTO’s GP
  surrogate (`method = "gp"`) and kernR’s conditional mean embedding
  (`method = "cme"`) are available to `fit_surrogate()` when the member
  package is installed, registered through the existing surrogate
  registry. PESTO’s GP already reports a per-output predictive mean and
  variance; kernR’s CME is a kernel linear smoother, so the adapter
  supplies the standard linear-smoother predictive variance
  `sigma^2 (1 + ||w||^2)` to honour the distribution-preserving `<o>_sd`
  requirement. Both register opportunistically on load, or explicitly
  via `register_pesto_surrogates()` / `register_kernr_surrogates()`. The
  surrogate methods are consumed through the contract, never copied into
  flexyBayes, and the built-in `rff` reference keeps the package
  standalone.
- **`fb_surrogate_ies()` calibrates a parameter ensemble through a
  surrogate.** An Ensemble Smoother with Multiple Data Assimilation
  (ES-MDA) drives a prior parameter ensemble toward observed data using
  the surrogate in place of the forward model; the predictive `<o>_sd`
  is folded into the data covariance, so emulator uncertainty inflates
  the assimilation noise.
- **`fb_surrogate_ppc()` checks whether a surrogate distorted the output
  distribution.** It draws from the surrogate’s predictive distribution
  at held-out parameter sets and runs kernR’s `mmd_ppc()` against the
  true simulator outputs there – a distribution-level check, not a
  mean-only accuracy measure.

### Bug fixes

- **Registered surrogate conformers now predict on the standardised
  design.** `predict.fb_surrogate()` previously passed raw `newdata` to
  a registered conformer’s
  [`predict()`](https://rdrr.io/r/stats/predict.html) while training it
  on the standardised design, so train and predict disagreed on scale.
  The standardised design is now passed consistently.

------------------------------------------------------------------------

Older release history (0.6.0 and earlier) lives in `NEWS_ARCHIVE.md`.
