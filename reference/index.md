# Package index

## Verbs

The universal entry, the single-engine pins, the feasibility planner,
cross-engine triangulation, and approximation validation.

- [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  : Bayesian Mixed Models with ASReml Syntax
- [`flexybayes_stream()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes_stream.md)
  : Fit a mixed model to an out-of-core dataset by streaming aggregation
- [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
  : Fit a flexyBayes model via the greta engine
- [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
  : Fit a flexyBayes model via the INLA engine
- [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
  : Fit a flexyBayes model via the brms (Stan) engine
- [`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md)
  : Plan a flexyBayes fit without firing the backend
- [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  : Cross-engine posterior triangulation
- [`validate_approximation()`](https://aagi-aus.github.io/flexyBayes/reference/validate_approximation.md)
  : Validate an approximate model fit against its bias bound
- [`proceed()`](https://aagi-aus.github.io/flexyBayes/reference/proceed.md)
  : Advance a deferred-execution object into its fit

## Ingest adapters

Turn a model specification into a backend-agnostic representation.

- [`fb_from_asreml()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_asreml.md)
  : Ingest an ASReml-format model specification into the flexyBayes IR
- [`fb_from_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_brms.md)
  : Ingest a brms-format formula into the flexyBayes IR
- [`fb_from_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_from_greta.md)
  : Ingest a user-built greta model into the flexyBayes IR

## Constructor nouns

The four classed input constructors: priors, covariance carriers,
approximation schemes, and engine selections.

- [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  : Specify priors via the PC-canonical hybrid DSL

- [`fb_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_cov.md)
  : Construct a structured-covariance carrier

- [`fb_approx()`](https://aagi-aus.github.io/flexyBayes/reference/fb_approx.md)
  : Construct an approximation-scheme specification

- [`fb_engine()`](https://aagi-aus.github.io/flexyBayes/reference/fb_engine.md)
  : Construct an inference-engine specification

- [`is_fb_cov()`](https://aagi-aus.github.io/flexyBayes/reference/is_fb_cov.md)
  :

  Test whether an object is an `fb_cov` carrier

- [`is_fb_approx()`](https://aagi-aus.github.io/flexyBayes/reference/is_fb_approx.md)
  :

  Test whether an object is an `fb_approx` specification

- [`is_fb_engine()`](https://aagi-aus.github.io/flexyBayes/reference/is_fb_engine.md)
  :

  Test whether an object is an `fb_engine` specification

## Accessors and discovery

Routing trace, refusal vocabulary, canonical names, prior summary,
status.

- [`backend_decision()`](https://aagi-aus.github.io/flexyBayes/reference/backend_decision.md)
  : Backend dispatch trace for a flexyBayes fit
- [`canonical_names()`](https://aagi-aus.github.io/flexyBayes/reference/canonical_names.md)
  : Canonical parameter-name view for a flexyBayes fit
- [`fb_backend_status()`](https://aagi-aus.github.io/flexyBayes/reference/fb_backend_status.md)
  : Report inference-backend readiness
- [`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md)
  : Identified covariance for factor-analytic structured-covariance
  terms
- [`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md)
  : List flexyBayes refusal reasons
- [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
  : Resolved-prior summary for a flexyBayes fit
- [`gretaR_status()`](https://aagi-aus.github.io/flexyBayes/reference/gretaR_status.md)
  : Introspect the gretaR backend slot
- [`cat_code()`](https://aagi-aus.github.io/flexyBayes/reference/cat_code.md)
  : Emit the generated backend code for a deferred review object

## Cross-package interop

Tidiers, emmeans, marginaleffects, and draws conversion.

- [`tidy(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes.md)
  : Tidy a flexyBayes fit into a one-row-per-term data frame
- [`tidy(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/tidy.flexybayes_inla.md)
  : Tidy a per-row INLA fit into a one-row-per-term data frame
- [`glance(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/glance.flexybayes.md)
  [`glance(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/glance.flexybayes.md)
  : Glance at a flexyBayes fit
- [`augment(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/augment.flexybayes.md)
  [`augment(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/augment.flexybayes.md)
  : Augment a flexyBayes fit with fitted values and residuals
- [`emm_basis(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/emm_basis.flexybayes.md)
  [`emm_basis(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/emm_basis.flexybayes.md)
  : emmeans support: estimation basis (greta backend)
- [`recover_data(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/recover_data.flexybayes.md)
  [`recover_data(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/recover_data.flexybayes.md)
  : emmeans support: recover model data (greta backend)
- [`get_coef(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_coef.flexybayes.md)
  [`get_coef(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_coef.flexybayes.md)
  : marginaleffects support: fixed-effect coefficients (greta backend)
- [`get_vcov(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_vcov.flexybayes.md)
  [`get_vcov(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_vcov.flexybayes.md)
  : marginaleffects support: covariance (greta backend)
- [`get_predict(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_predict.flexybayes.md)
  [`get_predict(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_predict.flexybayes.md)
  : marginaleffects support: population-level predictions (greta
  backend)
- [`get_data(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_data.flexybayes.md)
  [`get_data(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_data.flexybayes.md)
  : Model data accessor (greta backend)
- [`set_coef(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/set_coef.flexybayes.md)
  [`set_coef(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/set_coef.flexybayes.md)
  : marginaleffects support: set coefficients (greta backend)
- [`fb_as_draws_simple()`](https://aagi-aus.github.io/flexyBayes/reference/fb_as_draws_simple.md)
  : Extract per-parameter posterior draws from a model fit

## Methods

S3 methods on the fit objects. Most are reached through their generic
(for example [`coef()`](https://rdrr.io/r/stats/coef.html),
[`predict()`](https://rdrr.io/r/stats/predict.html)); they are listed
here for completeness.

- [`summary(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md)
  : Summarise a flexybayes object

- [`summary(`*`<flexybayes_aggregated>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes_aggregated.md)
  : Summarise a flexybayes_aggregated object

- [`summary(`*`<flexybayes_glm>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes_glm.md)
  : Summary for flexybayes GLM-compatible object

- [`summary(`*`<fb_plan>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.fb_plan.md)
  :

  Summarise an `<fb_plan>` — verbose form

- [`print(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes.md)
  : Print a flexybayes object

- [`print(`*`<flexybayes_aggregated>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes_aggregated.md)
  : Print a flexybayes_aggregated object

- [`print(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes_brms.md)
  : Print method for the brms-passthrough flexybayes subclass

- [`print(`*`<flexybayes_direct_greta>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes_direct_greta.md)
  : Print a flexybayes object built via fb_greta() (direct greta entry)

- [`print(`*`<fb_approx>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_approx.md)
  :

  Print an `fb_approx` specification

- [`print(`*`<fb_cov>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_cov.md)
  :

  Print an `fb_cov` carrier

- [`print(`*`<fb_engine>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_engine.md)
  :

  Print an `fb_engine` specification

- [`print(`*`<fb_plan>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_plan.md)
  :

  Print an `<fb_plan>` — flight-checklist form

- [`coef(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes.md)
  : Extract fixed effect coefficients

- [`coef(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes_inla.md)
  : Fixed-effect coefficients of a per-row INLA fit

- [`coef(`*`<flexybayes_direct_greta>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes_direct_greta.md)
  : Extract canonical-named posterior means from an fb_greta() fit

- [`confint(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes.md)
  : Credible intervals for fixed effects

- [`confint(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes_brms.md)
  : Credible intervals on the brms path

- [`confint(`*`<flexybayes_glm>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes_glm.md)
  : Credible intervals for flexybayes_glm

- [`vcov(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/vcov.flexybayes.md)
  : Extract variance-covariance matrix of fixed effects

- [`vcov(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/vcov.flexybayes_inla.md)
  : Posterior covariance of a per-row INLA fit's fixed effects

- [`predict(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes.md)
  : Predict from a flexybayes model

- [`predict(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_brms.md)
  : Predict from a brms-passthrough flexybayes fit

- [`predict(`*`<flexybayes_direct_greta>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_direct_greta.md)
  : Predict from a flexybayes_direct_greta fit

- [`predict(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_inla.md)
  : Population-level predictions from a per-row INLA fit

- [`fitted(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/fitted.flexybayes.md)
  : Extract fitted values

- [`fitted(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/fitted.flexybayes_inla.md)
  : In-sample fitted values from a per-row INLA fit

- [`residuals(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/residuals.flexybayes.md)
  : Extract residuals

- [`residuals(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/residuals.flexybayes_inla.md)
  : Response residuals from a per-row INLA fit

- [`plot(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
  [`plot(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
  [`plot(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
  [`plot(`*`<flexybayes_aggregated>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
  [`plot(`*`<flexybayes_direct_greta>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
  [`plot(`*`<flexybayes_glm>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
  : Plot diagnostics for a flexybayes model

- [`anova(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/anova.flexybayes.md)
  : Compare flexybayes models

- [`logLik(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes.md)
  : Log-likelihood (approximate)

- [`logLik(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes_inla.md)
  : Log-likelihood of a per-row INLA fit (not computed)

- [`logLik(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes_brms.md)
  : Log-likelihood on the brms path

- [`nobs(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/nobs.flexybayes.md)
  : Number of observations

- [`family(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/family.flexybayes.md)
  : Extract model family

- [`family(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/family.flexybayes_inla.md)
  : Response family of a per-row INLA fit

- [`formula(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/formula.flexybayes.md)
  : Extract model formula

- [`formula(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/formula.flexybayes_inla.md)
  : Fixed-effect model formula of a per-row INLA fit

- [`model.matrix(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/model.matrix.flexybayes.md)
  : Extract model matrix

- [`update(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/update.flexybayes.md)
  : Update a flexybayes model

- [`as.data.frame(`*`<fb_plan>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/as.data.frame.fb_plan.md)
  :

  Coerce an `<fb_plan>` to data.frame — one row, stable columns

## Specialised models

Extreme-value and Dirichlet regression, genome-wide association, and
genomic selection, with their fit methods and the cross-engine
log-posterior producer.

- [`fb_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gev.md)
  : Fit a generalised extreme value distribution to block maxima
- [`print(`*`<fb_gev_fit>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_gev_fit.md)
  : Print a GEV fit
- [`tidy(`*`<fb_gev_fit>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/tidy.fb_gev_fit.md)
  : Tidy a GEV fit
- [`fb_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_dirichlet.md)
  : Fit a Dirichlet distribution to compositional data
- [`print(`*`<fb_dirichlet_fit>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_dirichlet_fit.md)
  : Print a Dirichlet fit
- [`tidy(`*`<fb_dirichlet_fit>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/tidy.fb_dirichlet_fit.md)
  : Tidy a Dirichlet fit
- [`fb_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gwas.md)
  : Genome-wide association scan (EMMAX / P3D)
- [`triangulate_gwas()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_gwas.md)
  : Triangulate two genome-wide association scans
- [`fb_gblup_cv()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gblup_cv.md)
  : Genomic-prediction accuracy by cross-validation
- [`genomic_summary()`](https://aagi-aus.github.io/flexyBayes/reference/genomic_summary.md)
  : Genomic summary of a fitted relationship model
- [`triangulate_genomic()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_genomic.md)
  : Triangulate genomic model outputs across engines or against a field
  lens
- [`fb_met_summary()`](https://aagi-aus.github.io/flexyBayes/reference/fb_met_summary.md)
  : Breeder summary of a factor-analytic multi-environment-trial fit
- [`fb_log_posterior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_log_posterior.md)
  : Emit a flexyBayes posterior as a log-density producer

## Distribution helpers

Custom response families and their random-number generators.

- [`fb_family_gev()`](https://aagi-aus.github.io/flexyBayes/reference/fb_family_gev.md)
  : Generalised extreme value (GEV) family object
- [`rgev()`](https://aagi-aus.github.io/flexyBayes/reference/rgev.md) :
  Simulate from a generalised extreme value distribution
- [`fb_family_dirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/fb_family_dirichlet.md)
  : Dirichlet family object
- [`rdirichlet()`](https://aagi-aus.github.io/flexyBayes/reference/rdirichlet.md)
  : Simulate from a Dirichlet distribution

## Datasets

- [`met_example`](https://aagi-aus.github.io/flexyBayes/reference/met_example.md)
  : Example Multi-Environment Trial (MET) Dataset

## Package

- [`flexyBayes`](https://aagi-aus.github.io/flexyBayes/reference/flexyBayes-package.md)
  [`flexyBayes-package`](https://aagi-aus.github.io/flexyBayes/reference/flexyBayes-package.md)
  : flexyBayes: Bayesian Mixed Models with ASReml Syntax via greta,
  INLA, brms
