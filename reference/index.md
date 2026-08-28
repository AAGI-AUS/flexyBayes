# Package index

## Verbs

The universal entry, the single-engine pins, the feasibility planner,
cross-engine triangulation, and approximation validation.

- [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  : Bayesian Mixed Models with ASReml Syntax
- [`flexybayes_stream()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes_stream.md)
  : Fit a mixed model to an out-of-core dataset by streaming aggregation
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
- [`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md)
  : List flexyBayes refusal reasons
- [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
  : Resolved-prior summary for a flexyBayes fit
- [`ranef()`](https://aagi-aus.github.io/flexyBayes/reference/ranef.md)
  : Random-effect predictions from a flexyBayes fit
- [`fb_complete_grid()`](https://aagi-aus.github.io/flexyBayes/reference/fb_complete_grid.md)
  : Complete a design grid before fitting
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
  : emmeans support: estimation basis (default method)
- [`recover_data(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/recover_data.flexybayes.md)
  [`recover_data(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/recover_data.flexybayes.md)
  : emmeans support: recover model data (default method)
- [`get_coef(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_coef.flexybayes.md)
  [`get_coef(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_coef.flexybayes.md)
  : marginaleffects support: fixed-effect coefficients (default method)
- [`get_vcov(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_vcov.flexybayes.md)
  [`get_vcov(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_vcov.flexybayes.md)
  : marginaleffects support: covariance (default method)
- [`get_predict(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_predict.flexybayes.md)
  [`get_predict(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_predict.flexybayes.md)
  : marginaleffects support: population-level predictions (default
  method)
- [`get_data(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_data.flexybayes.md)
  [`get_data(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/get_data.flexybayes.md)
  : Model data accessor (default method)
- [`set_coef(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/set_coef.flexybayes.md)
  [`set_coef(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/set_coef.flexybayes.md)
  : marginaleffects support: set coefficients (default method)
- [`fb_as_draws_simple()`](https://aagi-aus.github.io/flexyBayes/reference/fb_as_draws_simple.md)
  : Extract per-parameter posterior draws from a model fit
- [`as_draws_df(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/as_draws_df.flexybayes.md)
  [`as_draws(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/as_draws_df.flexybayes.md)
  [`as_draws_matrix(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/as_draws_df.flexybayes.md)
  : Posterior draws from a flexyBayes fit
- [`loo(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/loo.flexybayes.md)
  : Approximate leave-one-out cross-validation for a flexyBayes fit
- [`pp_check(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/pp_check.flexybayes.md)
  : Posterior predictive check for a flexyBayes fit

## Methods

S3 methods on the fit objects. Most are reached through their generic
(for example [`coef()`](https://rdrr.io/r/stats/coef.html),
[`predict()`](https://rdrr.io/r/stats/predict.html)); they are listed
here for completeness.

- [`summary(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md)
  : Summarise a flexyBayes fit

- [`print(`*`<summary.flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.summary.flexybayes.md)
  : Print a flexyBayes summary

- [`summary(`*`<flexybayes_aggregated>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes_aggregated.md)
  : Summarise a fit run on the aggregated representation

- [`summary(`*`<flexybayes_glm>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes_glm.md)
  : Summary for flexybayes GLM-compatible object

- [`summary(`*`<fb_plan>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/summary.fb_plan.md)
  :

  Summarise an `<fb_plan>` — verbose form

- [`print(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes.md)
  : Print a compact description of a flexyBayes fit

- [`print(`*`<flexybayes_aggregated>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes_aggregated.md)
  : Print a flexybayes_aggregated object

- [`print(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.flexybayes_brms.md)
  : Print method for the brms-passthrough flexybayes subclass

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
  : Extract coefficients from a flexyBayes fit

- [`coef(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/coef.flexybayes_inla.md)
  : Coefficients of a per-row INLA fit

- [`confint(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes.md)
  : Credible intervals for the fixed effects of a flexyBayes fit

- [`confint(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes_brms.md)
  : Credible intervals on the brms path

- [`confint(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes_inla.md)
  : Credible intervals for the fixed effects of a per-row INLA fit

- [`confint(`*`<flexybayes_glm>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/confint.flexybayes_glm.md)
  : Credible intervals for flexybayes_glm

- [`vcov(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/vcov.flexybayes.md)
  : Extract variance-covariance matrix of fixed effects

- [`vcov(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/vcov.flexybayes_inla.md)
  : Posterior covariance of a per-row INLA fit's fixed effects

- [`predict(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_brms.md)
  : Predict from a brms-passthrough flexybayes fit

- [`predict(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_inla.md)
  : Population-level predictions from a per-row INLA fit

- [`print(`*`<fb_predict_classify>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/print.fb_predict_classify.md)
  : Print a classify means table

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
  [`plot(`*`<flexybayes_glm>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/plot.flexybayes.md)
  : Plot diagnostics for a flexyBayes model

- [`anova(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/anova.flexybayes.md)
  : Compare flexyBayes models on a plug-in information criterion

- [`logLik(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes.md)
  : Plug-in conditional log-likelihood of a flexyBayes fit

- [`logLik(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes_inla.md)
  : Log-likelihood of a per-row INLA fit (not computed)

- [`logLik(`*`<flexybayes_brms>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes_brms.md)
  : Log-likelihood on the brms path

- [`nobs(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/nobs.flexybayes.md)
  : Number of observations a flexyBayes fit was fitted to

- [`family(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/family.flexybayes.md)
  : Extract model family

- [`family(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/family.flexybayes_inla.md)
  : Response family of a per-row INLA fit

- [`formula(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/formula.flexybayes.md)
  : Extract model formula

- [`formula(`*`<flexybayes_inla>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/formula.flexybayes_inla.md)
  : Fixed-effect model formula of a per-row INLA fit

- [`model.matrix(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/model.matrix.flexybayes.md)
  : Fixed-effect model matrix of a flexyBayes fit

- [`update(`*`<flexybayes>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/update.flexybayes.md)
  : Re-fit a flexyBayes model with modified arguments

- [`as.data.frame(`*`<fb_plan>`*`)`](https://aagi-aus.github.io/flexyBayes/reference/as.data.frame.fb_plan.md)
  :

  Coerce an `<fb_plan>` to data.frame — one row, stable columns

## Genomics

Genomic selection and the cross-engine genomic accessors.

- [`fb_gblup_cv()`](https://aagi-aus.github.io/flexyBayes/reference/fb_gblup_cv.md)
  : Genomic-prediction accuracy by cross-validation
- [`genomic_summary()`](https://aagi-aus.github.io/flexyBayes/reference/genomic_summary.md)
  : Genomic summary of a fitted relationship model
- [`triangulate_genomic()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate_genomic.md)
  : Triangulate genomic model outputs across engines or against a field
  lens

## Datasets

- [`met_example`](https://aagi-aus.github.io/flexyBayes/reference/met_example.md)
  : Example Multi-Environment Trial (MET) Dataset

## Package

- [`flexyBayes`](https://aagi-aus.github.io/flexyBayes/reference/flexyBayes-package.md)
  [`flexyBayes-package`](https://aagi-aus.github.io/flexyBayes/reference/flexyBayes-package.md)
  : flexyBayes: Bayesian Mixed Models with ASReml Syntax via INLA and
  brms
