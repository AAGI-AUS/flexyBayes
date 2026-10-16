# Bayesian Mixed Models with ASReml Syntax

Specify mixed models using ASReml formula syntax and estimate them on
one of two active engines: INLA (integrated nested Laplace
approximation) or brms (the Stan passthrough). Returns a three-part
result: a GLM-compatible object for use with standard R packages
(emmeans, marginaleffects, etc.), the native backend object for
engine-specific diagnostics, and extras for secondary analyses.

## Usage

``` r
flexybayes(
  fixed,
  random = NULL,
  residual = NULL,
  data,
  family = "gaussian",
  link = NULL,
  known_matrices = list(),
  weights = NULL,
  na_action = c("augment", "omit", "fail"),
  n_samples = 1000,
  warmup = 500,
  chains = 4,
  seed = NULL,
  control = NULL,
  prior = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  review_code = FALSE,
  backend = c("auto", "inla", "brms"),
  aggregate = "auto",
  plan = FALSE,
  syntax = c("auto", "asreml", "brms"),
  known_variance = NULL,
  fixed_alias = c("refuse", "drop"),
  rcov
)

fb(
  fixed,
  random = NULL,
  residual = NULL,
  data,
  family = "gaussian",
  link = NULL,
  known_matrices = list(),
  weights = NULL,
  na_action = c("augment", "omit", "fail"),
  n_samples = 1000,
  warmup = 500,
  chains = 4,
  seed = NULL,
  control = NULL,
  prior = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  review_code = FALSE,
  backend = c("auto", "inla", "brms"),
  aggregate = "auto",
  plan = FALSE,
  syntax = c("auto", "asreml", "brms"),
  known_variance = NULL,
  fixed_alias = c("refuse", "drop"),
  rcov
)
```

## Arguments

- fixed:

  Two-sided formula `response ~ fixed_effects`. This is the universal
  entry's model slot: it accepts the ASReml `fixed` form (paired with
  `random` / `residual`) **or** a brms / lme4-style bar-grouped formula
  such as `response ~ x + (1 | g)` (in which case the grouping lives in
  the formula and `random` / `residual` must be left `NULL`). The
  grammar is detected from the call shape. Use `syntax` to force it.

  **Offsets.** [`offset()`](https://rdrr.io/r/stats/offset.html) terms
  in the fixed formula enter the linear predictor with a fixed
  coefficient of one on both engines, and several add. Each is evaluated
  against `data` and has to be a finite number on every row
  (`fixed_offset_invalid` otherwise). brms receives the terms in its
  formula; INLA receives the evaluated offset through `inla()`'s
  `offset` argument, and a Gaussian identity-link fit reads its centre,
  scale and default priors from the response minus the offset. The
  aggregated route carries an offset exactly on the Gaussian family
  (subtracted before the cell statistics are taken) and the Poisson
  family (as exposure), and refuses one on the binomial family
  (`fixed_offset_unsupported`), where `aggregate = "auto"` takes the
  per-row route. [`predict()`](https://rdrr.io/r/stats/predict.html) on
  new data evaluates the offset there and refuses when a column it reads
  is absent (`newdata_offset_missing`). An
  [`offset()`](https://rdrr.io/r/stats/offset.html) in `random` or
  `residual` is refused.

- random:

  One-sided formula: `~ random_terms` using ASReml syntax. "Parsed" and
  "fitted" are not the same claim, so the term catalogue is two lists.

  **Parsed and fitted** (emitted by at least one active engine): `vm()`,
  `at()`, [`diag()`](https://rdrr.io/r/base/diag.html), `idh()`, `us()`,
  `ar1()`, `spl()`, `ped()`, `dsum()`, `id()`, and nested colon terms,
  including the per-trial separable AR1 field described under "Spatial
  and temporal autoregressive fields" below.

  **Parsed and refused by name.** `fa(f, k):g`, the factor-analytic
  genotype-by-environment covariance, is refused on both engines
  (`factor_analytic_not_repeatable`). On simulated trials with a known
  covariance, the genetic-correlation intervals of INLA fits covered the
  true value less often than required (134 and 125 of 150 for k = 1 and
  2, against 136), and flexyBayes has no factor-analytic lowering for
  brms. Write `diag(f):g`, which fits on both engines, or `us(f):g` on
  brms. The bare spelling `fa(v, k)`, with no grouping factor for the
  covariance to apply to, keeps its own refusal
  (`fa_not_representable`). `corh(f):g` is likewise ASReml syntax that
  flexyBayes parses and refuses by name. It is described just below.

  **Heterogeneous variances.** `diag(f):g`, `idh(f):g` and `at(f):g` all
  fit one variance per level of `f` with no covariance between levels:
  the usual multi-environment structure, where the genotype variance
  differs by site. ASReml treats
  [`diag()`](https://rdrr.io/r/base/diag.html) and `idh()` as the same
  structure, and so does flexyBayes. The three spellings emit identical
  code. `us(f):g` estimates every pairwise correlation instead, at
  `k(k+1)/2` parameters against
  [`diag()`](https://rdrr.io/r/base/diag.html)'s `k`.

  `corh(f):g`, heterogeneous variances with a single shared correlation,
  is refused on both engines: flexyBayes has no lowering for it, and
  approximating it with either neighbour would change the parameter
  count. With one variance for every level and a correlation of zero or
  more, `random = ~ g + f:g` is the same covariance and fits on both
  engines.

  **Selected levels.** `at(f, level):g` gives the random effect `g` only
  to the rows where `f` takes one of the named levels, with one variance
  per named level, as in ASReml. `at(design, "split"):mainplot` is a
  main-plot effect in the split-plot trials alone. Levels are written as
  values, `at(f, c("a", "b"))`, or as positions in `levels(f)`,
  `at(f, 2)`. A level the data do not have is refused
  (`at_level_not_in_data`). It fits on both engines and reports
  `sd_<g>__<f><level>`, the name `at(f):g` gives the same level.

  `diag(f):g` fits on both engines. `us(f):g` fits on brms, and INLA
  refuses it, so `auto` routes such a model to brms.

  **Spatial and temporal autoregressive fields.** `ar1(t)` and the
  separable `ar1(row):ar1(col)` are written here, on the random side,
  and emitted on INLA as a latent autoregressive field plus the Gaussian
  observation nugget. Four hyperparameters on a field grid: the row
  correlation, the column correlation, the field SD and the nugget SD,
  all four printed by
  [`summary()`](https://rdrr.io/r/base/summary.html). Each correlation
  is named after the factor it runs along, `rho_<factor>`: `rho_row` and
  `rho_col` for `ar1(row):ar1(col)`, and the same two names for
  `ar1(col):ar1(row)`. The field is faithful only with one observation
  per grid node, so an incomplete or replicated grid refuses (keep the
  unobserved nodes as design cells with the default
  `na_action = "augment"`). flexyBayes has no brms lowering for the
  field, so brms refuses.

  That four-parameter model is **not** ASReml's separable residual,
  which is one correlated residual with no independent plot error: three
  parameters, nested inside this one as the nugget goes to zero. On data
  with real plot-to-plot noise the two return different correlations,
  because a nugget-free model must absorb independent noise into the
  correlated structure and is pulled towards zero by it. The gap grows
  with the share of the total variance the nugget holds, and it falls on
  the correlations a reader acts on. Neither model is wrong. The
  residual form is written on the residual side,
  `residual = ~ ar1(row):ar1(col)`, described under `residual` below.

  **Per-trial nested field.** `id(trial):ar1(row):ar1(col)` fits one
  separable AR1 field per level of `trial` on INLA, via the
  `replicate =` mechanism: one complete grid per trial, each with its
  own realised field. The row correlation, column correlation and field
  SD are **shared across every level of `trial`**, not estimated per
  trial. This is the model ASReml fits under the same spelling.

  **A field for each trial.** `at(trial):ar1(row):ar1(col)` fits a
  separate field on each trial, with its own SD and correlations, as
  ASReml does: three parameters per trial, reported as
  `sd_spatial__trial<level>`, `rho_row__trial<level>` and
  `rho_col__trial<level>` for `ar1(row):ar1(col)`, plus the shared
  Gaussian nugget. `at(trial, level):ar1(row):ar1(col)` fits the field
  on the named trials only. Each trial's plots must fill its own grid.
  brms refuses both.

  **Units followed through time.** `ar1(t):g`, or the same term written
  `g:ar1(t)`, fits one AR1 series along `t` within each level of `g`,
  independent across levels, with one correlation and one field SD
  shared by every level and the Gaussian nugget on top, three
  hyperparameters in all. This is the units-by-occasions
  repeated-measures covariance of Verbyla and Cullis (1992) written on
  the random side: `g` indexes the plots, animals or series and `t` the
  occasions, whose levels must be in time order. An occasion a unit
  missed stays as a design row with the response `NA` under the default
  `na_action = "augment"`, and the `t` by `g` grid must be complete in
  that sense, one row per node. brms refuses, as for the other fields.

- residual:

  One-sided formula: `~ residual_structure`. Default `~ units` (iid
  residuals). `~ dsum(~ units | env)`, and the equivalent
  `~ at(env):units`, give a separate residual variance per level of
  `env`. Only families that carry a residual scale parameter can have
  one made heterogeneous. A Poisson's dispersion is a function of its
  mean, so there is nothing there to vary, and the request is refused.

  `~ ar1(row):ar1(col)` is ASReml's separable residual: one correlated
  residual with no nugget, three parameters (`sd_spatial`, `rho_row`,
  `rho_col`), and no `sigma` row. `~ dsum(~ ar1(row):ar1(col) | trial)`
  fits it separately in each trial, with names suffixed
  `__trial<level>`. Both fit on INLA for a Gaussian response with the
  identity link and one plot per grid cell, where the observation error
  is fixed at a variance of `exp(-12)` on the standardised response;
  brms refuses them. They cannot be combined with a random-side field.
  `~ ar1(t)` and `~ ar1(row):col` are refused on both engines, with a
  pointer to the random-side field described under `random` above.

- data:

  A data.frame containing all variables referenced in the formulas.

- family:

  Character: `"gaussian"`, `"binomial"`, `"poisson"`,
  `"negative_binomial"`, `"gamma"`, `"beta"`, `"hurdle_gamma"` or
  `"lognormal"`. The last two are brms only: INLA's likelihood roster
  has no hurdle gamma, and no INLA lognormal fit is validated, so the
  family gate refuses both there and `backend = "auto"` routes them to
  brms. A family object
  ([`binomial()`](https://rdrr.io/r/stats/family.html),
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  `brms::lognormal()`) is also accepted, goes through the same gate as
  its name, and supplies its own link, so
  [`Gamma()`](https://rdrr.io/r/stats/family.html) means the inverse
  link it names, not the log link `family = "gamma"` defaults to.
  Passing a family object together with a contradicting `link` is
  refused, not resolved. A family outside this set is refused at the
  door. Where the boundary is flexyBayes's own and not the engines', the
  refusal says so and `inst/KNOWN_ISSUES.md` records it. A response
  outside the family's domain (0 or negative for `"gamma"` or
  `"lognormal"`, negative for `"hurdle_gamma"`, at or outside \[0, 1\]
  for `"beta"`) is refused by name (`family_response_domain_violation`)
  before any engine runs. A missing (`NA`) response is `na_action`'s
  concern, not a domain violation.

- link:

  Character or NULL: override the default link function (e.g.,
  `"probit"` for binomial). Leave it unset when `family` is a family
  object, which carries its own link. The links fitted, default first:
  gaussian `identity`, `log`, `inverse`, `logit`, `probit`, `cloglog`,
  `cauchit`; binomial `logit`, `probit`, `probit_approx`, `cloglog`,
  `cauchit`, `softit`, `identity`, `log`; poisson `log`, `identity`,
  `sqrt`; negative_binomial `log`, `identity`, `sqrt`, `softplus`,
  `squareplus`; gamma `log`, `identity`, `inverse`; hurdle_gamma `log`,
  `identity`, `inverse`, `softplus`, `squareplus`; beta `logit`;
  lognormal `identity`. `backend = "inla"` takes the default link only,
  refusing any other by name (`inla_gate_refused`), and
  `backend = "auto"` fits another link on brms. A link outside the list
  is refused by name (`unsupported_family`) on every backend.

- known_matrices:

  Named list of matrices referred to in the random formula (e.g.,
  `list(Gmat = G_mat, Amat = A_mat)`). The carrier is declared with the
  [`fb_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_cov.md)
  constructor inside the random term: a dense covariance
  (`vm(group, cov = fb_cov(G, type = "dense"))`; brms), a user-supplied
  lower-triangular Cholesky factor
  (`vm(group, cov = fb_cov(L, type = "chol"))`; brms), or a sparse
  precision matrix `Matrix::dgCMatrix`
  (`vm(group, cov = fb_cov(Q, type = "precision"))`; INLA and brms).
  INLA takes the precision and pedigree-precision carriers only, and
  refuses by name a dense or Cholesky carrier, pointing at `solve(V)` or
  at brms, and a block carrier, which INLA would fit with one standard
  deviation per block, pointing at the precision carrier of the
  block-diagonal inverse. The bare dense forms `vm(group, V = ...)` /
  `ped(group, A = ...)` remain the default. The legacy v0.3.7 keyword
  carriers (`chol = `, `precision = `, `blocks = `,
  `low_rank_factor = `) are deprecated and emit a migration warning. See
  the formula-surface vignette for per-type worked examples.

- weights:

  Optional numeric weight vector (length N). Lowered for
  `family = "gaussian"` with the identity link only, in the ASReml /
  lme4 / `glm(weights=)` precision sense, `Var(y_i) = sigma^2 / w_i`:
  INLA's `scale = w` per-observation precision multiplier, and on brms a
  known offset on the log-link sigma distributional parameter (not
  brms's own [`weights()`](https://rdrr.io/r/stats/weights.html)
  addition term, which implements a different, likelihood-power
  quantity). A constant vector (including `rep(1, N)`) is the unweighted
  model under a different spelling and passes through unchanged. A
  non-constant vector on any other family, or a non-identity link on
  Gaussian, is refused by name (`weights_requires_gaussian`) instead of
  silently producing the unweighted posterior. Likelihood-power,
  frequency, trials and exposure are different models and cannot share
  this argument silently. `aggregate = TRUE` together with weights is
  also refused by name (`weights_not_aggregatable`): the aggregated
  route's closed-form sufficient statistics do not carry a
  per-observation weight through the compression.

- na_action:

  How to treat observations whose response is missing, and, through the
  list form below, observations whose predictors are missing.

  `"augment"` (default) retains missing-response rows, carrying the
  missing response as a latent quantity the engine marginalises, and
  completes the design grid where the absent cells are determinable.
  This is ASReml's `na.method(y = "include")` behaviour, and it matters
  whenever the model carries a covariance indexed by the design:
  deleting the row of a lost plot changes the index set a separable AR1
  field is built over, so the fitted model is no longer the model that
  was written down. `"omit"` drops the rows (complete-case), and a
  structured covariance over the resulting broken grid then refuses
  downstream. `"fail"` refuses if any response is missing, and leaves
  the design grid alone.

  A response value that is `NaN` or infinite, as `log(0)` gives, is not
  a missing observation. It refuses by name under every setting
  (`response_not_finite`). Set it to `NA` to treat it as missing.

  The argument also accepts the object an ASReml user already writes,
  `asreml::na.method(y = , x = )`, and the bare list of the same shape
  for readers without an asreml licence:

      flexybayes(..., na_action = list(y = "include", x = "fail"))

  The response words map `include` to `"augment"`, `omit` to `"omit"`,
  `fail` to `"fail"`. The covariate words are ASReml's own. `x = "fail"`
  (the default, and ASReml's) refuses a missing predictor. `x = "omit"`
  drops the affected rows with a warning naming the count and the
  columns. `x = "include"`, ASReml's zero-fill (Reference Manual 4.2,
  section 3.11), is refused by name, because a zero is a value the plot
  did not have. A value arriving as an unreduced `na.method()` default
  vector is read as its first element, which is the policy ASReml itself
  would use.

  Under ignorability, and only when both representations put the same
  prior on the same model parameters, the posterior for those parameters
  is the same whether a missing response is augmented or omitted –
  augmentation preserves the representation, not information. A backend
  default prior can depend on the realised design matrix, so it does not
  establish that common-prior condition automatically. The identity is a
  statement about the posterior and not a promise about what an
  optimiser returns: the two settings hand the engine an intact design
  and a ragged one, and on a weakly identified mode they can land in
  different places. Where missingness depends on the unobserved response
  itself, both are biased and neither this argument nor any diagnostic
  here will tell you so.

- n_samples:

  Integer: number of posterior samples per chain.

- warmup:

  Integer: number of warmup (burn-in) iterations per chain.

- chains:

  Integer: number of MCMC chains.

- seed:

  Integer or `NULL`: the seed for the sampler, on the engine that
  samples.

  On brms it is forwarded to `brms::brm(seed = )` and fixes Stan's
  draws: two fits of the same model to the same data under the same
  seed, chains and backend return identical draws. With `NULL` (the
  default) brms passes no seed and Stan runs with one drawn for it. On
  the cmdstanr backend that draw comes from R's random-number stream, so
  it follows [`set.seed()`](https://rdrr.io/r/base/Random.html) only
  when nothing else draws from the stream first. Either way the seed
  Stan ran with is kept at `fit$extras$reproducibility$seed$as_run`, and
  a refit under that value returns the same draws. Draws taken later
  from the fitted object, such as `brms::posterior_epred()` on new
  levels, follow R's random-number stream and not this seed.

  On INLA the seed has no effect: `INLA::inla()` takes no seed, and the
  argument says once that it was not used. A refit in a fresh R session
  can reach a different mode whatever the seed (see the section *INLA
  threads*). The thread setting an INLA fit ran with is recorded beside
  the engine versions, in the same record.

- control:

  A named list of sampler control settings or `NULL`, forwarded to
  `brms::brm(control = )`. This is the route to `adapt_delta` and
  `max_treedepth`: `control = list(adapt_delta = 0.95)` is the standard
  first response to divergent transitions. `NULL` (the default) uses
  brms's own settings. A no-op on the INLA path, as `seed` is.

- prior:

  An optional
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  object specifying priors via the PC-canonical hybrid DSL (preferred).
  When supplied it overrides `prior_vc_sd` for the variance components
  it covers. See
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md).
  On INLA, a random-effect standard deviation the prior does not name
  takes the package's default, `half_normal(scale = s)` with the
  family-aware `s` described under `prior_vc_sd` (divided by the
  covariate's standard deviation for an uncorrelated slope);
  [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
  lists what was filled. A standard deviation the prior does name is
  left exactly as written.

- prior_fixed_sd:

  Numeric: SD for fixed-effect normal priors, applied uniformly to the
  intercept, factor contrasts, continuous slopes, factor x continuous
  interactions, and [`I()`](https://rdrr.io/r/base/AsIs.html)-expression
  terms. Default `100`: weakly informative on the natural response scale
  for the vast majority of agricultural / clinical responses (covers
  responses with central tendency up to several hundred without crushing
  the posterior toward zero), while still regularising at sample sizes
  below ~ 30 per coefficient. Set wider (e.g. `1000`) for responses on a
  larger natural scale, or narrower for explicit shrinkage. A
  weakly-informative normal prior on the data scale, in the spirit of
  the weakly-informative-prior literature for regression coefficients
  (e.g. Gelman et al. 2008, *Annals of Applied Statistics*
  2(4):1360-1383). **Applied when you supply it**, on every backend:
  brms receives one `normal(0, sd)` row per fixed-effect class, and INLA
  receives `control.fixed` with `prec = 1 / sd^2` on the slopes and the
  intercept alike. Left unsupplied, brms keeps its own fixed-effect
  default (flat coefficients and a response-centred `student_t`
  intercept), because replacing a response-centred intercept prior with
  `normal(0, 100)` is a different model for any response far from zero.
  INLA gives each column built from a numeric variable the
  covariate-scaled default `normal(0, 1e5 * sd(y) / sd(x))` on an
  identity link and `normal(0, sqrt(1000) / sd(x))` on a log or logit
  link, with `sd(x)` read over the rows whose response is observed, so a
  covariate's units do not change the fit, and keeps its own
  `prec = 0.001` on factor contrasts and its flat intercept.
  [`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
  names what a fit ran under. A supplied value is never rescaled by a
  covariate: it is a standard deviation in response units per covariate
  unit, so on a covariate recorded in small units (a large coefficient)
  it can bind. Measure the covariate in natural units, or leave the
  argument unsupplied on INLA.

- prior_vc_sd:

  Numeric: hyperparameter for the legacy `lognormal(0, prior_vc_sd)`
  variance-component prior. **Note:** when both `prior` and
  `prior_vc_sd` are left at their defaults, the package applies a
  half-normal default on the SD scale (`half_normal(scale = s)` with
  family-aware `s`: `5 * sd(y)` for Gaussian, `5` on the logit scale for
  binomial / beta, `3` on the log scale for Poisson / negative-binomial
  / gamma) to `sigma` and every named random-effect group, on both
  engines, with the scale on the standard deviation of an uncorrelated
  random slope divided by the standard deviation of its numeric
  covariate. It is weakly informative and has no upper bound. A half-t
  or half-Cauchy (Gelman 2006, *Bayesian Analysis* 1(3):515-534) is
  stated through
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md).
  Before 0.13.0 the default was `uniform(0, U)` with `U` equal to the
  scale above. The legacy `lognormal(0, 1)` default fires only when
  `prior_vc_sd` is passed explicitly. Silence the one-time announcement
  message via `options(flexyBayes.silence_default_prior_note = TRUE)`.
  When it is passed explicitly it reaches both active engines as the
  same density on the standard-deviation scale: brms as a
  `lognormal(0, prior_vc_sd)` prior row, INLA as the expression prior
  that writes that density in INLA's log-precision parameterisation. The
  two engines therefore carry the same variance-component prior.

- verbose:

  Logical: print the generated backend code to the console.

- mcmc_verbose:

  Logical: show the sampler's progress bar. Has no effect on INLA, which
  does not sample.

- return_code:

  Logical: if TRUE, return the generated backend code without fitting
  the model: the Stan program on `backend = "brms"`, the INLA formula,
  family and hyperparameter list on `backend = "inla"`.

- review_code:

  Logical: if TRUE, do not fit the model immediately; instead return a
  `<flexybayes_review>` deferred-execution object carrying the generated
  Stan code, the resolved prior, the parsed intermediate representation
  (IR), the captured call, and a snapshot of `.Random.seed`. Inspect the
  code with
  [`cat_code()`](https://aagi-aus.github.io/flexyBayes/reference/cat_code.md)
  and run the deferred fit with
  [`proceed()`](https://aagi-aus.github.io/flexyBayes/reference/proceed.md).
  A second
  [`proceed()`](https://aagi-aus.github.io/flexyBayes/reference/proceed.md)
  call returns the cached fit. Useful as a teaching / auditing surface
  before a long MCMC run. Default `FALSE` preserves the existing
  run-immediately semantics. A session-level override is available via
  `options(flexyBayes.review_code_default = TRUE)`. The argument value
  at call time always wins. Closest published precedent:
  `brms::stancode()` plus the `chains = 0` idiom for "do everything
  except sampling". `review_code = TRUE` and `return_code = TRUE` are
  mutually exclusive. brms is the only engine with a code slot, so
  `review_code = TRUE` is available under `backend = "brms"` and under
  `"auto"` (which resolves the code modes to brms). Under
  `backend = "inla"` it refuses with `review_code_backend_unsupported`:
  the deferred-execution token would need an INLA-side code slot, which
  is queued for a later release.

- backend:

  Character: one of `"auto"` (**default**), `"inla"`, or `"brms"`. INLA
  fits the latent-Gaussian model class by nested Laplace approximation,
  without sampling. brms is the Stan passthrough. Any other value,
  including a formerly-registered engine name this package withdrew
  entirely (see `NEWS.md`), refuses with `unknown_backend` naming the
  two active engines.

  **What `"auto"` does.** The call runs `lgm_gate()` and routes to INLA
  when the model is latent-Gaussian feasible and INLA is installed.
  Otherwise it routes to brms when brms is installed *and* its
  capability predicate accepts the model. When neither can represent the
  model, the call refuses with `auto_no_active_route` instead of fitting
  something else. Nothing is translated silently. `auto` does change
  engine in one case after starting: an INLA program failure is caught
  and re-routed to brms, and a message reports it. The fit's class and
  `backend` field always record which engine produced it.
  [`backend_decision()`](https://aagi-aus.github.io/flexyBayes/reference/backend_decision.md)
  surfaces the full dispatch trace (including `rejected_routes`)
  post-fit, and
  [`fb_backend_status()`](https://aagi-aus.github.io/flexyBayes/reference/fb_backend_status.md)
  reports which engines are usable. The one-time notes on the auto path
  are silenceable via
  `options(flexyBayes.silence_auto_fallback_note = TRUE)` (gate refusal
  or an INLA numerical failure) and
  `options(flexyBayes.silence_auto_inla_missing_note = TRUE)` (INLA not
  installed).

  **What each engine represents.** INLA fits: Gaussian LMM, simple
  random intercept; GLMM (binomial, Poisson, negative binomial, gamma,
  beta), simple random effect; Uncorrelated random slope;
  Factor-by-numeric fixed interaction; Nested / interaction random
  effects, multi-stratum; Heterogeneous variance by factor level;
  Heterogeneous residual by factor level; Baseline heteroscedastic MET:
  iid genotype, iid genotype-by-environment and a sectioned residual;
  Known-covariance genomic / pedigree random effect; Separable AR1
  spatial field; Per-trial separable AR1 field; Grouped temporal AR1
  field; Univariate P-spline; Observation weights (Gaussian, identity
  link); Exact sufficient-statistic aggregation (binomial, Poisson).
  INLA refuses: Hurdle gamma (zero mass plus a positive gamma part);
  Correlated random slope; Unstructured genotype-by-environment
  covariance; Heterogeneous variances with one shared correlation;
  Factor-analytic genotype-by-environment covariance; Multi-trait
  covariance; Exact sufficient-statistic aggregation (Gaussian). brms
  fits: Gaussian LMM, simple random intercept; GLMM (binomial, Poisson,
  negative binomial, gamma, beta), simple random effect; Hurdle gamma
  (zero mass plus a positive gamma part); Uncorrelated random slope;
  Factor-by-numeric fixed interaction; Correlated random slope; Nested /
  interaction random effects, multi-stratum; Heterogeneous variance by
  factor level; Unstructured genotype-by-environment covariance;
  Heterogeneous residual by factor level; Baseline heteroscedastic MET:
  iid genotype, iid genotype-by-environment and a sectioned residual;
  Known-covariance genomic / pedigree random effect; Observation weights
  (Gaussian, identity link). brms refuses: Heterogeneous variances with
  one shared correlation; Factor-analytic genotype-by-environment
  covariance; Multi-trait covariance; Separable AR1 spatial field;
  Per-trial separable AR1 field; Grouped temporal AR1 field; Univariate
  P-spline. The per-class table is in `README.md`, and with a note on
  each row in `system.file("KNOWN_ISSUES.md", package = "flexyBayes")`.

  **Caution (small-group random effects).** Both engines apply the same
  default prior, the half-normal on each standard deviation, which the
  INLA path represents exactly as an expression prior on the
  log-precision, so they agree on the variance component far more
  closely than in earlier versions. A model with very few groups
  nonetheless carries a weakly-identified variance component (the data
  say little about the between-group spread), and INLA's Laplace
  approximation is less accurate there than full MCMC. For a flagship
  random-intercept model with few groups, pin `backend = "brms"` or
  supply an explicit informative prior (a half-Cauchy, per Gelman 2006,
  as in the *priors and regularisation* vignette).

  **Convergence.** The brms fit emits a warning when the sampler may not
  have converged (a parameter with Rhat at or above 1.1, or a low
  effective sample size). Treat such a posterior with caution: increase
  `warmup` / `n_samples`, simplify the model, or supply a more
  informative prior, and inspect the full diagnostics with
  [`summary()`](https://rdrr.io/r/base/summary.html). Silence the
  warning (for intentionally short fits) via
  `options(flexyBayes.silence_convergence_warning = TRUE)`. An INLA fit
  does not sample and carries no such warning.

  **Code inspection.** brms is the only engine with a code slot, so
  `return_code = TRUE` / `review_code = TRUE` under `"auto"` resolve to
  brms, but only when brms can represent the model. Where it cannot, the
  call refuses with `auto_no_active_route` instead of returning the code
  for a different model. `backend = "inla"` with `return_code = TRUE`
  returns the INLA formula, family and hyperparameter list. With
  `review_code = TRUE` it refuses until the INLA-side code slot lands.

- aggregate:

  One of `"auto"` (default), `TRUE`, or `FALSE`. Exact
  sufficient-statistics aggregation gate. INLA is the only engine with
  an aggregated emit. `"auto"` consults the aggregation plan and routes
  a binomial-logit or poisson-log model through the per-cell path when
  the IR is in scope (fixed plus random intercept; productive
  compression) and the model reaches INLA. A Gaussian model is fitted
  per row: the aggregated Gaussian route cannot yet carry the residual
  prior, so `TRUE` on a Gaussian model refuses by name
  (`aggregation_residual_prior_not_composed`). `TRUE` also raises a
  structured refusal when the plan declares ineligibility, or when no
  active engine offers an aggregated route. `FALSE` skips the gate
  entirely. Aggregated fits carry `$exactness == "aggregated_exact"`,
  the dispatch trace's `path` slot reads `"aggregated_count"`, and their
  summary reports no marginal likelihood
  (`aggregated_evidence_unavailable`). For out-of-core datasets that do
  not fit in memory, see
  [`flexybayes_stream()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes_stream.md).

- plan:

  Logical: if `TRUE`, short-circuit after intermediate representation
  (IR) build and return a `<fb_plan>` object carrying the IR, the
  routing decision, the representation plan, the aggregation plan, and
  the prediction plan, *without* emitting backend code or fitting.
  Equivalent to calling
  [`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md)
  with the same formula triple, exposed inline as
  `flexybayes(plan = TRUE)` so users can reach the planning object
  without re-typing the call. Default `FALSE` preserves the run-
  immediately semantics. Mutually exclusive with `return_code` and
  `review_code`.

- syntax:

  One of `"auto"` (default), `"asreml"`, or `"brms"`. Selects how
  `fixed` is interpreted. `"auto"` detects the grammar from the call
  shape (a bar-grouped formula is read as brms, otherwise ASReml). The
  other values force a grammar.

- known_variance:

  Optional numeric vector of length N, or an N x N covariance matrix
  (below). A vector is a known variance for each response, added to the
  estimated residual variance, `Var(y_i) = sigma^2 + v_i`. This is the
  second stage of a two-stage analysis, where each response is a
  first-stage estimate (a variety mean in one trial, say) and `v_i` is
  the variance of that estimate. `weights` is a different model,
  `Var(y_i) = sigma^2 / w_i`. Given both, the residual variance is
  `sigma^2 / w_i + v_i`. Lowered for `family = "gaussian"` with the
  identity link: brms fits `y | se(sqrt(v), sigma = TRUE)`, and INLA a
  row effect with known variance `v_i`. The summaries report `sigma` as
  the estimated part, and fitted values on both engines leave the known
  part out. A vector of zeros fits the model without it. Any other
  family refuses with `known_variance_not_supported`, and a vector that
  is not one finite, non-negative value per row with
  `known_variance_invalid`. A symmetric positive-definite N x N matrix
  `K` is a known covariance, `Var(y) = sigma^2 I + K`, as when the first
  stage returns the full covariance of its estimates. It is fitted on
  INLA as a latent effect with precision `K^-1`. brms has no term for it
  and refuses with `known_variance_not_supported`. A row and column of
  zeros marks a row without a known part. A matrix that is not square,
  symmetric and positive definite over the rows of `data` refuses with
  `known_variance_invalid`.

- fixed_alias:

  `"refuse"` (the default) or `"drop"`: what to do with a fixed-effect
  column that is a linear combination of the columns before it on the
  rows with a response, as in a factorial with empty cells. INLA and
  brms would report such a column at its prior, so by default the fit
  refuses with `fixed_effects_aliased` and names the columns. `"drop"`
  drops them as asreml 4.2 does, on INLA: terms are taken in model
  order, and within a term by descending factor levels, the term's first
  factor varying slowest. A column that is a linear combination of the
  columns before it is fixed at zero.
  [`summary()`](https://rdrr.io/r/base/summary.html) reports it with
  estimate 0 and no standard error, the columns are listed in
  `fit$extras$fixed_alias$dropped`, and a prediction whose row is not a
  combination of the fitting rows is `NA` with status `"Aliased"` (see
  [`predict.flexybayes_inla()`](https://aagi-aus.github.io/flexyBayes/reference/predict.flexybayes_inla.md)).

- rcov:

  Defunct in flexyBayes 0.9.0. This was the ASReml 3 name for the
  residual-structure argument. ASReml 4 renamed it to `residual`, and
  flexyBayes follows the ASReml 4 name. Supplying `rcov` now raises an
  error. Use `residual` instead.

## Value

An object of class `"flexybayes"`, a list with three components.

- `$glm`:

  A GLM-compatible object (class `c("flexybayes_glm", "glm", "lm")`)
  with posterior mean coefficients, vcov, residuals, fitted values, etc.
  Works with [`summary()`](https://rdrr.io/r/base/summary.html),
  `emmeans()` and `marginaleffects()`.

- `$brms` or `$inla`:

  The native backend object: a live `brmsfit` on the brms path, or
  INLA's own fitted object on the INLA path. Use with `bayesplot`,
  `posterior::as_draws()`, and each engine's own accessors.

- `$extras`:

  Additional outputs: posterior `summary`, `convergence` diagnostics,
  `variance_comps`, `blups`, `predictions`, generated `code`,
  `param_names`, `parse_info`, `call_info`, `run_time`, `model_info`,
  and `reproducibility`: the R, flexyBayes and engine versions, the
  thread settings and the seed the fit ran with, which
  [`print()`](https://rdrr.io/r/base/print.html) shows.

If `return_code = TRUE`, returns the generated backend code instead: a
character string holding the Stan program on the brms path, or a list of
the formula, family and hyperparameter specification on the INLA path.

## INLA threads

An INLA fit runs on one thread, `num.threads = "1:1"`, whatever
`INLA::inla.setOption()` holds. With more threads, identical calls can
return different variance components where the data inform a variance
only weakly. On a 48-trial barley model with 22 localities, eight
identical fits at `"4:1"` gave the locality variance from 0.58 to 42.3,
and at `"1:1"` from 0.75 to 0.79, with the fixed effects unchanged. Set
`options(flexyBayes.inla_threads = "4:1")` for another plan, or `"inla"`
for INLA's own default. The plan a fit ran with is in
`$extras$reproducibility`.

One thread does not guarantee the same result in a new R session.
Identical calls in fresh R processes can reach different modes. On a
48-trial barley series, four processes reached two modes, 0.10 apart in
log marginal likelihood, whose variance component standard deviations
differed by at most 3.4%. On six trials of the same series, the modes
differed by up to 26%.

Each INLA fit is therefore run three times: once in this session and
twice in fresh R processes started through callr. The fit with the
highest log marginal likelihood is returned (the first of equal ones).
`options(flexyBayes.inla_best_of = k)` sets the number of fits, a whole
number of one or more; `1` fits once in this session, as 0.16.0 did, and
any other value refuses with `inla_best_of_invalid`. Above 1 without
callr installed, the fit refuses with `isolation_requires_callr`. A fit
that fails in a separate process stops the call with that process's
message. `$extras$reproducibility$best_of` holds each fit's log marginal
likelihood and mode status, the fit kept, the range of the marginal
likelihoods, and each standard deviation per fit with its largest
difference across the fits, relative to the smallest.
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) add a line when a
standard deviation differs by more than 5%. The restart probe of the
numerical-confirm gate examines the fit kept. The fresh-process fits run
beside the session fit, at most
`options(flexyBayes.inla_best_of_workers)` fits at once (default: all of
them) and fewer, with a message, when the free memory holds fewer. The
record gives each fit's time in seconds (`wall_s`) and whether they ran
at once (`parallel`).

## Fitting in a separate process

`options(flexyBayes.isolation = "process")` runs the engine fit in a
separate R process started through callr, which is in Suggests. Without
callr the request refuses with `isolation_requires_callr` before
anything is fitted. The option applies to `flexybayes()`,
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md),
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
and
[`proceed()`](https://aagi-aus.github.io/flexyBayes/reference/proceed.md),
and
[`fb_isolated_lapply()`](https://aagi-aus.github.io/flexyBayes/reference/fb_isolated_lapply.md)
runs any loop of calls the same way. The default, `"none"`, fits in this
session.

If the operating system kills that process, or compiled code crashes
inside it, the call refuses with `isolated_process_died`, naming the
engine, the function called, and the exit status or signal, and this
session carries on. A failure inside the fit, a refusal for instance,
arrives exactly as it would without the option. The call is not retried.

The separate process loads the flexyBayes code this session loaded and
starts from this session's options, INLA settings, environment variables
and random-number state, so it fits the same model with the same seed,
and the random-number stream moves on as it would have. Messages,
warnings and printed output come back in the order they were raised, and
the fit records the process in `$extras$isolation`. A brms fit with a
`seed` returns the same draws either way. INLA does not return identical
numbers for two identical calls in one session, and an isolated INLA fit
agrees with an in-process one as closely as two in-process fits agree.
The separate process starts with an empty global environment, so the
formula finds its variables in `data`. Starting the process and moving
the data and the fit costs time on every call.

## Deprecated

`fb()` is deprecated in flexyBayes 0.13.0 and will be removed in a later
release. Call `flexybayes()`, which takes the same arguments. A call to
`fb()` warns with a condition of class `flexybayes_deprecated`.

## Examples

``` r
if (FALSE) { # \dontrun{
# live fit, needs an active backend (INLA or brms/Stan)
data(met_example, package = "flexyBayes")
# Simple random intercept model (small budget for example purposes)
fit <- flexybayes(
  fixed  = yield ~ env,
  random = ~ geno,
  data   = met_example$dat,
  n_samples = 100, warmup = 100, chains = 1, verbose = FALSE
)
summary(fit)
coef(fit)
} # }
```
