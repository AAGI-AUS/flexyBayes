#' Bayesian Mixed Models with ASReml Syntax
#'
#' Specify mixed models using ASReml formula syntax and estimate them on one
#' of two active engines: INLA (integrated nested Laplace approximation) or
#' brms (the Stan passthrough). Returns a three-part result: a
#' GLM-compatible object for use with standard R packages (emmeans,
#' marginaleffects, etc.), the native backend object for engine-specific
#' diagnostics, and extras for secondary analyses.
#'
#' @param fixed Two-sided formula `response ~ fixed_effects`. This is the
#'   universal entry's model slot: it accepts the ASReml `fixed` form
#'   (paired with `random` / `residual`) **or** a brms / lme4-style
#'   bar-grouped formula such as `response ~ x + (1 | g)` (in which case
#'   the grouping lives in the formula and `random` / `residual` must be
#'   left `NULL`). The grammar is detected from the call shape; use
#'   `syntax` to force it.
#' @param random One-sided formula: `~ random_terms` using ASReml syntax.
#'   "Parsed" and "fitted" are not the same claim, so the term catalogue is
#'   two lists.
#'
#'   **Parsed and fitted** (emitted by at least one active engine): `vm()`,
#'   `at()`, `diag()`, `idh()`, `us()`, `ar1()`, `spl()`, `ped()`, `dsum()`,
#'   `id()`, and nested colon terms, including the per-trial separable AR1
#'   field described under "Spatial and temporal autoregressive fields"
#'   below.
#'
#'   **Parsed and refused by name** (recognised syntax, no active engine has
#'   a lowering for it): `fa()`, the factor-analytic genotype-by-environment
#'   covariance, is refused on both engines (`fa_not_representable`) -- it
#'   parses into the formula catalogue and stops there. `corh(f):g` and
#'   `at(f, level):g` are likewise real syntax with no active emit rather
#'   than unsupported syntax; both are described in full just below.
#'
#'   **Heterogeneous variances.** `diag(f):g`, `idh(f):g` and `at(f):g` all fit
#'   one variance per level of `f` with no covariance between levels -- the
#'   usual multi-environment structure, where the genotype variance differs by
#'   site. ASReml treats `diag()` and `idh()` as the same structure and so does
#'   flexyBayes; the three spellings emit identical code. `us(f):g` estimates
#'   every pairwise correlation instead, at `k(k+1)/2` parameters against
#'   `diag()`'s `k`.
#'
#'   `corh(f):g` -- heterogeneous variances with a single shared correlation --
#'   is refused: no active backend has an equicorrelation group-level
#'   structure, and approximating it with either neighbour would change the
#'   parameter count. `at(f, level):g` is also refused, because conditioning
#'   on selected levels is a different model from varying the variance across
#'   all of them.
#'
#'   These structures fit on brms. INLA refuses them structurally, so `auto`
#'   routes such a model to brms.
#'
#'   **Spatial and temporal autoregressive fields.** `ar1(t)` and the
#'   separable `ar1(row):ar1(col)` are written here, on the random side, and
#'   emitted on INLA as a latent autoregressive field plus the Gaussian
#'   observation nugget. Four hyperparameters on a field grid: the row
#'   correlation, the column correlation, the field SD and the nugget SD, all
#'   four printed by `summary()`. The field is faithful only with one
#'   observation per grid node, so an incomplete or replicated grid refuses
#'   (keep the unobserved nodes as design cells with the default
#'   `na_action = "augment"`). brms has no lowering for a Kronecker
#'   autoregressive precision and refuses.
#'
#'   That four-parameter model is **not** ASReml's separable residual, which
#'   is one correlated residual with no independent plot error -- three
#'   parameters, nested inside this one as the nugget goes to zero. On data
#'   with real plot-to-plot noise the two return different correlations,
#'   because a nugget-free model must absorb independent noise into the
#'   correlated structure and is pulled towards zero by it. The gap grows
#'   with the share of the total variance the nugget holds, and it falls on
#'   the correlations a reader acts on. Neither model is wrong. They are
#'   different models, and writing the field on the random side is what
#'   keeps them distinguishable by name.
#'
#'   **Per-trial nested field.** `at(trial):ar1(row):ar1(col)`, with no
#'   `level` argument on `at()`, fits one separable AR1 field per level of
#'   `trial` on INLA, via the `replicate =` mechanism -- one complete grid
#'   per trial, each with its own realised field, but the row correlation,
#'   column correlation and field SD are **shared across every level of
#'   `trial`**, not estimated per trial. `at(trial, level):ar1(row):ar1(col)`
#'   (a `level` argument) asks for either a single conditioned trial or
#'   per-trial (unshared) hyperparameters, and both are refused
#'   (`at_field_per_level_hyper_not_representable`): there is no lowering
#'   for a level-conditioned field or for hyperparameters that vary by
#'   level. brms has no lowering for this field at all and refuses.
#' @param residual One-sided formula: `~ residual_structure`. Default `~ units`
#'   (iid residuals). `~ dsum(~ units | env)`, and the equivalent
#'   `~ at(env):units`, give a separate residual variance per level of `env`.
#'   Only families that carry a residual scale parameter can have one made
#'   heterogeneous -- a Poisson's dispersion is a function of its mean, so
#'   there is nothing there to vary, and the request is refused.
#'
#'   `~ ar1(row):ar1(col)` and `~ ar1(t)` are **refused** here: they name
#'   ASReml's nugget-free separable residual, which neither active engine
#'   fits. The refusal points at the random-side field described under
#'   `random` above, and at ASReml for the residual formulation itself.
#' @param rcov Defunct in flexyBayes 0.9.0. This was the ASReml 3 name for
#'   the residual-structure argument; ASReml 4 renamed it to `residual` and
#'   flexyBayes follows the ASReml 4 name. Supplying `rcov` now raises an
#'   error -- use `residual` instead.
#' @param data A data.frame containing all variables referenced in the formulas.
#' @param family Character: `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"negative_binomial"`, `"gamma"`, `"beta"`, or `"hurdle_gamma"`
#'   (brms only -- INLA's likelihood roster has no counterpart, so the
#'   family gate refuses it there and `backend = "auto"` routes to brms).
#'   A `stats` family object (`binomial()`, `gaussian()`) is also
#'   accepted and supplies its own link, so `Gamma()` means the inverse
#'   link it names rather than the log link `family = "gamma"` defaults
#'   to. Passing a family object together with a contradicting `link` is
#'   refused rather than resolved. A family outside this set is refused
#'   at the door; where the boundary is flexyBayes's own rather than the
#'   engines', the refusal says so and `inst/KNOWN_ISSUES.md` records it.
#' @param link Character or NULL: override the default link function
#'   (e.g., `"probit"` for binomial). Leave it unset when `family` is a
#'   family object, which carries its own link.
#' @param known_matrices Named list of matrices referred to in the random
#'   formula (e.g., `list(Gmat = G_mat, Amat = A_mat)`). The carrier is
#'   declared with the [fb_cov()] constructor inside the random term:
#'   a dense covariance (`vm(group, cov = fb_cov(G, type = "dense"))`;
#'   brms), a user-supplied lower-triangular Cholesky factor
#'   (`vm(group, cov = fb_cov(L, type = "chol"))`; brms),
#'   or a sparse precision matrix `Matrix::dgCMatrix`
#'   (`vm(group, cov = fb_cov(Q, type = "precision"))`; INLA and brms).
#'   INLA takes the precision, pedigree-precision and block carriers only,
#'   and refuses a dense or Cholesky carrier by name, pointing at
#'   `solve(V)` or at brms. The bare dense forms `vm(group, V = ...)` /
#'   `ped(group, A = ...)` remain the default; the legacy v0.3.7
#'   keyword carriers (`chol = `, `precision = `, `blocks = `,
#'   `low_rank_factor = `) are deprecated and emit a migration warning.
#'   See the formula-surface vignette for per-type worked examples.
#' @param weights Optional numeric weight vector (length N). Lowered for
#'   `family = "gaussian"` with the identity link only, in the ASReml /
#'   lme4 / `glm(weights=)` precision sense, `Var(y_i) = sigma^2 / w_i`:
#'   INLA's `scale = w` per-observation precision multiplier, and on
#'   brms a known offset on the log-link sigma distributional parameter
#'   (not brms's own `weights()` addition term, which implements a
#'   different, likelihood-power quantity). A constant vector (including
#'   `rep(1, N)`) is the unweighted model under a different spelling and
#'   passes through unchanged. A non-constant vector on any other family,
#'   or a non-identity link on Gaussian, is refused by name
#'   (`weights_requires_gaussian`) rather than silently producing the
#'   unweighted posterior -- likelihood-power, frequency, trials, and
#'   exposure are different models and cannot share this argument
#'   silently. `aggregate = TRUE` together with weights is also refused
#'   by name (`weights_not_aggregatable`): the aggregated route's
#'   closed-form sufficient statistics do not carry a per-observation
#'   weight through the compression.
#' @param na_action How to treat observations whose response is missing,
#'   and -- through the list form below -- observations whose predictors
#'   are missing.
#'
#'   `"augment"` (default) retains missing-response rows, carrying the
#'   missing response as a latent quantity the engine marginalises, and
#'   completes the design grid where the absent cells are determinable.
#'   This is ASReml's `na.method(y = "include")` behaviour, and it
#'   matters whenever the model carries a covariance indexed by the
#'   design: deleting the row of a lost plot changes the index set a
#'   separable AR1 field is built over, so the fitted model is no longer
#'   the model that was written down. `"omit"` drops the rows
#'   (complete-case); a structured covariance over the resulting broken
#'   grid then refuses downstream. `"fail"` refuses if any response is
#'   missing, and leaves the design grid alone.
#'
#'   The argument also accepts the object an ASReml user already writes,
#'   `asreml::na.method(y = , x = )`, and the bare list of the same shape
#'   for readers without an asreml licence:
#'
#'   ```
#'   flexybayes(..., na_action = list(y = "include", x = "fail"))
#'   ```
#'
#'   The response words map `include` to `"augment"`, `omit` to
#'   `"omit"`, `fail` to `"fail"`. The covariate words are ASReml's own:
#'   `x = "fail"` (the default, and ASReml's) refuses a missing
#'   predictor; `x = "omit"` drops the affected rows with a warning
#'   naming the count and the columns; `x = "include"` -- ASReml's
#'   zero-fill (Reference Manual 4.2, section 3.11) -- is refused by
#'   name, because a zero is a value the plot did not have. A value
#'   arriving as an unreduced `na.method()` default vector is read as
#'   its first element, which is the policy ASReml itself would use.
#'
#'   Under ignorability the posterior for the model parameters is the
#'   same whether a missing response is augmented or omitted --
#'   augmentation preserves the representation, not information. Where
#'   missingness depends on the unobserved response itself, both are
#'   biased and neither this argument nor any diagnostic here will tell
#'   you so.
#' @param n_samples Integer: number of posterior samples per chain.
#' @param warmup Integer: number of warmup (burn-in) iterations per chain.
#' @param chains Integer: number of MCMC chains.
#' @param seed Integer or `NULL`: the sampler's random seed, forwarded to
#'   `brms::brm(seed = )`. Two brms fits of the same model to the same data
#'   under the same seed return identical draws. `set.seed()` before the
#'   call does not achieve this -- Stan draws from its own stream -- so a
#'   posterior quoted to more than about two significant figures needs this
#'   argument. `NULL` (the default) leaves the seed to brms and the run is
#'   not reproducible. The INLA path is a deterministic Laplace
#'   approximation and has no random stream, so the argument is a no-op
#'   there and says so once.
#' @param control A named list of sampler control settings or `NULL`,
#'   forwarded to `brms::brm(control = )`. This is the route to
#'   `adapt_delta` and `max_treedepth`: `control = list(adapt_delta = 0.95)`
#'   is the standard first response to divergent transitions. `NULL` (the
#'   default) uses brms's own settings. A no-op on the INLA path, as
#'   `seed` is.
#' @param prior An optional `fb_prior()` object specifying priors via the
#'   PC-canonical hybrid DSL (preferred). When supplied it overrides
#'   `prior_vc_sd` for the variance components it covers. See
#'   [fb_prior()].
#' @param prior_fixed_sd Numeric: SD for fixed-effect normal priors,
#'   applied uniformly to the intercept, factor contrasts, continuous
#'   slopes, factor x continuous interactions, and `I()`-expression
#'   terms. Default `100` -- weakly informative on the natural response
#'   scale for the vast majority of agricultural / clinical responses
#'   (covers responses with central tendency up to several hundred
#'   without crushing the posterior toward zero), while still
#'   regularising at sample sizes below ~ 30 per coefficient. Set
#'   wider (e.g. `1000`) for responses on a larger natural scale, or
#'   narrower for explicit shrinkage. A weakly-informative normal prior
#'   on the data scale, in the spirit of the weakly-informative-prior
#'   literature for regression coefficients (e.g. Gelman et al. 2008,
#'   *Annals of Applied Statistics* 2(4):1360-1383). **Applied when you
#'   supply it**, on every backend: brms receives one `normal(0, sd)`
#'   row per fixed-effect class, and INLA receives `control.fixed` with
#'   `prec = 1 / sd^2` on the slopes and the intercept alike. Left
#'   unsupplied, each engine keeps its own
#'   fixed-effect default -- brms's flat coefficients and
#'   response-centred `student_t` intercept, INLA's `prec = 0.001` and
#'   flat intercept -- because replacing a response-centred intercept
#'   prior with `normal(0, 100)` is a different model for any response
#'   far from zero. [prior_summary()] names which of the two a fit ran
#'   under.
#' @param prior_vc_sd Numeric: hyperparameter for the legacy
#'   `lognormal(0, prior_vc_sd)` variance-component prior. **Note:**
#'   when both `prior` and `prior_vc_sd` are left at their defaults,
#'   v0.1 activates a bounded-uniform default on the SD scale
#'   (`uniform(0, U)` with family-aware `U`: `5 * sd(y)` for Gaussian;
#'   `5` on the logit scale for binomial / beta; `3` on the log scale
#'   for Poisson / negative-binomial / gamma) for `sigma` and every
#'   named random-effect group -- a weakly-informative choice for
#'   moderate group counts; for very small `J`, Gelman (2006), *Bayesian
#'   Analysis* 1(3):515-534, recommends a half-t / half-Cauchy instead
#'   (see the priors vignette). The legacy `lognormal(0, 1)` default
#'   fires only when `prior_vc_sd` is passed explicitly. Silence the
#'   one-time announcement message via
#'   `options(flexyBayes.silence_default_prior_note = TRUE)`. When it is
#'   passed explicitly it reaches both active engines as the same
#'   density on the standard-deviation scale: brms as a
#'   `lognormal(0, prior_vc_sd)` prior row, INLA as the expression prior
#'   that writes that density in INLA's log-precision parameterisation.
#'   The two engines therefore carry the same variance-component prior,
#'   which is what [triangulate()] requires before it compares them.
#' @param verbose Logical: print the generated backend code to the console.
#' @param mcmc_verbose Logical: show the sampler's progress bar. Has no
#'   effect on the deterministic INLA path.
#' @param return_code Logical: if TRUE, return the generated backend code
#'   without fitting the model -- the Stan program on `backend = "brms"`,
#'   the INLA formula, family and hyperparameter list on
#'   `backend = "inla"`.
#' @param review_code Logical: if TRUE, do not fit the model immediately;
#'   instead return a `<flexybayes_review>` deferred-execution object
#'   carrying the generated Stan code, the resolved prior, the
#'   parsed intermediate representation (IR), the captured call, and
#'   a snapshot of `.Random.seed`. Inspect the code with [cat_code()];
#'   run the deferred fit with [proceed()]; a second [proceed()] call
#'   returns the cached fit. Useful as a teaching / auditing surface
#'   before a long MCMC run. Default `FALSE` preserves the existing
#'   run-immediately semantics. A session-level override is available
#'   via `options(flexyBayes.review_code_default = TRUE)`; the
#'   argument value at call time always wins. Closest published
#'   precedent: [brms::make_stancode()] plus the `chains = 0` idiom
#'   for "do everything except sampling". `review_code = TRUE` and
#'   `return_code = TRUE` are mutually exclusive. brms is the only engine
#'   with a code slot, so `review_code = TRUE` is available under
#'   `backend = "brms"` and under `"auto"` (which resolves the code modes
#'   to brms). Under `backend = "inla"` it refuses with
#'   `review_code_backend_unsupported`: the deferred-execution token would
#'   need an INLA-side code slot, which is queued for a later release.
#' @param plan Logical: if `TRUE`, short-circuit after intermediate
#'   representation (IR) build and return a `<fb_plan>` object
#'   carrying the IR, the routing decision, the representation plan,
#'   the aggregation plan, and the prediction plan, *without* emitting
#'   backend code or fitting.  Equivalent to calling [fb_plan()] with
#'   the same formula triple --- exposed inline as `flexybayes(plan =
#'   TRUE)` so users can reach the planning object without re-typing
#'   the call.  Default `FALSE` preserves the run-
#'   immediately semantics.  Mutually exclusive with `return_code` and
#'   `review_code`.
#' @param backend Character: one of `"auto"` (**default**), `"inla"`, or
#'   `"brms"`. INLA is the deterministic Laplace path over the
#'   latent-Gaussian model class; brms is the Stan passthrough. Any other
#'   value -- including a formerly-registered engine name this package
#'   withdrew entirely (see `NEWS.md`) -- refuses with `unknown_backend`
#'   naming the two active engines.
#'
#'   **What `"auto"` does.** The call runs `lgm_gate()` and routes to INLA
#'   when the model is latent-Gaussian feasible and INLA is installed.
#'   Otherwise it routes to brms when brms is installed *and* its
#'   capability predicate accepts the model. When neither can represent
#'   the model, the call refuses with `auto_no_active_route` rather than
#'   fitting something else. There is no silent fallback.
#'   [backend_decision()] surfaces the full dispatch trace (including
#'   `rejected_routes`) post-fit, and [fb_backend_status()] reports which
#'   engines are usable. The one-time notes on the auto path are
#'   silenceable via `options(flexyBayes.silence_auto_fallback_note =
#'   TRUE)` (gate refusal or an INLA numerical failure) and
#'   `options(flexyBayes.silence_auto_inla_missing_note = TRUE)` (INLA not
#'   installed).
#'
#'   **What each engine represents.** INLA takes simple iid random
#'   effects, P-splines, the sparse-precision and pedigree carriers of
#'   `vm()` / `ped()`, and separable AR1 fields; it refuses
#'   heterogeneous variances (`diag` / `idh` / `at`), unstructured `us()`,
#'   interaction random effects, and heterogeneous residuals. brms takes
#'   interaction random effects `(1 | a:b)`, uncorrelated random slopes
#'   `(x || g)`, heterogeneous genotype variances (`diag` / `idh` / `at`),
#'   unstructured `us(f):g`, per-level residual variances
#'   (`dsum(~ units | f)`), and the dense / Cholesky / precision carriers
#'   of `vm()` / `ped()`; it refuses correlated random slopes `(x | g)`,
#'   AR1 residual structures, `corh()`, factor-analytic GxE, and splines.
#'   The generated per-class table is in `README.md` and
#'   `system.file("KNOWN_ISSUES.md", package = "flexyBayes")`.
#'
#'   **Caution (small-group random effects).** Both engines apply the same
#'   default prior --- the exact uniform-on-SD, which the INLA path
#'   represents via an expression-prior on the log-precision rather than
#'   the former PC approximation --- so they agree on the variance
#'   component far more closely than in earlier versions. A model with
#'   very few groups nonetheless carries a weakly-identified variance
#'   component (the data say little about the between-group spread), and
#'   INLA's Laplace approximation is less accurate there than full MCMC.
#'   For a flagship random-intercept model with few groups, pin
#'   `backend = "brms"` or supply an explicit informative prior (a
#'   half-Cauchy, per Gelman 2006; see the *priors and regularisation*
#'   vignette).
#'
#'   **Convergence.** The brms fit emits a warning when the sampler may
#'   not have converged (a parameter with Rhat at or above 1.1, or a low
#'   effective sample size). Treat such a posterior with caution ---
#'   increase `warmup` / `n_samples`, simplify the model, or supply a more
#'   informative prior --- and inspect the full diagnostics with
#'   [summary()]. Silence the warning (for intentionally short fits) via
#'   `options(flexyBayes.silence_convergence_warning = TRUE)`. The INLA
#'   path is deterministic and carries no such warning.
#'
#'   **Code inspection.** brms is the only engine with a code slot, so
#'   `return_code = TRUE` / `review_code = TRUE` under `"auto"` resolve to
#'   brms --- but only when brms can represent the model. Where it cannot,
#'   the call refuses with `auto_no_active_route` rather than returning
#'   the code for a different model. `backend = "inla"` with
#'   `return_code = TRUE` returns the INLA formula, family and
#'   hyperparameter list; with `review_code = TRUE` it refuses until the
#'   INLA-side code slot lands.
#'
#' @param aggregate One of `"auto"` (default), `TRUE`, or `FALSE`.
#'   Exact sufficient-statistics aggregation gate. INLA is the only
#'   engine with an aggregated emit. `"auto"` consults the aggregation
#'   plan and routes through the per-cell path when the IR is in scope
#'   (gaussian-identity, binomial-logit, or poisson-log; fixed plus
#'   random intercept; productive compression) and the model reaches
#'   INLA. `TRUE` forces aggregation and raises a structured refusal when
#'   the plan declares ineligibility, or when no active engine offers an
#'   aggregated route. `FALSE` skips the gate entirely. Aggregated fits
#'   carry `$exactness == "aggregated_exact"` and the dispatch trace's
#'   `path` slot reads `"aggregated_gaussian"` (gaussian) or
#'   `"aggregated_count"` (binomial / poisson). For out-of-core datasets
#'   that do not fit in memory, see [flexybayes_stream()].
#' @param syntax One of `"auto"` (default), `"asreml"`, or `"brms"`.
#'   Selects how `fixed` is interpreted. `"auto"` detects the grammar
#'   from the call shape (a bar-grouped formula is read as brms,
#'   otherwise ASReml); the other values force a grammar.
#'
#' @returns An object of class `"flexybayes"`, a list with three
#'   components.
#'
#' \describe{
#'   \item{`$glm`}{A GLM-compatible object (class `c("flexybayes_glm", "glm",
#'     "lm")`) with posterior mean coefficients, vcov, residuals, fitted values,
#'     etc. Works with `summary()`, `emmeans()`, `marginaleffects()`,
#'     `effectsize()`.}
#'   \item{`$brms` or `$inla`}{The native backend object --- a live `brmsfit`
#'     on the brms path, or INLA's own fitted object on the INLA path. Use
#'     with `bayesplot`, `posterior::as_draws()`, and each engine's own
#'     accessors.}
#'   \item{`$extras`}{Additional outputs: posterior `summary`, `convergence`
#'     diagnostics, `variance_comps`, `blups`, `predictions`, generated `code`,
#'     `param_names`, `parse_info`, `call_info`, `run_time`, `model_info`.}
#' }
#'
#' If `return_code = TRUE`, returns the generated backend code instead: a
#' character string holding the Stan program on the brms path, or a list of
#' the formula, family and hyperparameter specification on the INLA path.
#'
#' @examples
#' \dontrun{
#' # live fit -- needs an active backend (INLA or brms/Stan)
#' data(met_example, package = "flexyBayes")
#' # Simple random intercept model (small budget for example purposes)
#' fit <- flexybayes(
#'   fixed  = yield ~ env,
#'   random = ~ geno,
#'   data   = met_example$dat,
#'   n_samples = 100, warmup = 100, chains = 1, verbose = FALSE
#' )
#' summary(fit)
#' coef(fit)
#' }
#'
#' @importFrom stats terms model.frame model.matrix model.response
#'   formula family gaussian binomial poisson Gamma
#'   fitted residuals predict confint coef vcov logLik
#'   nobs update anova na.omit var cov quantile median
#'   qnorm pnorm dnorm setNames df.residual deviance
#'   .getXlevels dbinom density dpois lowess
#'   printCoefmat qqline qqnorm runif
#' @importFrom graphics abline axis hist legend lines par segments
#' @importFrom methods is
#' @importFrom coda effectiveSize gelman.diag
#' @importFrom splines bs
#' @export
flexybayes <- function(
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
  rcov = lifecycle::deprecated()
) {
  # `rcov` was the ASReml 3 name for the residual-structure argument;
  # ASReml 4 renamed it to `residual`, and flexyBayes followed suit in
  # 0.9.0. The formal is kept only as a tripwire so a stray `rcov =` --
  # including one forwarded through an engine pin -- raises a guiding
  # error instead of an opaque "unused argument".
  if (lifecycle::is_present(rcov)) {
    lifecycle::deprecate_stop(
      when = "0.9.0",
      what = "flexybayes(rcov)",
      with = "flexybayes(residual)",
      details = paste0(
        "ASReml-R renamed this argument from `rcov` (ASReml 3) to ",
        "`residual` (ASReml 4)."
      )
    )
  }

  # Refuse approximate-scheme requests
  # before match.arg fires its generic "should be one of" error.
  # The structured refusal points to the future
  # approximation registry and the available exact-route
  # alternatives. Runs first so the user sees the architectural
  # rationale rather than match.arg's surface message.
  # Accept backend = fb_engine(...) directly.
  # Resolve to the engine-name string before the approximate-scheme
  # check + match.arg, then apply the engine's sampler-control opts.
  engine_in <- backend
  backend <- .resolve_engine_string(engine_in)
  .check_approximate_scheme(backend)
  .check_known_backend_name(backend, allowed = c("auto", "inla", "brms"))
  backend <- match.arg(backend)
  eng_opts <- .fb_engine_opts(engine_in)
  if (!is.null(eng_opts)) {
    if (!is.null(eng_opts$n_samples)) {
      n_samples <- eng_opts$n_samples
    }
    if (!is.null(eng_opts$warmup)) {
      warmup <- eng_opts$warmup
    }
    if (!is.null(eng_opts$chains)) chains <- eng_opts$chains
  }
  aggregate <- .normalise_aggregate(aggregate)
  syntax <- match.arg(syntax)

  # Resolve session-level review-mode default BEFORE the unsupported-
  # backend guard, otherwise options(flexyBayes.review_code_default =
  # TRUE) would slip past the refusal on backend = "inla" / "auto".
  # Argument at call time wins.
  if (missing(review_code)) {
    review_code <- isTRUE(getOption("flexyBayes.review_code_default", FALSE))
  }

  # The code-inspection modes return generated backend code. brms is the
  # only active code-producing engine, so when backend resolves to
  # "auto" the code modes pick brms rather than a non-code INLA object.
  # The pick is CONDITIONAL on brms being able to represent the model --
  # see the capability gate below, which runs once the IR exists.
  # Rewriting unconditionally here would hand back code for a model the
  # engine cannot express (an AR1xAR1 residual lowered to an
  # intercept-only iid Gaussian, say), which is the silent-substitution
  # failure this guards against.
  code_mode_auto <- identical(backend, "auto") &&
    (isTRUE(return_code) || isTRUE(review_code))

  # review_code = TRUE is scoped to the code-emitting engine: brms
  # (Stan source via brms::make_stancode()). Under backend = "inla"
  # review_code is deferred to a future release -- the deferred-execution
  # token would need
  # an INLA-side `code` slot (the inla() formula + family + hyper list)
  # and a different proceed() target. Refuse cleanly rather than
  # silently emitting code for an engine that did not author it.
  # (brms support folded in here from the recast
  # fb_brms() pin, which now routes through this shared review branch.)
  if (
    isTRUE(review_code) &&
      !identical(backend, "brms") &&
      !code_mode_auto
  ) {
    stop(.fb_refusal_condition(
      reason_code = "review_code_backend_unsupported",
      message = paste0(
        "`review_code = TRUE` is supported with backend = \"brms\" (Stan ",
        "code via brms::make_stancode()). Under backend = \"",
        backend,
        "\" ",
        "the inspect-then-fit token would need an INLA-side code slot ",
        "(queued for a subsequent release). Pass ",
        "backend = \"brms\", or drop review_code."
      )
    ))
  }

  if (isTRUE(review_code) && isTRUE(return_code)) {
    stop(.fb_refusal_condition(
      reason_code = "code_flags_mutually_exclusive",
      message = paste0(
        "`return_code` and `review_code` are mutually exclusive. ",
        "Use `review_code = TRUE` for inspect-then-fit (returns a ",
        "<flexybayes_review> object); use `return_code = TRUE` for ",
        "the code string only."
      )
    ))
  }

  the_call <- match.call()
  data_name <- deparse(substitute(data))

  # Detect "all defaults" -- user supplied neither an fb_prior() nor
  # an explicit prior_vc_sd. In that case the v0.1 default fires:
  # build a bounded-uniform default keyed to the response scale (per
  # family / link) for sigma + every random group surfaced by the IR,
  # and emit the one-time announcement message. The uniform default
  # supersedes the earlier PC default.
  default_prior_active <- is.null(prior) && missing(prior_vc_sd)

  # Which of the two scalar prior arguments the caller actually wrote.
  # Both are documented as applied and neither was: `prior_fixed_sd` was
  # absent from the condition above, so passing it alone left the
  # auto-default in charge and the only consumer of the scalar was never
  # reached, and the INLA route never consumed `prior_vc_sd` at all.
  # Supplied-ness rather than value is what the emits need, because
  # applying the documented default unconditionally is not the fix:
  # brms centres its own intercept prior on the response, and replacing
  # that with `normal(0, 100)` would be informative in the wrong
  # direction for any response far from zero. Supplied means honoured;
  # not supplied means the engine's own default, which prior_summary()
  # now says out loud.
  prior_scalars <- list(
    fixed_sd = prior_fixed_sd,
    vc_sd = prior_vc_sd,
    supplied = c(
      fixed_sd = !missing(prior_fixed_sd),
      vc_sd = !missing(prior_vc_sd)
    )
  )

  # Build the flexyBayes intermediate representation (IR). The universal
  # entry detects the grammar from the call shape -- ASReml
  # fixed/random/residual, or a brms-style bar-grouped formula -- and
  # routes to the matching ingest adapter. ASReml
  # ingest is byte-identical to the historical direct fb_from_asreml()
  # call; `syntax = ` forces a grammar. See .build_ir_polymorphic() in
  # R/fb.R.
  fb <- .build_ir_polymorphic(
    fixed = fixed,
    random = random,
    residual = residual,
    data = data,
    family = family,
    link = link,
    weights = weights,
    known_matrices = known_matrices,
    prior = prior,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    syntax = syntax
  )
  fb$prior_scalars <- prior_scalars

  # Missing responses, resolved once the IR exists so the layer can see
  # which variables index a structured covariance. Under the default,
  # a row whose response is missing is retained and carried as a latent
  # quantity, keeping the design index set intact; see R/na_action.R.
  #
  # The policy is normalised here rather than inside the layer so the
  # recorded call carries the native word for whatever spelling the user
  # wrote -- an asreml na.method() value, a bare list, or a string.
  na_policy <- .fb_normalise_na_action(na_action)
  na_meta <- .fb_apply_na_action(fb, data, na_policy)
  data <- na_meta$data
  # The IR caches the row count for the preflight size estimates, so it
  # has to follow the augmented data rather than the data as supplied.
  if (!is.null(fb$data_summary)) {
    fb$data_summary$n <- nrow(data)
  }

  # Capability gate for the code-inspection modes on backend = "auto"
  # (flagged above, resolved here now that the IR exists). brms is the
  # only active code-producing engine, so auto picks it -- but only for a
  # model brms can actually represent. Where it cannot, returning its
  # code would present a DIFFERENT model as the answer to the user's
  # request; a structured residual, for instance, has no brms lowering
  # and would come back as an intercept-only iid Gaussian program. Refuse
  # with the capability predicate's own reason code instead.
  if (code_mode_auto) {
    # Named structures first. The generic capability message below says only
    # that brms cannot represent the model, which is true and unhelpful when
    # the reason is that the structure has a name and an alternative -- a
    # residual-side AR1 spelling above all. The guard is the same one the
    # dispatch choke point runs, so the two routes give the same answer.
    .refuse_unrepresentable_structures(fb)
    cap <- .backend_can_fit("brms", fb)
    if (!isTRUE(cap$ok)) {
      stop(.fb_refusal_condition(
        reason_code = "auto_no_active_route",
        message = paste0(
          "backend = \"auto\" with code inspection resolves to brms (the ",
          "only active code-producing engine), but ",
          "brms cannot represent this model (",
          cap$reason_code,
          "). ",
          "Returning its code would show a different model from the one ",
          "requested. Pass backend = \"inla\" to fit it, or reformulate ",
          "for an engine that can express it."
        ),
        family_class = "flexybayes_auto_no_route_refusal",
        backend = "auto"
      ))
    }
    backend <- "brms"
  }

  if (default_prior_active) {
    # Which terms the default reaches is decided in one place,
    # .fb_default_prior_targets(), because the model fingerprint has to
    # record the same split: a variance component the walker skipped
    # carries whatever the engine chose, and triangulate() must exclude it
    # from a cross-engine comparison rather than compare two answers to
    # two different questions.
    prior_targets <- .fb_default_prior_targets(fb)
    unif_default <- .default_uniform_prior(
      data = data,
      response = fb$response,
      family = family,
      link = link,
      random_groups = prior_targets$shared,
      vm_ped_groups = prior_targets$vm_ped
    )
    fb$priors <- unif_default
    prior <- unif_default
    .default_prior_note_once(
      scale = attr(unif_default, "fb_prior_default_scale"),
      basis = attr(unif_default, "fb_prior_default_basis")
    )
  }

  # plan = TRUE: short-circuit after IR build. The plan
  # surface lives in fb_plan(); flexybayes(plan = TRUE) is the
  # courtesy alternative invocation that lets asreml-style callers
  # reach the same planning object without re-typing the formula in
  # brms shape.
  if (isTRUE(plan)) {
    return(.fb_plan_from_ir(
      fb = fb,
      data = data,
      backend = backend,
      known_matrices = known_matrices,
      aggregate = aggregate,
      memory_ceiling_gb = NULL,
      predict_plan = NULL,
      the_call = the_call,
      data_name = data_name
    ))
  }

  # Review-mode branch. Build the deferred-execution token instead
  # of firing the backend. Code generation does not consume RNG; the
  # .Random.seed snapshot is captured before any RNG-touching step so
  # that proceed(rev) reproduces the chain a direct call would have
  # produced at the same outer seed. verbose printing is suppressed
  # because the review object owns the code surface (cat_code(rev)).
  if (isTRUE(review_code)) {
    # The emit engine for the review code is always brms (Stan source via
    # brms::make_stancode()): the guard above (review_code_backend_unsupported)
    # ensures `backend` is already "brms" by the time this branch runs,
    # whether explicitly requested or resolved here via code_mode_auto.
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop(
        "Package 'brms' is required to generate Stan review code. ",
        "Install with: install.packages('brms'). A working C++ ",
        "toolchain (rstan or cmdstanr) is required for the ",
        "downstream fit.",
        call. = FALSE
      )
    }
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      runif(1) # initialise RNG so the snapshot is reproducible
    }
    seed_snapshot <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

    # Build a `the_call` for the deferred fit that does not reach
    # back into the caller's frame for symbol resolution. `match.call()`
    # captured `data = <symbol>` (and similarly for any other
    # non-literal argument such as `weights`, `known_matrices`,
    # `prior`); when proceed(rev) fires later, those caller-frame
    # bindings may have gone out of scope. R's standard fit-object
    # pipeline (model.frame / terms / class dispatch on the glm
    # surface inside .build_glm) lazily evaluates the stored call in
    # certain paths, which raises "object '<symbol>' not found" if
    # the binding is gone. Replace the symbol slots with their
    # current values so the proceed-side call is self-contained.
    # The user-side captured call (stored on the review object as
    # `$call`) keeps the original symbol form for diagnostic clarity.
    the_call_proceed <- the_call
    the_call_proceed$data <- data
    if (!is.null(weights)) {
      the_call_proceed$weights <- weights
    }
    if (length(known_matrices)) {
      the_call_proceed$known_matrices <- known_matrices
    }
    if (!is.null(prior)) {
      the_call_proceed$prior <- prior
    }

    # Run preflight upstream of the code emit so the
    # review token carries the design-memory summary and a refusal
    # short-circuits before any code generation happens. Below the
    # 1e5-row threshold .maybe_preflight() returns NULL and the
    # review object's $preflight slot stays NULL (v0.2 behaviour).
    review_preflight <- .maybe_preflight(
      fb = fb,
      data = data,
      the_call = the_call
    )

    # Emit the review code (Stan source via make_stancode()). Derives the
    # model from the IR on return_code = TRUE, so the fixed / random /
    # residual args (display-only) pass through as-is.
    review_code_str <- emit_brms(
      fb = fb,
      data = data,
      known_matrices = known_matrices,
      weights = weights,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = FALSE,
      mcmc_verbose = mcmc_verbose,
      return_code = TRUE,
      the_call = the_call,
      fixed = fixed,
      random = random,
      residual = residual,
      family = family,
      link = link,
      data_name = data_name
    )

    return(.new_flexybayes_review(
      code = review_code_str,
      backend = "stan_via_brms",
      ir = fb,
      prior = prior,
      data_name = data_name,
      call = the_call,
      seed = seed_snapshot,
      preflight = review_preflight,
      proceed_args = list(
        fb = fb,
        data = data,
        known_matrices = known_matrices,
        weights = weights,
        n_samples = n_samples,
        warmup = warmup,
        chains = chains,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd,
        verbose = FALSE,
        mcmc_verbose = mcmc_verbose,
        return_code = FALSE,
        the_call = the_call_proceed,
        fixed = fixed,
        random = random,
        residual = residual,
        family = family,
        link = link,
        data_name = data_name
      )
    ))
  }

  # Backend dispatch lives in R/dispatch.R as the shared helper
  # `.dispatch_backend()` (lifted so fb_brms() drives
  # the same routing). Semantics: backend = "brms" skips the gate;
  # backend = "inla" calls lgm_gate() and raises the refusal on
  # non-LGM; backend = "auto" gate-then-route with silenceable
  # fall-back notes.
  fit <- .dispatch_backend(
    fb = fb,
    data = data,
    backend = backend,
    known_matrices = known_matrices,
    weights = weights,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = control,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    verbose = verbose,
    mcmc_verbose = mcmc_verbose,
    return_code = return_code,
    the_call = the_call,
    fixed = fixed,
    random = random,
    residual = residual,
    family = family,
    link = link,
    data_name = data_name,
    aggregate = aggregate,
    na_action = na_policy$y
  )
  # Record what the missing-response layer did. Whether a posterior was
  # computed on the design as laid out or on the plots that survived is
  # not recoverable from the fit otherwise, and the two are different
  # models whenever a covariance is indexed by the design.
  if (!isTRUE(return_code) && is.list(fit) && !is.null(fit$extras)) {
    fit$extras$na_action <- na_meta$meta
    # The active emits record the policy in the call record themselves.
    # The aggregated and worker-hosted routes reach their emit by another
    # path, so the record is completed here rather than left absent on
    # one route out of several.
    if (
      !is.null(fit$extras$call_info) &&
        is.null(fit$extras$call_info$na_action)
    ) {
      fit$extras$call_info$na_action <- na_policy$y
    }
  }
  .fb_warn_poor_convergence(fit)
  fit
}


# ----------------------------------------------------------------- #
# fb -- alias for flexybayes()                                      #
# ----------------------------------------------------------------- #
#
# `fb` is a literal alias for `flexybayes()`. One canonical
# asreml-format implementation; two exported names for typing
# economy. The brms-format ingest path (former fb() body) is
# deferred to v0.2 as `fb_brms()`; the internal helper
# `fb_from_brms()` in R/fb_from_brms.R remains unexported as v0.2
# work-continuity code.

#' @rdname flexybayes
#' @export
fb <- flexybayes
