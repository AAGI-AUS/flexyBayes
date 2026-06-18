# Bayesian Mixed Models with ASReml Syntax

Specify mixed models using ASReml formula syntax and estimate them via
Bayesian MCMC using greta. Returns a three-part result: a GLM-compatible
object for use with standard R packages (emmeans, marginaleffects,
etc.), native greta output for Bayesian diagnostics, and extras for
secondary analyses.

## Usage

``` r
flexybayes(
  fixed,
  random = NULL,
  rcov = NULL,
  data,
  family = "gaussian",
  link = NULL,
  known_matrices = list(),
  weights = NULL,
  n_samples = 1000,
  warmup = 500,
  chains = 4,
  prior = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  review_code = FALSE,
  backend = c("auto", "greta", "inla", "brms", "gretaR"),
  aggregate = "auto",
  plan = FALSE,
  syntax = c("auto", "asreml", "brms", "greta")
)

fb(
  fixed,
  random = NULL,
  rcov = NULL,
  data,
  family = "gaussian",
  link = NULL,
  known_matrices = list(),
  weights = NULL,
  n_samples = 1000,
  warmup = 500,
  chains = 4,
  prior = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  review_code = FALSE,
  backend = c("auto", "greta", "inla", "brms", "gretaR"),
  aggregate = "auto",
  plan = FALSE,
  syntax = c("auto", "asreml", "brms", "greta")
)
```

## Arguments

- fixed:

  Two-sided formula `response ~ fixed_effects`. This is the universal
  entry's model slot: it accepts the ASReml `fixed` form (paired with
  `random` / `rcov`) **or** a brms / lme4-style bar-grouped formula such
  as `response ~ x + (1 | g)` (in which case the grouping lives in the
  formula and `random` / `rcov` must be left `NULL`). The grammar is
  detected from the call shape; use `syntax` to force it.

- random:

  One-sided formula: `~ random_terms` using ASReml syntax. Supports
  `vm()`, `at()`, `us()`, `fa()`, `ar1()`, `spl()`, `ped()`, `dsum()`,
  `id()`, and nested colon terms.

- rcov:

  One-sided formula: `~ residual_structure`. Default `~ units` (iid
  residuals). Use `~ at(env):units` for heterogeneous variance.

- data:

  A data.frame containing all variables referenced in the formulas.

- family:

  Character: `"gaussian"`, `"binomial"`, `"poisson"`,
  `"negative_binomial"`, `"gamma"`, or `"beta"`.

- link:

  Character or NULL: override the default link function (e.g.,
  `"probit"` for binomial).

- known_matrices:

  Named list of matrices referred to in the random formula (e.g.,
  `list(Gmat = G_mat, Amat = A_mat)`). The carrier is declared with the
  [`fb_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_cov.md)
  constructor inside the random term: a dense covariance
  (`vm(group, cov = fb_cov(G, type = "dense"))`), a user-supplied
  lower-triangular Cholesky factor
  (`vm(group, cov = fb_cov(L, type = "chol"))`; greta backend only), or
  a sparse precision matrix `Matrix::dgCMatrix`
  (`vm(group, cov = fb_cov(Q, type = "precision"))`; greta and INLA
  backends). The bare dense forms `vm(group, V = ...)` /
  `ped(group, A = ...)` remain the default; the legacy v0.3.7 keyword
  carriers (`chol = `, `precision = `, `blocks = `,
  `low_rank_factor = `) are deprecated and emit a migration warning. See
  the structured-covariance vignette for per-type worked examples.

- weights:

  Optional numeric weight vector (length N); sets Var = sigma^2 / w.

- n_samples:

  Integer: number of posterior samples per chain.

- warmup:

  Integer: number of warmup (burn-in) iterations per chain.

- chains:

  Integer: number of MCMC chains.

- prior:

  An optional
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  object specifying priors via the PC-canonical hybrid DSL (preferred).
  When supplied it overrides `prior_vc_sd` for the variance components
  it covers. See
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md).

- prior_fixed_sd:

  Numeric: SD for fixed-effect normal priors, applied uniformly to the
  intercept, factor contrasts, continuous slopes, factor x continuous
  interactions, and [`I()`](https://rdrr.io/r/base/AsIs.html)-expression
  terms. Default `100` — weakly informative on the natural response
  scale for the vast majority of agricultural / clinical responses
  (covers responses with central tendency up to several hundred without
  crushing the posterior toward zero), while still regularising at
  sample sizes below ~ 30 per coefficient. Set wider (e.g. `1000`) for
  responses on a larger natural scale, or narrower for explicit
  shrinkage. A weakly-informative normal prior on the data scale, in the
  spirit of the weakly-informative-prior literature for regression
  coefficients (e.g. Gelman et al. 2008, *Annals of Applied Statistics*
  2(4):1360-1383).

- prior_vc_sd:

  Numeric: hyperparameter for the legacy `lognormal(0, prior_vc_sd)`
  variance-component prior. **Note:** when both `prior` and
  `prior_vc_sd` are left at their defaults, v0.1 activates a
  bounded-uniform default on the SD scale (`uniform(0, U)` with
  family-aware `U`: `5 * sd(y)` for Gaussian; `5` on the logit scale for
  binomial / beta; `3` on the log scale for Poisson / negative-binomial
  / gamma) for `sigma` and every named random-effect group – a
  weakly-informative choice for moderate group counts; for very small
  `J`, Gelman (2006), *Bayesian Analysis* 1(3):515-534, recommends a
  half-t / half-Cauchy instead (see the priors vignette). The legacy
  `lognormal(0, 1)` default fires only when `prior_vc_sd` is passed
  explicitly. Silence the one-time announcement message via
  `options(flexyBayes.silence_default_prior_note = TRUE)`.

- verbose:

  Logical: print generated greta code to console.

- mcmc_verbose:

  Logical: show MCMC progress bar from greta.

- return_code:

  Logical: if TRUE, return the generated code string without fitting the
  model.

- review_code:

  Logical: if TRUE, do not fit the model immediately; instead return a
  `<flexybayes_review>` deferred-execution object carrying the generated
  backend code, the resolved prior, the parsed intermediate
  representation (IR), the captured call, and a snapshot of
  `.Random.seed`. Inspect the code with
  [`cat_code()`](https://aagi-aus.github.io/flexyBayes/reference/cat_code.md);
  run the deferred fit with
  [`proceed()`](https://aagi-aus.github.io/flexyBayes/reference/proceed.md);
  a second
  [`proceed()`](https://aagi-aus.github.io/flexyBayes/reference/proceed.md)
  call returns the cached fit. Useful as a teaching / auditing surface
  before a long MCMC run. Default `FALSE` preserves the existing
  run-immediately semantics. A session-level override is available via
  `options(flexyBayes.review_code_default = TRUE)`; the argument value
  at call time always wins. Closest published precedent:
  `brms::make_stancode()` plus the `chains = 0` idiom for "do everything
  except sampling". `review_code = TRUE` and `return_code = TRUE` are
  mutually exclusive.

- backend:

  Character: one of `"auto"` (**default**), `"greta"`, `"inla"`, or
  `"brms"`. Under `"auto"` the call consults `lgm_gate()` and routes to
  INLA when the model is latent-Gaussian feasible (deterministic and
  faster on that certified class), and to greta (full MCMC, the
  universal fallback) otherwise — so a no-`backend` call reaches every
  available engine the model supports.

  **Caution (small-group random effects).** Both backends now apply the
  same default prior — the exact uniform-on-SD (the INLA path represents
  it faithfully via an expression-prior on the log-precision rather than
  the former PC approximation), so the two engines agree on the variance
  component far more closely than in earlier versions. A model with very
  few groups nonetheless carries a weakly-identified variance component
  (the data say little about the between-group spread), and INLA's
  Laplace approximation is less accurate there than full MCMC. For a
  flagship random-intercept model with few groups, prefer
  `backend = "greta"` or supply an explicit informative prior (a
  half-Cauchy, per Gelman 2006; see the *priors* vignette).
  [`fb_backend_status()`](https://aagi-aus.github.io/flexyBayes/reference/fb_backend_status.md)
  reports which engines are usable.

  **Convergence.** MCMC fits (`"greta"`, `"brms"`) emit a warning when
  the sampler may not have converged (a parameter with Rhat at or above
  1.1, or a low effective sample size). Treat such a posterior with
  caution — increase `warmup` / `n_samples`, simplify the model, or
  supply a more informative prior — and inspect the full diagnostics
  with [`summary()`](https://rdrr.io/r/base/summary.html). Silence the
  warning (for intentionally short fits) via
  `options(flexyBayes.silence_convergence_warning = TRUE)`. The INLA
  path is deterministic and carries no such warning.

  On the greta backend, generalised mixed models (a `poisson` or
  `binomial` family with random effects) adapt more slowly than the
  Gaussian case and need a larger `warmup` than the default — a few
  thousand iterations is typical, and the convergence warning above will
  tell you when more is needed. This is an adaptation-budget matter, not
  a parameterisation one: the random effects already use the non-centred
  form, which mixes at least as well as the centred form across the
  regimes tested. For a latent-Gaussian GLMM the `"auto"` route
  sidesteps the question entirely by fitting it with INLA.
  Factor-analytic / unstructured covariance terms (greta only) are
  harder still and may not mix at modest budgets; judge their
  convergence on the identified covariance via
  [`fb_structured_cov()`](https://aagi-aus.github.io/flexyBayes/reference/fb_structured_cov.md)
  rather than the rotation-/sign-ambiguous raw loadings.

  `"greta"` forces full MCMC; `"inla"` forces INLA and raises a
  structured refusal if the model is not LGM-feasible; `"brms"` is the
  Stan passthrough via brms (refuses asreml structured-covariance terms
  — `vm`/`ped`/`fa`/`us`/`ar1` — that have no lossless Stan
  translation). Under `"auto"`, falling back to greta (gate refusal,
  INLA not installed, or an INLA numerical failure) emits a one-time
  silenceable note;
  `options(flexyBayes.silence_auto_fallback_note = TRUE)` silences the
  gate-refusal / numerical-fallback note and
  `options(flexyBayes.silence_auto_inla_missing_note = TRUE)` the
  INLA-not-installed note.
  [`backend_decision()`](https://aagi-aus.github.io/flexyBayes/reference/backend_decision.md)
  surfaces the full dispatch trace (including `rejected_routes`)
  post-fit. When `return_code = TRUE` or `review_code = TRUE` is
  requested under `"auto"`, the call resolves to greta (the
  code-producing engine); an explicit `backend = "inla"` with
  `review_code = TRUE` raises a structured refusal until INLA-side
  `code`-slot support lands.

- aggregate:

  One of `"auto"` (default), `TRUE`, or `FALSE`. Exact
  sufficient-statistics aggregation gate. `"auto"` consults the
  aggregation plan and routes through the per-cell emit path when the IR
  is in scope (gaussian-identity, binomial-logit, or poisson-log;
  fixed + random-intercept; productive compression) **on
  `backend = "inla"`**; on `backend = "greta"`, `"auto"` falls through
  to the per-row path even when the plan is eligible (`TRUE` is required
  to activate the greta aggregated emit explicitly). `TRUE` forces
  aggregation on either backend – raises a structured refusal when the
  plan declares ineligibility. `FALSE` skips the gate entirely.
  Aggregated fits carry `$exactness == "aggregated_exact"` and the
  dispatch trace's `path` slot reads `"aggregated_gaussian"` (gaussian)
  or `"aggregated_count"` (binomial / poisson). For out-of-core datasets
  that do not fit in memory, see
  [`flexybayes_stream()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes_stream.md).
  The asymmetry is documented behaviour rather than a bug: greta uses
  dense-matrix per-row emit by default for predictable RAM, and
  switching to the aggregated path silently would change the posterior
  numerical profile without an explicit opt-in.

- plan:

  Logical: if `TRUE`, short-circuit after intermediate representation
  (IR) build and return a `<fb_plan>` object carrying the IR, the
  routing decision, the representation plan, the aggregation plan, and
  the prediction plan, *without* emitting backend code or fitting.
  Equivalent to calling
  [`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md)
  with the same formula triple — exposed inline as
  `flexybayes(plan = TRUE)` so users can reach the planning object
  without re-typing the call. Default `FALSE` preserves the run-
  immediately semantics. Mutually exclusive with `return_code` and
  `review_code`.

- syntax:

  One of `"auto"` (default), `"asreml"`, `"brms"`, or `"greta"`. Selects
  how `fixed` is interpreted. `"auto"` detects the grammar from the call
  shape (a bar-grouped formula is read as brms, otherwise ASReml); the
  other values force a grammar. `"greta"` (a native `greta_model`) is
  reserved: pass such models to
  [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
  in v0.4.x; direct ingest through the universal entry lands at v0.5.0.

## Value

An object of class `"flexybayes"` — a list with three components:

- `$glm`:

  A GLM-compatible object (class `c("flexybayes_glm", "glm", "lm")`)
  with posterior mean coefficients, vcov, residuals, fitted values, etc.
  Works with [`summary()`](https://rdrr.io/r/base/summary.html),
  `emmeans()`, `marginaleffects()`, `effectsize()`.

- `$greta`:

  Native greta output: `model`, `draws` (mcmc.list), `greta_arrays`, and
  `env`. Use with `bayesplot`, `greta::calculate()`,
  `posterior::as_draws()`.

- `$extras`:

  Additional outputs: posterior `summary`, `convergence` diagnostics,
  `variance_comps`, `blups`, `predictions`, generated `code`,
  `param_names`, `parse_info`, `call_info`, `run_time`, `model_info`.

If `return_code = TRUE`, returns a character string of greta code
instead.

## Examples

``` r
if (FALSE) { # \dontrun{
# live fit -- needs a backend (greta Python/TF, INLA, or brms/Stan)
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
