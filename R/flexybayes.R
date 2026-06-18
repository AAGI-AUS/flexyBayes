#' Bayesian Mixed Models with ASReml Syntax
#'
#' Specify mixed models using ASReml formula syntax and estimate them via
#' Bayesian MCMC using greta. Returns a three-part result: a GLM-compatible
#' object for use with standard R packages (emmeans, marginaleffects, etc.),
#' native greta output for Bayesian diagnostics, and extras for secondary
#' analyses.
#'
#' @param fixed Two-sided formula `response ~ fixed_effects`. This is the
#'   universal entry's model slot: it accepts the ASReml `fixed` form
#'   (paired with `random` / `rcov`) **or** a brms / lme4-style
#'   bar-grouped formula such as `response ~ x + (1 | g)` (in which case
#'   the grouping lives in the formula and `random` / `rcov` must be
#'   left `NULL`). The grammar is detected from the call shape; use
#'   `syntax` to force it.
#' @param random One-sided formula: `~ random_terms` using ASReml syntax.
#'   Supports `vm()`, `at()`, `us()`, `fa()`, `ar1()`, `spl()`, `ped()`,
#'   `dsum()`, `id()`, and nested colon terms.
#' @param rcov One-sided formula: `~ residual_structure`. Default `~ units`
#'   (iid residuals). Use `~ at(env):units` for heterogeneous variance.
#' @param data A data.frame containing all variables referenced in the formulas.
#' @param family Character: `"gaussian"`, `"binomial"`, `"poisson"`,
#'   `"negative_binomial"`, `"gamma"`, or `"beta"`.
#' @param link Character or NULL: override the default link function
#'   (e.g., `"probit"` for binomial).
#' @param known_matrices Named list of matrices referred to in the random
#'   formula (e.g., `list(Gmat = G_mat, Amat = A_mat)`). The carrier is
#'   declared with the [fb_cov()] constructor inside the random term:
#'   a dense covariance (`vm(group, cov = fb_cov(G, type = "dense"))`),
#'   a user-supplied lower-triangular Cholesky factor
#'   (`vm(group, cov = fb_cov(L, type = "chol"))`; greta backend only),
#'   or a sparse precision matrix `Matrix::dgCMatrix`
#'   (`vm(group, cov = fb_cov(Q, type = "precision"))`; greta and INLA
#'   backends). The bare dense forms `vm(group, V = ...)` /
#'   `ped(group, A = ...)` remain the default; the legacy v0.3.7
#'   keyword carriers (`chol = `, `precision = `, `blocks = `,
#'   `low_rank_factor = `) are deprecated and emit a migration warning.
#'   See the structured-covariance vignette for per-type worked examples.
#' @param weights Optional numeric weight vector (length N); sets
#'   Var = sigma^2 / w.
#' @param n_samples Integer: number of posterior samples per chain.
#' @param warmup Integer: number of warmup (burn-in) iterations per chain.
#' @param chains Integer: number of MCMC chains.
#' @param prior An optional `fb_prior()` object specifying priors via the
#'   PC-canonical hybrid DSL (preferred). When supplied it overrides
#'   `prior_vc_sd` for the variance components it covers. See
#'   [fb_prior()].
#' @param prior_fixed_sd Numeric: SD for fixed-effect normal priors,
#'   applied uniformly to the intercept, factor contrasts, continuous
#'   slopes, factor x continuous interactions, and `I()`-expression
#'   terms. Default `100` — weakly informative on the natural response
#'   scale for the vast majority of agricultural / clinical responses
#'   (covers responses with central tendency up to several hundred
#'   without crushing the posterior toward zero), while still
#'   regularising at sample sizes below ~ 30 per coefficient. Set
#'   wider (e.g. `1000`) for responses on a larger natural scale, or
#'   narrower for explicit shrinkage. A weakly-informative normal prior
#'   on the data scale, in the spirit of the weakly-informative-prior
#'   literature for regression coefficients (e.g. Gelman et al. 2008,
#'   *Annals of Applied Statistics* 2(4):1360-1383).
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
#'   `options(flexyBayes.silence_default_prior_note = TRUE)`.
#' @param verbose Logical: print generated greta code to console.
#' @param mcmc_verbose Logical: show MCMC progress bar from greta.
#' @param return_code Logical: if TRUE, return the generated code string
#'   without fitting the model.
#' @param review_code Logical: if TRUE, do not fit the model immediately;
#'   instead return a `<flexybayes_review>` deferred-execution object
#'   carrying the generated backend code, the resolved prior, the
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
#'   `return_code = TRUE` are mutually exclusive.
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
#' @param backend Character: one of `"auto"` (**default**), `"greta"`,
#'   `"inla"`, or `"brms"`. Under `"auto"` the call consults
#'   `lgm_gate()` and routes to INLA when the model is latent-Gaussian
#'   feasible (deterministic and faster on that certified class), and to
#'   greta (full MCMC, the universal fallback) otherwise --- so a
#'   no-`backend` call reaches every available engine the model supports.
#'
#'   **Caution (small-group random effects).** Both backends now apply
#'   the same default prior --- the exact uniform-on-SD (the INLA path
#'   represents it faithfully via an expression-prior on the log-precision
#'   rather than the former PC approximation), so the two engines agree on
#'   the variance component far more closely than in earlier versions. A
#'   model with very few groups nonetheless carries a weakly-identified
#'   variance component (the data say little about the between-group
#'   spread), and INLA's Laplace approximation is less accurate there than
#'   full MCMC. For a flagship random-intercept model with few groups,
#'   prefer `backend = "greta"` or supply an explicit informative prior
#'   (a half-Cauchy, per Gelman 2006; see the *priors* vignette).
#'   [fb_backend_status()] reports which engines are usable.
#'
#'   **Convergence.** MCMC fits (`"greta"`, `"brms"`) emit a warning when
#'   the sampler may not have converged (a parameter with Rhat at or above
#'   1.1, or a low effective sample size). Treat such a posterior with
#'   caution --- increase `warmup` / `n_samples`, simplify the model, or
#'   supply a more informative prior --- and inspect the full diagnostics
#'   with [summary()]. Silence the warning (for intentionally short fits)
#'   via `options(flexyBayes.silence_convergence_warning = TRUE)`. The INLA
#'   path is deterministic and carries no such warning.
#'
#'   On the greta backend, generalised mixed models (a `poisson` or
#'   `binomial` family with random effects) adapt more slowly than the
#'   Gaussian case and need a larger `warmup` than the default --- a few
#'   thousand iterations is typical, and the convergence warning above will
#'   tell you when more is needed. This is an adaptation-budget matter, not
#'   a parameterisation one: the random effects already use the
#'   non-centred form, which mixes at least as well as the centred form
#'   across the regimes tested. For a latent-Gaussian GLMM the `"auto"`
#'   route sidesteps the question entirely by fitting it with INLA.
#'   Factor-analytic / unstructured covariance terms (greta only) are
#'   harder still and may not mix at modest budgets; judge their
#'   convergence on the identified covariance via [fb_structured_cov()]
#'   rather than the rotation-/sign-ambiguous raw loadings.
#'
#'   `"greta"` forces full MCMC; `"inla"` forces INLA and raises a
#'   structured refusal if the model is not LGM-feasible; `"brms"` is the
#'   Stan passthrough via brms (refuses asreml structured-covariance
#'   terms --- `vm`/`ped`/`fa`/`us`/`ar1` --- that have no lossless Stan
#'   translation). Under `"auto"`, falling back to greta (gate refusal,
#'   INLA not installed, or an INLA numerical failure) emits a one-time
#'   silenceable note; `options(flexyBayes.silence_auto_fallback_note =
#'   TRUE)` silences the gate-refusal / numerical-fallback note and
#'   `options(flexyBayes.silence_auto_inla_missing_note = TRUE)` the
#'   INLA-not-installed note. [backend_decision()] surfaces the full
#'   dispatch trace (including `rejected_routes`) post-fit. When
#'   `return_code = TRUE` or `review_code = TRUE` is requested under
#'   `"auto"`, the call resolves to greta (the code-producing engine);
#'   an explicit `backend = "inla"` with `review_code = TRUE` raises a
#'   structured refusal until INLA-side `code`-slot support lands.
#'
#' @param aggregate One of `"auto"` (default), `TRUE`, or `FALSE`.
#'   Exact sufficient-statistics aggregation gate. `"auto"`
#'   consults the aggregation plan and routes through the per-cell
#'   emit path when the IR is in scope (gaussian-identity, binomial-logit,
#'   or poisson-log; fixed + random-intercept; productive compression)
#'   **on `backend = "inla"`**; on `backend = "greta"`, `"auto"` falls
#'   through to the per-row path even when the plan is eligible (`TRUE`
#'   is required to activate the greta aggregated emit explicitly).
#'   `TRUE` forces aggregation on either backend -- raises a structured
#'   refusal when the plan declares ineligibility. `FALSE` skips the gate
#'   entirely. Aggregated fits carry `$exactness == "aggregated_exact"`
#'   and the dispatch trace's `path` slot reads `"aggregated_gaussian"`
#'   (gaussian) or `"aggregated_count"` (binomial / poisson). For
#'   out-of-core datasets that do not fit in memory, see
#'   [flexybayes_stream()]. The asymmetry is documented behaviour rather
#'   than a bug: greta
#'   uses dense-matrix per-row emit by default for predictable RAM,
#'   and switching to the aggregated path silently would change the
#'   posterior numerical profile without an explicit opt-in.
#' @param syntax One of `"auto"` (default), `"asreml"`, `"brms"`, or
#'   `"greta"`. Selects how `fixed` is interpreted. `"auto"` detects the
#'   grammar from the call shape (a bar-grouped formula is read as brms,
#'   otherwise ASReml); the other values force a grammar. `"greta"` (a
#'   native `greta_model`) is reserved: pass such models to `fb_greta()`
#'   in v0.4.x; direct ingest through the universal entry lands at v0.5.0.
#'
#' @return An object of class `"flexybayes"` — a list with three components:
#' \describe{
#'   \item{`$glm`}{A GLM-compatible object (class `c("flexybayes_glm", "glm",
#'     "lm")`) with posterior mean coefficients, vcov, residuals, fitted values,
#'     etc. Works with `summary()`, `emmeans()`, `marginaleffects()`,
#'     `effectsize()`.}
#'   \item{`$greta`}{Native greta output: `model`, `draws` (mcmc.list),
#'     `greta_arrays`, and `env`. Use with `bayesplot`, `greta::calculate()`,
#'     `posterior::as_draws()`.}
#'   \item{`$extras`}{Additional outputs: posterior `summary`, `convergence`
#'     diagnostics, `variance_comps`, `blups`, `predictions`, generated `code`,
#'     `param_names`, `parse_info`, `call_info`, `run_time`, `model_info`.}
#' }
#'
#' If `return_code = TRUE`, returns a character string of greta code instead.
#'
#' @examples
#' \dontrun{
#' # live fit -- needs a backend (greta Python/TF, INLA, or brms/Stan)
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
) {
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

  # The universal entry accepts a native greta_model
  # graph or a prebuilt greta-source IR (from fb_from_greta()) on the
  # model slot. Those carry their own data into the greta graph, so
  # `data` is optional on that path; the formula grammars still require
  # it. Detect before `data` is forced below.
  spec_is_greta_native <- inherits(fixed, "greta_model") ||
    (inherits(fixed, "fb_terms") && identical(fixed$source, "greta"))
  if (spec_is_greta_native && missing(data)) {
    data <- NULL
  }

  # Defer the greta package check to the emit-greta branch. For
  # backend = "inla" we don't need greta installed; for backend =
  # "auto" we only require greta if the gate refuses (or INLA is
  # unavailable) and we fall through to the greta emit.

  # Resolve session-level review-mode default BEFORE the unsupported-
  # backend guard, otherwise options(flexyBayes.review_code_default =
  # TRUE) would slip past the refusal on backend = "inla" / "auto".
  # Argument at call time wins.
  if (missing(review_code)) {
    review_code <- isTRUE(getOption("flexyBayes.review_code_default", FALSE))
  }

  # The code-inspection modes return generated
  # backend code, which only greta produces under this entry. When
  # backend resolves to "auto", pick greta -- the engine that satisfies
  # the request -- rather than refusing or returning a non-code INLA
  # object. auto's contract is to resolve to a capable engine.
  if (
    identical(backend, "auto") &&
      (isTRUE(return_code) || isTRUE(review_code))
  ) {
    backend <- "greta"
  }

  # review_code = TRUE is scoped to the code-emitting engines: greta
  # (greta source via the codegen path) and brms (Stan source via
  # brms::make_stancode()). Under backend = "inla" review_code is
  # deferred to a future release -- the deferred-execution token would need
  # an INLA-side `code` slot (the inla() formula + family + hyper list)
  # and a different proceed() target. Refuse cleanly rather than
  # silently emitting code for an engine that did not author it.
  # (brms support folded in here from the recast
  # fb_brms() pin, which now routes through this shared review branch.)
  if (isTRUE(review_code) && !backend %in% c("greta", "brms")) {
    stop(.fb_refusal_condition(
      reason_code = "review_code_backend_unsupported",
      message = paste0(
        "`review_code = TRUE` is supported with backend = \"greta\" ",
        "(greta code) or backend = \"brms\" (Stan code via ",
        "brms::make_stancode()). Under backend = \"",
        backend,
        "\" ",
        "the inspect-then-fit token would need an INLA-side code ",
        "slot, queued for a subsequent release. Pass ",
        "backend = \"greta\" / \"brms\", or drop review_code."
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
  data_name <- if (spec_is_greta_native) {
    "<greta-native>"
  } else {
    deparse(substitute(data))
  }

  # Detect "all defaults" -- user supplied neither an fb_prior() nor
  # an explicit prior_vc_sd. In that case the v0.1 default fires:
  # build a bounded-uniform default keyed to the response scale (per
  # family / link) for sigma + every random group surfaced by the IR,
  # and emit the one-time announcement message. The uniform default
  # supersedes the earlier PC default.
  default_prior_active <- is.null(prior) && missing(prior_vc_sd)

  # Build the flexyBayes intermediate representation (IR). The universal
  # entry detects the grammar from the call shape -- ASReml
  # fixed/random/rcov, a brms-style bar-grouped formula, or (reserved) a
  # greta_model -- and routes to the matching ingest adapter. ASReml
  # ingest is byte-identical to the historical direct fb_from_asreml()
  # call; `syntax = ` forces a grammar. See .build_ir_polymorphic() in
  # R/fb.R.
  fb <- .build_ir_polymorphic(
    fixed = fixed,
    random = random,
    rcov = rcov,
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

  # A native greta model graph is fit directly by
  # greta::mcmc() (no shared emit path, no formula-specific machinery).
  # .dispatch_native_greta() pins the backend (greta only), refuses the
  # code-inspection / plan modes that do not apply to a user-built graph,
  # and assembles the flexybayes_direct_greta result. Everything below
  # (default-prior expansion, plan, review, the formula emit dispatch) is
  # formula-specific and bypassed.
  if (identical(fb$source, "greta")) {
    fit <- .dispatch_native_greta(
      fb = fb,
      backend = backend,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      return_code = return_code,
      review_code = review_code,
      plan = plan,
      the_call = the_call
    )
    .fb_warn_poor_convergence(fit)
    return(fit)
  }

  if (default_prior_active) {
    # simple_slope_uncor contributes BOTH the standard
    # intercept-variance group name AND a slope-variance group name
    # (paste0(slope_var, "_", grouping_factor)) so the existing
    # per-group uniform-on-SD machinery covers both hyperparameters.
    grp_names <- character(0)
    for (t in fb$random_terms) {
      if (is.null(t$var)) {
        next
      }
      if (t$type %in% c("simple", "ide", "id")) {
        grp_names <- c(grp_names, t$var)
      }
      if (identical(t$type, "simple_slope_uncor")) {
        if (isTRUE(t$with_intercept)) {
          grp_names <- c(grp_names, t$var)
        }
        if (!is.null(t$slope_var) && nzchar(t$slope_var)) {
          grp_names <- c(grp_names, paste0(t$slope_var, "_", t$var))
        }
      }
    }
    grp_names <- unique(grp_names)
    # vm() + ped() structured-cov groups
    # also receive the uniform-on-SD default.
    vm_ped_names <- vapply(
      fb$random_terms,
      function(t) {
        if (!is.null(t$var) && t$type %in% c("vm", "ped")) {
          t$var
        } else {
          NA_character_
        }
      },
      character(1)
    )
    vm_ped_names <- vm_ped_names[!is.na(vm_ped_names)]
    unif_default <- .default_uniform_prior(
      data = data,
      response = fb$response,
      family = family,
      link = link,
      random_groups = grp_names,
      vm_ped_groups = vm_ped_names
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
    # The emit engine for the review code: greta source on the greta
    # path; Stan source (brms::make_stancode()) on the brms path.
    review_emit_backend <- if (identical(backend, "brms")) "brms" else "greta"
    if (
      review_emit_backend == "greta" &&
        !requireNamespace("greta", quietly = TRUE)
    ) {
      stop(
        "Package 'greta' is required to generate the review-code ",
        "string. Install with:\n",
        "  install.packages('greta'); greta::install_greta_deps()",
        call. = FALSE
      )
    }
    if (
      review_emit_backend == "brms" &&
        !requireNamespace("brms", quietly = TRUE)
    ) {
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

    # Emit the review code from the engine the user
    # pinned -- greta source (codegen) or Stan source (make_stancode()).
    # Both emit paths derive the model from the IR on return_code = TRUE,
    # so the fixed / random / rcov args (display-only) pass through as-is.
    emit_review_fn <- if (review_emit_backend == "brms") {
      emit_brms
    } else {
      emit_greta
    }
    review_code_str <- emit_review_fn(
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
      rcov = rcov,
      family = family,
      link = link,
      data_name = data_name
    )

    return(.new_flexybayes_review(
      code = review_code_str,
      backend = if (review_emit_backend == "brms") {
        "stan_via_brms"
      } else {
        "greta"
      },
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
        rcov = rcov,
        family = family,
        link = link,
        data_name = data_name
      )
    ))
  }

  # Backend dispatch lives in R/dispatch.R as the shared helper
  # `.dispatch_backend()` (lifted so fb_brms() drives
  # the same routing). Semantics: backend = "greta"
  # skips the gate; backend = "inla" calls lgm_gate() and raises
  # the refusal on non-LGM; backend = "auto" gate-then-route with
  # silenceable fall-back notes.
  fit <- .dispatch_backend(
    fb = fb,
    data = data,
    backend = backend,
    known_matrices = known_matrices,
    weights = weights,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    verbose = verbose,
    mcmc_verbose = mcmc_verbose,
    return_code = return_code,
    the_call = the_call,
    fixed = fixed,
    random = random,
    rcov = rcov,
    family = family,
    link = link,
    data_name = data_name,
    aggregate = aggregate
  )
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
