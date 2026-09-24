# Package load hooks

.onLoad <- function(libname, pkgname) {
  # Register emmeans methods if emmeans is available -- for both the
  # brms (`flexybayes`) and INLA (`flexybayes_inla`) fit classes.
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register(c("flexybayes", "flexybayes_inla"), pkgname)
  }

  # marginaleffects gates on a class allow-list before dispatching its
  # get_predict / get_coef / get_vcov methods, but exposes a custom-class
  # hook via this option (see
  # marginaleffects:::sanity_model_supported_class). Append the flexyBayes
  # classes without clobbering a user's existing setting.
  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    .register_marginaleffects_classes()
  }

  # The broom-style tidy / glance / augment methods are registered
  # statically in NAMESPACE: `generics` is now an Imports dependency, so
  # roxygen emits the `S3method()` entries directly and no runtime
  # registerS3method() shim is needed. See R/tidiers.R.

  # Initialise the emit-once message latches into their fresh state.
  # The store lives in R/emit_state.R inside a package-private env so
  # it cannot be read or mutated via getOption() / options(); test
  # code reaches it through .reset_emit_state_for_test() exposed via
  # helper-emit-state.R.
  .init_emit_state()

  # Populate the representation
  # registry with the fifteen v0.4.0-open classes (the ten plan-named
  # entries plus indexed_random_intercept + dense_baseline discovered
  # in the v0.3.10 source plus indexed_fixed_numeric /
  # indexed_fixed_factor / indexed_fixed_factor_numeric discovered
  # during the .preflight_fixed_term() producer-site sweep).
  # The registry is the single source of truth for legitimate
  # representation_class strings at preflight / planning / display
  # sites; .lookup_representation() / .representation_class() refuse
  # unknown names so typos cannot silently propagate. low_rank
  # registers alongside the approximation registry.
  .populate_representation_registry_v0400()

  # Lock the representation registry. Bindings remain readable
  # indefinitely; further .register_representation() calls after
  # this point raise rather than mutate the locked environment.
  .lock_representation_registry()

  # Register the first three
  # canonical refusal-reason codes into the v0.3.8 refusal-registry
  # scaffold --- block_partition_incomplete,
  # block_not_positive_definite, and the
  # upgraded approximate_route_not_yet_registered. This is the first
  # migration; the v0.4.0 bulk migration follows the
  # same .register_refusal() shape for the ~28 remaining sites.
  .populate_refusal_registry_v0310()

  # Bulk-register the remaining
  # user-facing refusal codes (the v0.3.10 scaffold above seeds the
  # first three). With the family gate and the parse-time spec
  # refusals and the entry-point argument guards folded in, the
  # registry held 39 user-facing
  # codes; the fb_cov() work adds three
  # (fb_cov_type_unknown, fb_cov_missing_matrix, cov_arg_not_fb_cov)
  # for a v0.4.0 total of 42. Routing-
  # decision reasons (backend_decision trace) and internal aggregate-
  # out-of-scope control-flow signals are deliberately not registered
  # --- they are not refusals a user can hit. Every refusal site now
  # routes through .fb_refusal_condition(), which gates on this
  # registry.
  .populate_refusal_registry_v0400()

  # 0.9.0 design-fidelity: refuse a structured dsum() inner rather than
  # silently reducing it to a per-region heteroscedastic variance.
  .populate_refusal_registry_v0900()

  # 0.9.0 fidelity of names: the separable AR1 respelling, the closed
  # parser vocabulary, the typed brms term refusal, and the numeric-in-a-
  # random-interaction guard.
  .populate_refusal_registry_v0900_names()

  # Field-hardening: the prior mini-language refuses at construction, the
  # scalar prior arguments arrive or are refused by name, and the INLA
  # family gate speaks the reconciled spelling.
  .populate_refusal_registry_field_hardening()

  # 0.9.3 scale strategy: the cell-count integer-cast sibling, typed
  # engine-death wraps for INLA and brms, and the per-level AR1-field
  # boundary refusal.
  .populate_refusal_registry_v093()

  # 0.10.0.9000: the streaming route's transformed-response refusal.
  .populate_refusal_registry_v01000_9000()

  # 0.10.0.9000: the special-form argument contract -- arity, named
  # arguments, literal counts, and the malformed-dsum fall-through.
  .populate_refusal_registry_grammar_contract()

  # 0.10.0.9000: the selection module's covariance-builder refusals.
  .populate_refusal_registry_selection()

  # 0.10.0.9000: the selection module's design-side refusals -- the
  # protected reduction's rank deficiency and the coordinate-semantics gate.
  .populate_refusal_registry_selection_design()

  # 0.10.0.9000: the selection module's weight, offset and missing-response
  # absorption refusals.
  .populate_refusal_registry_selection_weights()

  # 0.10.0.9000: the selection module's resampling, prediction and scoring
  # refusals, including the leakage gate.
  .populate_refusal_registry_selection_resamples()

  # 0.11.0: the joint-prior layer's grammar, resolution, targeting and
  # scaling refusals.
  .populate_refusal_registry_penalty()

  # 0.11.0: the reference model's adequacy gates and the projection
  # provider's structure and family gates.
  .populate_refusal_registry_selection_projection()

  # 0.11.0: the Bayesian selection route's absent-provider and runtime-budget
  # refusals.
  .populate_refusal_registry_selection_bayes_cv()

  # 0.11.0: the penalised-optimiser provider -- an absent provider, a
  # withdrawn penalty family, an uncertified proximal mode, and penalty
  # blocks no single glmnet problem expresses.
  .populate_refusal_registry_selection_glmnet()

  # 0.11.0: the public verb's own refusals -- an unauthorised fallback, a
  # joint shrinkage prior routed at INLA, and an incomplete decision rule.
  .populate_refusal_registry_selection_verb()

  # 2026-09-13: INLA containment -- estimated correlation matrices go to
  # brms.
  .populate_refusal_registry_inla_containment()

  # 2026-09-15: Gaussian aggregation containment -- Gaussian models are
  # fitted per row, and aggregated fits report no marginal likelihood.
  .populate_refusal_registry_aggregation_containment()

  # 2026-09-15: the compiled model manifest -- a residual prior on a family
  # without a residual scale, a superseded prior_vc_sd, and evidence an
  # engine did not compute.
  .populate_refusal_registry_model_manifest()

  # 2026-09-15: the working response -- NaN and infinite response values.
  .populate_refusal_registry_working_response()

  # 2026-09-15: accessors on aggregated and streamed count fits -- a streamed
  # fit's head sample and re-fit, and the refusals an INLA fit gives for
  # augment(), a log-density and a prediction on a non-identity link.
  .populate_refusal_registry_aggregated_accessors()

  # 2026-09-16: an INLA fit outside the bound of a stated prior.
  .populate_refusal_registry_prior_support()

  # 2026-09-16: flexybayes_stream()'s own model-shape and data-availability
  # refusals -- a model column absent from the source, and a source with
  # no rows to build from or fold.
  .populate_refusal_registry_stream_shape()

  # 2026-09-16: the response-domain gate -- a family whose likelihood is
  # defined over a restricted response domain (gamma / lognormal
  # positive, hurdle_gamma non-negative, beta an open unit interval) on a
  # response that violates it.
  .populate_refusal_registry_family_support()
  # 2026-09-16: the factor-analytic covariance is withdrawn from the fitting
  # route, because repeated identical calls did not agree within the
  # tolerance declared before they were measured.
  .populate_refusal_registry_factor_analytic()

  # 2026-09-16: a fold worker that neither completed nor refused, which is a
  # child the operating system killed rather than a fold that failed.
  .populate_refusal_registry_selection_parallel()

  # 2026-09-19: a fit asked to run in a separate R process -- callr absent,
  # and a process that ended before it returned.
  .populate_refusal_registry_process_isolation()

  # 2026-09-19: a known variance per row that is malformed, or asked of a
  # model or an fb_select() route that cannot carry it.
  .populate_refusal_registry_known_variance()

  # 2026-09-24: fixed-effect columns the rows with a response do not
  # estimate.
  .populate_refusal_registry_fixed_alias()

  # Lock the refusal-reason registry (v0.3.8 scaffold).
  # The lock makes the registry immutable to user code once the
  # package is loaded; .register_refusal() calls beyond this point
  # raise rather than mutate the locked environment.
  .lock_refusal_registry()

  # Register the
  # first approximation scheme, low_rank_smooth, with its full
  # five-field schema and its validation_fn (.validate_low_rank_smooth,
  # R/emit_smooth_low_rank.R) and the exported validate_approximation()
  # surface --- the registry invariant ("every registered scheme
  # carries a working validation procedure") is honoured because the
  # validation_fn reads the truncation metadata a real low-rank fit
  # records on fit$extras$parse_info$approx. No active engine currently
  # has an emit path for this scheme (see the
  # low_rank_smooth_unsupported refusal).
  .populate_approximation_registry_v0400()

  # Lock the approximation-scheme registry. Further
  # .register_approximation() calls beyond this point raise rather
  # than mutate the locked environment; further schemes (inla_laplace,
  # predictive_process, vb) register here in subsequent releases.
  .lock_approximation_registry()

  # Populate the backend-independence
  # registry with the pair between the two active backends
  # (inla, brms) and its independence axes
  # (algorithmic / implementation / specification).
  .populate_backend_independence_registry_v0400()

  # Lock the backend-independence registry. The closed axis vocabulary
  # and the per-pair records are immutable to user code once loaded.
  .lock_backend_independence_registry()

  # Populate the backend
  # registry --- the fifth closed-vocabulary registry, modelling backend
  # as a first-class axis.
  # Registers the two active backends (inla, brms). The registry is the
  # single source of truth for backend facts and dispatch CONSUMES them:
  # availability via .available_backend_names(), capability via
  # .backend_can_fit(), the emit entry-point via .backend_emit_fn(). The
  # per-paradigm routing order / fallback policy stays explicit in
  # dispatch by design. test-backend-registry.R guards consistency
  # against the match.arg vocabularies and .routing_policy_table().
  .populate_backend_registry_v050()

  # Lock the backend registry. Further .register_backend() calls beyond
  # this point raise rather than mutate the locked environment; future
  # engines (NIMBLE) register here in later releases.
  .lock_backend_registry()
}

# Append the flexyBayes fit classes to marginaleffects' custom-class
# allow-list. Safe to call repeatedly and safe on any pre-existing user
# value: the new option is the set union of whatever is already registered
# with the flexyBayes classes, so it never drops a user's own classes
# (non-clobbering) and a second call adds nothing (idempotent). Factored
# out of .onLoad() so the contract is unit-testable
# (test-onload-options.R).
.register_marginaleffects_classes <- function(
  classes = c("flexybayes", "flexybayes_inla")
) {
  cur <- getOption("marginaleffects_model_classes", default = NULL)
  options(
    marginaleffects_model_classes = union(as.character(cur), classes)
  )
  invisible(getOption("marginaleffects_model_classes"))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "flexyBayes ",
    utils::packageVersion(pkgname),
    " -- multi-backend Bayesian mixed models (INLA / brms)\n",
    "  Development release: all exports experimental; not on CRAN. ",
    "See system.file(\"KNOWN_ISSUES.md\", package = \"flexyBayes\")."
  )

  # Backend-readiness note: only when no inference engine is available, so the
  # message is actionable rather than noisy on a normally-configured machine.
  have_engine <- requireNamespace("INLA", quietly = TRUE) ||
    requireNamespace("brms", quietly = TRUE)
  if (!have_engine) {
    packageStartupMessage(
      "  Note: no inference backend is installed. Install at least one of:\n",
      "    install.packages('INLA', repos = c(getOption('repos'),\n",
      "      INLA = 'https://inla.r-inla-download.org/R/stable'))\n",
      "    install.packages('brms')"
    )
  }
}
