# Package load hooks

.onLoad <- function(libname, pkgname) {
  # Register emmeans methods if emmeans is available, for both the
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
  # it cannot be read or mutated via getOption() / options(). Test
  # code reaches it through .reset_emit_state_for_test() exposed via
  # helper-emit-state.R.
  .init_emit_state()

  # Populate the representation registry. It is the single source of
  # truth for representation_class strings at preflight, planning and
  # display sites. .lookup_representation() and .representation_class()
  # refuse unknown names, so typos cannot propagate silently. low_rank
  # registers alongside the approximation registry.
  .populate_representation_registry_v0400()

  # Lock the representation registry. Bindings remain readable
  # indefinitely. Further .register_representation() calls after
  # this point raise, not mutate the locked environment.
  .lock_representation_registry()

  # Register the first three refusal-reason codes:
  # block_partition_incomplete, block_not_positive_definite and
  # approximate_route_not_yet_registered.
  .populate_refusal_registry_v0310()

  # Bulk-register the remaining user-facing refusal codes, including the
  # fb_cov() codes (fb_cov_type_unknown, fb_cov_missing_matrix,
  # cov_arg_not_fb_cov). Routing-decision reasons (backend_decision
  # trace) and internal aggregate-out-of-scope control-flow signals are
  # not registered, because a user cannot hit them as refusals. Every
  # refusal site routes through .fb_refusal_condition(), which gates on
  # this registry.
  .populate_refusal_registry_v0400()

  # 0.9.0 design-fidelity: refuse a structured dsum() inner, not
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

  # The streaming route's transformed-response refusal.
  .populate_refusal_registry_v01000_9000()

  # The special-form argument contract, arity, named
  # arguments, literal counts, and the malformed-dsum fall-through.
  .populate_refusal_registry_grammar_contract()

  # INLA containment, estimated correlation matrices go to
  # brms.
  .populate_refusal_registry_inla_containment()

  # Gaussian aggregation containment, Gaussian models are
  # fitted per row, and aggregated fits report no marginal likelihood.
  .populate_refusal_registry_aggregation_containment()

  # The compiled model manifest, a residual prior on a family
  # without a residual scale, a superseded prior_vc_sd, and evidence an
  # engine did not compute.
  .populate_refusal_registry_model_manifest()

  # The working response, NaN and infinite response values.
  .populate_refusal_registry_working_response()

  # Accessors on aggregated and streamed count fits, a streamed
  # fit's head sample and re-fit, and the refusals an INLA fit gives for
  # augment(), a log-density and a prediction on a non-identity link.
  .populate_refusal_registry_aggregated_accessors()

  # An INLA fit outside the bound of a stated prior.
  .populate_refusal_registry_prior_support()

  # flexybayes_stream()'s own model-shape and data-availability
  # refusals: a model column absent from the source, and a source with
  # no rows to build from or fold.
  .populate_refusal_registry_stream_shape()

  # The response-domain gate, a family whose likelihood is
  # defined over a restricted response domain (gamma / lognormal
  # positive, hurdle_gamma non-negative, beta an open unit interval) on a
  # response that violates it.
  .populate_refusal_registry_family_support()
  # The factor-analytic covariance is withdrawn from the fitting
  # route, because repeated identical calls did not agree within the
  # tolerance declared before they were measured.
  .populate_refusal_registry_factor_analytic()

  # A fit asked to run in a separate R process, callr absent,
  # and a process that ended before it returned.
  .populate_refusal_registry_process_isolation()

  # A known variance per row that is malformed, or asked of a
  # model or a route that cannot carry it.
  .populate_refusal_registry_known_variance()

  # Fixed-effect columns the rows with a response do not
  # estimate.
  .populate_refusal_registry_fixed_alias()

  # fb_refit() on a fit it cannot rebuild, and brms objects
  # handed to a refit that is not fitted with brms.
  .populate_refusal_registry_refit()

  # A model whose marginal covariance is not built.
  .populate_refusal_registry_marginal_cov()

  # A call to a penalised-selection function that left the package.
  .populate_refusal_registry_selection_removed()

  # An argument fb_plan() does not take.
  .populate_refusal_registry_plan_arguments()

  # predict(allow_new_levels = TRUE) on a structure, a level or a
  # family the new-level draw is not built for.
  .populate_refusal_registry_predict_new_levels()
  # An argument of ASReml's predict() that predict(classify = )
  # does not implement.
  .populate_refusal_registry_classify_random()
  # An INLA configuration precision that no constraint makes
  # positive definite.
  .populate_refusal_registry_inla_config()
  # Entry arguments and model variables the routes stopped on
  # untyped, and the aggregated routes' scope reasons.
  .populate_refusal_registry_entry_arguments()
  # Option flexyBayes.inla_best_of that is not a count.
  .populate_refusal_registry_inla_best_of()
  # A pooled() prior sent to INLA.
  .populate_refusal_registry_pooled_prior()

  # Lock the refusal-reason registry (v0.3.8 scaffold).
  # The lock makes the registry immutable to user code once the
  # package is loaded. .register_refusal() calls beyond this point
  # raise, not mutate the locked environment.
  .lock_refusal_registry()

  # Register the
  # first approximation scheme, low_rank_smooth, with its full
  # five-field schema and its validation_fn (.validate_low_rank_smooth,
  # R/emit_smooth_low_rank.R) and the exported validate_approximation()
  # surface: the registry invariant ("every registered scheme
  # carries a working validation procedure") is honoured because the
  # validation_fn reads the truncation metadata a real low-rank fit
  # records on fit$extras$parse_info$approx. No active engine currently
  # has an emit path for this scheme (see the
  # low_rank_smooth_unsupported refusal).
  .populate_approximation_registry_v0400()

  # Lock the approximation-scheme registry. Further
  # .register_approximation() calls beyond this point raise,
  # not mutate the locked environment. Further schemes (inla_laplace,
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
  # registry: the fifth closed-vocabulary registry, modelling backend
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
  # this point raise, not mutate the locked environment. Future
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
    ": multi-backend Bayesian mixed models (INLA / brms)\n",
    "  Development release: all exports are experimental, ",
    "and it is not on CRAN. ",
    "See system.file(\"KNOWN_ISSUES.md\", package = \"flexyBayes\")."
  )

  # Backend-readiness note: only when no inference engine is available, so the
  # message is actionable, not noisy on a normally-configured machine.
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
