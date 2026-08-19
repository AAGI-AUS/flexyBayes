# Package load hooks

.onLoad <- function(libname, pkgname) {
  # Register emmeans methods if emmeans is available -- for both the
  # greta (`flexybayes`) and INLA (`flexybayes_inla`) fit classes.
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

  # Borrow greta's as.greta_array methods into this namespace, for the
  # reasons set out at .register_greta_shim() below. Two registration
  # points, because loadedness alone is order-dependent:
  # isNamespaceLoaded() answers for the instant .onLoad() runs, and a
  # user who loads greta afterwards -- the usual order, since flexyBayes
  # attaches first -- would never get the shim. The hook covers that
  # direction, the direct call covers greta already being loaded, and
  # registerS3method() is idempotent so a session hitting both is fine.
  if (isNamespaceLoaded("greta")) {
    .register_greta_shim(pkgname)
  }
  setHook(
    packageEvent("greta", "onLoad"),
    function(...) .register_greta_shim(pkgname),
    action = "append"
  )

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

  # gretaR backend activation: gretaR-specific refusal codes.
  .populate_refusal_registry_gretaR()

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

  # Lock the refusal-reason registry (v0.3.8 scaffold).
  # The lock makes the registry immutable to user code once the
  # package is loaded; .register_refusal() calls beyond this point
  # raise rather than mutate the locked environment.
  .lock_refusal_registry()

  # Register the
  # first approximation scheme, low_rank_smooth, with its full
  # five-field schema and its validation_fn (.validate_low_rank_smooth,
  # R/emit_smooth_low_rank.R). The scheme lands together with its emit
  # path (R/codegen.R low-rank branch), its predict-side projection
  # (R/predict_kernel.R), and the exported validate_approximation()
  # surface --- the registry invariant ("every registered scheme
  # carries a working validation procedure") is honoured because the
  # validation_fn reads the truncation metadata a real low-rank fit
  # records on fit$extras$parse_info$approx.
  .populate_approximation_registry_v0400()

  # Lock the approximation-scheme registry. Further
  # .register_approximation() calls beyond this point raise rather
  # than mutate the locked environment; further schemes (inla_laplace,
  # predictive_process, vb) register here in subsequent releases.
  .lock_approximation_registry()

  # Populate the backend-independence
  # registry with the three pairs among the v0.3.x triangulatable
  # backends (greta, inla, brms). triangulate() consumes this to label
  # each cross-engine comparison with its independence axis
  # (algorithmic / implementation / specification). The stan_brms
  # backend registers its own pairs when it lands at v0.4.1.
  .populate_backend_independence_registry_v0400()

  # Lock the backend-independence registry. The closed axis vocabulary
  # and the per-pair records are immutable to user code once loaded.
  .lock_backend_independence_registry()

  # Populate the backend
  # registry --- the fifth closed-vocabulary registry, modelling backend
  # as a first-class axis.
  # Registers the three active backends (greta, inla, brms) plus the
  # dormant gretaR slot. The registry is the single source of
  # truth for backend facts and dispatch CONSUMES them: availability via
  # .available_backend_names(), capability via .backend_can_fit(), the
  # emit entry-point via .backend_emit_fn(). The per-paradigm routing
  # order / fallback policy stays explicit in dispatch by design.
  # test-backend-registry.R guards consistency
  # against the match.arg vocabularies and .routing_policy_table().
  .populate_backend_registry_v050()

  # Lock the backend registry. Further .register_backend() calls beyond
  # this point raise rather than mutate the locked environment; future
  # engines (gretaR activation, NIMBLE) register here in later releases.
  .lock_backend_registry()
}

# Append the flexyBayes fit classes to marginaleffects' custom-class
# allow-list. Safe to call repeatedly and safe on any pre-existing user
# value: the new option is the set union of whatever is already registered
# with the flexyBayes classes, so it never drops a user's own classes
# (non-clobbering) and a second call adds nothing (idempotent). Factored
# out of .onLoad() so the contract is unit-testable
# (test-onload-options.R).
# Borrow greta's own as.greta_array methods into the flexyBayes namespace.
#
# The aggregated-greta emit calls into greta from a package-namespace
# context, where greta's internal check_dims() -> as.greta_array() dispatch
# routes a greta_array operand to as.greta_array.matrix rather than to
# as.greta_array.greta_array, and reports a spurious missing/infinite-values
# error on a valid array. Registering greta's methods locally fixes the
# dispatch.
#
# The guard on the call sites is loadedness, never requireNamespace():
# requireNamespace("greta") LOADS greta, and with it reticulate, tensorflow
# and the Python bridge, on every library(flexyBayes) -- paying that cost,
# and putting greta in every user's sessionInfo(), for a path dispatch
# cannot reach while greta is quarantined. Registering only for a user who
# has greta loaded keeps the shim available for the re-entry candidate
# without pulling the stack in.
#
# Idempotent, and silent on any failure: an absent or restructured greta
# is not a reason for library(flexyBayes) to fail.
.register_greta_shim <- function(pkgname = "flexyBayes") {
  # asNamespace() loads the namespace when it is not already loaded, so
  # this guard is what keeps the helper from doing the very thing the
  # registration policy forbids. It holds for the hook route too: an
  # onLoad user hook fires after loadNamespace() has registered the
  # namespace, so greta is loaded by the time the hook runs.
  if (!isNamespaceLoaded("greta")) {
    return(invisible(NULL))
  }
  tryCatch(
    {
      greta_ns <- asNamespace("greta")
      flx_ns <- asNamespace(pkgname)
      for (cls in c(
        "greta_array",
        "matrix",
        "numeric",
        "logical",
        "integer",
        "data.frame",
        "array"
      )) {
        m_name <- paste0("as.greta_array.", cls)
        if (exists(m_name, envir = greta_ns, inherits = FALSE)) {
          registerS3method(
            "as.greta_array",
            cls,
            get(m_name, envir = greta_ns, inherits = FALSE),
            envir = flx_ns
          )
        }
      }
    },
    error = function(e) NULL
  )
  invisible(NULL)
}

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
    " -- multi-backend Bayesian mixed models (INLA / brms)",
    " with cross-engine triangulation\n",
    "  Development release: all exports experimental; not on CRAN. ",
    "See system.file(\"KNOWN_ISSUES.md\", package = \"flexyBayes\")."
  )

  # Backend-readiness note: only when no inference engine is available, so the
  # message is actionable rather than noisy on a normally-configured machine.
  # greta is excluded from this check: it is registered but quarantined
  # (lifecycle state, not an installation gap), so its presence does not
  # make an engine available.
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
