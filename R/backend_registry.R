# Backend registry (v0.5.0 backend-axis recovery).
#
# The fifth closed-vocabulary registry, alongside representation
# (R/representation_registry.R), approximation (R/approximation_registry.R),
# backend-independence (R/backend_independence_registry.R), and refusal
# (R/refusal_taxonomy.R). It models *backend* as a first-class axis
# rather than a value hard-coded into each verb's match.arg vocabulary
# and into .routing_policy_table(). The package's whole architecture is
# built on "extend by registration, not by API growth" but backend was
# the one concept still hard-coded; this registry supplies the
# extension point that philosophy always implied.
#
# Same shape as the other four registries: an environment allocated at
# namespace load, filled by .register_backend() during .onLoad(), locked
# immediately after by .lock_backend_registry(). parent = emptyenv() so a
# key miss cannot fall through to the package namespace.
#
# DISPATCH ROLE (v0.5.0).
# The registry is the single source of truth for backend *facts*:
# legitimate names, per-backend capability_predicate, auto-membership
# (default_in_auto), availability (available_pkg), and the emit
# entry-point (engine). Dispatch (R/dispatch.R) CONSUMES these facts ---
# backend availability via .available_backend_names(); per-backend
# capability via .backend_can_fit(); the emit function via the `engine`
# field. The brms structured-cov gate and the INLA lgm_gate reconcile with
# the registry predicates (the predicates delegate to the same authorities,
# so refusal semantics cannot drift).
#
# What stays EXPLICIT, by design: the per-paradigm routing ORDER and
# fallback policy (INLA Laplace fast-path; brms/Stan compile-latency
# opt-out from auto; greta universal HMC fallback) is genuine policy, not
# mechanical iteration, so dispatch keeps it as code rather than deriving
# it. Consequence (the honest extensibility claim): a new same-paradigm
# formula engine is added by registration plus an emit hookup; a new
# inference paradigm additionally needs an orchestration step.
#
# NAMING. The proposed user-facing "brms" -> "stan" engine rename was
# REVERSED (2026-05-31): brms is retained as the engine label so the
# backend axis stays consistent (greta / inla / brms are all front-ends)
# and gretaR is a natural front-end sibling; the Stan/HMC sampler is
# recorded as the paradigm attribute (paradigm = "hmc_nuts") so
# triangulate() / the backend-independence registry grade pairs honestly.
# `rename_to` is therefore NA for every backend.
#
# Internal -- not exported.

# --- closed status vocabulary ------------------------------------- #

# A backend entry is one of three lifecycle states. `active` backends are
# reachable today; `dormant` backends have a provisioned but inactive
# slot (gretaR -- see R/gretaR_slot.R) and refuse at dispatch until
# activated; `reserved` is documentation-only and not registered here
# (NIMBLE is the reserved next slot; it has no slot yet, so
# registering it would be a stub -- it is named in this comment
# instead). The vocabulary is closed; a new state requires a
# deliberate amendment.
.BACKEND_STATUS_VOCABULARY <- c("active", "dormant")

# The grammars a backend can ingest-and-fit. A formula model (asreml or
# brms/lme4 dialect) lowers to the shared fb_terms IR and can target any
# formula-capable engine subject to capability; a native greta model
# graph (the "greta" grammar) is greta-only by construction. Closed
# vocabulary mirrored from the ingest-adapter family (fb_from_asreml /
# fb_from_brms / fb_from_greta).
.BACKEND_GRAMMAR_VOCABULARY <- c("asreml", "brms", "greta")

# --- container ---------------------------------------------------- #

.backend_registry <- new.env(parent = emptyenv())

# --- registration helper ------------------------------------------ #

# .register_backend() --- the one-shot registration call. Validates the
# status + grammar vocabularies (an unknown value is a hard error naming
# the offender), refuses a duplicate, refuses once the registry is
# locked. Field schema:
#   name            chr(1)  canonical backend name (the routing label).
#   status          chr(1)  one of .BACKEND_STATUS_VOCABULARY.
#   engine          chr(1)  emit entry-point function name, or NA for a
#                           dormant backend (resolved lazily at dispatch;
#                           stored as a name string, not a closure, so a
#                           load-time reference cannot fail).
#   grammars        chr     subset of .BACKEND_GRAMMAR_VOCABULARY.
#   paradigm        chr(1)  inference-paradigm label (aligns with the
#                           backend-independence registry's paradigms).
#   available_pkg   chr(1)  the R package whose presence makes the
#                           backend usable, or NA (dormant / always-on).
#   default_in_auto lgl(1)  whether `auto` considers it by default.
#   capability_predicate function(fb) -> TRUE (capable) or a single
#                           reason-code string naming why the engine
#                           cannot represent the model. Default permissive.
#   rename_to       chr(1)  reserved for a future user-facing rename; NA
#                           for every backend today (the proposed
#                           brms -> stan rename was reversed 2026-05-31).
#   registered_in_adr chr(1)
.register_backend <- function(
  name,
  status,
  engine,
  grammars,
  paradigm,
  available_pkg,
  default_in_auto,
  capability_predicate = function(fb) TRUE,
  rename_to = NA_character_,
  registered_in_adr = "0031"
) {
  # Validate inputs before checking mutability, so a caller passing an
  # unknown status / grammar sees the actionable vocabulary error rather
  # than a lock error even on a loaded (locked) package.
  if (
    !is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)
  ) {
    stop(
      ".register_backend(): `name` must be a single non-empty string.",
      call. = FALSE
    )
  }
  if (
    !identical(length(status), 1L) ||
      !status %in% .BACKEND_STATUS_VOCABULARY
  ) {
    stop(
      ".register_backend(): unknown status '",
      status,
      "'. The closed vocabulary is ",
      paste(.BACKEND_STATUS_VOCABULARY, collapse = ", "),
      "; expanding it requires an ADR 0031 amendment.",
      call. = FALSE
    )
  }
  unknown_g <- setdiff(grammars, .BACKEND_GRAMMAR_VOCABULARY)
  if (length(unknown_g) > 0L) {
    stop(
      ".register_backend(): unknown grammar ",
      paste0("'", unknown_g, "'", collapse = ", "),
      ". The closed vocabulary is ",
      paste(.BACKEND_GRAMMAR_VOCABULARY, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (
    !is.logical(default_in_auto) ||
      length(default_in_auto) != 1L ||
      is.na(default_in_auto)
  ) {
    stop(
      ".register_backend(): `default_in_auto` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (environmentIsLocked(.backend_registry)) {
    stop(
      ".register_backend(): the backend registry is locked; ",
      "backends must be registered before .lock_backend_registry() ",
      "fires at end of .onLoad().",
      call. = FALSE
    )
  }

  if (exists(name, envir = .backend_registry, inherits = FALSE)) {
    stop(
      ".register_backend(): backend '",
      name,
      "' is already ",
      "registered. Backend vocabulary is append-only.",
      call. = FALSE
    )
  }

  assign(
    name,
    list(
      name = name,
      status = status,
      engine = engine,
      grammars = grammars,
      paradigm = paradigm,
      available_pkg = available_pkg,
      default_in_auto = default_in_auto,
      capability_predicate = capability_predicate,
      rename_to = rename_to,
      registered_in_adr = registered_in_adr
    ),
    envir = .backend_registry
  )

  invisible(NULL)
}

# --- accessors ---------------------------------------------------- #

# .lookup_backend() --- internal accessor returning the full entry, or
# NULL when the name is not registered.
.lookup_backend <- function(name) {
  if (
    !is.character(name) ||
      length(name) != 1L ||
      is.na(name) ||
      !exists(name, envir = .backend_registry, inherits = FALSE)
  ) {
    return(NULL)
  }
  get(name, envir = .backend_registry, inherits = FALSE)
}

# .registered_backend_names() --- every registered backend name (any
# status), sorted for determinism.
.registered_backend_names <- function() {
  sort(ls(envir = .backend_registry, all.names = FALSE))
}

# .available_backend_names() --- the active backends whose required
# package is installed (NA available_pkg means always-on). Dormant
# backends are excluded. This is the candidate set the universal entry
# defaults to; it is also read by the consistency test.
.available_backend_names <- function() {
  out <- character(0)
  for (nm in .registered_backend_names()) {
    e <- .lookup_backend(nm)
    if (!identical(e$status, "active")) {
      next
    }
    ok <- is.na(e$available_pkg) ||
      requireNamespace(e$available_pkg, quietly = TRUE)
    if (ok) out <- c(out, nm)
  }
  out
}

# .auto_default_backend_names() --- registered backends flagged for
# auto's default candidate set (any status; dormant ones are filtered by
# availability at dispatch). Mirrors the current dispatch candidate list.
.auto_default_backend_names <- function() {
  out <- character(0)
  for (nm in .registered_backend_names()) {
    e <- .lookup_backend(nm)
    if (isTRUE(e$default_in_auto)) out <- c(out, nm)
  }
  sort(out)
}

# --- lock --------------------------------------------------------- #

.lock_backend_registry <- function() {
  if (!environmentIsLocked(.backend_registry)) {
    lockEnvironment(.backend_registry, bindings = TRUE)
  }
  invisible(NULL)
}

# --- capability predicates ---------------------------------------- #

# Each predicate takes the fb_terms IR and returns TRUE (the engine can
# represent the model) or a single reason-code string naming why it
# cannot. The registry stores them; .backend_can_fit() is the
# dispatch-facing accessor. The predicates are the systematic
# replacement for the special-cased gates (low_rank_requires_greta,
# lgm_gate) -- they delegate to the existing authorities rather than
# duplicate them, so refusal semantics do not drift.

# greta fits every model currently in scope -- the universal fallback.
.capability_greta <- function(fb) TRUE

# inla is capable iff lgm_gate() accepts the model. The 11-rule gate is
# INLA's capability predicate; the closure delegates to it (single
# authority, no duplication).
.capability_inla <- function(fb) {
  gated <- tryCatch(lgm_gate(fb), error = function(e) NULL)
  if (is.null(gated) || is_lgm_refusal(gated)) {
    return("inla_not_lgm_feasible")
  }
  TRUE
}

# The asreml structured-covariance term types with no lossless brms /
# Stan translation (dispatch.R names fa / us / ar1; vm / ped add
# known-matrix and pedigree carriers). Closed set; extending it is a
# deliberate amendment.
.STRUCTURED_COV_TYPES <- c("vm", "ped", "fa", "us", "ar1")

# brms / Stan reaches the vm() / ped() relationship random effects via
# its native known-covariance group term, (1 | gr(var, cov = K)) -- brms
# Cholesky-factors the supplied covariance internally (the K = L L'
# decorrelation Stan fits directly), so GBLUP / pedigree BLUP become
# three-engine triangulatable. This holds only for an exact dense-able
# carrier (dense / chol / precision / pedigree sparse precision); the
# remaining asreml structured-covariance terms (fa / us / ar1), a
# block-diagonal or low-rank vm() carrier, and a low_rank_smooth
# approximation have no lossless brms / Stan translation.
.BRMS_VM_DENSEABLE_CARRIERS <- c(
  "dense", "chol", "precision", "pedigree_sparse_precision"
)

.capability_brms <- function(fb) {
  rt <- fb$random_terms %||% list()
  for (t in rt) {
    ty <- t$type %||% ""
    if (!ty %in% .STRUCTURED_COV_TYPES) {
      next
    }
    if (ty %in% c("vm", "ped")) {
      fmt <- (t$cov_representation$format %||% "dense")
      if (fmt %in% .BRMS_VM_DENSEABLE_CARRIERS) {
        next
      }
    }
    return("stan_cannot_represent_structured_cov")
  }
  if (length(.collect_approx(rt)) > 0L) {
    return("stan_cannot_represent_low_rank_approx")
  }
  TRUE
}

# .backend_can_fit() --- dispatch-facing capability check. Returns
# list(ok = TRUE) or list(ok = FALSE, reason_code = <chr>). An
# unregistered backend or one without a predicate returns ok = TRUE
# (its own dispatch-side handling owns the outcome -- e.g. the gretaR
# dormant refusal).
.backend_can_fit <- function(backend, fb) {
  e <- .lookup_backend(backend)
  pred <- if (!is.null(e)) e$capability_predicate else NULL
  if (is.null(pred)) {
    return(list(ok = TRUE))
  }
  res <- pred(fb)
  if (isTRUE(res)) {
    list(ok = TRUE)
  } else {
    list(ok = FALSE, reason_code = res)
  }
}

# .backend_emit_fn() --- dispatch-facing resolver for a backend's emit
# entry-point. The registry stores the emit function by NAME (a string,
# e.g. "emit_greta") rather than as a closure, so a load-time reference
# cannot fail; this resolves the name to the function within the package
# namespace at dispatch time. Errors loudly on an unregistered backend or
# one whose engine is NA (a dormant slot reached in error) -- the caller
# is expected to have handled dormancy upstream (.gretaR_dormant_refusal).
# This is the seam that lets a newly-registered same-paradigm engine be
# dispatched without hard-coding its emit symbol at the call site.
.backend_emit_fn <- function(name) {
  e <- .lookup_backend(name)
  if (is.null(e)) {
    stop(
      ".backend_emit_fn(): backend '",
      name,
      "' is not registered.",
      call. = FALSE
    )
  }
  if (is.na(e$engine)) {
    stop(
      ".backend_emit_fn(): backend '",
      name,
      "' has no emit engine ",
      "(status '",
      e$status,
      "'); it cannot be dispatched.",
      call. = FALSE
    )
  }
  get(e$engine, envir = asNamespace("flexyBayes"), inherits = FALSE)
}

# --- v0.5.0 population --------------------------------------------- #

# .populate_backend_registry_v050() --- registers the three active
# backends reachable today plus the dormant gretaR slot. brms is retained
# as the engine label (the brms -> stan rename was reversed 2026-05-31);
# the Stan/HMC sampler is the paradigm attribute, not the name.
# `default_in_auto` encodes the auto candidate set
# c("greta", "gretaR", "inla"): greta + inla + gretaR are auto-considered,
# brms is opt-in only (its 30--60 s Stan compile would break the auto
# fast-path promise). Dispatch reads this set via
# .available_backend_names() / .auto_default_backend_names().
.populate_backend_registry_v050 <- function() {
  .register_backend(
    name = "greta",
    status = "active",
    engine = "emit_greta",
    grammars = c("asreml", "brms", "greta"),
    paradigm = "hmc_nuts",
    available_pkg = "greta",
    default_in_auto = TRUE,
    capability_predicate = .capability_greta
  )
  .register_backend(
    name = "inla",
    status = "active",
    engine = "emit_inla",
    grammars = c("asreml", "brms"),
    paradigm = "laplace_approximation",
    available_pkg = "INLA",
    default_in_auto = TRUE,
    capability_predicate = .capability_inla
  )
  .register_backend(
    name = "brms",
    status = "active",
    engine = "emit_brms",
    grammars = c("asreml", "brms"),
    paradigm = "hmc_nuts",
    available_pkg = "brms",
    default_in_auto = FALSE,
    capability_predicate = .capability_brms
  )
  .register_backend(
    name = "gretaR",
    status = "active",
    engine = "emit_gretaR",
    grammars = c("asreml", "brms", "greta"),
    paradigm = "hmc_nuts",
    available_pkg = "gretaR",
    default_in_auto = TRUE, # an auto CANDIDATE; actual auto-selection
    # stays gated on options(flexyBayes.gretaR_
    # activated) + the lgm flag (both dormant by
    # default), so auto still picks greta/INLA
    capability_predicate = .capability_gretaR,
    registered_in_adr = "0013/0031"
  )
  # (koine, the dormant synthesis fourth-opinion slot, moved to
  # flexyBayesOrchestra in the lean-core split, 2026-06-06. The backend
  # registry locks at .onLoad and resolves engines from the flexyBayes
  # namespace, so a companion-hosted dispatchable backend is out of scope;
  # koine now ships as an informational surface (koine_status()) there.)
  invisible(NULL)
}
