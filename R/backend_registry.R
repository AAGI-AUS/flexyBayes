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
# it. Consequence, and the limit of the extensibility claim: a new
# same-paradigm
# formula engine is added by registration plus an emit hookup; a new
# inference paradigm additionally needs an orchestration step.
#
# NAMING. The proposed user-facing "brms" -> "stan" engine rename was
# REVERSED (2026-05-31): brms is retained as the engine label so the
# backend axis stays consistent (greta / inla / brms are all front-ends)
# and gretaR is a natural front-end sibling; the Stan/HMC sampler is
# recorded as the paradigm attribute (paradigm = "hmc_nuts") so
# triangulate() / the backend-independence registry grade pairs on what
# the engines actually share.
# `rename_to` is therefore NA for every backend.
#
# Internal -- not exported.

# --- closed status vocabulary ------------------------------------- #

# A backend entry is one of these lifecycle states. `active` backends are
# reachable today; `dormant` backends have a provisioned but inactive
# slot and refuse at dispatch until activated; `quarantined` backends were
# active and have been deliberately withdrawn as fitting engines (greta /
# gretaR, 2026-07-24 reshape R1) -- their descriptor + emit code are
# RETAINED as a re-entry candidate but dispatch refuses them and `auto`
# never selects them, and re-entry is repair + conform to the backend
# conformance battery, never a bare re-add (§4.1). `reserved` is
# documentation-only and not registered here (NIMBLE is the reserved next
# slot; it has no slot yet, so registering it would be a stub -- it is
# named in this comment instead). The vocabulary is closed; a new state
# requires a deliberate amendment (this file adds `quarantined`, 2026-07-24).
.BACKEND_STATUS_VOCABULARY <- c("active", "dormant", "quarantined")

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

# .backend_is_quarantined() --- TRUE iff `name` is a registered backend
# whose lifecycle status is `quarantined` (a deliberately-withdrawn fitting
# engine; dispatch refuses it and `auto` never selects it -- §4.1).
# Consulted by dispatch (the explicit-request refusal) and the test
# skip-guards (a greta-fitting test skips as a re-entry guard, not fails).
.backend_is_quarantined <- function(name) {
  e <- .lookup_backend(name)
  !is.null(e) && identical(e$status, "quarantined")
}

# .auto_default_backend_names() --- registered backends flagged for
# auto's default candidate set AND currently `active`. A non-active
# (dormant / quarantined) backend can never be an auto candidate, so the
# invariant "auto never selects a quarantined engine" is structural here,
# not merely a consequence of each descriptor's default_in_auto flag.
.auto_default_backend_names <- function() {
  out <- character(0)
  for (nm in .registered_backend_names()) {
    e <- .lookup_backend(nm)
    if (isTRUE(e$default_in_auto) && identical(e$status, "active")) {
      out <- c(out, nm)
    }
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

# The residual-term types .fb_to_brms_formula() actually lowers:
#
#   units      the iid residual, which brms carries implicitly as the family
#              scale parameter (sigma).
#   at_units   a heterogeneous residual sectioned by a factor -- ASReml's
#              dsum(~ units | f) and the at(f):units spelling -- lowered to
#              distributional regression on sigma, `sigma ~ 0 + f`. Validated
#              against ASReml's dsum on the same data to within 4.8% on the
#              per-level variances (posterior mean against REML), largest on
#              the smallest of them.
#
# Every structured residual form -- ar1, separable ar1_spatial, and any
# future correlated residual -- has no emit path at all: the reconstructed
# brms formula rebuilds the fixed and random blocks and would drop the
# residual block silently. Naming the supported set here (rather than the
# unsupported one) keeps a newly parsed residual type refused by default
# instead of silently dropped.
.BRMS_RESIDUAL_TYPE_ALLOWLIST <- c("units", "at_units")

.capability_brms <- function(fb) {
  rt <- fb$random_terms %||% list()
  for (t in rt) {
    ty <- t$type %||% ""
    # A random-side AR1 field -- ar1(t) or the separable ar1(row):ar1(col)
    # -- is a Kronecker autoregressive precision that only the INLA emit
    # builds. Naming it here rather than leaving it to the generic
    # structured-covariance code matters on the auto path: without this
    # branch the separable form was invisible to the predicate (only the
    # 1-D `ar1` is in .STRUCTURED_COV_TYPES) and auto would route a failed
    # INLA spatial fit to brms, which reconstructs the fixed and random
    # blocks and would drop the field.
    if (ty %in% c("ar1", "ar1_spatial")) {
      return("stan_cannot_represent_ar1_field")
    }
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
  for (t in fb$residual_terms %||% list()) {
    ty <- t$type %||% ""
    if (ty %in% .BRMS_RESIDUAL_TYPE_ALLOWLIST) {
      next
    }
    if (ty %in% c("ar1", "ar1_spatial")) {
      return("stan_cannot_represent_ar1_residual")
    }
    return("stan_cannot_represent_structured_residual")
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
  # greta: QUARANTINED 2026-07-24 (reshape R1). Descriptor + emit code are
  # RETAINED as a re-entry candidate; dispatch refuses it and `auto` never
  # selects it. Re-entry is repair + conform (§4.1), never a bare re-add.
  .register_backend(
    name = "greta",
    status = "quarantined",
    engine = "emit_greta",
    grammars = c("asreml", "brms", "greta"),
    paradigm = "hmc_nuts",
    available_pkg = "greta",
    default_in_auto = FALSE,
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
  # gretaR: QUARANTINED 2026-07-24 (reshape R1). Registered `active` with a
  # runtime option-gate before the reshape; now a retained re-entry
  # descriptor like greta -- dispatch refuses it, `auto` never selects it.
  .register_backend(
    name = "gretaR",
    status = "quarantined",
    engine = "emit_gretaR",
    grammars = c("asreml", "brms", "greta"),
    paradigm = "hmc_nuts",
    available_pkg = "gretaR",
    default_in_auto = FALSE,
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
