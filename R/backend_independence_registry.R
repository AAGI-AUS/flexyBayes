# Backend-independence registry.
#
# A small controlled vocabulary plus a per-pair registry that lets
# triangulate() label *what kind* of convergence claim a cross-engine
# agreement underwrites, not just the numerical metric. A small
# Wasserstein distance between two HMC fits ("implementation"
# independence) underwrites a weaker claim than the same distance
# between an HMC fit and a Laplace-approximation fit ("algorithmic"
# independence). The label is the report's honesty machinery -- the
# same posture R/approximation_registry.R takes for
# approximation labels and R/refusal_taxonomy.R takes for
# refusal codes: a closed vocabulary, an internal registry the backends
# populate, surfaced on every report, snapshot-gated against silent
# vocabulary drift.
#
# Independence is a *pair* property, not a *backend* property -- "greta
# is HMC" tells you nothing about which other backend it is HMC-against
# -- so the registry is keyed by the canonicalised (sorted) backend
# pair. The vocabulary is closed; adding a fourth axis or splitting one
# of the three requires a deliberate vocabulary amendment.

# --- closed axis vocabulary --------------------------------------- #

# The three controlled axis values. Closed;
# expansion requires a deliberate vocabulary amendment. Glossary:
#   algorithmic    -- different inference paradigms (HMC vs Laplace vs
#                     VI vs Gibbs). Catches algorithm-specific failure
#                     modes (mixing pathologies, tail-approximation
#                     error, mode-finding vs tail-coverage tradeoffs).
#   implementation -- same paradigm, different code base (different AD
#                     framework, codegen, numerical-precision regime).
#                     Catches numerical bugs, prior-translation bugs,
#                     indexing errors at codegen boundaries.
#   specification  -- same likelihood, different parameterisation
#                     (centred vs non-centred RE, marginal vs
#                     conditional). No current backend pair exercises
#                     this axis alone; it is the placeholder for future
#                     within-backend reparameterisation cross-checks.
.AXIS_VOCABULARY <- c("algorithmic", "implementation", "specification")

# --- container ---------------------------------------------------- #

# `.backend_independence_registry` is allocated at namespace load,
# filled by .register_pair_independence() during .onLoad(), and locked
# immediately after by .lock_backend_independence_registry(). Keyed by
# the canonical sorted-pair string so triangulate(a, b) and
# triangulate(b, a) resolve to the same record. parent = emptyenv() so
# a key miss cannot fall through to the package namespace.
.backend_independence_registry <- new.env(parent = emptyenv())

# Canonical registry key for a backend pair: sorted, "||"-joined. The
# sort enforces the symmetry contract (order of triangulate() arguments
# does not change the lookup).
.pair_key <- function(pair) {
  if (
    !is.character(pair) ||
      length(pair) != 2L ||
      anyNA(pair) ||
      any(!nzchar(pair))
  ) {
    stop(
      ".pair_key(): `pair` must be two non-empty backend-name ",
      "strings.",
      call. = FALSE
    )
  }
  paste(sort(pair), collapse = "||")
}

# --- registration helper ------------------------------------------ #

# .register_pair_independence() --- the one-shot registration call.
# Validates `axes` against the closed vocabulary (an unknown axis is a
# hard error naming the offending value); canonicalises the pair;
# refuses a duplicate; refuses once the registry is locked. Backend
# additions register their independence axes against every existing
# backend at the same call site that introduces the new backend label.
.register_pair_independence <- function(
  pair,
  axes,
  paradigms,
  justification,
  registered_in_adr
) {
  # Validate inputs (vocabulary, shape) before checking the registry's
  # mutability state, so a caller passing an unknown axis sees the
  # vocabulary error -- the actionable problem -- rather than a lock
  # error, even on a loaded (locked) package.
  if (!is.character(axes) || length(axes) == 0L || anyNA(axes)) {
    stop(
      ".register_pair_independence(): `axes` must be a non-empty ",
      "character vector.",
      call. = FALSE
    )
  }
  unknown <- setdiff(axes, .AXIS_VOCABULARY)
  if (length(unknown) > 0L) {
    stop(
      ".register_pair_independence(): unknown axis ",
      paste0("'", unknown, "'", collapse = ", "),
      ". The closed vocabulary is ",
      paste(.AXIS_VOCABULARY, collapse = ", "),
      "; expanding it requires an ADR 0029 amendment.",
      call. = FALSE
    )
  }

  if (!is.list(paradigms) || is.null(names(paradigms))) {
    stop(
      ".register_pair_independence(): `paradigms` must be a named ",
      "list (one inference-paradigm string per backend in the pair).",
      call. = FALSE
    )
  }

  key <- .pair_key(pair)

  if (environmentIsLocked(.backend_independence_registry)) {
    stop(
      ".register_pair_independence(): the backend-independence ",
      "registry is locked; pairs must be registered before ",
      ".lock_backend_independence_registry() fires at end of ",
      ".onLoad().",
      call. = FALSE
    )
  }

  if (exists(key, envir = .backend_independence_registry, inherits = FALSE)) {
    stop(
      ".register_pair_independence(): pair '",
      key,
      "' is already ",
      "registered. Independence vocabulary is append-only.",
      call. = FALSE
    )
  }

  assign(
    key,
    list(
      pair = sort(pair),
      axes = axes,
      paradigms = paradigms,
      justification = justification,
      registered_in_adr = registered_in_adr
    ),
    envir = .backend_independence_registry
  )

  invisible(NULL)
}

# --- accessor ----------------------------------------------------- #

# .lookup_pair_independence() --- internal accessor returning the full
# pair-record, or NULL when the pair is not registered (a same-backend
# pair, or a pair involving a backend whose independence claims have
# not been registered). triangulate() consumes the NULL gracefully
# rather than refusing, so a report on an unregistered pair still
# renders (without the axis label).
.lookup_pair_independence <- function(pair) {
  key <- .pair_key(pair)
  if (!exists(key, envir = .backend_independence_registry, inherits = FALSE)) {
    return(NULL)
  }
  get(key, envir = .backend_independence_registry, inherits = FALSE)
}

# --- lock --------------------------------------------------------- #

.lock_backend_independence_registry <- function() {
  if (!environmentIsLocked(.backend_independence_registry)) {
    lockEnvironment(.backend_independence_registry, bindings = TRUE)
  }
  invisible(NULL)
}

# --- v0.4.0 population --------------------------------------------- #

# .populate_backend_independence_registry_v0400() --- registers the
# three pairs among the v0.3.x triangulatable backends (greta = HMC on
# TensorFlow; inla = Laplace approximation on C; brms = HMC on Stan via
# the fb_brms() surface). The axis assignments follow the authoritative
# table: Laplace-vs-HMC pairs differ on BOTH paradigm and
# code base (algorithmic + implementation); the two HMC backends differ
# only on code base (implementation). The stan_brms backend
# registers its own pairs when it lands at v0.4.1.
.populate_backend_independence_registry_v0400 <- function() {
  .register_pair_independence(
    pair = c("greta", "inla"),
    axes = c("algorithmic", "implementation"),
    paradigms = list(greta = "hmc_nuts", inla = "laplace_approximation"),
    justification = paste0(
      "HMC (greta on TensorFlow) versus Laplace approximation (INLA on ",
      "C): different inference paradigms and different code bases."
    ),
    registered_in_adr = "0029"
  )
  .register_pair_independence(
    pair = c("greta", "brms"),
    axes = "implementation",
    paradigms = list(greta = "hmc_nuts", brms = "hmc_nuts"),
    justification = paste0(
      "HMC (greta on TensorFlow) versus HMC (brms on Stan): the same ",
      "inference paradigm through different code bases (different AD ",
      "framework, codegen, and numerical regime)."
    ),
    registered_in_adr = "0029"
  )
  .register_pair_independence(
    pair = c("brms", "inla"),
    axes = c("algorithmic", "implementation"),
    paradigms = list(brms = "hmc_nuts", inla = "laplace_approximation"),
    justification = paste0(
      "HMC (brms on Stan) versus Laplace approximation (INLA on C): ",
      "different inference paradigms and different code bases."
    ),
    registered_in_adr = "0029"
  )
  invisible(NULL)
}
