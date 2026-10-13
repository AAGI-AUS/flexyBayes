# Backend-independence registry.
#
# A small controlled vocabulary plus a per-pair registry that records
# *what kind* of independence separates two engines: agreement between
# two HMC fits ("implementation" independence) underwrites a weaker claim
# than agreement between an HMC fit and a Laplace-approximation fit
# ("algorithmic" independence). The same posture
# R/approximation_registry.R takes for approximation labels and
# R/refusal_taxonomy.R takes for refusal codes: a closed vocabulary, an
# internal registry the backends populate, snapshot-gated against silent
# vocabulary drift.
#
# Independence is a property of a pair of backends, not of one backend,
# so the registry is keyed by the canonicalised (sorted) backend pair.

# --- closed axis vocabulary --------------------------------------- #

# The three axis values. Glossary:
#   algorithmic:    different inference paradigms (HMC vs Laplace vs
#                     VI vs Gibbs). Catches algorithm-specific failure
#                     modes (mixing pathologies, tail-approximation
#                     error, mode-finding vs tail-coverage tradeoffs).
#   implementation: same paradigm, different code base (different AD
#                     framework, codegen, numerical-precision regime).
#                     Catches numerical bugs, prior-translation bugs,
#                     indexing errors at codegen boundaries.
#   specification:  same likelihood, different parameterisation
#                     (centred vs non-centred RE, marginal vs
#                     conditional). No current backend pair exercises
#                     this axis alone.
.AXIS_VOCABULARY <- c("algorithmic", "implementation", "specification")

# --- container ---------------------------------------------------- #

# `.backend_independence_registry` is allocated at namespace load,
# filled by .register_pair_independence() during .onLoad(), and locked
# immediately after by .lock_backend_independence_registry(). Keyed by
# the canonical sorted-pair string so the pairs (a, b) and (b, a)
# resolve to the same record. parent = emptyenv() so a key miss cannot
# fall through to the package namespace.
.backend_independence_registry <- new.env(parent = emptyenv())

# Canonical registry key for a backend pair: sorted, "||"-joined. The
# sort enforces the symmetry contract (the order of the two backends
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
# hard error naming the offending value), canonicalises the pair, and
# refuses a duplicate or a locked registry.
.register_pair_independence <- function(
  pair,
  axes,
  paradigms,
  justification,
  registered_in_adr
) {
  # Inputs are checked before the lock, so an unknown axis reports the
  # vocabulary error even on a loaded (locked) package.
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
      ". Expanding it requires an ADR 0029 amendment.",
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
      "registry is locked. Pairs must be registered before ",
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
# not been registered).
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

# .populate_backend_independence_registry_v0400() --- registers the pair
# between the two active backends (inla, Laplace approximation
# on C, and brms, HMC on Stan via the fb_brms() surface). The axis
# assignment follows the authoritative table: a Laplace-vs-HMC pair
# differs on BOTH paradigm and code base (algorithmic + implementation).
.populate_backend_independence_registry_v0400 <- function() {
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
