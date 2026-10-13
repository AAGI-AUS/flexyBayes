# Representation registry, v0.4.0
#
# Container + registration helper + lookup accessor + lock helper
# for the canonical representation-class vocabulary that classifies
# every term entering the preflight / planning / dispatch stack by
# the shape of the design or covariance machinery it carries.
#
# Same shape as the R/refusal_taxonomy.R scaffold: an
# environment (parent = emptyenv()) populated by .register_*() at
# .onLoad() time and locked immediately after via .lock_*(). The
# registry is the single source of truth for "is this string a
# legitimate representation class?". Inline string literals at
# consumer sites in fb_preflight.R / fb_plan.R / methods_truth_display.R
# route through .representation_class() so typos refuse at runtime
# instead of silently propagating an unknown class.
#
# Vocabulary (16 entries):
#   - block_diagonal, indexed_structured_known, indexed_structured_estimate
#   - dense_cov, chol_cov, sparse_precision, pedigree_sparse_precision
#   - dense_smooth, sparse_smooth, banded_smooth
#   - indexed_random_intercept, dense_baseline (both used by the per-term
#     INLA memory estimator)
#   - indexed_fixed_numeric, indexed_fixed_factor,
#     indexed_fixed_factor_numeric (every fixed-effect preflight entry,
#     in the `indexed_fixed` category)
#   - low_rank, which registers alongside the approximation registry (the
#     low_rank_smooth scheme's truncated smooth-basis representation)
# axis_* entries are reserved for the backend-independence axis
# vocabulary.
#
# The string "unknown" at fb_preflight.R:608 is a sentinel for the
# unknown_representation = TRUE meta-state and is intentionally
# outside the registry (a missing class flag, not a class).

# --- container ---------------------------------------------------- #

# `.representation_registry` is allocated at namespace load. The
# parent environment is `emptyenv()` so symbol lookup cannot fall
# through to the package namespace and pick up an unrelated binding.
# Filled by .register_representation() during .onLoad(). Locked
# immediately after by .lock_representation_registry().
.representation_registry <- new.env(parent = emptyenv())


# --- registration helper ------------------------------------------ #

# .register_representation(): the one-shot registration call.
# Idempotent only in the sense that a duplicate `name` raises,
# not silently overwrites. The vocabulary is append-only at
# .onLoad() time and the registry then locks.
#
# Arguments
#   name               canonical machine-readable identifier. The
#                      same string consumers see at preflight /
#                      planning / display sites. Must be a single
#                      non-empty string.
#   description        one-line human-readable description. What the
#                      class names in plain language.
#   category           coarse-grained grouping. One of
#                      "dense_covariance", "factored_covariance",
#                      "sparse_covariance", "smooth_basis", or
#                      "indexed_structured". Free-form for forward
#                      categories. The v0.4.0 vocabulary uses the
#                      five listed.
#   registered_in_adr  the ADR number the entry traces to (e.g.,
#                      "ADR 0030"). Free-form provenance. Does not gate
#                      registration.
#   since_version      the flexyBayes version the entry was
#                      introduced. Useful for backward-compat audits.
#
# Returns invisible NULL. Side effect: assign to
# .representation_registry.
.register_representation <- function(
  name,
  description,
  category,
  registered_in_adr,
  since_version
) {
  .check_string_scalar(
    name,
    ".register_representation(): `name` must be a non-empty ",
    "single string."
  )
  if (!is.character(description) || length(description) != 1L) {
    stop(
      ".register_representation(): `description` must be a single ",
      "string.",
      call. = FALSE
    )
  }
  if (!is.character(category) || length(category) != 1L || !nzchar(category)) {
    stop(
      ".register_representation(): `category` must be a non-empty ",
      "single string.",
      call. = FALSE
    )
  }
  if (
    !is.character(registered_in_adr) ||
      length(registered_in_adr) != 1L
  ) {
    stop(
      ".register_representation(): `registered_in_adr` must be a ",
      "single string.",
      call. = FALSE
    )
  }
  if (
    !is.character(since_version) ||
      length(since_version) != 1L ||
      !nzchar(since_version)
  ) {
    stop(
      ".register_representation(): `since_version` must be a ",
      "non-empty single string.",
      call. = FALSE
    )
  }

  if (environmentIsLocked(.representation_registry)) {
    stop(
      ".register_representation(): representation registry is ",
      "locked. New entries must be registered before ",
      ".lock_representation_registry() fires at end of .onLoad().",
      call. = FALSE
    )
  }

  if (exists(name, envir = .representation_registry, inherits = FALSE)) {
    stop(
      ".register_representation(): name '",
      name,
      "' is already ",
      "registered. Representation vocabulary is append-only. Use a ",
      "distinct name for new classes.",
      call. = FALSE
    )
  }

  assign(
    name,
    list(
      name = name,
      description = description,
      category = category,
      registered_in_adr = registered_in_adr,
      since_version = since_version
    ),
    envir = .representation_registry
  )

  invisible(NULL)
}


# --- accessor: full entry ----------------------------------------- #

# .lookup_representation(): internal accessor returning the full
# registered entry list. Refuses with a structured error if `name`
# is not registered: typos at consumer sites surface here and not
# silently propagate an unknown class through preflight /
# planning / display.
.lookup_representation <- function(name) {
  .check_string_scalar(
    name,
    ".lookup_representation(): `name` must be a non-empty ",
    "single string."
  )
  if (!exists(name, envir = .representation_registry, inherits = FALSE)) {
    stop(
      ".lookup_representation(): '",
      name,
      "' is not a registered ",
      "representation class. Known classes: ",
      paste(
        sort(ls(envir = .representation_registry, all.names = FALSE)),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
  get(name, envir = .representation_registry, inherits = FALSE)
}


# --- accessor: validated name string ----------------------------- #

# .representation_class(): thin convenience returning the
# validated name back. The edit-light replacement for inline string
# literals at consumer sites: `representation_class = "block_diagonal"`
# becomes `representation_class = .representation_class("block_diagonal")`,
# trading a single character of verbosity for vocabulary-lock
# enforcement.
.representation_class <- function(name) {
  .lookup_representation(name)$name
}


# --- v0.4.0 registry population ---------------------------------- #

# .populate_representation_registry_v0400(), registers the 16
# v0.4.0 representation classes (15 at open plus low_rank added
# later). Called from .onLoad() *before*
# .lock_representation_registry() so the entries are in place when
# user code first observes the registry.
.populate_representation_registry_v0400 <- function() {
  .register_representation(
    name = "dense_cov",
    description = paste0(
      "Dense covariance matrix carrier. V supplied directly (full ",
      "K x K dense matrix)."
    ),
    category = "dense_covariance",
    registered_in_adr = "ADR 0030",
    since_version = "0.4.0"
  )
  .register_representation(
    name = "chol_cov",
    description = paste0(
      "Cholesky-factored covariance carrier. User supplies L with ",
      "V = L L^T."
    ),
    category = "factored_covariance",
    registered_in_adr = "ADR 0025 + ADR 0030",
    since_version = "0.4.0"
  )
  .register_representation(
    name = "sparse_precision",
    description = paste0(
      "Sparse precision matrix Q carrier. V = Q^{-1} via Cholesky ",
      "back-solve in the linear-Gaussian-model path."
    ),
    category = "sparse_covariance",
    registered_in_adr = "ADR 0025 + ADR 0030",
    since_version = "0.4.0"
  )
  .register_representation(
    name = "pedigree_sparse_precision",
    description = paste0(
      "Pedigree-derived sparse precision. ped(group, A_inv, ",
      "use_sparse_precision = TRUE) routes through the sparse-Q ",
      "engine path."
    ),
    category = "sparse_covariance",
    registered_in_adr = "ADR 0025 + ADR 0030",
    since_version = "0.3.7"
  )
  .register_representation(
    name = "block_diagonal",
    description = paste0(
      "Block-diagonal vm/ped via cov = list(V_1, ..., V_K). Each ",
      "V_k is a within-block covariance and the cross-block off-",
      "diagonal is zero."
    ),
    category = "factored_covariance",
    registered_in_adr = "ADR 0025 + ADR 0030",
    since_version = "0.3.10"
  )
  .register_representation(
    name = "dense_smooth",
    description = paste0(
      "Dense smooth basis B (typical s(x) at low-to-moderate k). ",
      "B is materialised as a dense matrix at preflight."
    ),
    category = "smooth_basis",
    registered_in_adr = "ADR 0030",
    since_version = "0.4.0"
  )
  .register_representation(
    name = "sparse_smooth",
    description = paste0(
      "Sparse smooth basis (B-spline / penalised-regression bases ",
      "with sparse structure), and B is carried as a sparseMatrix."
    ),
    category = "smooth_basis",
    registered_in_adr = "ADR 0030",
    since_version = "0.4.0"
  )
  .register_representation(
    name = "banded_smooth",
    description = paste0(
      "Banded smooth basis (e.g., random-walk / AR-style precision ",
      "smooths), and B is carried as a banded matrix structure."
    ),
    category = "smooth_basis",
    registered_in_adr = "ADR 0030",
    since_version = "0.4.0"
  )
  .register_representation(
    name = "indexed_structured_known",
    description = paste0(
      "Indexed structured term with a known cov/precision matrix. ",
      "The design carries a single integer index column and the ",
      "engine receives the matrix via known_matrices."
    ),
    category = "indexed_structured",
    registered_in_adr = "ADR 0024 + ADR 0030",
    since_version = "0.3.6"
  )
  .register_representation(
    name = "indexed_structured_estimate",
    description = paste0(
      "Indexed structured term with cov/precision parameters ",
      "estimated from data (reserved for future extensions, no ",
      "active emit path yet)."
    ),
    category = "indexed_structured",
    registered_in_adr = "ADR 0030",
    since_version = "0.4.0"
  )
  .register_representation(
    name = "indexed_random_intercept",
    description = paste0(
      "Random-intercept (1 | g). INLA memory estimator categorises ",
      "this separately from the structured family because the ",
      "design block is a single integer index and the precision is ",
      "a scalar variance hyperparameter."
    ),
    category = "indexed_structured",
    registered_in_adr = "ADR 0024 + ADR 0030",
    since_version = "0.3.6"
  )
  .register_representation(
    name = "dense_baseline",
    description = paste0(
      "Dense baseline term (default fallback for smooth / spline ",
      "blocks whose design is already dense at preflight), and the INLA ",
      "memory estimator counts dense baseline blocks via the ",
      "per-term design_memory_bytes surrogate, not a ",
      "structure-specific formula."
    ),
    category = "dense_covariance",
    registered_in_adr = "ADR 0024 + ADR 0030",
    since_version = "0.3.6"
  )
  .register_representation(
    name = "indexed_fixed_numeric",
    description = paste0(
      "Indexed fixed-effect numeric term: one double design ",
      "column plus one slope coefficient. Covers bare numeric ",
      "predictors, numeric:numeric interactions (single product ",
      "column), and expression terms (I() / arithmetic / function ",
      "calls evaluated against data)."
    ),
    category = "indexed_fixed",
    registered_in_adr = "ADR 0030",
    since_version = "0.3.6"
  )
  .register_representation(
    name = "indexed_fixed_factor",
    description = paste0(
      "Indexed fixed-effect factor term: integer level index ",
      "column plus (L - 1) dummy-contrast beta coefficients under ",
      "the treatment-coding default. NA level cardinality escalates ",
      "to the `unknown` sentinel, not silently sizing the ",
      "dummy block at zero."
    ),
    category = "indexed_fixed",
    registered_in_adr = "ADR 0030",
    since_version = "0.3.6"
  )
  .register_representation(
    name = "indexed_fixed_factor_numeric",
    description = paste0(
      "Indexed fixed-effect factor:continuous interaction: ",
      "integer factor index column + continuous column ",
      "+ (L - 1) per-level slope coefficients (reference level's ",
      "slope pinned at zero, no coefficient slot). NA level ",
      "cardinality escalates to the `unknown` sentinel."
    ),
    category = "indexed_fixed",
    registered_in_adr = "ADR 0019 + ADR 0030",
    since_version = "0.3.6"
  )
  # The deferred low-rank reservation lands here, alongside the
  # approximation registry.
  # Classifies an s() smooth whose dense mgcv basis is replaced by its
  # rank-K principal-component truncation B_K = B V_K (the
  # low_rank_smooth approximation scheme). A consuming emit path would
  # carry a design block of n x K, not n x k (no active engine
  # consumes it today, see R/emit_smooth_low_rank.R).
  .register_representation(
    name = "low_rank",
    description = paste0(
      "Rank-K low-rank truncation of a smooth basis (low_rank_smooth ",
      "scheme): the n x k dense mgcv basis B is replaced by ",
      "B_K = B V_K with V_K the top-K right singular vectors, so the ",
      "model carries K basis coefficients and prediction projects the ",
      "newdata basis through the same V_K."
    ),
    category = "smooth_basis",
    registered_in_adr = "ADR 0027 + ADR 0030",
    since_version = "0.4.0"
  )
  invisible(NULL)
}


# --- lock helper -------------------------------------------------- #

# .lock_representation_registry(): locks the environment so no
# further .register_representation() calls or user-side assign()s
# can mutate it. Called once at the end of .onLoad(). The lock is
# enforced by R's environmentIsLocked() machinery. Bindings
# registered before the lock remain readable indefinitely.
.lock_representation_registry <- function() {
  if (!environmentIsLocked(.representation_registry)) {
    lockEnvironment(.representation_registry, bindings = TRUE)
  }
  invisible(NULL)
}
