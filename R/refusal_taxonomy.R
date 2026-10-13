# Refusal registry: the canonical refusal-reason vocabulary
#
# Container, registration helper, and the user-facing fb_refusals()
# accessor. Populated: 136 codes registered across versioned
# .populate_refusal_registry_*() calls in .onLoad(), against 139 distinct
# reason_code slugs thrown in R/.
#
# The header here described this file as "scaffold only at v0.3.8, empty,
# full population lands at v0.4.0" until 2026-09-12, six minor versions
# after that had happened, in a file 2,300 lines long. A registry that
# understates itself is the same claim-surface drift the registry exists
# to prevent elsewhere, so the count above is stated and the test at
# test-refusal-registry-scaffold.R asserts it.
#
# Rationale: the deliberate choice here is to ship the container
# without the migration risk of touching ~30 existing refusal sites
# in one release.
#
# The registry is an environment (not a list) so .register_refusal()
# can write idempotently without copy-on-modify cost. Locked at the
# end of .onLoad() so user-side assign() raises a contract violation.

# --- container ---------------------------------------------------- #

# `.refusal_registry` is allocated at namespace load. The parent
# environment is `emptyenv()` so symbol lookup cannot accidentally
# fall through to the package namespace and pick up an unrelated
# binding. Filled by .register_refusal() during .onLoad() (currently
# zero calls: the scaffold-only posture). Locked immediately
# after by .lock_refusal_registry().
.refusal_registry <- new.env(parent = emptyenv())

# Plug-in partitions: one environment per registering package, filled by
# fb_register_refusals() at any time after load. Each entry's binding is
# locked as it is written.
.refusal_plugin_registry <- new.env(parent = emptyenv())


# --- registration helper ------------------------------------------ #

# .register_refusal(): the one-shot registration call. Idempotent
# only in the sense that a duplicate reason_code raises instead of
# silently overwriting. Intentional refusals (the v0.4.0 migration)
# fire once per reason at .onLoad() time and the registry then
# locks.
#
# Arguments
#   reason_code        canonical machine-readable identifier. The
#                      same string that appears on the structured
#                      condition object thrown by the refusal site.
#                      Must be a single non-empty string.
#   description        one-line human-readable description. What the
#                      refusal means in plain language.
#   message_template   the template string used to render the
#                      message, conventionally uses base::sprintf()
#                      placeholders so the refusal site can fill in
#                      the data-specific slots.
#   registered_in_adr  the ADR number the refusal traces to (e.g.,
#                      "ADR 0024" for the routing-policy refusals,
#                      "ADR 0025" for the structured-cov refusals).
#                      DESIGN_DECISIONS.md indexes every ADR number.
#                      Free-form. This field documents provenance,
#                      it does not gate registration.
#   plan_field         optional name of the <fb_plan> slot that
#                      surfaces this refusal at plan time (e.g.,
#                      "rejected_routes"). NA_character_ if the
#                      refusal is not surfaced by fb_plan() (e.g.,
#                      runtime-only refusals).
#   since_version      the flexyBayes version the refusal was
#                      introduced, useful for backward-compat
#                      audits.
#
# Returns invisible NULL. Side effect: assign to .refusal_registry.
.register_refusal <- function(
  reason_code,
  description,
  message_template,
  registered_in_adr,
  plan_field = NA_character_,
  since_version
) {
  if (
    !is.character(reason_code) ||
      length(reason_code) != 1L ||
      !nzchar(reason_code)
  ) {
    stop(
      ".register_refusal(): `reason_code` must be a non-empty ",
      "single string.",
      call. = FALSE
    )
  }
  if (!is.character(description) || length(description) != 1L) {
    stop(
      ".register_refusal(): `description` must be a single ",
      "string.",
      call. = FALSE
    )
  }
  if (!is.character(message_template) || length(message_template) != 1L) {
    stop(
      ".register_refusal(): `message_template` must be a single ",
      "string.",
      call. = FALSE
    )
  }
  if (!is.character(registered_in_adr) || length(registered_in_adr) != 1L) {
    stop(
      ".register_refusal(): `registered_in_adr` must be a single ",
      "string.",
      call. = FALSE
    )
  }
  if (!is.character(plan_field) || length(plan_field) != 1L) {
    stop(
      ".register_refusal(): `plan_field` must be a single string ",
      "(use NA_character_ for unsurfaced refusals).",
      call. = FALSE
    )
  }
  if (
    !is.character(since_version) ||
      length(since_version) != 1L ||
      !nzchar(since_version)
  ) {
    stop(
      ".register_refusal(): `since_version` must be a non-empty ",
      "single string.",
      call. = FALSE
    )
  }

  if (environmentIsLocked(.refusal_registry)) {
    stop(
      ".register_refusal(): refusal registry is locked. New ",
      "entries must be registered before .lock_refusal_registry() ",
      "fires at end of .onLoad().",
      call. = FALSE
    )
  }

  if (exists(reason_code, envir = .refusal_registry, inherits = FALSE)) {
    stop(
      ".register_refusal(): reason_code '",
      reason_code,
      "' is ",
      "already registered. Refusal vocabulary is append-only. ",
      "Use a distinct reason_code for new refusals.",
      call. = FALSE
    )
  }

  assign(
    reason_code,
    list(
      reason_code = reason_code,
      description = description,
      message_template = message_template,
      registered_in_adr = registered_in_adr,
      plan_field = plan_field,
      since_version = since_version,
      package = "flexyBayes"
    ),
    envir = .refusal_registry
  )

  invisible(NULL)
}


# --- accessor ----------------------------------------------------- #

# .lookup_refusal(): internal accessor for refusal sites that
# want to consult the registry for a description or message template.
# Returns the registered entry list, or NULL if the reason_code is
# not yet registered (the v0.3.8 scaffold-only state returns NULL
# for every code). Callers must tolerate NULL, the legacy
# inline-message refusal sites continue to function untouched.
.lookup_refusal <- function(reason_code) {
  if (
    !is.character(reason_code) ||
      length(reason_code) != 1L ||
      !nzchar(reason_code)
  ) {
    return(NULL)
  }
  if (exists(reason_code, envir = .refusal_registry, inherits = FALSE)) {
    return(get(reason_code, envir = .refusal_registry, inherits = FALSE))
  }
  for (pkg in ls(.refusal_plugin_registry)) {
    partition <- get(pkg, envir = .refusal_plugin_registry, inherits = FALSE)
    if (exists(reason_code, envir = partition, inherits = FALSE)) {
      return(get(reason_code, envir = partition, inherits = FALSE))
    }
  }
  NULL
}

# Every registered entry, the core partition first, then each plug-in
# partition in package-name order.
.refusal_entries <- function() {
  core <- mget(ls(.refusal_registry), envir = .refusal_registry)
  plugins <- lapply(ls(.refusal_plugin_registry), function(pkg) {
    partition <- get(pkg, envir = .refusal_plugin_registry, inherits = FALSE)
    mget(ls(partition), envir = partition)
  })
  unname(c(core, unlist(plugins, recursive = FALSE)))
}


# --- unified refusal-condition constructor ----------------------- #

# .fb_refusal_condition(): the single constructor every
# user-facing refusal site routes through (v0.4.0).
# It assembles one canonical condition-class taxonomy so the whole
# refusal surface is assertable the same way, and it gates on the
# registry so the refusal vocabulary cannot drift silently.
#
# Class vector (most specific first):
#
#   c("flexybayes_refusal_<reason_code>",   per-code assertion handle
#     <family_class...>,                     retained legacy class(es)
#     "flexybayes_refusal",                  umbrella
#     "error", "condition")
#
# The per-code class makes
#   expect_error(class = "flexybayes_refusal_<reason_code>")
# the preferred test contract. The retained family
# class keeps every pre-existing class-based handler and assertion
# working unchanged (additive migration).
#
# Registry gate: `reason_code` MUST be registered (via
# .register_refusal() at .onLoad()). An unregistered code is a
# flexyBayes internal contract violation: the refusal does not
# silently surface under an unknown vocabulary. This is the
# "register or the gate refuses it" closing posture.
#
# The fully-rendered, site-specific `message` is supplied by the
# call site and carried verbatim: one reason_code may be raised at
# several sites with genuinely different messages, so the message
# body is not forced through a single template (the registry's
# message_template is documentation surfaced by fb_refusals(), not a
# runtime straitjacket).
#
# Arguments
#   reason_code   canonical machine-readable identifier. Must be in
#                 the locked .refusal_registry or a plug-in partition.
#   message       the fully-rendered, site-specific message string.
#   family_class  zero or more legacy family-class strings to retain
#                 below the per-code class (e.g.
#                 "flexybayes_structured_cov_refusal"). Character
#                 vector, default character(0).
#   call          the condition call carried on the object. Defaults
#                 to NULL (matching every site except the preflight
#                 wrapper, which threads the originating call through).
#   ...           extra named fields carried on the condition object
#                 (e.g. supplied, format, scheme, backend, factor)
#                 exactly as the pre-migration sites carried them.
#
# Returns the condition object, which also carries `package`, the
# partition the code is registered in. The caller raises it via stop().
.fb_refusal_condition <- function(
  reason_code,
  message,
  family_class = character(0),
  call = NULL,
  ...
) {
  if (
    !is.character(reason_code) ||
      length(reason_code) != 1L ||
      !nzchar(reason_code)
  ) {
    stop(
      ".fb_refusal_condition(): `reason_code` must be a non-empty ",
      "single string.",
      call. = FALSE
    )
  }

  entry <- .lookup_refusal(reason_code)
  if (is.null(entry)) {
    stop(
      ".fb_refusal_condition(): reason_code '",
      reason_code,
      "' is not registered in the refusal registry. Every ",
      "user-facing refusal must be registered via ",
      ".register_refusal() at .onLoad(). This is a ",
      "flexyBayes internal contract violation, not a user error.",
      call. = FALSE
    )
  }

  structure(
    c(
      list(
        message = message,
        reason_code = reason_code,
        package = entry$package,
        call = call
      ),
      list(...)
    ),
    class = c(
      paste0("flexybayes_refusal_", reason_code),
      family_class,
      "flexybayes_refusal",
      "error",
      "condition"
    )
  )
}


# --- v0.3.10 registry population ------------------------------- #

# .populate_refusal_registry_v0310(): the first migration of
# existing refusal sites into the v0.3.8 scaffold. Called from
# .onLoad() *before* .lock_refusal_registry() so the entries are
# in place when user code first observes the registry.
#
# Three reason codes register at v0.3.10:
#
#   block_partition_incomplete: blocks structural refusal.
#                                          sum of block sizes !=
#                                          grouping factor level count.
#   block_not_positive_definite: blocks structural refusal.
#                                          at least one V_k fails the
#                                          PD probe. The offending
#                                          block index is named in
#                                          the message text.
#   approximate_route_not_yet_registered: low-rank refusal-stub
#                                            upgrade naming
#                                            validate_approximation()
#                                            as the v0.4.0 forward
#                                            export.
#
# Forward pattern: the v0.4.0 bulk migration moves the ~28 remaining
# inline-message refusal sites into this scaffold using the same
# .register_refusal() call shape.
.populate_refusal_registry_v0310 <- function() {
  .register_refusal(
    reason_code = "block_partition_incomplete",
    description = paste0(
      "Block-diagonal vm/ped: the block sizes do not partition ",
      "the grouping factor's level count."
    ),
    message_template = paste0(
      "vm(%s, blocks = %s): block sizes (%s = %d) do not match ",
      "levels(`%s`) (%d levels)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.3.10"
  )
  .register_refusal(
    reason_code = "block_not_positive_definite",
    description = paste0(
      "Block-diagonal vm/ped: at least one V_k failed the ",
      "positive-definite probe."
    ),
    message_template = paste0(
      "vm(%s, blocks = %s): block %d failed the PD probe (size ",
      "%d x %d)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.3.10"
  )
  .register_refusal(
    reason_code = "approximate_route_not_yet_registered",
    description = paste0(
      "Approximate covariance / dispatch carriers refuse until an ",
      "approximation scheme is registered ",
      "(validate_approximation())."
    ),
    message_template = paste0(
      "vm(%s, low_rank_factor = %s, low_rank_scheme = '%s'): ",
      "low-rank carriers require validate_approximation()."
    ),
    registered_in_adr = "ADR 0025+0027",
    plan_field = "representation_plan",
    since_version = "0.3.10"
  )
  invisible(NULL)
}


# --- v0.4.0 bulk migration --------------------------------------- #

# .populate_refusal_registry_v0400(): the bulk migration of the
# remaining 28 user-facing refusal codes into the registry.
# Together with the 3 v0.3.10 codes the registry now
# holds the complete user-facing refusal vocabulary (31 codes). The
# routing-decision reasons (backend_decision trace) and the internal
# aggregate-out-of-scope control-flow signals are deliberately NOT
# registered: they are not refusals the user can hit, and listing
# them in fb_refusals() would misrepresent the surface.
#
# message_template here is documentation: a single representative
# rendering surfaced by fb_refusals(). It is intentionally NOT the
# runtime message: five codes are raised at multiple sites with
# different bodies, and each site keeps its verbatim message via
# .fb_refusal_condition(). plan_field follows the v0.3.10 precedent:
# the representation-level (structured covariance) refusals carry
# "representation_plan", runtime-only refusals carry NA_character_.
.populate_refusal_registry_v0400 <- function() {
  # ---- structured-covariance carriers (structured_cov.R,
  #    emit_inla.R) --------------------------------------------------
  .register_refusal(
    reason_code = "vm_redundant_specification",
    description = paste0(
      "vm()/ped(): more than one covariance carrier supplied. ",
      "Exactly one of V / chol / precision / blocks / ",
      "low_rank_factor is allowed."
    ),
    message_template = paste0(
      "%s(): exactly one covariance carrier may be supplied. Got ",
      "%s. Pick one of V / chol / precision / blocks / ",
      "low_rank_factor."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "low_rank_scheme_required",
    description = paste0(
      "vm()/ped(): low_rank_factor supplied without an explicit ",
      "low_rank_scheme naming a registered approximation."
    ),
    message_template = paste0(
      "%s(): low_rank_factor requires an explicit low_rank_scheme ",
      "naming a registered approximation. See ?validate_approximation ",
      "for the available schemes."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "chol_not_in_known_matrices",
    description = paste0(
      "vm(..., chol = ): the named Cholesky factor is absent from ",
      "known_matrices."
    ),
    message_template = paste0(
      "vm(%s, chol = %s): the Cholesky factor '%s' is not in ",
      "known_matrices. Pass it via known_matrices = list(%s = ",
      "<your Cholesky factor>)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "chol_not_square",
    description = "vm(..., chol = ): the Cholesky factor is not square.",
    message_template = paste0(
      "vm(%s, chol = %s): the Cholesky factor must be square. Got ",
      "dim %s."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "chol_not_triangular",
    description = paste0(
      "vm(..., chol = ): the Cholesky factor is not ",
      "lower-triangular."
    ),
    message_template = paste0(
      "vm(%s, chol = %s): the Cholesky factor must be ",
      "lower-triangular. If you have the upper factor U with V = ",
      "U^T U, pass t(U) or chol = t(U)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "precision_not_in_known_matrices",
    description = paste0(
      "vm(..., precision = ): the named precision matrix is absent ",
      "from known_matrices."
    ),
    message_template = paste0(
      "vm(%s, precision = %s): the precision matrix '%s' is not in ",
      "known_matrices. Pass it via known_matrices = list(%s = ",
      "<your precision matrix>)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "precision_not_square",
    description = "vm(..., precision = ): the precision matrix is not square.",
    message_template = paste0(
      "vm(%s, precision = %s): the precision matrix must be square. ",
      "Got dim %s."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "precision_not_symmetric",
    description = paste0(
      "vm(..., precision = ): the precision matrix ",
      "is not symmetric."
    ),
    message_template = paste0(
      "vm(%s, precision = %s): the precision matrix must be ",
      "symmetric. Matrix::isSymmetric(Q) returned FALSE on values ",
      "(ignoring dimnames; the level-alignment check is separate)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "precision_not_positive_definite",
    description = paste0(
      "vm(..., precision = ): the precision matrix failed the ",
      "positive-definite probe."
    ),
    message_template = paste0(
      "vm(%s, precision = %s): the precision matrix must be ",
      "positive-definite. chol(as.matrix(Q)) probe failed. Pass ",
      "options(flexyBayes.trust_pd = TRUE) if you have externally ",
      "verified Q's positive-definiteness."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "blocks_not_in_known_matrices",
    description = paste0(
      "vm(..., blocks = ): the named block list is absent from ",
      "known_matrices."
    ),
    message_template = paste0(
      "vm(%s, blocks = %s): the block list '%s' is not in ",
      "known_matrices. Pass it via known_matrices = list(%s = ",
      "list(V_1, ..., V_K))."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "blocks_not_a_list",
    description = paste0(
      "vm(..., blocks = ): the block carrier is not a base-R list ",
      "of covariance matrices."
    ),
    message_template = paste0(
      "vm(%s, blocks = %s): the block carrier must be a base-R list ",
      "of K covariance matrices. Got class '%s'."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "blocks_empty_list",
    description = "vm(..., blocks = ): the block list is empty.",
    message_template = paste0(
      "vm(%s, blocks = %s): the block list is empty. Expected at ",
      "least one covariance matrix."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "known_matrix_dim_mismatch",
    description = paste0(
      "vm(): the known matrix dimension does not match the grouping ",
      "factor's level count."
    ),
    message_template = paste0(
      "vm(%s): the known matrix '%s' has dimension %s, but the ",
      "grouping factor `%s` carries %d levels. The matrix must be ",
      "%d x %d (one row + column per level)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "known_matrix_dimnames_mismatch",
    description = "vm(): the known matrix has differing row and column names.",
    message_template = paste0(
      "vm(%s): the known matrix '%s' has different rownames and ",
      "colnames. The matrix represents the per-level covariance / ",
      "precision structure. Rownames and colnames must be identical ",
      "and must equal levels(`%s`)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "known_matrix_level_mismatch",
    description = paste0(
      "vm(): the known matrix dimnames do not match (or are ",
      "mis-ordered relative to) the grouping factor levels."
    ),
    message_template = paste0(
      "vm(%s): the known matrix '%s' carries dimnames that do not ",
      "match (or are ordered differently from) levels(`%s`)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "known_matrices_data_name_collision",
    description = paste0(
      "INLA emit: a known-matrices / blocks carrier name collides ",
      "with a data column name."
    ),
    message_template = paste0(
      "known_matrices entry/entries '%s' collide with data column ",
      "names, rename to disambiguate (INLA's data-list lookup is ",
      "name-keyed)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )

  # ---- approximation surface (parse_formula.R,
  #    validate_approximation.R, emit_smooth_low_rank.R, dispatch.R) -
  .register_refusal(
    reason_code = "low_rank_smooth_unsupported",
    description = paste0(
      "A smooth requesting the low_rank_smooth approximation. flexyBayes ",
      "has no route for the rank-K basis truncation on either engine, so ",
      "the request refuses under `auto` as well as under an explicit ",
      "`inla` / `brms`."
    ),
    message_template = paste0(
      "flexyBayes cannot fit a smooth with the low_rank_smooth ",
      "approximation: it has no route for the approximation on either ",
      "engine. Drop the representation = ... argument and write the ",
      "smooth as spl(x), which fits on INLA."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "approximation_spec_invalid",
    description = paste0(
      "s(..., representation = ): the representation spec is not a ",
      "list / fb_approx() carrying a single-string scheme."
    ),
    message_template = paste0(
      "s(%s, representation = ...): the representation spec must be ",
      "a list (or fb_approx() object) carrying a single-string ",
      "`scheme`, e.g. representation = list(scheme = ",
      "\"low_rank_smooth\", rank = 5L)."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "approximation_no_smooth_path",
    description = paste0(
      "s(..., representation = ): the named scheme is registered ",
      "but has no smooth-basis emit path at this release."
    ),
    message_template = paste0(
      "s(%s, representation = ...): scheme '%s' is registered but ",
      "has no smooth-basis emit path at this release. The only ",
      "smooth approximation scheme with an emit path is ",
      "low_rank_smooth."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "approximation_scheme_unknown",
    description = paste0(
      "validate_approximation(): the fit carries no recognised ",
      "approximation to validate (it is exact)."
    ),
    message_template = paste0(
      "validate_approximation(): this fit carries no recognised ",
      "approximation to validate (exactness = '%s'). ",
      "validate_approximation() applies only to approximate fits ",
      "(those routed through a registered approximation scheme such ",
      "as low_rank_smooth). An exact fit has no bias bound to ",
      "report."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "approximation_absent",
    description = paste0(
      "validate_approximation(): the low_rank_smooth scheme is ",
      "registered but no smooth term was routed through the ",
      "truncation path on this fit."
    ),
    message_template = paste0(
      "validate_approximation(): the fit carries no low_rank_smooth ",
      "smooth to validate. The scheme registered but no smooth term ",
      "was routed through the truncation path."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "low_rank_rank_invalid",
    description = paste0(
      "low_rank_smooth: the requested rank is not a single positive ",
      "integer."
    ),
    message_template = paste0(
      "low_rank_smooth rank%s must be a single positive integer. ",
      "Got %s."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "low_rank_rank_exceeds_basis",
    description = paste0(
      "low_rank_smooth: the requested rank meets or exceeds the ",
      "truncation ceiling min(basis dimension k, n) and so is not ",
      "an approximation."
    ),
    message_template = paste0(
      "low_rank_smooth rank%s (%d) exceeds the truncation ceiling ",
      "min(basis dimension k = %d, n = %d) = %d. A rank at or above ",
      "the basis dimension reproduces the exact basis. Request the ",
      "exact smooth (drop the approximation) instead, or choose ",
      "rank <= %d."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  # (The dispatch-side approximate-route refusal reuses the existing
  # v0.3.10 code approximate_route_not_yet_registered, one code,
  # two site-specific messages. Nothing new to register here.)

  # ---- aggregated-Gaussian emit (emit_gaussian_aggregated.R)
  .register_refusal(
    reason_code = "heterogeneous_residual_factor_not_in_cell_key",
    description = paste0(
      "Aggregated Gaussian emit: an at(f):units heterogeneous ",
      "residual factor is not in the cell key, so the ",
      "cell-constant sigma property does not hold."
    ),
    message_template = paste0(
      "emit_gaussian_aggregated(): heterogeneous residual ",
      "at(%s):units refused, '%s' is not in the cell key {%s}. ",
      "The cell-constant sigma property does not hold, so the ",
      "per-row / per-cell algebraic identity breaks. Pass ",
      "aggregate = FALSE for the per-row path."
    ),
    registered_in_adr = "ADR 0022",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "residual_type_unsupported_for_aggregation",
    description = paste0(
      "Aggregated Gaussian emit: the residual term type is outside the ",
      "supported aggregation scope."
    ),
    message_template = paste0(
      "emit_gaussian_aggregated(): residual term type '%s' is not ",
      "supported by the aggregated path. Only homogeneous (units / ",
      "id) and at_units heterogeneous residual are supported. Pass ",
      "aggregate = FALSE for the per-row path."
    ),
    registered_in_adr = "ADR 0022",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # ---- preflight refusals (fb_preflight.R, dispatch.R) -----
  .register_refusal(
    reason_code = "design_memory_exceeds_ceiling",
    description = paste0(
      "Preflight: the design matrix is estimated to exceed the ",
      "active memory ceiling. Dispatch is short-circuited before ",
      "any backend code runs."
    ),
    message_template = paste0(
      "flexyBayes preflight refused: the design exceeds the active ",
      "memory ceiling. The dispatch was short-circuited before any ",
      "backend code ran."
    ),
    registered_in_adr = "ADR 0024",
    plan_field = "preflight",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "representation_unknown_for_preflight",
    description = paste0(
      "Preflight: the design representation is not characterised by ",
      "the preflight memory estimator."
    ),
    message_template = paste0(
      "flexyBayes preflight refused: the design representation is ",
      "not characterised by the preflight estimator. The dispatch ",
      "was short-circuited before any backend code ran, raising ",
      "the memory ceiling will not help."
    ),
    registered_in_adr = "ADR 0024",
    plan_field = "preflight",
    since_version = "0.4.0"
  )

  # ---- family support (utils.R::.resolve_family) ---------------------
  # .resolve_family() is the single family gate every user entry
  # passes through (asreml via fb.R, brms via fb_from_brms.R), for a
  # family name and a family object alike. It
  # admits only the families flexyBayes can emit. Any other family is
  # refused up front, never silently fitted, including those INLA's
  # likelihood roster recognises but flexyBayes has no emit path for
  # (e.g. survival / time-to-event). The same holds for a link: each
  # family is admitted with the links an engine receives as written
  # (.fb_family_links(), R/utils.R), and any other link is refused here.
  .register_refusal(
    reason_code = "unsupported_family",
    description = paste0(
      "The requested family, or the link requested with it, is outside ",
      "the set flexyBayes can emit. Refused at the family gate ",
      "(.resolve_family) before any backend code runs."
    ),
    message_template = paste0(
      "Unsupported family \"%s\". flexyBayes supports: gaussian, ",
      "binomial, binary, poisson, negative_binomial, negbinom, ",
      "gamma, beta, hurdle_gamma, lognormal. Other families (including ",
      "survival / time-to-event and multivariate responses) are planned ",
      "future work. See fb_refusals()."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # ---- formula / spec parse-time refusals (parse_formula.R) ----------
  # Early model-spec and data-binding refusals raised by the formula
  # parser. Migrated from raw stop() to the structured registry so
  # downstream tooling pattern-matches on the condition class, not
  # on the free-text message. Missing-suggested-package errors
  # (mgcv) and argument-combination guards stay as plain
  # stop(): they are environment / call-shape errors, not
  # model-scope refusals.
  .register_refusal(
    reason_code = "formula_not_two_sided",
    description = paste0(
      "The model formula must be two-sided (response ~ predictors). ",
      "A formula carrying no left-hand-side response was supplied."
    ),
    message_template = paste0(
      "The model formula must be two-sided: response ~ predictors."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "response_not_in_data",
    description = paste0(
      "The response variable named on the formula's left-hand side ",
      "is not a column of `data`."
    ),
    message_template = paste0(
      "Response variable \"%s\" not found in data."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "fa_rank_invalid",
    description = paste0(
      "A factor-analytic term fa(x, k) was given a rank that is not a ",
      "positive whole number written in the formula. Registered at 0.4.0 ",
      "as the k < 1 floor and widened by the 0.10.0.9000 grammar contract ",
      "to every non-count rank, a symbol, a string, a logical, a ",
      "fraction, each of which as.integer() used to coerce into a rank ",
      "the user did not write."
    ),
    message_template = paste0(
      "fa() requires k >= 1. Got k = %s."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "fa_rank_exceeds_dim",
    description = paste0(
      "A factor-analytic term fa(x, k) was given a rank k that is not ",
      "strictly below the number of levels of the outer factor. A ",
      "factor-analytic covariance is identifiable only for k < n_outer: ",
      "at k = n_outer the loadings and specific variances form an ",
      "over-parameterised reparameterisation of the unstructured form, ",
      "and at k > n_outer the lower-triangular loadings carry empty ",
      "columns. This is a data-aware preflight (n_outer is known only ",
      "after the term is matched against the data), complementing the ",
      "data-free fa_rank_invalid (k < 1) check."
    ),
    message_template = paste0(
      "fa() requires k < the number of levels of the outer factor, ",
      "got k = %s with %s level(s)."
    ),
    registered_in_adr = "06-06 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.7.0"
  )
  .register_refusal(
    reason_code = "smooth_variable_not_in_data",
    description = paste0(
      "The variable inside a smooth term s(x) is not a column of ",
      "`data`."
    ),
    message_template = paste0(
      "Smooth s() variable \"%s\" not found in data."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "tensor_smooth_unsupported",
    description = paste0(
      "A tensor-product or multivariate smooth (te(), ti(), t2()) ",
      "was supplied. flexyBayes fits univariate penalised splines ",
      "(s(), spl()) only."
    ),
    message_template = paste0(
      "Tensor-product / multivariate smooth %s is not supported. ",
      "flexyBayes fits univariate penalised splines only, use ",
      "s(x) or spl(x) per smooth dimension."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # ---- entry-point argument guards (flexybayes.R) -------------------
  # Call-shape guards on the review_code / return_code flags. They are
  # user-reachable refusals (the registry's inclusion criterion), so
  # they carry catchable classes even though they are argument-
  # combination errors instead of model-scope refusals.
  .register_refusal(
    reason_code = "review_code_backend_unsupported",
    description = paste0(
      "review_code = TRUE was requested on a backend with no code to ",
      "review. The inspect-then-fit token is scoped to the ",
      "code-emitting engines, which among the active pair is brms."
    ),
    message_template = paste0(
      "`review_code = TRUE` is supported with backend = \"brms\" ",
      "(Stan code). Pass backend = \"brms\", or drop review_code."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "code_flags_mutually_exclusive",
    description = paste0(
      "return_code = TRUE and review_code = TRUE were both supplied. ",
      "The two code-return modes are mutually exclusive."
    ),
    message_template = paste0(
      "`return_code` and `review_code` are mutually exclusive."
    ),
    registered_in_adr = "05-30 parse-time refusals",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # ---- fb_cov() constructor ----------
  # Carrier-construction refusals reachable from both the standalone
  # fb_cov() constructor and the inline `cov = fb_cov(...)` formula
  # front door.
  .register_refusal(
    reason_code = "fb_cov_type_unknown",
    description = paste0(
      "fb_cov(): the requested carrier `type` is not one of the five ",
      "known types (dense / chol / precision / blocks / low_rank)."
    ),
    message_template = paste0(
      "fb_cov(): '%s' is not a known carrier type. Supported types: ",
      "dense, chol, precision, blocks, low_rank."
    ),
    registered_in_adr = "ADR 0030 (C3)",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "fb_cov_missing_matrix",
    description = paste0(
      "fb_cov(): the carrier matrix `M` (the first argument) was not ",
      "supplied."
    ),
    message_template = paste0(
      "fb_cov(type = \"%s\"): the carrier matrix `M` is required."
    ),
    registered_in_adr = "ADR 0030 (C3)",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "cov_arg_not_fb_cov",
    description = paste0(
      "vm() / ped(): the `cov` argument must be written inline as an ",
      "fb_cov() carrier."
    ),
    message_template = paste0(
      "%s(%s, cov = ...): the `cov` argument must be an inline fb_cov() ",
      "carrier, e.g. cov = fb_cov(G, type = \"chol\")."
    ),
    registered_in_adr = "ADR 0030 (C3)",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )

  # ---- backend capability gates (dispatch.R) -------
  # Raised when flexybayes(backend = "brms") meets a structured-covariance
  # term the brms lowering does not take: a vm() / ped() term with a blocks
  # or low-rank carrier, or a structured term with no grouping factor. The
  # brms capability predicate (R/backend_registry.R, .capability_brms) is
  # the declarative source. The dispatch brms branch raises this code and
  # uses the template below as the message.
  .register_refusal(
    reason_code = "stan_cannot_represent_structured_cov",
    description = paste0(
      "flexyBayes has no brms lowering for a structured-covariance term in ",
      "the model: a vm() / ped() term with a blocks or low-rank carrier, ",
      "or a structured term written without a grouping factor."
    ),
    message_template = paste0(
      "flexyBayes cannot fit this model on brms: it has no brms lowering ",
      "for one of its structured-covariance terms: a vm() or ped() term ",
      "with a blocks or low-rank carrier, or a us() term written without ",
      "a grouping factor. Give vm() a dense, Cholesky or precision carrier, ",
      "or write the structure as us(f):g."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.5.0"
  )

  # Raised when a model carrying a STRUCTURED RESIDUAL reaches brms,
  # whether requested explicitly or reached through the auto fallback
  # after an INLA runtime failure. brms reconstructs the fixed and
  # random blocks only. A residual covariance has no lowering, so
  # emitting would silently return an iid model under the requested
  # model's name. Refusing is the only faithful outcome until a
  # residual lowering exists and is validated against an independent
  # oracle.
  .register_refusal(
    reason_code = "stan_cannot_represent_ar1_residual",
    description = paste0(
      "flexyBayes has no brms lowering for an AR1 or separable AR1xAR1 ",
      "residual covariance."
    ),
    message_template = paste0(
      "flexyBayes cannot fit this model on brms: it has no brms lowering ",
      "for an AR1 or separable AR1xAR1 residual, and fitting without it ",
      "would return an independent-residual model. Fit it with ",
      "backend = \"inla\", on a complete grid with one observation per ",
      "node."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # Raised when a missing response reaches brms on a family whose
  # augmentation is not admitted. brms drops NA-response rows silently, so
  # without this the augmentation layer would quietly become complete-case
  # deletion, its `mi()` addition term, which fixes that, is admitted only
  # for the gaussian and lognormal families, the two validated by live fits.
  .register_refusal(
    reason_code = "brms_cannot_augment_nongaussian",
    description = paste0(
      "A missing response on brms for a family other than gaussian and ",
      "lognormal, the two flexyBayes has checked with brms's mi() term. ",
      "Without mi() brms drops the rows."
    ),
    message_template = paste0(
      "flexyBayes carries a missing response on brms only for the ",
      "gaussian and lognormal families, the two checked with brms's `mi()` ",
      "term, for any other family brms would drop the rows whose response ",
      "is NA, deleting the incomplete cells under the name of ",
      "augmentation. Use backend = \"inla\", which carries a missing ",
      "response as a prediction target for the families it fits, or pass ",
      "na_action = \"omit\" to delete them deliberately."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # ---- missing responses (R/na_action.R) -----------
  # The augmentation layer keeps a missing response as a latent quantity
  # so the design index set stays intact. These three mark its edges.
  .register_refusal(
    reason_code = "missing_response_refused",
    description = paste0(
      "The response has missing values and na_action = \"fail\" was ",
      "requested."
    ),
    message_template = paste0(
      "The response has missing values and na_action = \"fail\". Use ",
      "na_action = \"augment\" to carry them as latent quantities and ",
      "keep the design intact, or na_action = \"omit\" to analyse the ",
      "observed cells only."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  .register_refusal(
    reason_code = "missing_covariate_not_supported",
    description = paste0(
      "A predictor has missing values under the covariate policy `fail`, ",
      "which is both this package's default and ASReml's."
    ),
    message_template = paste0(
      "A predictor has missing values and the covariate policy is `fail`, ",
      "ASReml's default too. A missing predictor is not a missing ",
      "observation: there is no latent quantity for a value nobody ",
      "recorded. Supply the values, drop the affected rows yourself, or ",
      "pass na_action = list(y = \"include\", x = \"omit\") to drop them ",
      "with a reported count."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # The one ASReml covariate behaviour this package will not reproduce.
  # Sibling of missing_covariate_not_supported: same input, different
  # requested policy, and the two messages point at each other's remedy.
  .register_refusal(
    reason_code = "covariate_zero_fill_not_supported",
    description = paste0(
      "A predictor has missing values under the covariate policy ",
      "`include`, which in ASReml means treating the missing covariate ",
      "as zero (Reference Manual 4.2, section 3.11)."
    ),
    message_template = paste0(
      "A predictor has missing values and the covariate policy is ",
      "`include`. ASReml treats a missing covariate as zero under that ",
      "setting (ASReml-R Reference Manual 4.2, section 3.11), which puts ",
      "a value in the model that the plot did not have. Supply the ",
      "values, or pass na_action = list(y = \"include\", x = \"omit\") to ",
      "drop the affected rows and see the count."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.1"
  )

  .register_refusal(
    reason_code = "augment_cell_not_determinable",
    description = paste0(
      "A design cell is absent from the data and the model depends on a ",
      "variable whose value at that cell was never recorded."
    ),
    message_template = paste0(
      "A cell of the design grid is absent from the data altogether, and ",
      "a model variable is not constant across the trial, so completing ",
      "the grid would mean inventing a level of it for a plot nobody ",
      "recorded. ASReml's nin89 example codes such plots as LANCER. ",
      "flexyBayes will not write a level you did not. Supply the ",
      "field-book rows (one row per sown plot, response set to NA) ",
      "pad the array with fb_complete_grid(..., unused_level = ) to name ",
      "the level the empty cells carry, or use na_action = \"omit\"."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # A suggested-package guard that refuses by name instead of through a
  # bare requireNamespace() error, so a caller can catch it by class and
  # the message can name the alternative the fit already answers.
  .register_refusal(
    reason_code = "classify_requires_emmeans",
    description = paste0(
      "predict(classify = ) builds its marginal-means table through the ",
      "emmeans estimation seam, and emmeans is a suggested package, not ",
      "a hard dependency."
    ),
    message_template = paste0(
      "predict(classify = ) builds the marginal-means table through ",
      "emmeans, which is not installed. Install it with ",
      "install.packages(\"emmeans\"), or read the fixed effects from ",
      "summary(fit)$fixed and the random-effect predictions from ",
      "coef(fit, what = \"random\")."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.9.1"
  )

  # A residual variogram measures separation along a design index. A
  # model with no such array has no lag to plot against, which is a
  # property of the model instead of a gap in the display.
  .register_refusal(
    reason_code = "variogram_requires_design_index",
    description = paste0(
      "plot(type = \"variogram\") needs the row / column array a ",
      "structured covariance is indexed by, and the fit's terms name ",
      "none."
    ),
    message_template = paste0(
      "plot(type = \"variogram\") needs a design index to measure ",
      "separation along, and this fit carries none. A residual variogram ",
      "is a picture of how far apart two plots are. Fit a spatial term (",
      "random = ~ ar1(row):ar1(col)) or use plot(fit, type = ",
      "\"residuals\") for the fitted-value diagnostics."
    ),
    registered_in_adr = "0.9.1 figures",
    plan_field = NA_character_,
    since_version = "0.9.1"
  )

  # The spatial field display draws the posterior mean of a separable
  # autoregressive field. A fit whose terms name no row-by-column array has
  # no field, which is a property of the model instead of a gap in the
  # display.
  .register_refusal(
    reason_code = "spatial_field_requires_ar1_field",
    description = paste0(
      "plot(type = \"spatial_field\") draws the posterior mean of the ",
      "separable autoregressive field random = ~ ar1(row):ar1(col), and the ",
      "fit's random terms name no such field."
    ),
    message_template = paste0(
      "plot(type = \"spatial_field\") draws the posterior mean of a ",
      "separable spatial field, and this fit carries none. Fit one with ",
      "random = ~ ar1(row):ar1(col)."
    ),
    registered_in_adr = "0.14.0 figures",
    plan_field = NA_character_,
    since_version = "0.14.0"
  )

  # Raised when a model carrying observation weights reaches dispatch.
  # The ingest adapters parse `weights` into an addition term, but no
  # active emitter consumes it: brms and INLA both generate identical
  # code with and without weights, so a weighted fit silently returns
  # the unweighted posterior. "Weights" also has no single meaning:
  # inverse-residual-variance (the documented ASReml sense), frequency,
  # likelihood-power, trials and exposure are different models, and for
  # a Gaussian with unknown sigma raising each likelihood contribution
  # to w_i is not the normalised model Var(y_i) = sigma^2 / w_i.
  # Refusing until each mapping is implemented and validated against an
  # analytic oracle is safer than a no-op that looks like a fit.
  .register_refusal(
    reason_code = "weights_not_supported",
    description = paste0(
      "Observation weights are parsed but not consumed by any active ",
      "emitter, so a weighted fit is refused instead of silently ",
      "returning the unweighted posterior."
    ),
    message_template = paste0(
      "`weights` are not yet supported by the active backends. They are ",
      "parsed and recorded, but neither the brms nor the INLA emitter ",
      "consumes them, a weighted fit would return the UNWEIGHTED ",
      "posterior with no indication. Refusing instead of fitting the ",
      "wrong model. Drop `weights`, or pre-scale the response if the ",
      "intended semantics is exactly Var(y_i) = sigma^2 / w_i."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # at(f, level):units sections the residual on selected levels only, a
  # different model from dsum(~ units | f), which sections it on all of
  # them. The random-side at(f, level):g fits from 0.16.0.
  .register_refusal(
    reason_code = "at_level_conditioning_unsupported",
    description = paste0(
      "at(f, level):units gives selected levels of f a residual variance ",
      "of their own, and flexyBayes has no lowering for it on either ",
      "engine. It is a different model from dsum(~ units | f), and is ",
      "refused by name."
    ),
    message_template = paste0(
      "at(<factor>, <level>):units gives a residual variance to selected ",
      "levels only, which is a different model from a residual variance ",
      "per level, and flexyBayes has no lowering for it. Drop the level to ",
      "fit dsum(~ units | <factor>)."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # corh() is the structure between diag() and us(): k variances plus ONE
  # shared correlation. flexyBayes lowers only the two ends (uncorrelated,
  # or fully unstructured) so the middle has no emit. Approximating it in
  # either direction changes the parameter count and the meaning, so it is
  # named and refused instead. (INLA has an exchangeable group model, with
  # one variance and one shared correlation, that flexyBayes does not use.
  # Tests/testthat/test-engine-claims.R.)
  .register_refusal(
    reason_code = "corh_no_equicorrelation_representation",
    description = paste0(
      "corh(f):g, or cor(f), asks for one correlation shared by every pair ",
      "of levels. flexyBayes has no lowering for it on either engine."
    ),
    message_template = paste0(
      "corh(<factor>):<term> asks for separate variances with one shared ",
      "correlation. flexyBayes has no lowering for it on either engine. Use ",
      "diag(<factor>):<term> for separate variances without correlation, ",
      "us(<factor>):<term> on brms to estimate every correlation, or ",
      "random = ~ <term> + <factor>:<term> for one variance with a shared ",
      "correlation of zero or more."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  .register_refusal(
    reason_code = "heterogeneous_residual_multiple_factors",
    description = paste0(
      "More than one heterogeneous-residual term was supplied. A residual ",
      "is sectioned by one factor. Nesting several would require their ",
      "interaction, which the user must state explicitly."
    ),
    message_template = paste0(
      "More than one heterogeneous-residual term was supplied. A residual ",
      "is sectioned by ONE factor, if the intent is a variance per ",
      "combination, form the interaction and section by that."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # A residual variance can only vary where a residual variance exists.
  # Poisson and binomial have no free scale parameter, their dispersion is
  # a function of the mean, so `sigma ~ f` would put a predictor on a
  # parameter the model does not have.
  .register_refusal(
    reason_code = "heterogeneous_residual_family_has_no_sigma",
    description = paste0(
      "A heterogeneous residual was requested for a family with no residual ",
      "scale parameter (Poisson, binomial and relatives), whose dispersion ",
      "is determined by the mean."
    ),
    message_template = paste0(
      "A heterogeneous residual variance was requested for a family that ",
      "has no residual scale parameter to vary, its dispersion is a ",
      "function of the mean. Model the dispersion directly for that family, ",
      "or fit a Gaussian on a transformed response."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  .register_refusal(
    reason_code = "stan_cannot_represent_structured_residual",
    description = paste0(
      "flexyBayes has no brms lowering for a residual structure other ",
      "than the iid `units` residual and a residual per section."
    ),
    message_template = paste0(
      "flexyBayes cannot fit this model on brms: its residual term is ",
      "neither the iid `units` residual nor a residual per section, and ",
      "flexyBayes has no brms lowering for it, so fitting would return an ",
      "independent-residual model. Write the residual as ~ units or ",
      "~ dsum(~ units | f)."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )

  # ---- grammar polymorphism on the universal entry ----
  # fb() / flexybayes() detect the input grammar from the call shape.
  # These two guard the brms-grammar branch of
  # .build_ir_polymorphic() (R/fb.R).
  .register_refusal(
    reason_code = "grammar_brms_with_asreml_terms",
    description = paste0(
      "A brms-style bar-grouped formula was combined with ASReml ",
      "`random` / `residual` arguments on the universal entry."
    ),
    message_template = paste0(
      "A brms-style formula (with `(... | g)` grouping) cannot be ",
      "combined with `random` / `residual` (ASReml grammar). Put every ",
      "grouping term inside the formula, or use the ASReml form ",
      "throughout."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.4.1"
  )
  .register_refusal(
    reason_code = "grammar_brms_known_matrices_unsupported",
    description = paste0(
      "`known_matrices` was supplied with brms-grammar ingest via the ",
      "universal entry, which has no known-matrix carrier."
    ),
    message_template = paste0(
      "`known_matrices` is not supported with brms-grammar ingest via ",
      "fb() / flexybayes(). Use the ASReml form for known-matrix ",
      "carriers."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.4.1"
  )
  .register_refusal(
    reason_code = "engine_pin_backend_conflict",
    description = paste0(
      "An engine pin (fb_inla / fb_brms) was given a `backend` ",
      "argument that conflicts with the engine it pins."
    ),
    message_template = paste0(
      "An engine pin fits via one engine only, so it cannot take a ",
      "conflicting `backend`. Drop `backend`, or use flexybayes() ",
      "to choose a backend."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.5.0"
  )

  invisible(NULL)
}


# 0.9.0 design-fidelity: a dsum() residual with a structured inner term
# would otherwise be silently reduced to a per-region heteroscedastic
# variance, dropping the spatial (ar1/us/...) structure the user wrote.
.populate_refusal_registry_v0900 <- function() {
  .register_refusal(
    reason_code = "dsum_structured_inner_unsupported",
    description = paste0(
      "Residual dsum() with a structured inner term (ar1, us, ...). ",
      "flexyBayes represents dsum() only as a per-region heteroscedastic ",
      "variance, so a structured inner is refused instead of silently ",
      "reduced to a heteroscedastic residual."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "auto_no_active_route",
    description = paste0(
      "backend = \"auto\" found neither engine able to fit the model: the ",
      "INLA route refused it and flexyBayes has no brms lowering for it. ",
      "There is no silent fallback, so auto refuses."
    ),
    message_template = paste0(
      "backend = \"auto\": flexyBayes cannot fit this model on INLA or on ",
      "brms. Pass backend = \"inla\" or backend = \"brms\" to see which ",
      "term each route refuses, and change that term."
    ),
    registered_in_adr = "Design specification",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )
  # WP16 (2026-07-25): INLA separable-spatial latent field.
  .register_refusal(
    reason_code = "ar1_spatial_requires_complete_grid",
    description = paste0(
      "A separable ar1(row):ar1(col) spatial term on INLA is emitted as an ",
      "AR1xAR1 latent field, faithful only with one observation per grid ",
      "node. The data are an incomplete or replicated grid."
    ),
    message_template = "%s",
    registered_in_adr = "0.9.0 refusal registry",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  invisible(NULL)
}


# 0.9.0 fidelity-of-names pass. Two families of refusal register here.
#
# The first separates the random-side AR1 field, a latent field plus the
# Gaussian observation nugget, from a residual-side AR1 the package does not
# fit.
#
# The second closes the parser: an unrecognised call, and every IR type that
# parses for the formula catalogue without an emit path, now reaches a named
# refusal instead of an untyped stop or, worse, a silent reinterpretation
# as a plain variable.
.populate_refusal_registry_v0900_names <- function() {
  .register_refusal(
    reason_code = "ar1_residual_not_representable",
    description = paste0(
      "An AR1 process written on the residual side that flexyBayes does not ",
      "fit: the 1-D and grouped forms, and the separable form for a ",
      "response other than Gaussian with the identity link. The separable ",
      "residual on a Gaussian response fits on INLA without a nugget."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "stan_cannot_represent_ar1_field",
    description = paste0(
      "A random-side AR1 or separable AR1xAR1 latent field on brms. ",
      "flexyBayes builds the field on INLA only and has no brms lowering ",
      "for it."
    ),
    message_template = paste0(
      "flexyBayes fits a random-side AR1 field (ar1(t), or the separable ",
      "ar1(row):ar1(col)) on INLA only, as a latent autoregressive field ",
      "plus the Gaussian observation nugget. It has no brms lowering for ",
      "the field, and fitting without it would return an ",
      "independent-effect model under the same name. Fit it with ",
      "backend = \"inla\" on a complete grid with one observation per node."
    ),
    registered_in_adr = "Design specification",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "asreml_function_not_recognised",
    description = paste0(
      "A call in a random or residual formula whose function is outside the ",
      "recognised ASReml vocabulary. Previously such a call was read as a ",
      "plain variable named after its own source text."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "ar2_not_representable",
    description = paste0(
      "ar2() parses for the formula catalogue and has no emit path on ",
      "either active engine."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "str_not_representable",
    description = paste0(
      "str() declares a covariance structure over a group of terms. It ",
      "parses for the formula catalogue and has no emit path on either ",
      "active engine."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "fa_not_representable",
    description = paste0(
      "A factor-analytic covariance written without a grouping factor, ",
      "fa(f, k). It parses for the formula catalogue, and flexyBayes does ",
      "not fit a factor-analytic covariance in any form."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "interaction_not_representable",
    description = paste0(
      "An interaction of two structured terms that matches none of the ",
      "recognised interaction patterns, us(trait):vm(gen), say. It parsed ",
      "into the catalogue's generic-interaction node, which flexyBayes has ",
      "no lowering for on either engine."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "brms_cannot_represent_term",
    description = paste0(
      "A random-effect term type reached the brms formula reconstruction ",
      "with no branch to lower it. It is refused by name instead of ",
      "dropped from the emitted model."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "numeric_variable_in_random_interaction",
    description = paste0(
      "A numeric variable inside a nested or multi-way random interaction. ",
      "The spelling reads as a random regression to an ASReml user and as ",
      "one independent deviation per cell to the emitted model, flexyBayes ",
      "refuses instead of choosing one of the two readings."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "row_count_exceeds_integer",
    description = paste0(
      "A single-file or in-memory source reports more rows than R's ",
      "integer type holds. Recording the count would coerce it to NA, so ",
      "the aggregation plan and the fitted object would carry no row ",
      "count. Partitioned .fst input keeps its total as a double and has ",
      "no such ceiling."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "family_argument_not_recognised",
    description = paste0(
      "The `family` argument was neither a single family name nor a ",
      "usable `stats` family object, or a family object and an explicit ",
      "`link =` named different links. A family object supplies its own ",
      "link, so resolving the disagreement either way would fit a model ",
      "the call did not ask for."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "met_summary_not_available",
    description = paste0(
      "fb_met_summary() was called on a fit it cannot summarise: the ",
      "breeder summary reads realised factor-analytic effects, and ",
      "flexyBayes fits no fa() term on either engine, so this refusal is ",
      "unconditional."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "inla_gate_refused",
    description = paste0(
      "An explicit backend = \"inla\" request that lgm_gate() refused. ",
      "The gate's structural check-list is the message and rides on the ",
      "condition in its `lgm_refusal` slot. The failing rule supplies a ",
      "second condition class so one rule can be caught on its own."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = "rejected_routes",
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "fit_lacks_posterior_draws",
    description = paste0(
      "A method needing the posterior itself (draws on a sampling ",
      "engine, marginal densities on INLA) was called on a fit that ",
      "carries neither. The refusal replaces an empty result that would ",
      "read as an answer. Its sibling is loo_requires_sampler_draws, ",
      "raised when the fit does carry a posterior but not the pointwise ",
      "log-likelihood draws that PSIS-LOO reads."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  # Sibling of fit_lacks_posterior_draws, and the two cross-reference each
  # other so the taxonomy does not read as one reason registered twice.
  # The distinction is what the fit is missing: that code is raised when
  # there is no posterior to read at all, this one when the posterior is
  # there but the quantity PSIS-LOO integrates, the log-likelihood of
  # each observation at each draw, was never stored.
  .register_refusal(
    reason_code = "loo_requires_sampler_draws",
    description = paste0(
      "loo() was called on a fit whose engine does not store the ",
      "pointwise log-likelihood draws PSIS-LOO reads. Sibling of ",
      "fit_lacks_posterior_draws, which is raised when the fit carries no ",
      "posterior at all."
    ),
    message_template = paste0(
      "loo() estimates the expected log pointwise predictive density by ",
      "importance-weighting the log-likelihood of each observation at ",
      "each posterior draw, and this fit does not store that quantity: ",
      "the engine returned an approximation to the posterior instead of ",
      "draws of the likelihood, so there is nothing to leave one ",
      "observation out of. An INLA fit does carry the information ",
      "criteria it computed at fit time (WAIC at fit$inla$waic$waic and ",
      "DIC at fit$inla$dic$dic) and a fit on the brms engine answers ",
      "loo() directly."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.9.1"
  )
  # A posterior predictive check overlays datasets DRAWN FROM the fitted
  # model on the observed response. A fit whose engine returns no
  # predictive draws has nothing to overlay, and an observed-versus-fitted
  # panel shown under the same name would be a different display wearing
  # the check's title.
  .register_refusal(
    reason_code = "pp_check_requires_predictive_draws",
    description = paste0(
      "pp_check() was called on a fit whose engine returns no posterior ",
      "predictive draws, so there are no replicated datasets to overlay ",
      "on the observed response."
    ),
    message_template = paste0(
      "pp_check() draws replicated datasets from the posterior predictive ",
      "distribution and overlays them on the observed response, and this ",
      "fit carries no predictive draws to replicate from. The residual ",
      "diagnostics it does answer are plot(fit, type = \"residuals\") and, ",
      "on a fit carrying a design index, plot(fit, type = \"variogram\"). ",
      "A fit on the brms engine answers pp_check() directly."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.9.1"
  )
  # A factor written on the fixed side AND on the random side is aliased
  # with itself: the fixed part already spends one parameter per level, so
  # the random copy's deviations are identified only by their own prior.
  # ASReml accepts the spelling, so a translated script reaches it, and
  # nothing downstream says the two term sets are estimating the same
  # thing. The failure it produces is not a clean one either, the INLA
  # marginal solve goes intermittently singular on it.
  .register_refusal(
    reason_code = "term_in_fixed_and_random",
    description = paste0(
      "The same factor appears as a fixed main effect and as a random ",
      "main effect, which aliases the random deviations against the fixed ",
      "level means: they are identified only by their own prior."
    ),
    message_template = paste0(
      "The term `%1$s` is on both the fixed and the random side. The ",
      "fixed part already estimates one mean per level of %1$s, so the ",
      "random copy has no information left to describe, its deviations ",
      "are held up by their prior alone, and the variance component ",
      "reported for %1$s would be a reading of that prior, not of ",
      "the data. Two ways to write what was probably meant: drop `%1$s` ",
      "from `random` to keep population-level means for it, or drop it ",
      "from the fixed part to keep shrunk level effects (a variance ",
      "component and BLUPs). ASReml accepts this spelling and fits it, so ",
      "a script translated from ASReml can arrive here unchanged."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.9.1"
  )
  # predict(classify = ) on a random factor whose random terms the table
  # cannot use: a structured covariance, a non-identity link, an offset, or
  # every term with the factor naming a variable the rules leave out.
  .register_refusal(
    reason_code = "classify_random_factor_not_supported",
    description = paste0(
      "predict(classify = ) was asked for a factor that enters the model as ",
      "a random term the predicted-value table cannot use."
    ),
    message_template = paste0(
      "predict(classify = \"%1$s\") is not available: %2$s. The effects ",
      "of %1$s with their intervals are in ranef(fit) and ",
      "coef(fit, what = \"random\")."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.9.1"
  )
  .register_refusal(
    reason_code = "conditional_loglik_not_available",
    description = paste0(
      "No plug-in conditional log-likelihood exists for this fit, the ",
      "family has no implemented form, the response or family record is ",
      "absent, or the engine reports a marginal likelihood instead. ",
      "Returning NA would let AIC() and anova() present a comparison that ",
      "was never computed."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "model_matrix_not_recoverable",
    description = paste0(
      "model.matrix() was called on a fit carrying no fixed-effect ",
      "formula or no fitted data, so the population-level design cannot ",
      "be rebuilt. Base R would silently return an intercept-only matrix ",
      "built from the calling frame."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  .register_refusal(
    reason_code = "update_call_not_reconstructable",
    description = paste0(
      "update() was called on a fit whose recorded argument set is ",
      "incomplete. Re-fitting would substitute package defaults for the ",
      "unrecorded arguments, a relationship matrix or a prior scale ",
      "among them, and return a different model under the same name."
    ),
    message_template = "%s",
    registered_in_adr = "Design specification",
    plan_field = NA_character_,
    since_version = "0.9.0"
  )
  invisible(NULL)
}


# --- field-hardening reasons -------------------------------------- #
#
# The 0.9.1 field-coverage sweep found three refusal-shaped holes that
# were not refusals: the prior mini-language accepted malformed
# specifications and substituted defaults, two scalar prior arguments
# were accepted and then discarded, and the INLA family gate compared a
# flexyBayes spelling against INLA's own roster. The codes below give the
# constructor, the arrival check, and the reconciled gate a vocabulary,
# so a caller can distinguish "this specification is not representable"
# from "R fell over".
#
# The prior-DSL codes cover the whole constructor surface instead of
# only the newly-caught classes: a mini-language is a public API, and a
# refusal a wrapper cannot pattern-match on is not a contract.
.populate_refusal_registry_field_hardening <- function() {
  .register_refusal(
    reason_code = "prior_spec_empty",
    description = paste0(
      "fb_prior() was called with no specification. A prior object with ",
      "no targets would silently leave every parameter on its engine ",
      "default."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_spec_not_formula",
    description = paste0(
      "An fb_prior() argument is not a formula. The DSL's only argument ",
      "shape is `target ~ distribution(...)`."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_spec_not_two_sided",
    description = paste0(
      "An fb_prior() formula is one-sided, so it names a distribution ",
      "without naming the parameter it applies to."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_target_unsupported",
    description = paste0(
      "The left-hand side of a prior formula is not one of the supported ",
      "targets. A bare name that is not `sigma` used to be carried as an ",
      "opaque target and dropped at emit time."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_target_not_a_single_name",
    description = paste0(
      "A b() prior target was given something other than one coefficient ",
      "name. as.character() on the unevaluated argument deparses a call ",
      "head-first, so b(c(\"x\", \"z\")) became the three names c(\"c\", ",
      "\"x\", \"z\") and survived fb_plan() before failing inside the ",
      "emitter with an untyped base error."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 prior-target contract",
    plan_field = NA_character_,
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "prior_target_argument_missing",
    description = paste0(
      "A prior target call is missing the argument that identifies which ",
      "parameter it targets, sd() and cor() need `group`, b() and ",
      "smooth() need a name."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_distribution_not_a_call",
    description = paste0(
      "The right-hand side of a prior formula is not a distribution call, ",
      "so it carries no hyperparameters."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_distribution_unknown",
    description = paste0(
      "The distribution named on the right-hand side of a prior formula ",
      "is outside the DSL's supported set."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_argument_unknown",
    description = paste0(
      "A prior distribution call carries an argument name the family does ",
      "not have. Before 0.9.2 the call fell through unmatched and every ",
      "hyperparameter then took its emit-time default, so a misspelt ",
      "argument became a silently substituted prior."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_argument_duplicated",
    description = paste0(
      "A prior distribution call binds the same hyperparameter twice, so ",
      "which of the two values applies is not decidable from the call."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_argument_missing",
    description = paste0(
      "A prior distribution call omits a hyperparameter the family ",
      "requires. A half-specified PC prior is the sharp case: half a ",
      "probability statement is not a probability statement, and the ",
      "completed value reached the engine unannounced."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_hyperparameter_not_scalar",
    description = paste0(
      "A prior hyperparameter is not a finite numeric scalar. A zero-length ",
      "or length-two value used to surface as one of R's own interpreter ",
      "messages several layers downstream."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_hyperparameter_out_of_domain",
    description = paste0(
      "A prior hyperparameter lies outside the domain its family defines: ",
      "a non-positive scale, rate, shape or degrees of freedom, a tail ",
      "probability outside (0, 1), or reversed uniform bounds. A negative ",
      "standard deviation used to reach the sampler and surface as a ",
      "missing-draws error."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_not_translatable_for_backend",
    description = paste0(
      "flexyBayes does not translate a supplied prior specification for the ",
      "selected backend. Before 0.9.2 the INLA route discarded such a row ",
      "silently while prior_summary() printed it as applied, and the brms ",
      "route refused it with an untyped error."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "prior_target_not_in_model",
    description = paste0(
      "A prior names a coefficient or a variance component the model does ",
      "not have. The check reads the engine's own parameter list, so the ",
      "refusal names the model's actual terms instead of a synthesised ",
      "Stan parameter name."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "fixed_smoother_not_supported",
    description = paste0(
      "An mgcv-style smoother (s(), te(), ti(), t2(), gp(), sos()) appears ",
      "in the fixed part of the formula. flexyBayes spells a univariate ",
      "smooth `random = ~ spl(x)` on the INLA backend. Before 0.9.2 the ",
      "term reached the engine and surfaced as `could not find function ",
      "\"s\"`, on brms only after a completed sampling run."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "inla_variable_used_twice",
    description = paste0(
      "The same variable indexes both a fixed-effect term and a latent ",
      "f() term in the emitted INLA formula. INLA refuses a duplicated ",
      "key, before 0.9.2 the refusal surfaced as the generic ",
      "inla-program-exited message."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "brms_factor_random_slope_unsupported",
    description = paste0(
      "A bar-grammar random slope on a factor ((f | g), (0 + f | g), ",
      "(f || g)) is outside the brms ingest. The ASReml surface of the ",
      "same model fits, so the refusal names it verbatim. Before 0.9.2 the ",
      "message listed other bar spellings, none of which expresses the ",
      "model."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "brms_random_effect_form_unsupported",
    description = paste0(
      "A bar-grammar random-effect specification outside the supported ",
      "set, multi-variable slopes and group-of-coefficients shorthand. ",
      "Untyped before 0.9.2."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "brms_ingest_feature_unsupported",
    description = paste0(
      "A brms formula feature outside the ingest layer's scope: an ",
      "extended term the walker collected but cannot lower, or a part of a ",
      "brms::bf() object the ingest does not carry to an engine (a ",
      "non-linear declaration, a formula or fixed value for a distributional ",
      "parameter, or a family that differs from the `family` argument). ",
      "The extended terms were untyped before 0.9.2, and the brms::bf() ",
      "parts were dropped with no message before 0.14.0."
    ),
    message_template = "%s",
    registered_in_adr = "Field sweep",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )

  # ---- aggregated-route entry gate (dispatch.R) ---------------------
  #
  # The execution grid's capability-matrix section reached three
  # aggregation refusals that the v0.4.0 bulk migration did not cover,
  # because they sit above the plan gate instead of inside it: the two
  # refusals below the plan gate already carry a condition class, these
  # three raised a bare stop(). They split on two boundaries, not one.
  # The response boundary is a property of the data (no sufficient
  # statistic exists), the route boundary is a property of the engine
  # roster (no aggregated emit is reachable), a caller handling one
  # has no business catching the other, so they get separate codes.
  .register_refusal(
    reason_code = "aggregation_response_incomplete",
    description = paste0(
      "aggregate = TRUE with missing values in the response. The ",
      "aggregated path compresses rows into per-cell sufficient ",
      "statistics, which are not defined for a missing response, so ",
      "the request has no aggregated reading. Untyped before 0.9.2."
    ),
    message_template = "%s",
    registered_in_adr = "2026-08-19 execution grid",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  .register_refusal(
    reason_code = "aggregation_route_unavailable",
    description = paste0(
      "aggregate = TRUE with no aggregated emit route available, ",
      "either the named backend carries no aggregated path at all, or ",
      "no aggregated backend resolves for the model, including a family ",
      "or link flexybayes_stream() does not aggregate. Raised with ",
      "site-specific messages. Untyped before 0.9.2."
    ),
    message_template = "%s",
    registered_in_adr = "2026-08-19 execution grid",
    plan_field = NA_character_,
    since_version = "0.9.2"
  )
  invisible(NULL)
}


# --- v0.9.3 registry population ---------------------------------- #

# .populate_refusal_registry_v093(): the 0.9.3 scale-strategy
# refusals: the cell-count integer-cast sibling of the 0.9.0
# row-count guard, typed engine-death wraps for INLA and brms, and the
# per-level AR1-field boundary this release still refuses, plus
# the unknown-backend refusal that replaces the two engine names this
# release withdraws entirely (see NEWS.md). Called from .onLoad() before
# .lock_refusal_registry().
.populate_refusal_registry_v093 <- function() {
  .register_refusal(
    reason_code = "unknown_backend",
    description = paste0(
      "An explicit `backend` request named a value that is not a ",
      "registered engine, including a name this release withdrew ",
      "entirely (see NEWS.md). The active engines are inla and brms."
    ),
    message_template = paste0(
      "backend = \"%s\" is not a recognised flexyBayes engine. The ",
      "active engines are \"inla\" and \"brms\" (or \"auto\" to let ",
      "flexyBayes choose)."
    ),
    registered_in_adr = "0.9.3 backend withdrawal",
    plan_field = "rejected_routes",
    since_version = "0.9.3"
  )
  .register_refusal(
    reason_code = "cell_count_exceeds_integer",
    description = paste0(
      "The aggregation planner's estimated cell count (the product of ",
      "the cell key's factor level counts) exceeds R's integer limit. ",
      "Recording it would coerce it to NA, so the plan's cell count and ",
      "compression estimate would both silently disappear. The row-count ",
      "sibling (row_count_exceeds_integer) fixed the identical defect for ",
      "N, this is the cell-count member of the same class."
    ),
    message_template = "%s",
    registered_in_adr = "0.9.3 scale strategy",
    plan_field = NA_character_,
    since_version = "0.9.3"
  )
  .register_refusal(
    reason_code = "inla_program_failed",
    description = paste0(
      "The INLA engine died, a non-zero exit from the inla-program ",
      "subprocess, an inla.core.safe failure, or an empty/NULL result ",
      "with no R-level error at all, past the design-memory preflight, ",
      "which models per-term storage but not the sparse-Cholesky solver ",
      "cost of a large latent field. Before 0.9.3 this reached the user ",
      "as a raw pass-through of the engine's own message after tens of ",
      "minutes of runtime, with no reason code and no remedy. ",
      "Carries the design size, the largest design this package has ",
      "verified an INLA fit to complete (read from ",
      "inst/validation/benchmark_scaling.csv), the binding random-effect ",
      "term where determinable, and the remedies."
    ),
    message_template = "%s",
    registered_in_adr = "0.9.3 scale strategy",
    plan_field = NA_character_,
    since_version = "0.9.3"
  )
  .register_refusal(
    reason_code = "weights_requires_gaussian",
    description = paste0(
      "`weights` is lowered for family = \"gaussian\" with the identity ",
      "link only (brms's y | weights(w) addition term; INLA's ",
      "scale = w per-observation precision multiplier, both give ",
      "Var(y_i) = sigma^2 / w_i on the Gaussian likelihood). Superseded ",
      "the all-families weights_not_supported refusal for this one ",
      "well-defined mapping. Every other family, and a non-identity ",
      "link on gaussian, still refuse, weights_not_supported remains ",
      "registered for the historical all-families refusal text."
    ),
    message_template = "%s",
    registered_in_adr = "0.9.3 weights refusals",
    plan_field = NA_character_,
    since_version = "0.9.3"
  )
  .register_refusal(
    reason_code = "weights_not_aggregatable",
    description = paste0(
      "aggregate = TRUE with non-constant observation weights. Weights ",
      "are lowered for the per-row Gaussian route only. The aggregated ",
      "route's sufficient statistics (per-cell sums and counts) do not ",
      "fold in a per-observation weight, so a weighted aggregated fit is ",
      "refused instead of silently ignoring the weight or answering a ",
      "different question under the same name."
    ),
    message_template = "%s",
    registered_in_adr = "0.9.3 weights refusals",
    plan_field = NA_character_,
    since_version = "0.9.3"
  )
  .register_refusal(
    reason_code = "at_level_not_in_data",
    description = paste0(
      "at(f, level) names a level that is not a level of f in the data. ",
      "Levels are literal values or positions in levels(f), as ASReml ",
      "reads them."
    ),
    message_template = "%s",
    registered_in_adr = "0.16.0 at() and AR1 residual checks",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
  .register_refusal(
    reason_code = "ar1_residual_with_random_field",
    description = paste0(
      "The nugget-free autoregressive residual together with an ",
      "autoregressive field on the random side. flexyBayes fits one of ",
      "the two per model."
    ),
    message_template = "%s",
    registered_in_adr = "0.16.0 at() and AR1 residual checks",
    plan_field = NA_character_,
    since_version = "0.16.0"
  )
  .register_refusal(
    reason_code = "binomial_response_not_binary",
    description = paste0(
      "family = \"binomial\" on the main (non-streaming) entry requires ",
      "the response to be a numeric vector of only 0 and 1: there is no ",
      "cbind(success, failure) or trials = spelling on flexybayes() / ",
      "fb_inla() / fb_brms() (that exists only on flexybayes_stream()). ",
      "Before this refusal, a response outside {0, 1} reached the engine ",
      "unrefused and failed there: INLA with a raw subprocess error ",
      "('Binomial data ... is void', wrapped into inla_program_failed ",
      "with an unrelated design-memory diagnosis), brms with an untyped ",
      "condition ('Family ... requires responses to contain only two ",
      "different values'), and the aggregated route with an untyped ",
      "stop(). This refusal fires once, before any of the three, naming ",
      "the column, the offending values, and the remedy: trials = on ",
      "flexybayes_stream() for pre-aggregated counts, or recode to a ",
      "0/1 numeric vector for a per-row binary outcome."
    ),
    message_template = "%s",
    registered_in_adr = "0.9.3 argument checks",
    plan_field = NA_character_,
    since_version = "0.9.3"
  )
  .register_refusal(
    reason_code = "update_unnamed_argument_not_supported",
    description = paste0(
      "update() re-issues every recorded flexybayes() argument by name, ",
      "matching an override in ... to it by name too. An unnamed ",
      "argument, stats::update()'s classic replacement-formula idiom, ",
      "update(fit, . ~ . + z), has no name to match and was, before ",
      "this refusal, silently discarded: the re-fit ran with the ",
      "unchanged recorded call, a valid object answering a question ",
      "nobody asked instead of a re-fit or a refusal. Name the slot ",
      "instead: update(fit, fixed = y ~ x + z) (brms-style), or ",
      "random = ~ ... / residual = ~ ... on the ASReml grammar."
    ),
    message_template = "%s",
    registered_in_adr = "0.9.3 argument checks",
    plan_field = NA_character_,
    since_version = "0.9.3"
  )
  invisible(NULL)
}


# --- lock helper -------------------------------------------------- #

# .lock_refusal_registry(): locks the environment so no further
# .register_refusal() calls or user-side assign()s can mutate it.
# Called once at the end of .onLoad(). The lock is enforced by R's
# environmentIsLocked() machinery. Bindings registered before the
# lock remain readable indefinitely.
.lock_refusal_registry <- function() {
  if (!environmentIsLocked(.refusal_registry)) {
    lockEnvironment(.refusal_registry, bindings = TRUE)
  }
  invisible(NULL)
}


# --- user-facing discovery surface ------------------------------- #

#' List flexyBayes refusal reasons
#'
#' `fb_refusals()` exposes the locked refusal-reason registry as a
#' browsable table: the canonical vocabulary of conditions under which
#' flexyBayes declines to fit, route, or validate a model, each
#' with a one-line description and the release it was introduced in.
#' It is the discovery surface for the structured refusals the package
#' raises. Every such refusal carries
#' a condition class `flexybayes_refusal_<reason_code>`, so a reason
#' listed here can be caught precisely, for example with
#' `tryCatch(fit, flexybayes_refusal_precision_not_symmetric =
#' handler)`.
#'
#' Two optional filters narrow the listing. `reason_code` selects rows
#' by exact reason-code match (a single code or a vector). The
#' `since_version` filter selects reasons introduced in a matching
#' release by version-string prefix: `since_version = "0.4"`
#' returns every reason added in the 0.4 series.
#'
#' Routing-decision reasons (surfaced by [fb_plan()] and
#' [backend_decision()]) and internal control-flow signals are
#' deliberately excluded: this table lists only refusals a user can
#' actually encounter.
#'
#' Packages that extend flexyBayes add their own reasons with
#' [fb_register_refusals()]. The `package` column names the package
#' each reason belongs to: `"flexyBayes"` for the core, and the
#' plug-in's name for the rest.
#'
#' @param reason_code Optional character vector of exact reason codes
#'   to filter to. `NULL` (default) returns all registered reasons.
#' @param since_version Optional single version-string prefix to
#'   filter to. `NULL` (default) returns all.
#'
#' @return A data frame of subclass `fb_refusals_table`, one row per
#'   matching refusal reason, with columns `reason_code`, `package`,
#'   `description`, `since_version`, and `plan_field`. The print
#'   method renders it as a compact checklist.
#'
#' @examples
#' fb_refusals()
#' fb_refusals(reason_code = "precision_not_symmetric")
#' fb_refusals(since_version = "0.4")
#'
#' @export
fb_refusals <- function(reason_code = NULL, since_version = NULL) {
  if (!is.null(reason_code) && !is.character(reason_code)) {
    stop("`reason_code` must be NULL or a character vector.", call. = FALSE)
  }
  if (
    !is.null(since_version) &&
      (!is.character(since_version) || length(since_version) != 1L)
  ) {
    stop(
      "`since_version` must be NULL or a single version string.",
      call. = FALSE
    )
  }

  entries <- .refusal_entries()

  df <- data.frame(
    reason_code = vapply(entries, `[[`, character(1L), "reason_code"),
    package = vapply(entries, `[[`, character(1L), "package"),
    description = vapply(entries, `[[`, character(1L), "description"),
    since_version = vapply(entries, `[[`, character(1L), "since_version"),
    plan_field = vapply(entries, `[[`, character(1L), "plan_field"),
    stringsAsFactors = FALSE
  )

  if (!is.null(reason_code)) {
    df <- df[df$reason_code %in% reason_code, , drop = FALSE]
  }
  if (!is.null(since_version)) {
    df <- df[startsWith(df$since_version, since_version), , drop = FALSE]
  }

  df <- df[order(df$reason_code), , drop = FALSE]
  rownames(df) <- NULL
  structure(
    df,
    class = c("fb_refusals_table", "data.frame"),
    filter = list(reason_code = reason_code, since_version = since_version)
  )
}

#' @export
print.fb_refusals_table <- function(x, ...) {
  n <- nrow(x)
  flt <- attr(x, "filter")

  filt_txt <- character(0L)
  if (!is.null(flt$reason_code)) {
    filt_txt <- c(
      filt_txt,
      paste0("reason_code in {", paste(flt$reason_code, collapse = ", "), "}")
    )
  }
  if (!is.null(flt$since_version)) {
    filt_txt <- c(filt_txt, paste0("since_version ~ '", flt$since_version, "'"))
  }

  header <- paste0(
    "flexyBayes refusal registry: ",
    n,
    " reason",
    if (n == 1L) "" else "s"
  )
  if (length(filt_txt)) {
    header <- paste0(
      header,
      "  (filter: ",
      paste(filt_txt, collapse = "; "),
      ")"
    )
  }
  cat(header, "\n", sep = "")

  if (n == 0L) {
    cat("  (no matching refusal reasons)\n")
    return(invisible(x))
  }

  cat("\n")
  for (i in seq_len(n)) {
    owner <- if (identical(x$package[[i]], "flexyBayes")) {
      ""
    } else {
      paste0("  (", x$package[[i]], ")")
    }
    cat(sprintf(
      "  [since %s] %s%s\n",
      x$since_version[[i]],
      x$reason_code[[i]],
      owner
    ))
    cat(strwrap(x$description[[i]], width = 72L, prefix = "      "), sep = "\n")
    cat("\n")
  }
  invisible(x)
}


# --- plug-in partitions ------------------------------------------ #

#' Register a package's refusal reasons
#'
#' A package that extends flexyBayes calls `fb_register_refusals()` from its
#' own `.onLoad()` to add the reasons it refuses with. Its users then catch
#' those refusals with the same `flexybayes_refusal` condition class as the
#' core's, find them in [fb_refusals()], and the package raises them with
#' [fb_refuse()].
#'
#' Each package's reasons go into a partition named after the package.
#' A reason code is unique across all partitions: a code that flexyBayes or
#' another package already registered is refused, and so is any attempt to
#' register under the name `"flexyBayes"`. Registering the same code again
#' from the same package with identical fields does nothing, so a package
#' can be loaded twice in one session. Registering it with different fields
#' is refused. A table that fails any check registers none of its rows.
#'
#' @param refusals A data frame with one row per reason and these
#'   character columns:
#'   * `reason_code` (required): lower-case letters, digits and
#'     underscores, starting with a letter. Prefix codes with something
#'     particular to the package so they do not collide with others.
#'   * `description` (required): one sentence on what the refusal means,
#'     shown by [fb_refusals()] and used by [fb_refuse()] when no message
#'     is given.
#'   * `since_version` (required): the package version that introduced
#'     the reason.
#'   * `message_template` (optional, default `"%s"`): the message shape,
#'     for documentation.
#'   * `registered_in_adr` (optional, default `NA`): the design record the
#'     reason traces to.
#'   * `plan_field` (optional, default `NA`): the [fb_plan()] slot that
#'     reports the refusal before a fit, if any.
#' @param package The package that owns the reasons. `NULL` (default)
#'   takes the package whose code makes the call, which is the right value
#'   inside `.onLoad()`. Outside a package it must be given.
#'
#' @return The registered reasons as an `fb_refusals_table`, invisibly.
#'
#' @seealso [fb_refuse()] to raise a registered reason, [fb_refusals()] to
#'   list them.
#'
#' @examples
#' \dontrun{
#' # In the plug-in's R/zzz.R:
#' .onLoad <- function(libname, pkgname) {
#'   flexyBayes::fb_register_refusals(
#'     data.frame(
#'       reason_code = "mypkg_threshold_negative",
#'       description = "The threshold must be zero or positive.",
#'       since_version = "0.1.0"
#'     )
#'   )
#' }
#' }
#'
#' @export
fb_register_refusals <- function(refusals, package = NULL) {
  if (is.null(package)) {
    package <- .fb_calling_package(parent.frame())
  }
  if (is.null(package)) {
    stop(
      "fb_register_refusals() was called from outside a package, so ",
      "`package` must name the package that owns these refusal codes.",
      call. = FALSE
    )
  }
  if (
    !is.character(package) ||
      length(package) != 1L ||
      is.na(package) ||
      !grepl("^[A-Za-z][A-Za-z0-9.]*[A-Za-z0-9]$", package)
  ) {
    stop("`package` must be a single package name.", call. = FALSE)
  }
  if (identical(package, "flexyBayes")) {
    stop(
      "The flexyBayes partition of the refusal registry is locked. ",
      "Register a plug-in's codes under the plug-in's own package name.",
      call. = FALSE
    )
  }
  if (!is.data.frame(refusals)) {
    stop(
      "`refusals` must be a data frame with one row per refusal code.",
      call. = FALSE
    )
  }
  required <- c("reason_code", "description", "since_version")
  absent <- setdiff(required, names(refusals))
  if (length(absent) > 0L) {
    stop(
      "`refusals` has no column ",
      paste0("`", absent, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  n <- nrow(refusals)
  column <- function(name, default) {
    x <- if (name %in% names(refusals)) refusals[[name]] else rep(default, n)
    if (is.logical(x) && all(is.na(x))) {
      x <- as.character(x)
    }
    if (!is.character(x)) {
      stop(
        "Column `", name, "` of `refusals` must be character.",
        call. = FALSE
      )
    }
    x
  }
  fields <- list(
    reason_code = column("reason_code"),
    description = column("description"),
    message_template = column("message_template", "%s"),
    registered_in_adr = column("registered_in_adr", NA_character_),
    plan_field = column("plan_field", NA_character_),
    since_version = column("since_version")
  )
  required <- c(
    "reason_code", "description", "message_template", "since_version"
  )
  for (name in required) {
    if (any(is.na(fields[[name]]) | !nzchar(fields[[name]]))) {
      stop(
        "Column `", name, "` of `refusals` has an empty or missing value.",
        call. = FALSE
      )
    }
  }

  codes <- fields$reason_code
  malformed <- codes[!grepl("^[a-z][a-z0-9_]*$", codes)]
  if (length(malformed) > 0L) {
    stop(
      "Each reason_code must be lower-case letters, digits and ",
      "underscores, starting with a letter: ",
      paste0("'", malformed, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  repeated <- unique(codes[duplicated(codes)])
  if (length(repeated) > 0L) {
    stop(
      "`refusals` gives reason_code ",
      paste0("'", repeated, "'", collapse = ", "),
      " more than once.",
      call. = FALSE
    )
  }

  entries <- lapply(seq_len(n), function(i1) {
    list(
      reason_code = fields$reason_code[[i1]],
      description = fields$description[[i1]],
      message_template = fields$message_template[[i1]],
      registered_in_adr = fields$registered_in_adr[[i1]],
      plan_field = fields$plan_field[[i1]],
      since_version = fields$since_version[[i1]],
      package = package
    )
  })
  is_new <- logical(n)
  for (i1 in seq_len(n)) {
    held <- .lookup_refusal(codes[[i1]])
    if (is.null(held)) {
      is_new[[i1]] <- TRUE
    } else if (!identical(held$package, package)) {
      stop(
        "reason_code '", codes[[i1]], "' is already registered by ",
        held$package, ". Choose a code particular to ", package, ".",
        call. = FALSE
      )
    } else if (!identical(held, entries[[i1]])) {
      stop(
        "reason_code '", codes[[i1]], "' is already registered by ",
        package, " with different fields. Restart R to load the ",
        "changed registration.",
        call. = FALSE
      )
    }
  }

  if (any(is_new)) {
    if (!exists(package, envir = .refusal_plugin_registry, inherits = FALSE)) {
      assign(
        package, new.env(parent = emptyenv()),
        envir = .refusal_plugin_registry
      )
    }
    partition <- get(
      package, envir = .refusal_plugin_registry, inherits = FALSE
    )
    for (i1 in which(is_new)) {
      assign(codes[[i1]], entries[[i1]], envir = partition)
      lockBinding(codes[[i1]], partition)
    }
  }

  invisible(fb_refusals(reason_code = codes))
}


#' Raise a registered refusal
#'
#' `fb_refuse()` stops with a refusal condition for a reason registered by
#' flexyBayes or by a package through [fb_register_refusals()]. The
#' condition has the classes `flexybayes_refusal_<reason_code>`, any
#' classes given in `class`, and `flexybayes_refusal`, so users catch a
#' plug-in's refusals exactly as they catch the core's.
#'
#' @param reason_code A registered reason code (see [fb_refusals()]).
#' @param message The message shown to the user. Defaults to the reason's
#'   registered description. A message that names the user's own
#'   variables or values is usually more helpful.
#' @param ... Further named fields to carry on the condition object, for
#'   handlers that need the offending value. The names `message`,
#'   `reason_code`, `package` and `call` are set by `fb_refuse()` and
#'   cannot be used.
#' @param class Additional condition classes, placed between the per-code
#'   class and `flexybayes_refusal`.
#' @param call The call to record on the condition. Defaults to `NULL`,
#'   as for the core's refusals.
#'
#' @return Does not return. It signals an error. The condition carries
#'   `message`, `reason_code`, `package` (the package that registered the
#'   code), `call` and the fields given in `...`.
#'
#' @seealso [fb_register_refusals()], [fb_refusals()].
#'
#' @examples
#' res <- tryCatch(
#'   fb_refuse("unsupported_family"),
#'   flexybayes_refusal = function(e) e
#' )
#' res$reason_code
#' res$package
#'
#' @export
fb_refuse <- function(
  reason_code,
  message = NULL,
  ...,
  class = character(0),
  call = NULL
) {
  if (
    !is.character(reason_code) ||
      length(reason_code) != 1L ||
      is.na(reason_code) ||
      !nzchar(reason_code)
  ) {
    stop("`reason_code` must be a single non-empty string.", call. = FALSE)
  }
  entry <- .lookup_refusal(reason_code)
  if (is.null(entry)) {
    stop(
      "reason_code '", reason_code, "' is not registered. A package ",
      "registers its codes with fb_register_refusals() in its .onLoad() ",
      "before raising them.",
      call. = FALSE
    )
  }
  message <- message %||% entry$description
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    stop("`message` must be a single string.", call. = FALSE)
  }
  if (!is.character(class) || anyNA(class)) {
    stop("`class` must be a character vector.", call. = FALSE)
  }
  fields <- list(...)
  if (
    length(fields) > 0L &&
      (is.null(names(fields)) || any(!nzchar(names(fields))))
  ) {
    stop("Every field passed in `...` must be named.", call. = FALSE)
  }
  reserved <- intersect(
    names(fields),
    c("message", "reason_code", "package", "call")
  )
  if (length(reserved) > 0L) {
    stop(
      "fb_refuse() sets ",
      paste0("`", reserved, "`", collapse = ", "),
      " itself, give the extra field another name.",
      call. = FALSE
    )
  }
  cnd <- do.call(
    .fb_refusal_condition,
    c(
      list(
        reason_code = reason_code,
        message = message,
        family_class = class,
        call = call
      ),
      fields
    )
  )
  stop(cnd)
}


# The package whose namespace `env` belongs to, or NULL outside one.
.fb_calling_package <- function(env) {
  top <- topenv(env)
  if (!isNamespace(top) || isBaseNamespace(top)) {
    return(NULL)
  }
  unname(getNamespaceName(top))
}

# Removes a plug-in partition. Used by tests to clean up fake plug-ins.
.unregister_refusal_partition <- function(package) {
  if (exists(package, envir = .refusal_plugin_registry, inherits = FALSE)) {
    rm(list = package, envir = .refusal_plugin_registry)
  }
  invisible(NULL)
}


# 0.10.0.9000, the transformed response. A left-hand side that is an
# expression over columns instead of a column name is now accepted on the
# per-row route (issue #14, fit1), and the streaming route, which reads named
# columns off disk, refuses it by name instead of asking a reader for a column
# that does not exist.
.populate_refusal_registry_v01000_9000 <- function() {
  .register_refusal(
    reason_code = "stream_response_must_be_a_column",
    description = paste0(
      "Streaming aggregation reads named columns from the source chunk by ",
      "chunk, so a transformed response, log(yield), has no column to ",
      "read. The per-row route evaluates the expression against the data. ",
      "The streaming route cannot, and says so instead of failing inside ",
      "the reader."
    ),
    message_template = paste0(
      "Streaming aggregation needs the response to be a column on disk, ",
      "and \"%s\" is an expression. Materialise it as a column in the ",
      "source before streaming, or fit the per-row route, which evaluates ",
      "the expression."
    ),
    registered_in_adr = "0.10.0.9000 transformed response",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "nested_variable_not_in_data",
    description = paste0(
      "An interaction random effect names a variable that is not a column ",
      "of the data, so the interaction index cannot be built."
    ),
    message_template = paste0(
      "Interaction random effect `%s` names a variable not in data: %s."
    ),
    registered_in_adr = "0.10.0.9000 nested on INLA",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "us_gxe_prior_not_translatable_inla",
    description = paste0(
      "us(f):g on INLA is an iidKd field carrying a joint Wishart prior over ",
      "the whole covariance. It has no standard-deviation slot, so a ",
      "supplied SD prior cannot be translated. Refused instead of dropped."
    ),
    message_template = paste0(
      "A prior was supplied for `%s`, and us() on INLA carries a joint ",
      "Wishart prior over the whole covariance. Fit with backend = \"brms\", ",
      "or drop the prior to accept INLA's default."
    ),
    registered_in_adr = "0.10.0.9000 us() on INLA",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "nested_index_name_collision",
    description = paste0(
      "Fitting an interaction random effect on INLA generates an index ",
      "column, and the data already carries a column of that name. ",
      "Overwriting it would index the field by whatever that column holds ",
      "with nothing downstream saying so, which is why this refuses ",
      "instead of renaming silently."
    ),
    message_template = paste0(
      "Fitting the interaction random effect `%s` on INLA needs an index ",
      "column named `%s`, and `data` already has one. Rename that column."
    ),
    registered_in_adr = "0.10.0.9000 nested on INLA",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "cor_slope_prior_not_translatable_inla",
    description = paste0(
      "(x | g) on INLA is an iid2d field whose hyperparameters are a single ",
      "joint Wishart over the 2-by-2 covariance of the intercept and the ",
      "slope. It has no standard-deviation slot, so a supplied SD prior ",
      "cannot be translated. Refused instead of dropped."
    ),
    message_template = paste0(
      "A prior was supplied for `%s`, and (x | g) on INLA carries a joint ",
      "Wishart prior over the whole 2-by-2 covariance. Fit with ",
      "backend = \"brms\", or write the model as `(x || g)`, whose two ",
      "independent precisions do take an SD prior."
    ),
    registered_in_adr = "0.10.0.9000 correlated slope on INLA",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "vpredict_no_functional",
    description = paste0(
      "vpredict() was called with no functional, or with an argument that ",
      "is not a two-sided `name ~ expression` formula."
    ),
    message_template = paste0(
      "vpredict() needs at least one functional, written as a two-sided ",
      "formula: vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2))."
    ),
    registered_in_adr = "0.10.0.9000 posterior functionals",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "vpredict_component_not_in_posterior",
    description = paste0(
      "A functional names a variance component the fit does not carry. ",
      "Resolving the name against the caller's workspace instead would ",
      "return a number computed from something other than this posterior, ",
      "so it refuses and lists the components that are available."
    ),
    message_template = paste0(
      "`%s` names a component this fit does not carry. It carries: %s."
    ),
    registered_in_adr = "0.10.0.9000 posterior functionals",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "vpredict_functional_not_one_value_per_draw",
    description = paste0(
      "A functional returned something other than one number per posterior ",
      "draw, so it is not a transform of the posterior and has no ",
      "posterior of its own to summarise."
    ),
    message_template = paste0(
      "`%s` must evaluate to one number per posterior draw."
    ),
    registered_in_adr = "0.10.0.9000 posterior functionals",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "brms_stan_toolchain_compile_failure",
    description = paste0(
      "The brms path compiles a Stan program, and a broken C++ toolchain ",
      "reaches R as a downstream symptom, `invalid connection`, not ",
      "as the compiler's own message. On a failed fit the toolchain ",
      "is probed directly, and this code is raised only when the probe ",
      "confirms it cannot compile."
    ),
    message_template = paste0(
      "brms::brm() failed, and this machine cannot compile a three-line ",
      "C++ file either, so the fault is in the C++ toolchain, not in the ",
      "model. backend = \"inla\" compiles nothing and needs no toolchain."
    ),
    registered_in_adr = "0.10.0.9000 Stan toolchain diagnostic",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  invisible(NULL)
}


# 0.10.0.9000 grammar contract, the special-form argument shape.
#
# A grammar fuzz corpus over all 22 ASReml special forms found four ways a
# malformed call escaped this taxonomy: an untyped base error from a branch's
# own indexing, a silent lowering to a plausible neighbour, a term carrying no
# variable, and a silently discarded surplus argument. The contract in
# R/grammar_contract.R closes all four at one locus. These are the codes it
# raises. See also .fb_check_grammar_arity() and .fb_literal_count().
# 2026-09-13: INLA containment. Estimated correlation matrices are fitted
# on brms only. An offset() in the fixed formula,
# which was dropped without a message on both engines, refuses at ingest.
.populate_refusal_registry_inla_containment <- function() {
  .register_refusal(
    reason_code = "correlated_covariance_requires_brms",
    description = paste0(
      "An estimated correlation matrix (a correlated random slope (x | g) ",
      "or an unstructured us(f):g) was sent to INLA. INLA carries it as an ",
      "iid2d or iidKd field under a joint Wishart prior that couples the ",
      "variances and the correlations, and on that route the fitted ",
      "covariance did not match REML's and was not invariant to rescaling ",
      "the response. flexyBayes fits these structures on brms."
    ),
    message_template = paste0(
      "%s estimates a correlation matrix, and flexyBayes fits estimated ",
      "correlation matrices on brms only. Pass backend = \"brms\", or leave ",
      "backend = \"auto\", which routes there."
    ),
    registered_in_adr = "2026-09-13 estimated correlation matrices",
    plan_field = "rejected_routes",
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "fixed_offset_unsupported",
    description = paste0(
      "An offset() term in the fixed formula on a route that cannot carry ",
      "it. The parser records every offset on the model description, and a ",
      "route with no exact lowering for it refuses by name instead of ",
      "fitting the model without the offset, which is what happened before ",
      "offsets were recorded."
    ),
    message_template = paste0(
      "The fixed formula carries %s, and %s does not carry an offset."
    ),
    registered_in_adr = "2026-09-13 fixed-formula offsets refused",
    plan_field = "rejected_routes",
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "fixed_offset_invalid",
    description = paste0(
      "An offset() term that does not evaluate against the data to a ",
      "finite number on every row: a name that is not a column, a factor or ",
      "character value, a length that is not the row count, or a missing or ",
      "infinite value. An offset enters the linear predictor with a fixed ",
      "coefficient of one, so an NA would be read as zero by one engine and ",
      "drop the row on the other."
    ),
    message_template = paste0(
      "The fixed formula carries %s, and %s. An offset has to be a finite ",
      "number on each row of the data."
    ),
    registered_in_adr = "0.11.0.9000 fixed-formula offsets",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "newdata_offset_missing",
    description = paste0(
      "predict() on new data from a model with an offset, where the new data ",
      "lack a column the offset reads. The offset is part of every ",
      "prediction and is never taken to be zero."
    ),
    message_template = paste0(
      "The model carries %s, and `newdata` has no column %s. Add the column ",
      "to `newdata`."
    ),
    registered_in_adr = "0.11.0.9000 fixed-formula offsets",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  invisible(NULL)
}

# 2026-09-15: Gaussian aggregation containment. The aggregated Gaussian
# route cannot carry a residual prior, so Gaussian models are fitted per row,
# and no fit on the aggregated representation reports a marginal likelihood.
.populate_refusal_registry_model_manifest <- function() {
  .register_refusal(
    reason_code = "prior_family_has_no_residual_scale",
    description = paste0(
      "A residual prior, stated or the package default, on a family whose ",
      "likelihood has no residual scale parameter: Poisson and binomial, ",
      "and on brms every family whose dispersion is a function of the ",
      "mean. No engine slot reads it, and the model manifest records the ",
      "input as not consumed under this reason."
    ),
    message_template = "%s",
    registered_in_adr = "Compiled model manifest (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "prior_vc_sd_superseded_by_prior",
    description = paste0(
      "`prior_vc_sd` supplied beside an fb_prior(). No route reads the ",
      "scalar once an fb_prior() is present: INLA gives a standard ",
      "deviation the prior leaves unstated the package's bounded default, ",
      "and brms keeps its own. The model manifest records the scalar as ",
      "not consumed under this reason."
    ),
    message_template = "%s",
    registered_in_adr = "Compiled model manifest (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "evidence_not_computed",
    description = paste0(
      "The marginal likelihood of a fit whose engine computed none: a brms ",
      "fit, where Stan sampling runs no bridge sampling. summary() and the ",
      "model manifest report the evidence as not available under this ",
      "reason."
    ),
    message_template = "%s",
    registered_in_adr = "Compiled model manifest (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
}

.populate_refusal_registry_prior_support <- function() {
  .register_refusal(
    reason_code = "prior_support_violated",
    description = paste0(
      "An INLA fit whose mode, or whose reported 2.5 or 97.5 per cent ",
      "quantile, for a hyperparameter under a bounded prior lies outside ",
      "that bound. The uniform prior on a standard deviation reaches INLA ",
      "as a density of -1e10 outside its support, and INLA can still stop ",
      "there. Every INLA fit is checked against INLA's own record of the ",
      "priors it received, and a fit outside its prior's support is refused."
    ),
    message_template = "%s",
    registered_in_adr = "2026-09-16 universal support audit",
    plan_field = NA_character_,
    since_version = "0.12.1"
  )
  invisible(NULL)
}

# .populate_refusal_registry_stream_shape(), flexybayes_stream()'s own
# model-shape and data-availability refusals, untyped before 0.13.0
# (2026-09-16). The route-availability boundary reuses
# aggregation_route_unavailable at its call site (R/stream_aggregate.R,
# .fb_stream_aggregate()): an intercept-only model has no factor cell-key
# term to aggregate by, which is the same "nothing routes into the
# aggregated representation" fact as a family or link the streamed path
# does not carry, and the remedy is identical, fit the model per row
# with flexybayes(). The two codes below cover the two boundaries that
# fact does not reach: a named model column absent from the source, and
# a source that returns no rows to build from or fold.
.populate_refusal_registry_stream_shape <- function() {
  .register_refusal(
    reason_code = "stream_model_column_not_in_source",
    description = paste0(
      "flexybayes_stream() reads named columns off the source chunk by ",
      "chunk, so every column the model names (the response, a fixed ",
      "or random term, `trials`, or `exposure`) has to be present in ",
      "the source's head sample. Untyped before 0.13.0."
    ),
    message_template = "%s",
    registered_in_adr = "2026-09-16 streaming model shape",
    plan_field = NA_character_,
    since_version = "0.13.0"
  )
  .register_refusal(
    reason_code = "stream_source_no_rows",
    description = paste0(
      "flexybayes_stream() found no rows to work with: either the head ",
      "sample used to build the model, or the full read across every ",
      "chunk of the source, returned zero rows. There is nothing to ",
      "aggregate. Untyped before 0.13.0."
    ),
    message_template = "%s",
    registered_in_adr = "2026-09-16 streaming model shape",
    plan_field = NA_character_,
    since_version = "0.13.0"
  )
  invisible(NULL)
}

# .populate_refusal_registry_family_support(), the response-domain gate
# (R/na_action.R, .fb_refuse_family_response_domain()). Admitting a family
# onto .fb_family_links() (R/utils.R) was necessary and sufficient to
# reach an engine. Nothing checked the response against the domain the
# family's likelihood is defined over. lognormal (added b979c5d) on a
# non-positive response reached brms::brm() and stopped there untyped
# ("Family 'lognormal' requires response greater than 0."). Gamma,
# hurdle_gamma and beta reached the same untyped brms death, and (bar
# lognormal, brms-only) an unrelated inla_program_failed diagnosis on
# INLA that named neither the family nor the offending row
# (2026-09-16).
.populate_refusal_registry_family_support <- function() {
  .register_refusal(
    reason_code = "family_response_domain_violation",
    description = paste0(
      "The requested family's likelihood is defined over a response ",
      "domain (gamma / lognormal strictly positive, hurdle_gamma non-",
      "negative, beta an open unit interval) that the observed response ",
      "violates. Refused before any backend code runs. A missing (NA) ",
      "response is exempt, it is na_action's concern, not a domain ",
      "violation. Untyped before 0.13.0."
    ),
    message_template = "%s",
    registered_in_adr = "2026-09-16 family support",
    plan_field = NA_character_,
    since_version = "0.13.0"
  )
  invisible(NULL)
}

.populate_refusal_registry_aggregation_containment <- function() {
  .register_refusal(
    reason_code = "aggregation_residual_prior_not_composed",
    description = paste0(
      "A Gaussian model sent to the aggregated representation, by ",
      "aggregate = TRUE or by flexybayes_stream(). The aggregated Gaussian ",
      "route composes the within-cell sum of squares with a fixed gamma ",
      "prior on the residual precision, so a stated or default residual ",
      "prior would not be the prior fitted while prior_summary() reported ",
      "it. Gaussian models are fitted per row. Binomial and Poisson models ",
      "still aggregate."
    ),
    message_template = "%s",
    registered_in_adr = "0.14.0 strict Gaussian containment",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "aggregated_evidence_unavailable",
    description = paste0(
      "The marginal likelihood of a fit on the aggregated representation. ",
      "The engine's value omits constants that depend on the data and, ",
      "under response standardisation, scales by the number of cells ",
      "instead of the number of rows, so it is not the marginal likelihood ",
      "of the rows. summary() reports it as not available, naming this ",
      "reason, instead of printing a number that cannot be compared."
    ),
    message_template = "%s",
    registered_in_adr = "0.14.0 strict Gaussian containment",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  invisible(NULL)
}

.populate_refusal_registry_working_response <- function() {
  .register_refusal(
    reason_code = "response_not_finite",
    description = paste0(
      "A response with NaN or infinite values, usually from a transformation ",
      "such as log(0) or from a data-entry error. They are not missing ",
      "observations, so na_action does not apply to them, and an engine ",
      "handed them either fitted the model without saying so or failed ",
      "without naming the row. The fit refuses before any engine is called."
    ),
    message_template = "%s",
    registered_in_adr = "Working-response resolution (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  invisible(NULL)
}

.populate_refusal_registry_grammar_contract <- function() {
  .register_refusal(
    reason_code = "grammar_arity_too_few",
    description = paste0(
      "An ASReml special form was written with fewer positional arguments ",
      "than it takes. Before the contract these parsed into a term whose ",
      "variable was NA, which reaches the engine as a model nobody wrote."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "grammar_arity_too_many",
    description = paste0(
      "An ASReml special form was written with more positional arguments ",
      "than it reads. The surplus used to be dropped in silence, so a ",
      "structure written into the call never reached the fitted model."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "grammar_named_argument_unknown",
    description = paste0(
      "A special form was given a named argument it does not answer to. An ",
      "unrecognised name used to be ignored, which reads as accepted."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "grammar_argument_duplicated",
    description = paste0(
      "A special form was given the same named argument twice. The contract ",
      "validated the set of names and not their multiplicity, so a doubled ",
      "name was either silently resolved to whichever copy the branch ",
      "happened to index, or escaped as an untyped matcher error."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "dsum_malformed",
    description = paste0(
      "A dsum() body carrying no grouping bar. dsum(~ units | env) requests ",
      "one residual variance per section. Written without the bar there is ",
      "nothing to section on, and the parser used to fall through to a ",
      "homogeneous `units` residual, a plausible neighbour that survives ",
      "convergence checks and prints credible output."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "fa_rank_missing",
    description = paste0(
      "fa() written without its rank. The rank is a number the user chooses, ",
      "so it is named separately from the general arity refusal."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  # fa()'s invalid-rank code has existed since 0.4.0 as the k < 1 floor. The
  # grammar contract widens what reaches it (a symbol, a string, a logical,
  # a fraction) instead of minting a second code for the same condition.
  # The rank-versus-levels bound has a code already, fa_rank_exceeds_dim,
  # registered at 0.4.0 and raised from .enrich() where the level count is
  # known. It was reachable only from the fa(f, k):g interaction form. The
  # bare fa(f, k) branch now raises the same code instead of a synonym.
  .register_refusal(
    reason_code = "grammar_argument_empty",
    description = paste0(
      "A special form written with a trailing or doubled comma, leaving ",
      "R's empty symbol where a value belongs. Touching it raised ",
      "\"argument is missing, with no default\" from wherever it was first ",
      "read, a base error about flexyBayes internals instead of a ",
      "statement about the formula."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "pol_degree_invalid",
    description = paste0(
      "pol()'s degree is not a positive whole number written in the formula. ",
      "Same coercion defect as fa()'s rank."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  .register_refusal(
    reason_code = "multivariate_smooth_unsupported",
    description = paste0(
      "s() written over more than one variable. flexyBayes fits univariate ",
      "penalised splines only. This used to raise an untyped names-attribute ",
      "error from mgcv's call handling."
    ),
    message_template = "%s",
    registered_in_adr = "0.10.0.9000 grammar contract",
    plan_field = "rejected_routes",
    since_version = "0.10.0.9000"
  )
  invisible(NULL)
}

.populate_refusal_registry_aggregated_accessors <- function() {
  .register_refusal(
    reason_code = "streamed_rows_not_retained",
    description = paste0(
      "A row-level accessor (fitted(), residuals(), model.matrix(), ",
      "insight::get_data(), emmeans::recover_data()) on a flexybayes_stream() ",
      "fit that keeps fewer rows than it was fitted to. The stream keeps only ",
      "the first rows it read, so values computed from them are not the ",
      "values of the data."
    ),
    message_template = "%s",
    registered_in_adr = "Aggregated count accessors (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "update_streamed_fit_not_supported",
    description = paste0(
      "update() on a flexybayes_stream() fit. The fit keeps neither the ",
      "source it read nor every row, so a re-fit through flexybayes() would ",
      "run on different data."
    ),
    message_template = "%s",
    registered_in_adr = "Aggregated count accessors (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "augment_not_available_inla",
    description = paste0(
      "augment() on an INLA fit, which keeps no per-row model frame for the ",
      "fitted-value and residual columns. tidy() and predict() are the ",
      "alternatives."
    ),
    message_template = "%s",
    registered_in_adr = "Aggregated count accessors (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "log_density_producer_unavailable",
    description = paste0(
      "fb_log_posterior() on a fit from either engine: flexyBayes has no ",
      "log-density producer for them. INLA's posterior is a Laplace ",
      "approximation, and brms's log-density is defined on the Stan ",
      "unconstrained scale. The condition also carries the class ",
      "fb_c4_unavailable."
    ),
    message_template = "%s",
    registered_in_adr = "Aggregated count accessors (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "predict_link_not_supported",
    description = paste0(
      "predict() without classify on an INLA fit whose link is not the ",
      "identity. Posterior draws from fb_as_draws_simple() are the ",
      "alternative."
    ),
    message_template = "%s",
    registered_in_adr = "Aggregated count accessors (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  .register_refusal(
    reason_code = "marginaleffects_link_not_supported",
    description = paste0(
      "A marginaleffects prediction on a fit whose link is not the identity. ",
      "Posterior draws from fb_as_draws_simple() are the alternative."
    ),
    message_template = "%s",
    registered_in_adr = "Aggregated count accessors (2026-09-15)",
    plan_field = NA_character_,
    since_version = "0.11.0.9000"
  )
  invisible(NULL)
}


# 2026-09-16: the factor-analytic covariance leaves the fitting route.
#
# The structure was admitted on INLA at 0.10.0.9000 on a simulated recovery
# at k = 1 and k = 2. What was never established is that a second identical
# call returns the same covariance. Measured before this refusal was written,
# under the 0.13.0 half-normal default and a tolerance declared in advance:
# on two simulated designs and on agridat::dasilva.maize the identified
# covariance repeated to better than 1.01, and on agridat::besag.met, a trial
# with weak genetic signal, eight identical calls spanned 1.23 in a variance
# and 0.12 in a genetic correlation against tolerances of 1.10 and 0.05.
# From 0.18.0 the refusal rests on simulated trials with a known G, where
# INLA's genetic-correlation intervals missed the coverage required.
#
# The emit is left in place, reachable by a direct emit_inla() call, as the
# us() branch is. What is withdrawn is the claim that the fitted numbers can
# be read as the covariance.
.populate_refusal_registry_factor_analytic <- function() {
  .register_refusal(
    reason_code = "factor_analytic_not_repeatable",
    description = paste0(
      "A factor-analytic genotype-by-environment covariance, fa(f, k):g, ",
      "sent to a fit. INLA represents it by entering each common factor as ",
      "scaled copies of a latent field, but on simulated trials with a ",
      "known covariance the genetic-correlation intervals of INLA fits ",
      "covered the true value less often than required. ",
      "flexyBayes has no factor-analytic lowering for brms, so neither ",
      "engine fits the structure here. The diagonal structure diag(f):g ",
      "fits on both engines, and us(f):g fits on brms."
    ),
    message_template = "%s",
    registered_in_adr = "2026-09-16 factor-analytic covariance",
    plan_field = "rejected_routes",
    since_version = "0.13.0"
  )
  invisible(NULL)
}
