# Refusal registry --- scaffold only at v0.3.8
#
# Container + registration helper for the canonical refusal-reason
# vocabulary. Empty at v0.3.8: this release ships the registry
# infrastructure without converting the ~30 existing refusal sites
# scattered through dispatch.R, lgm_gate.R, structured_cov.R, and
# emit_inla.R. Full population (one .register_refusal() call per
# existing site) + the user-facing fb_refusals() accessor land at
# v0.4.0.
#
# Rationale: the deliberate choice here is to ship the container
# without the migration risk of touching ~30 existing refusal sites
# in one release.
#
# The registry is an environment (not a list) so .register_refusal()
# can write idempotently without copy-on-modify cost; locked at the
# end of .onLoad() so user-side assign() raises a contract violation.

# --- container ---------------------------------------------------- #

# `.refusal_registry` is allocated at namespace load. The parent
# environment is `emptyenv()` so symbol lookup cannot accidentally
# fall through to the package namespace and pick up an unrelated
# binding. Filled by .register_refusal() during .onLoad() (currently
# zero calls --- the scaffold-only posture); locked immediately
# after by .lock_refusal_registry().
.refusal_registry <- new.env(parent = emptyenv())


# --- registration helper ------------------------------------------ #

# .register_refusal() --- the one-shot registration call. Idempotent
# only in the sense that a duplicate reason_code raises rather than
# silently overwriting; intentional refusals (the v0.4.0 migration)
# fire once per reason at .onLoad() time and the registry then
# locks.
#
# Arguments
#   reason_code        canonical machine-readable identifier; the
#                      same string that appears on the structured
#                      condition object thrown by the refusal site.
#                      Must be a single non-empty string.
#   description        one-line human-readable description; what the
#                      refusal means in plain language.
#   message_template   the template string used to render the
#                      message; conventionally uses base::sprintf()
#                      placeholders so the refusal site can fill in
#                      the data-specific slots.
#   registered_in_adr  the ADR number the refusal traces to (e.g.,
#                      "ADR 0024" for the routing-policy refusals,
#                      "ADR 0025" for the structured-cov refusals).
#                      DESIGN_DECISIONS.md indexes every ADR number.
#                      Free-form; this field documents provenance,
#                      it does not gate registration.
#   plan_field         optional name of the <fb_plan> slot that
#                      surfaces this refusal at plan time (e.g.,
#                      "rejected_routes"). NA_character_ if the
#                      refusal is not surfaced by fb_plan() (e.g.,
#                      runtime-only refusals).
#   since_version      the flexyBayes version the refusal was
#                      introduced; useful for backward-compat
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
      ".register_refusal(): refusal registry is locked; new ",
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
      "already registered. Refusal vocabulary is append-only; ",
      "use a distinct reason_code for new refusals.",
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
      since_version = since_version
    ),
    envir = .refusal_registry
  )

  invisible(NULL)
}


# --- accessor ----------------------------------------------------- #

# .lookup_refusal() --- internal accessor for refusal sites that
# want to consult the registry for a description or message template.
# Returns the registered entry list, or NULL if the reason_code is
# not yet registered (the v0.3.8 scaffold-only state returns NULL
# for every code). Callers must tolerate NULL --- the legacy
# inline-message refusal sites continue to function untouched.
.lookup_refusal <- function(reason_code) {
  if (
    !is.character(reason_code) ||
      length(reason_code) != 1L ||
      !nzchar(reason_code)
  ) {
    return(NULL)
  }
  if (!exists(reason_code, envir = .refusal_registry, inherits = FALSE)) {
    return(NULL)
  }
  get(reason_code, envir = .refusal_registry, inherits = FALSE)
}


# --- unified refusal-condition constructor ----------------------- #

# .fb_refusal_condition() --- the single constructor every
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
# the preferred test contract; the retained family
# class keeps every pre-existing class-based handler and assertion
# working unchanged (additive migration).
#
# Registry gate: `reason_code` MUST be registered (via
# .register_refusal() at .onLoad()). An unregistered code is a
# flexyBayes internal contract violation --- the refusal does not
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
#   reason_code   canonical machine-readable identifier; must be in
#                 the locked .refusal_registry.
#   message       the fully-rendered, site-specific message string.
#   family_class  zero or more legacy family-class strings to retain
#                 below the per-code class (e.g.
#                 "flexybayes_structured_cov_refusal"). Character
#                 vector; default character(0).
#   call          the condition call carried on the object; defaults
#                 to NULL (matching every site except the preflight
#                 wrapper, which threads the originating call through).
#   ...           extra named fields carried on the condition object
#                 (e.g. supplied, format, scheme, backend, factor)
#                 exactly as the pre-migration sites carried them.
#
# Returns the condition object; the caller raises it via stop().
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

  if (is.null(.lookup_refusal(reason_code))) {
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
      list(message = message, reason_code = reason_code, call = call),
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

# .populate_refusal_registry_v0310() --- the first migration of
# existing refusal sites into the v0.3.8 scaffold. Called from
# .onLoad() *before* .lock_refusal_registry() so the entries are
# in place when user code first observes the registry.
#
# Three reason codes register at v0.3.10:
#
#   block_partition_incomplete         --- blocks structural refusal:
#                                          sum of block sizes !=
#                                          grouping factor level count.
#   block_not_positive_definite        --- blocks structural refusal:
#                                          at least one V_k fails the
#                                          PD probe; the offending
#                                          block index is named in
#                                          the message text.
#   approximate_route_not_yet_registered --- low-rank refusal-stub
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

# .populate_refusal_registry_v0400() --- the bulk migration of the
# remaining 28 user-facing refusal codes into the registry.
# Together with the 3 v0.3.10 codes the registry now
# holds the complete user-facing refusal vocabulary (31 codes); the
# routing-decision reasons (backend_decision trace) and the internal
# aggregate-out-of-scope control-flow signals are deliberately NOT
# registered --- they are not refusals the user can hit, and listing
# them in fb_refusals() would misrepresent the surface.
#
# message_template here is documentation: a single representative
# rendering surfaced by fb_refusals(). It is intentionally NOT the
# runtime message --- five codes are raised at multiple sites with
# different bodies, and each site keeps its verbatim message via
# .fb_refusal_condition(). plan_field follows the v0.3.10 precedent:
# the representation-level (structured covariance) refusals carry
# "representation_plan"; runtime-only refusals carry NA_character_.
.populate_refusal_registry_v0400 <- function() {
  # -- structured-covariance carriers (structured_cov.R,
  #    emit_inla.R) --------------------------------------------------
  .register_refusal(
    reason_code = "vm_redundant_specification",
    description = paste0(
      "vm()/ped(): more than one covariance carrier supplied; ",
      "exactly one of V / chol / precision / blocks / ",
      "low_rank_factor is allowed."
    ),
    message_template = paste0(
      "%s(): exactly one covariance carrier may be supplied; got ",
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
      "vm(%s, chol = %s): the Cholesky factor must be square; got ",
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
      "vm(%s, precision = %s): the precision matrix must be square; ",
      "got dim %s."
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
      "of K covariance matrices; got class '%s'."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "blocks_empty_list",
    description = "vm(..., blocks = ): the block list is empty.",
    message_template = paste0(
      "vm(%s, blocks = %s): the block list is empty; expected at ",
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
      "precision structure; rownames and colnames must be identical ",
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
      "names; rename to disambiguate (INLA's data-list lookup is ",
      "name-keyed)."
    ),
    registered_in_adr = "ADR 0025",
    plan_field = "representation_plan",
    since_version = "0.4.0"
  )

  # -- approximation surface (parse_formula.R,
  #    validate_approximation.R, emit_smooth_low_rank.R, dispatch.R) -
  .register_refusal(
    reason_code = "low_rank_requires_greta",
    description = paste0(
      "A smooth requesting the low_rank_smooth approximation was ",
      "routed to a non-greta backend that cannot honour it."
    ),
    message_template = paste0(
      "A smooth requesting the low_rank_smooth approximation ",
      "requires the greta backend; the '%s' backend represents ",
      "smooths differently (INLA via rw2, brms via Stan spline ",
      "bases) and cannot honour the rank-K basis truncation. ",
      "Re-fit with backend = \"greta\", or drop the ",
      "representation = ... argument to fit the exact smooth."
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
      "`scheme`; e.g. representation = list(scheme = ",
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
      "has no smooth-basis emit path at this release; the only ",
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
      "as low_rank_smooth); an exact fit has no bias bound to ",
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
      "low_rank_smooth rank%s must be a single positive integer; ",
      "got %s."
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
      "the basis dimension reproduces the exact basis; request the ",
      "exact smooth (drop the approximation) instead, or choose ",
      "rank <= %d."
    ),
    registered_in_adr = "ADR 0027",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  # (The dispatch-side approximate-route refusal reuses the existing
  # v0.3.10 code approximate_route_not_yet_registered -- one code,
  # two site-specific messages; nothing new to register here.)

  # -- aggregated-Gaussian emit (emit_gaussian_aggregated.R)
  .register_refusal(
    reason_code = "heterogeneous_residual_factor_not_in_cell_key",
    description = paste0(
      "Aggregated Gaussian emit: an at(f):units heterogeneous ",
      "residual factor is not in the cell key, so the ",
      "cell-constant sigma property does not hold."
    ),
    message_template = paste0(
      "emit_gaussian_aggregated(): heterogeneous residual ",
      "at(%s):units refused -- '%s' is not in the cell key {%s}. ",
      "The cell-constant sigma property does not hold, so the ",
      "per-row / per-cell algebraic identity breaks. Pass ",
      "aggregate = FALSE for the per-row path."
    ),
    registered_in_adr = "ADR 0022",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "rcov_type_unsupported_for_aggregation",
    description = paste0(
      "Aggregated Gaussian emit: the rcov term type is outside the ",
      "supported aggregation scope."
    ),
    message_template = paste0(
      "emit_gaussian_aggregated(): rcov term type '%s' is not ",
      "supported by the aggregated path. Only homogeneous (units / ",
      "id) and at_units heterogeneous residual are supported. Pass ",
      "aggregate = FALSE for the per-row path."
    ),
    registered_in_adr = "ADR 0022",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # -- prediction kernel (predict_kernel.R) ---------------
  .register_refusal(
    reason_code = "predict_kernel_invalid_include",
    description = paste0(
      "predict(): `include` is empty or carries values outside the ",
      "prediction-kernel vocabulary."
    ),
    message_template = paste0(
      ".predict_linear_draws(): `include` must be a non-empty ",
      "character vector drawn from the kernel vocabulary."
    ),
    registered_in_adr = "ADR 0023",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # -- preflight refusals (fb_preflight.R, dispatch.R) -----
  .register_refusal(
    reason_code = "design_memory_exceeds_ceiling",
    description = paste0(
      "Preflight: the design matrix is estimated to exceed the ",
      "active memory ceiling; dispatch is short-circuited before ",
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
      "was short-circuited before any backend code ran -- raising ",
      "the memory ceiling will not help."
    ),
    registered_in_adr = "ADR 0024",
    plan_field = "preflight",
    since_version = "0.4.0"
  )

  # -- family support (utils.R::.resolve_family) ---------------------
  # .resolve_family() is the single family gate every user entry
  # passes through (asreml via fb.R, brms via fb_from_brms.R). It
  # admits only the families flexyBayes can emit; any other family --
  # including those INLA's likelihood roster recognises but flexyBayes
  # has no emit path for (e.g. survival / time-to-event) -- is refused
  # up front, never silently fitted.
  .register_refusal(
    reason_code = "unsupported_family",
    description = paste0(
      "The requested family is outside the set flexyBayes can emit. ",
      "Refused at the family gate (.resolve_family) before any ",
      "backend code runs."
    ),
    message_template = paste0(
      "Unsupported family \"%s\". flexyBayes supports: gaussian, ",
      "binomial, binary, poisson, negative_binomial, negbinom, ",
      "gamma, beta. Other families (including survival / ",
      "time-to-event and multivariate responses) are planned future ",
      "work; see fb_refusals()."
    ),
    registered_in_adr = "audit-2026-05-30",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # -- formula / spec parse-time refusals (parse_formula.R) ----------
  # Early model-spec and data-binding refusals raised by the formula
  # parser. Migrated from raw stop() to the structured registry so
  # downstream tooling pattern-matches on the condition class rather
  # than the free-text message. Missing-suggested-package errors
  # (mgcv, greta) and argument-combination guards stay as plain
  # stop() --- they are environment / call-shape errors, not
  # model-scope refusals.
  .register_refusal(
    reason_code = "formula_not_two_sided",
    description = paste0(
      "The model formula must be two-sided (response ~ predictors); ",
      "a formula carrying no left-hand-side response was supplied."
    ),
    message_template = paste0(
      "The model formula must be two-sided: response ~ predictors."
    ),
    registered_in_adr = "audit-2026-05-30",
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
    registered_in_adr = "audit-2026-05-30",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "fa_rank_invalid",
    description = paste0(
      "A factor-analytic term fa(x, k) was given a rank k below 1; ",
      "the factor-analytic rank must be a positive integer."
    ),
    message_template = paste0(
      "fa() requires k >= 1; got k = %s."
    ),
    registered_in_adr = "audit-2026-05-30",
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
      "fa() requires k < the number of levels of the outer factor; ",
      "got k = %s with %s level(s)."
    ),
    registered_in_adr = "audit-2026-06-06",
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
    registered_in_adr = "audit-2026-05-30",
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
      "flexyBayes fits univariate penalised splines only -- use ",
      "s(x) or spl(x) per smooth dimension."
    ),
    registered_in_adr = "audit-2026-05-30",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # -- entry-point argument guards (flexybayes.R) -------------------
  # Call-shape guards on the review_code / return_code flags. They are
  # user-reachable refusals (the registry's inclusion criterion), so
  # they carry catchable classes even though they are argument-
  # combination errors rather than model-scope refusals.
  .register_refusal(
    reason_code = "review_code_backend_unsupported",
    description = paste0(
      "review_code = TRUE was requested with a backend other than ",
      "greta; the inspect-then-fit token is currently greta-only."
    ),
    message_template = paste0(
      "`review_code = TRUE` is currently supported only with backend ",
      "= \"greta\". Pass backend = \"greta\" or drop review_code."
    ),
    registered_in_adr = "audit-2026-05-30",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )
  .register_refusal(
    reason_code = "code_flags_mutually_exclusive",
    description = paste0(
      "return_code = TRUE and review_code = TRUE were both supplied; ",
      "the two code-return modes are mutually exclusive."
    ),
    message_template = paste0(
      "`return_code` and `review_code` are mutually exclusive."
    ),
    registered_in_adr = "audit-2026-05-30",
    plan_field = NA_character_,
    since_version = "0.4.0"
  )

  # -- fb_cov() constructor ----------
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

  # -- backend capability gates (dispatch.R) -------
  # Raised when flexybayes(backend = "brms") meets an asreml structured-
  # covariance term (vm/ped/fa/us/ar1) that Stan cannot represent. The
  # brms capability predicate (R/backend_registry.R, .capability_brms)
  # is the declarative source; the dispatch brms branch raises this code.
  .register_refusal(
    reason_code = "stan_cannot_represent_structured_cov",
    description = paste0(
      "backend = \"brms\" (Stan) cannot represent an asreml structured-",
      "covariance term (vm/ped/fa/us/ar1)."
    ),
    message_template = paste0(
      "backend = \"brms\" (Stan) cannot represent this model: it ",
      "contains an asreml structured-covariance term (one of ",
      "vm/ped/fa/us/ar1) with no lossless brms/Stan translation. Re-fit ",
      "with backend = \"greta\", or backend = \"inla\" when latent-",
      "Gaussian feasible."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.5.0"
  )

  # -- grammar polymorphism on the universal entry --
  # fb() / flexybayes() detect the input grammar from the call shape.
  # These three guard the brms-grammar and reserved-greta branches of
  # .build_ir_polymorphic() (R/fb.R).
  .register_refusal(
    reason_code = "grammar_brms_with_asreml_terms",
    description = paste0(
      "A brms-style bar-grouped formula was combined with ASReml ",
      "`random` / `rcov` arguments on the universal entry."
    ),
    message_template = paste0(
      "A brms-style formula (with `(... | g)` grouping) cannot be ",
      "combined with `random` / `rcov` (ASReml grammar). Put every ",
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
  # The universal entry now FITS a native greta model graph (v0.5.0).
  # The deferral the v0.4.1 grammar_greta_via_fb_deferred code
  # described is gone -- it is removed from the vocabulary. Two
  # genuine refusals replace it: a native graph requested on a non-greta
  # engine, and an engine pin handed a conflicting `backend`.
  .register_refusal(
    reason_code = "native_greta_requires_greta_backend",
    description = paste0(
      "A native greta model graph was passed to the universal entry / the ",
      "greta pin with a non-greta backend. A native graph is greta-only ",
      "by construction."
    ),
    message_template = paste0(
      "A native greta model graph is fit by greta::mcmc() and is ",
      "greta-only by construction; the requested backend cannot fit it. ",
      "Use backend = \"greta\" (or the default), or rebuild the model in ",
      "the ASReml / brms formula grammar to reach another engine."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = "rejected_routes",
    since_version = "0.5.0"
  )
  .register_refusal(
    reason_code = "engine_pin_backend_conflict",
    description = paste0(
      "An engine pin (fb_greta / fb_inla / fb_brms) was given a `backend` ",
      "argument that conflicts with the engine it pins."
    ),
    message_template = paste0(
      "An engine pin fits via one engine only, so it cannot take a ",
      "conflicting `backend`. Drop `backend`, or use fb() / flexybayes() ",
      "to choose a backend."
    ),
    registered_in_adr = "ADR 0031",
    plan_field = NA_character_,
    since_version = "0.5.0"
  )

  invisible(NULL)
}

# gretaR backend activation: the refusal codes the gretaR
# backend (R/emit_gretaR.R + the dispatch branch) can raise. Registered
# alongside the others before the registry locks (see R/zzz.R).
.populate_refusal_registry_gretaR <- function() {
  reg <- function(code, desc, tmpl) {
    .register_refusal(
      reason_code = code,
      description = desc,
      message_template = tmpl,
      registered_in_adr = "ADR 0013/0031",
      plan_field = NA_character_,
      since_version = "0.6.0.9000"
    )
  }
  reg(
    "gretaR_cannot_represent_structured_cov",
    "gretaR backend: structured covariance (vm/ped/fa/us/ar1) unsupported.",
    "backend = \"gretaR\" cannot fit this model (%s)."
  )
  reg(
    "gretaR_random_term_type_unsupported",
    "gretaR backend: only random-intercept-class random terms are supported.",
    "backend = \"gretaR\" cannot fit this model (%s)."
  )
  reg(
    "gretaR_family_unsupported",
    "gretaR backend: family outside gaussian/binomial/poisson.",
    "%s"
  )
  reg(
    "gretaR_random_group_not_in_data",
    "gretaR backend: random-intercept grouping factor absent from data.",
    "%s"
  )
  reg(
    "gretaR_below_version_floor",
    "gretaR backend: the installed gretaR is older than the activation floor.",
    "%s"
  )
  reg(
    "gretaR_not_installed",
    "gretaR backend: gretaR not installed and no source home set.",
    "%s"
  )
  invisible(NULL)
}


# --- lock helper -------------------------------------------------- #

# .lock_refusal_registry() --- locks the environment so no further
# .register_refusal() calls or user-side assign()s can mutate it.
# Called once at the end of .onLoad(). The lock is enforced by R's
# environmentIsLocked() machinery; bindings registered before the
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
#' flexyBayes declines to fit, route, or validate a model --- each
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
#' release by version-string prefix --- `since_version = "0.4"`
#' returns every reason added in the 0.4 series.
#'
#' Routing-decision reasons (surfaced by [fb_plan()] and
#' [backend_decision()]) and internal control-flow signals are
#' deliberately excluded: this table lists only refusals a user can
#' actually encounter.
#'
#' @param reason_code Optional character vector of exact reason codes
#'   to filter to. `NULL` (default) returns all registered reasons.
#' @param since_version Optional single version-string prefix to
#'   filter to. `NULL` (default) returns all.
#'
#' @return A data frame of subclass `fb_refusals_table`, one row per
#'   matching refusal reason, with columns `reason_code`,
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

  codes <- ls(.refusal_registry)
  entries <- lapply(codes, .lookup_refusal)

  df <- data.frame(
    reason_code = vapply(entries, `[[`, character(1L), "reason_code"),
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
    cat(sprintf("  [since %s] %s\n", x$since_version[[i]], x$reason_code[[i]]))
    cat(strwrap(x$description[[i]], width = 72L, prefix = "      "), sep = "\n")
    cat("\n")
  }
  invisible(x)
}
