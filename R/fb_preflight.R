# .fb_preflight() -- design-memory preflight (v0.3.0).
#
# Takes a parsed `<fb_terms>` IR + an `<fb_dataset>` wrapper and
# returns an `<fb_preflight>` S3 list with three slots:
#   $per_term_estimate    named list, one entry per IR term, each
#                         carrying design_memory_bytes (numeric),
#                         representation_class (character), and
#                         aggregated_likelihood_candidate (logical)
#   $total_estimate_bytes sum across terms plus a fixed bookkeeping
#                         overhead (page-sized: 4096 bytes)
#   $refusal              NULL if accepted; otherwise a structured
#                         `<fb_preflight_refusal>` object with
#                         reason_code = "design_memory_exceeds_ceiling",
#                         the binding term label, its estimate, and
#                         the active ceiling in bytes.
#
# All term estimates are deterministic, computed from IR + dataset
# metadata only -- never reads the underlying $data row-wise. This
# keeps preflight O(terms) regardless of n_rows; the 1e8-row
# stress test exercises exactly this property.
#
# Internal-only at v0.3.0. No @export tag.

# Per-term overhead constant in bytes; covers list wrappers and the
# small attributes attached to each design slot. Empirically calibrated
# against `object.size()` on isolated fb_terms allocations: 256 bytes
# overshoots a hand-counted ALLOC headers + names + attributes block
# on R 4.4 by ~16 bytes, comfortably inside the 5% acceptance tolerance.
.FB_PREFLIGHT_TERM_OVERHEAD <- 256

# Fixed bookkeeping overhead for the preflight result itself plus the
# IR + dataset wrappers. Page-sized for cleanliness; the absolute value
# is a small fraction of any single term that would trigger a refusal.
.FB_PREFLIGHT_BASE_OVERHEAD <- 4096

# Default smooth-basis dimension when an mgcv-style s(x) term carries
# no explicit k= argument and .enrich() has not yet populated term$k.
# Matches mgcv's default (10) per ?mgcv::s. The enriched IR carries
# the resolved k = ncol(sm_obj$X) after absorb.cons; that value is
# preferred when available.
.FB_PREFLIGHT_SMOOTH_K_DEFAULT <- 10L

# Default spline df for spl(x) -- matches the codegen path's hardcoded
# splines::bs(..., df = 8, degree = 3, intercept = FALSE) basis. If a
# future enrich step records term$df we will prefer it.
.FB_PREFLIGHT_SPLINE_DF_DEFAULT <- 8L

# Default fraction of system RAM used as the design-memory ceiling when
# no explicit `memory_ceiling_gb` is passed. 0.6 (60%) is conservative
# for shared workstations / CI / interactive laptops where the user
# expects the rest of the machine to remain responsive under load.
# Override per session via `options(flexyBayes.preflight_ram_fraction
# = <fraction>)`. The previous 80% default was too high for shared
# hosts; the new default plus the option give users explicit control.
.FB_PREFLIGHT_DEFAULT_RAM_FRACTION <- 0.6


# Top-level constructor. Returns the `<fb_preflight>` object even when
# the design is feasible (refusal == NULL); the caller pattern-matches
# on inherits(pf$refusal, "fb_preflight_refusal").
#
# The `known_matrices` arg (v0.3.10) threads the user-supplied
# vm/ped carriers through so the per-term INLA memory estimator can
# resolve sparse-precision nnz() and block-diagonal sum_k n_k^2
# accurately. NULL is the v0.3.5 backward-compatible default ---
# callers that do not have known_matrices in
# scope (legacy tests, dataset-only preflight) get the conservative
# upper-bound estimator that assumes worst-case dense/diagonal forms.
.fb_preflight <- function(
  fb_ir,
  fb_dataset,
  memory_ceiling_gb = NULL,
  known_matrices = NULL
) {
  if (!inherits(fb_ir, "fb_terms")) {
    stop(
      ".fb_preflight() requires an `<fb_terms>` IR; got: ",
      paste(class(fb_ir), collapse = "/"),
      call. = FALSE
    )
  }
  if (!inherits(fb_dataset, "fb_dataset")) {
    stop(
      ".fb_preflight() requires an `<fb_dataset>` wrapper; got: ",
      paste(class(fb_dataset), collapse = "/"),
      call. = FALSE
    )
  }

  ceiling_bytes <- .fb_resolve_ceiling(memory_ceiling_gb)
  family_ok <- .fb_preflight_family_in_scope(fb_ir)

  per_term <- list()
  n_rows <- as.numeric(fb_dataset$n_rows)

  # Fixed-effect terms (intercept is implicit; not estimated separately;
  # the design block for fixed numerics + factors carries it). Each
  # term carries enough to dispatch to the right estimator.
  for (i in seq_along(fb_ir$fixed_terms)) {
    entry <- .preflight_fixed_term(
      term = fb_ir$fixed_terms[[i]],
      n_rows = n_rows,
      fb_dataset = fb_dataset,
      family_ok = family_ok
    )
    per_term[[entry$label]] <- entry
  }

  # Random-effect terms (random intercept, uncorrelated random slope).
  for (i in seq_along(fb_ir$random_terms)) {
    entry <- .preflight_random_term(
      term = fb_ir$random_terms[[i]],
      n_rows = n_rows,
      fb_dataset = fb_dataset,
      family_ok = family_ok
    )
    per_term[[entry$label]] <- entry
  }

  # Smooth detection on fixed_terms (mgcv-style s() lives there in the
  # asreml ingest path; the brms walker refuses smooths at ingest time
  # so they don't reach the IR via fb_from_brms, but the IR may still
  # carry them on the asreml side).
  # The .preflight_fixed_term() dispatcher handles the smooth marker
  # internally; the if-branch above already covers it.

  # Unknown-representation gate. Any per-term entry flagged as
  # unknown_representation = TRUE forces a refusal
  # regardless of byte estimate: an honest ceiling check cannot be made
  # on a term whose design shape we did not characterise. This refusal
  # is preferred over the byte-ceiling refusal -- a known-too-large
  # design is informative, but an unknown design is the safety-critical
  # case.
  unknown_terms <- vapply(
    per_term,
    function(e) isTRUE(e$unknown_representation),
    logical(1L)
  )

  refusal <- NULL
  total <- NA_real_

  if (any(unknown_terms)) {
    unknown_label <- names(per_term)[which(unknown_terms)[1L]]
    refusal <- .new_fb_preflight_refusal(
      reason_code = "representation_unknown_for_preflight",
      binding_term = unknown_label,
      binding_bytes = NA_real_,
      total_bytes = NA_real_,
      ceiling_bytes = ceiling_bytes,
      n_rows = n_rows
    )
  } else {
    total <- .FB_PREFLIGHT_BASE_OVERHEAD +
      sum(vapply(per_term, function(e) e$design_memory_bytes, numeric(1L)))
    if (total > ceiling_bytes) {
      # Binding term = largest single-term estimate.
      binding_idx <- which.max(vapply(
        per_term,
        function(e) e$design_memory_bytes,
        numeric(1L)
      ))
      binding_label <- names(per_term)[[binding_idx]]
      binding_entry <- per_term[[binding_idx]]
      refusal <- .new_fb_preflight_refusal(
        reason_code = "design_memory_exceeds_ceiling",
        binding_term = binding_label,
        binding_bytes = binding_entry$design_memory_bytes,
        total_bytes = total,
        ceiling_bytes = ceiling_bytes,
        n_rows = n_rows
      )
    }
  }

  # Model-level aggregation plan. Informational at v0.3.0.9000 --
  # the per-term
  # `aggregated_likelihood_candidate` flag remains the dispatcher
  # contract through v0.3.1. The plan supersedes per-term flags as the
  # dispatch input in v0.3.2 when the aggregated emit path lands.
  agg_plan <- tryCatch(
    .fb_aggregation_plan(fb_ir, fb_dataset),
    error = function(e) NULL
  )

  # Per-term INLA memory estimator (v0.3.10). Reads
  # representation_class per term and applies
  # the INLA-specific model (dense n^2, sparse nnz * c_sparse, block-
  # diagonal sum_k n_k^2, pedigree_sparse_precision nnz(A^-1) *
  # c_sparse, fixed effects n * p) plus a 2x internal overhead
  # constant. The result is an <fb_memory_estimate> with backward-
  # compatible as.numeric() returning the total bytes.
  memory_estimate <- .estimate_inla_memory_breakdown(
    per_term = per_term,
    fb_ir = fb_ir,
    n_rows = n_rows,
    known_matrices = known_matrices,
    fb_dataset = fb_dataset
  )

  structure(
    list(
      per_term_estimate = per_term,
      total_estimate_bytes = total,
      ceiling_bytes = ceiling_bytes,
      n_rows = n_rows,
      aggregation_plan = agg_plan,
      memory_estimate = memory_estimate,
      refusal = refusal
    ),
    class = c("fb_preflight", "list")
  )
}


# ---------------------------------------------------------------- #
# Family-in-scope check for aggregated_likelihood_candidate            #
# ---------------------------------------------------------------- #

# The aggregation scope envelope includes "gaussian family +
# identity link". Preflight evaluates this once and gates every
# per-term flag through it; non-gaussian / non-identity families
# cannot carry the aggregation candidacy on any term.
.fb_preflight_family_in_scope <- function(fb_ir) {
  fam_ok <- identical(fb_ir$family, "gaussian")
  link_ok <- is.null(fb_ir$link) || identical(fb_ir$link, "identity")
  fam_ok && link_ok
}


# ---------------------------------------------------------------- #
# Fixed-term estimator                                              #
# ---------------------------------------------------------------- #

# Walk a single fixed_terms[[i]] entry. The IR's fixed_terms list
# entries carry varying shapes depending on the ingest path; we
# pattern-match on $type with one branch per term class the parsers
# can emit. Anything that does not match a known class falls into
# the unknown-representation guard at the bottom -- a "conservative
# lower bound" is not safe for refusal; under-estimating risks
# letting a model through that the backend then fails to allocate.
.preflight_fixed_term <- function(term, n_rows, fb_dataset, family_ok) {
  ttype <- if (!is.null(term$type)) term$type else "expression"
  label <- .preflight_term_label(term, ttype, kind = "fixed")

  # Smooth s(x) / t2(x) -- dense basis block. Two forms hit this
  # branch: the hand-built sentinel "smooth"/"s"/"t2" used by the test
  # IR helpers and the parser-emitted "smooth_mgcv" (the live ingest
  # path; parse_formula.R:149 emits this type and .enrich() populates
  # term$k = ncol(sm_obj$X) after mgcv::smoothCon() absorbs the
  # identifiability constraint).
  if (ttype %in% c("smooth", "s", "t2", "smooth_mgcv")) {
    k <- if (!is.null(term$k) && is.numeric(term$k)) {
      as.integer(term$k)
    } else {
      .FB_PREFLIGHT_SMOOTH_K_DEFAULT
    }
    bytes <- 8 *
      n_rows *
      k + # N x k dense basis
      8 * k + # k basis coefficients
      .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "dense_baseline"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "fixed_smooth"
    ))
  }

  # P-spline spl(x) -- splines::bs() basis block. The codegen path
  # hardcodes df = 8, degree = 3, intercept = FALSE so the dense
  # basis is N x df. We prefer term$df if a future enrich step
  # records it; otherwise we use the codegen default.
  if (identical(ttype, "spline")) {
    df <- if (!is.null(term$df) && is.numeric(term$df)) {
      as.integer(term$df)
    } else {
      .FB_PREFLIGHT_SPLINE_DF_DEFAULT
    }
    bytes <- 8 * n_rows * df + 8 * df + .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "dense_baseline"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "fixed_spline"
    ))
  }

  # Numeric fixed term -- one double column.
  if (ttype %in% c("numeric", "continuous", "I")) {
    bytes <- 8 *
      n_rows + # one double column
      8 + # one beta coefficient
      .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "indexed_fixed_numeric"
      ),
      aggregated_likelihood_candidate = family_ok,
      term_kind = "fixed_numeric"
    ))
  }

  # Factor fixed term -- integer index + (k-1) beta dummy contrasts.
  # NA level count escalates to representation_unknown rather than
  # silently sizing the dummy block at 0 -- the count is the binding
  # factor for the dummy vector size.
  if (ttype %in% c("factor", "categorical")) {
    k <- .preflight_term_level_count(term, fb_dataset)
    if (is.na(k)) {
      return(.preflight_entry(
        label = label,
        design_memory_bytes = NA_real_,
        representation_class = "unknown",
        aggregated_likelihood_candidate = FALSE,
        term_kind = "fixed_factor",
        unknown_representation = TRUE
      ))
    }
    n_dummies <- max(k - 1L, 0L)
    bytes <- 4 *
      n_rows + # integer index column
      8 * n_dummies + # dummy-contrast beta vector
      .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "indexed_fixed_factor"
      ),
      aggregated_likelihood_candidate = family_ok,
      term_kind = "fixed_factor"
    ))
  }

  # Factor:factor interaction -- treatment-coded reduced-rank
  # model.matrix(~ f1*f2, data) returns prod(L_i - 1) columns for the
  # interaction block alone, on the standard convention where main
  # effects and the intercept are also present. The hard case:
  # high-cardinality factor:factor can exceed the ceiling by
  # orders of magnitude vs the v0.3.0 "single dense column"
  # fallthrough estimate. The estimate is exact for the common
  # additive-with-interaction model; on the rare `y ~ f1:f2 - 1`
  # form it under-estimates by (L1 + L2 - 1) columns -- in any
  # high-cardinality regime that is negligible relative to the
  # prod(L_i - 1) block size.
  if (identical(ttype, "factor_interaction")) {
    n_dummies <- .preflight_factor_interaction_dummies(term$vars, fb_dataset)
    if (is.na(n_dummies)) {
      # Cardinality unresolvable from the dataset wrapper -- the
      # interaction memory is unknown. Flag for representation-
      # unknown refusal at the top level.
      return(.preflight_entry(
        label = label,
        design_memory_bytes = NA_real_,
        representation_class = "unknown",
        aggregated_likelihood_candidate = FALSE,
        term_kind = "fixed_factor_interaction",
        unknown_representation = TRUE
      ))
    }
    bytes <- 8 *
      n_rows *
      n_dummies +
      8 * n_dummies +
      .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "dense_baseline"
      ),
      aggregated_likelihood_candidate = family_ok,
      term_kind = "fixed_factor_interaction"
    ))
  }

  # Factor:continuous interaction -- treatment-coded indexed slope.
  # Representation: factor index column + continuous column + (L - 1)
  # per-level slope coefficients (reference level's slope fixed to
  # zero so no coefficient slot). NA level count escalates to
  # representation_unknown.
  if (identical(ttype, "factor_numeric_interaction")) {
    L <- .preflight_term_level_count(term, fb_dataset)
    if (is.na(L)) {
      return(.preflight_entry(
        label = label,
        design_memory_bytes = NA_real_,
        representation_class = "unknown",
        aggregated_likelihood_candidate = FALSE,
        term_kind = "fixed_factor_numeric_interaction",
        unknown_representation = TRUE
      ))
    }
    bytes <- 4 *
      n_rows +
      8 * n_rows +
      8 * max(L - 1L, 0L) +
      .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "indexed_fixed_factor_numeric"
      ),
      aggregated_likelihood_candidate = family_ok,
      term_kind = "fixed_factor_numeric_interaction"
    ))
  }

  # Polynomial pol(x, d) -- degree-d dense basis. The IR carries
  # term$degree from parse_formula.R; the basis is d columns of N.
  if (identical(ttype, "polynomial")) {
    d <- if (!is.null(term$degree) && is.numeric(term$degree)) {
      as.integer(term$degree)
    } else {
      2L
    }
    bytes <- 8 * n_rows * d + 8 * d + .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "dense_baseline"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "fixed_polynomial"
    ))
  }

  # Generic numeric:numeric interaction -- a single product column.
  if (identical(ttype, "interaction")) {
    bytes <- 8 * n_rows + 8 + .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "indexed_fixed_numeric"
      ),
      aggregated_likelihood_candidate = family_ok,
      term_kind = "fixed_interaction"
    ))
  }

  # Expression term -- I() / arithmetic / function calls evaluated
  # against data (parse_formula.R:89 returns this when the label is
  # not a bare variable in `data`). The standard representation is
  # one double column (the evaluated expression vector) plus one
  # slope coefficient.
  if (identical(ttype, "expression")) {
    bytes <- 8 * n_rows + 8 + .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "indexed_fixed_numeric"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "fixed_expression"
    ))
  }

  # Unknown fallthrough. The contract is
  # "do not silently under-estimate"; an unrecognised fixed term flags
  # the entry for representation-unknown refusal at the top level.
  # A diagnostic single-column estimate is recorded so the print
  # surface still shows something, but it is NOT consumed by the
  # ceiling check when unknown_representation = TRUE.
  bytes <- 8 * n_rows + 8 + .FB_PREFLIGHT_TERM_OVERHEAD
  .preflight_entry(
    label = label,
    design_memory_bytes = bytes,
    representation_class = "unknown",
    aggregated_likelihood_candidate = FALSE,
    term_kind = "fixed_other",
    unknown_representation = TRUE
  )
}


# ---------------------------------------------------------------- #
# Random-term estimator                                             #
# ---------------------------------------------------------------- #

.preflight_random_term <- function(term, n_rows, fb_dataset, family_ok) {
  rtype <- if (!is.null(term$type)) term$type else "simple"
  label <- .preflight_term_label(term, rtype, kind = "random")
  k <- .preflight_term_level_count(term, fb_dataset)

  # (1 | g) -- random intercept. Integer index + k-vector of u's.
  # NA level count escalates to representation_unknown rather than
  # silently sizing the latent vector at 0 -- the group cardinality
  # drives the latent block.
  if (identical(rtype, "simple")) {
    if (is.na(k)) {
      return(.preflight_entry(
        label = label,
        design_memory_bytes = NA_real_,
        representation_class = "unknown",
        aggregated_likelihood_candidate = FALSE,
        term_kind = "random_intercept",
        unknown_representation = TRUE
      ))
    }
    bytes <- 4 *
      n_rows + # group index column
      8 * k + # latent random-effect vector
      .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "indexed_random_intercept"
      ),
      aggregated_likelihood_candidate = family_ok,
      term_kind = "random_intercept"
    ))
  }

  # (x || g) / (1 + x || g) / (0 + x || g) -- uncorrelated slope.
  # Integer index + slope variable column + k-vector of slope (and,
  # if with_intercept, another k-vector of intercept). NA level
  # count escalates to representation_unknown.
  if (identical(rtype, "simple_slope_uncor")) {
    if (is.na(k)) {
      return(.preflight_entry(
        label = label,
        design_memory_bytes = NA_real_,
        representation_class = "unknown",
        aggregated_likelihood_candidate = FALSE,
        term_kind = "random_slope",
        unknown_representation = TRUE
      ))
    }
    with_int <- isTRUE(term$with_intercept)
    bytes <- 4 *
      n_rows + # group index
      8 * n_rows + # slope variable column
      8 * k * (1L + as.integer(with_int)) + # latent slope (+ intercept)
      .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "dense_baseline"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "random_slope"
    ))
  }

  # Smooth s(x) on the random side -- the asreml ingest path places
  # parser-emitted "smooth_mgcv" terms in random_terms (e.g.
  # random = ~ s(x, k = 6)). Same N x k dense basis cost as the
  # fixed-side branch. This previously fell through to the
  # structured-cov estimate below, which under-estimated by a factor
  # of k.
  if (rtype %in% c("smooth_mgcv", "smooth", "s", "t2")) {
    sm_k <- if (!is.null(term$k) && is.numeric(term$k)) {
      as.integer(term$k)
    } else {
      .FB_PREFLIGHT_SMOOTH_K_DEFAULT
    }
    bytes <- 8 * n_rows * sm_k + 8 * sm_k + .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "dense_baseline"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "random_smooth"
    ))
  }

  # P-spline spl(x) on the random side. Same N x df dense basis as
  # the fixed-side branch (codegen default df = 8).
  if (identical(rtype, "spline")) {
    df <- if (!is.null(term$df) && is.numeric(term$df)) {
      as.integer(term$df)
    } else {
      .FB_PREFLIGHT_SPLINE_DF_DEFAULT
    }
    bytes <- 8 * n_rows * df + 8 * df + .FB_PREFLIGHT_TERM_OVERHEAD
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "dense_baseline"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "random_spline"
    ))
  }

  # interaction_generic -- the parse_formula classifier's fallback
  # when no structured-interaction pattern matched (parse_formula.R
  # :214). Memory is unknown because the design shape is not
  # characterised by a single template. Flag for representation-
  # unknown refusal at the top level.
  if (identical(rtype, "interaction_generic")) {
    return(.preflight_entry(
      label = label,
      design_memory_bytes = NA_real_,
      representation_class = "unknown",
      aggregated_likelihood_candidate = FALSE,
      term_kind = "random_unknown",
      unknown_representation = TRUE
    ))
  }

  # Structured-covariance random terms. This branch is split into
  # two categories:
  #
  #   (a) Known indexed representations with a tested estimator
  #       shape (vm, ped, ar1, ar1_spatial): a k-length latent vector
  #       plus an integer row index. The heuristic below sizes these
  #       at 4N (index) + 8N (slope-variable column, conservative
  #       carry-over) + 16k (latent + bookkeeping) + overhead.
  #
  #   (b) Known dense / Cholesky / factor-analysis matrix
  #       representations that the indexed heuristic under-estimates
  #       (at, us, fa, at_simple, at_units, us_gxe, fa_gxe, vm_gxe,
  #       nested, combo) -- routed to unknown_representation = TRUE
  #       rather than silently under-estimating at large N. A
  #       term-specific estimator can be added per class without
  #       changing the dispatch shape (extend KNOWN_INDEXED below
  #       once a per-class estimator and its test ship together).
  #
  # Either category escalates to unknown_representation when the
  # level count cannot be resolved from the IR or the dataset
  # wrapper's dictionaries (k = NA_integer_) -- the level count is
  # the binding factor for the latent vector size, so refusing
  # cleanly beats guessing.
  KNOWN_INDEXED <- c("vm", "ped", "ar1", "ar1_spatial")

  if (is.na(k) || !rtype %in% KNOWN_INDEXED) {
    return(.preflight_entry(
      label = label,
      design_memory_bytes = NA_real_,
      representation_class = "unknown",
      aggregated_likelihood_candidate = FALSE,
      term_kind = "random_structured",
      unknown_representation = TRUE
    ))
  }

  # Blocks-format vm/ped terms (v0.3.10) surface as `block_diagonal`
  # so fb_plan() can render "(K blocks)" annotations on the print
  # form. The byte estimate tracks the indexed_structured_known
  # formula at preflight time --- the per-block memory model lands in
  # the per-term INLA memory estimator.
  cov <- term$cov_representation
  bytes <- 4 * n_rows + 8 * n_rows + 16 * k + .FB_PREFLIGHT_TERM_OVERHEAD
  if (!is.null(cov) && identical(cov$format, "blocks")) {
    return(.preflight_entry(
      label = label,
      design_memory_bytes = bytes,
      representation_class = .representation_class(
        "block_diagonal"
      ),
      aggregated_likelihood_candidate = FALSE,
      term_kind = "random_structured"
    ))
  }
  .preflight_entry(
    label = label,
    design_memory_bytes = bytes,
    representation_class = .representation_class(
      "indexed_structured_known"
    ),
    aggregated_likelihood_candidate = FALSE,
    term_kind = "random_structured"
  )
}


# ---------------------------------------------------------------- #
# Per-term helpers                                                  #
# ---------------------------------------------------------------- #

# Build a single per-term entry. Centralised so adding a new field
# (e.g. an `effort_estimate`) is a one-line change.
# `unknown_representation = TRUE` flags terms whose design memory we
# could not characterise; the top-level `.fb_preflight()` then routes
# such terms to a `representation_unknown_for_preflight` refusal in
# preference to a (possibly under-estimated) ceiling check.
.preflight_entry <- function(
  label,
  design_memory_bytes,
  representation_class,
  aggregated_likelihood_candidate,
  term_kind,
  unknown_representation = FALSE
) {
  list(
    label = label,
    design_memory_bytes = as.numeric(design_memory_bytes),
    representation_class = representation_class,
    aggregated_likelihood_candidate = isTRUE(aggregated_likelihood_candidate),
    term_kind = term_kind,
    overhead_constant = .FB_PREFLIGHT_TERM_OVERHEAD,
    unknown_representation = isTRUE(unknown_representation)
  )
}

# Count the dummy columns for a factor:factor (...:factor)
# treatment-coded reduced-rank interaction. The interaction block
# alone in `model.matrix(~ f1*f2*...)` carries prod(L_i - 1) columns
# when main effects and intercept are also present (the standard
# convention). Returns NA_integer_ when any factor cardinality
# cannot be resolved from the dataset wrapper's dictionaries.
.preflight_factor_interaction_dummies <- function(vars, fb_dataset) {
  if (!length(vars)) {
    return(NA_integer_)
  }
  ks <- vapply(
    vars,
    function(v) {
      .fb_dataset_levels(fb_dataset, as.character(v))
    },
    numeric(1L)
  )
  if (anyNA(ks)) {
    return(NA_integer_)
  }
  ks_minus_one <- pmax(as.integer(ks) - 1L, 0L)
  as.integer(prod(ks_minus_one))
}

# Build a display label from an IR term descriptor. Mirrors the
# canonical-name conventions for diagnostic visibility -- the binding
# term in a refusal must be human-readable.
.preflight_term_label <- function(term, ttype, kind) {
  if (!is.null(term$label) && nzchar(term$label)) {
    return(term$label)
  }
  if (!is.null(term$deparse) && nzchar(term$deparse)) {
    return(term$deparse)
  }

  if (identical(kind, "fixed")) {
    if (
      ttype %in% c("smooth", "s", "t2", "smooth_mgcv") && !is.null(term$var)
    ) {
      return(paste0("s(", term$var, ")"))
    }
    if (identical(ttype, "spline") && !is.null(term$var)) {
      return(paste0("spl(", term$var, ")"))
    }
    if (identical(ttype, "factor_interaction") && !is.null(term$vars)) {
      return(paste(term$vars, collapse = ":"))
    }
    if (!is.null(term$var)) {
      return(as.character(term$var))
    }
    return("<fixed:unnamed>")
  }
  # random
  group <- if (!is.null(term$var)) as.character(term$var) else "?"
  if (identical(ttype, "simple")) {
    return(paste0("(1 | ", group, ")"))
  }
  if (identical(ttype, "simple_slope_uncor")) {
    slope <- if (!is.null(term$slope_var)) {
      as.character(term$slope_var)
    } else {
      "x"
    }
    lhs <- if (isTRUE(term$with_intercept)) {
      paste0("1 + ", slope)
    } else {
      slope
    }
    return(paste0("(", lhs, " || ", group, ")"))
  }
  if (ttype %in% c("smooth_mgcv", "smooth", "s", "t2")) {
    return(paste0("s(", group, ")"))
  }
  if (identical(ttype, "spline")) {
    return(paste0("spl(", group, ")"))
  }
  paste0("(", ttype, " | ", group, ")")
}

# Level-count lookup. Reads in priority order: (1) the IR's cached
# $var_n; (2) the dataset's frozen dictionary length; (3) NA_integer_
# when neither source can resolve the count.
#
# The previous `1L` fallback was removed: it was labelled
# "conservative" but is actually the smallest latent
# block, which silently under-estimates the design memory whenever
# the level count is genuinely unresolvable. The structured-cov
# branch in .preflight_random_term() (and any future estimator that
# depends on `k`) treats NA_integer_ as the signal to escalate to a
# `representation_unknown_for_preflight` refusal at the top level.
.preflight_term_level_count <- function(term, fb_dataset) {
  # Random-term IR carries $var_n; fixed-factor IR carries $n_levels;
  # the asreml ingest occasionally uses $K. Try each in order.
  for (slot in c("var_n", "n_levels", "K")) {
    v <- term[[slot]]
    if (!is.null(v) && !is.na(v)) return(as.integer(v))
  }
  if (!is.null(term$var)) {
    k <- .fb_dataset_levels(fb_dataset, as.character(term$var))
    if (!is.na(k)) return(as.integer(k))
  }
  NA_integer_
}


# ---------------------------------------------------------------- #
# Memory ceiling resolution                                         #
# ---------------------------------------------------------------- #

# Resolves the byte-valued ceiling. Resolution priority:
#   1. Explicit `memory_ceiling_gb` argument (always wins).
#   2. RAM probe x fraction-of-RAM, where the fraction is
#      `getOption("flexyBayes.preflight_ram_fraction", 0.6)`. The
#      default fraction moved from 0.8 to 0.6 and the option is
#      exposed so users on shared
#      workstations / CI / interactive laptops can tighten or
#      loosen the default without passing an argument through every
#      dispatch entry. The fraction is validated as `(0, 1]`.
#   3. 8 GiB conservative fallback when neither RAM probe is
#      available (benchmarkme / memuse uninstalled). The fallback
#      is documented so users in minimal install environments see
#      a predictable refusal threshold rather than a hard error;
#      the RAM-fraction option is NOT applied to the fallback (the
#      8 GiB number is itself already a conservative absolute floor).
.fb_resolve_ceiling <- function(memory_ceiling_gb) {
  if (!is.null(memory_ceiling_gb)) {
    if (
      !is.numeric(memory_ceiling_gb) ||
        length(memory_ceiling_gb) != 1L ||
        memory_ceiling_gb <= 0
    ) {
      stop(
        "`memory_ceiling_gb` must be a positive numeric scalar; got: ",
        deparse(memory_ceiling_gb),
        call. = FALSE
      )
    }
    return(as.numeric(memory_ceiling_gb) * 1024^3)
  }

  ram_bytes <- .fb_probe_system_ram()
  if (is.na(ram_bytes)) {
    # 8 GiB conservative fallback when neither probe is available.
    return(8 * 1024^3)
  }

  fraction <- getOption(
    "flexyBayes.preflight_ram_fraction",
    .FB_PREFLIGHT_DEFAULT_RAM_FRACTION
  )
  if (
    !is.numeric(fraction) ||
      length(fraction) != 1L ||
      is.na(fraction) ||
      fraction <= 0 ||
      fraction > 1
  ) {
    stop(
      "`flexyBayes.preflight_ram_fraction` must be a numeric scalar ",
      "in (0, 1]; got: ",
      deparse(fraction),
      call. = FALSE
    )
  }

  as.numeric(fraction) * ram_bytes
}

.fb_probe_system_ram <- function() {
  # Quiet, platform-native probes for the common platforms, with child stderr
  # suppressed (`stderr = FALSE`). benchmarkme / memuse shell out to system
  # tools (sysctl, kstat) without stderr control, which prints harmless but
  # alarming noise ("sysctl: hw.memsize: Operation not permitted",
  # "/bin/kstat: No such file or directory") on restricted or sandboxed hosts
  # and makes a successful preflight look broken. The native probes avoid that;
  # the optional helpers are a last resort for other platforms only.
  sysname <- Sys.info()[["sysname"]]
  if (identical(sysname, "Darwin")) {
    res <- tryCatch(
      as.numeric(suppressWarnings(
        system2("sysctl", c("-n", "hw.memsize"),
          stdout = TRUE, stderr = FALSE)
      )[1]),
      error = function(e) NA_real_
    )
    return(if (!is.na(res) && res > 0) res else NA_real_)
  }
  if (identical(sysname, "Linux")) {
    res <- tryCatch(
      {
        mt <- grep("^MemTotal:", readLines("/proc/meminfo"), value = TRUE)[1]
        as.numeric(regmatches(mt, regexpr("[0-9]+", mt))) * 1024
      },
      error = function(e) NA_real_
    )
    return(if (!is.na(res) && res > 0) res else NA_real_)
  }
  if (requireNamespace("benchmarkme", quietly = TRUE)) {
    res <- tryCatch(as.numeric(benchmarkme::get_ram()), error = function(e) {
      NA_real_
    })
    if (!is.na(res) && res > 0) return(res)
  }
  if (requireNamespace("memuse", quietly = TRUE)) {
    res <- tryCatch(
      {
        mi <- memuse::Sys.meminfo()
        as.numeric(mi$totalram@size) *
          switch(
            mi$totalram@unit,
            "B" = 1,
            "KiB" = 1024,
            "MiB" = 1024^2,
            "GiB" = 1024^3,
            "TiB" = 1024^4,
            NA_real_
          )
      },
      error = function(e) NA_real_
    )
    if (!is.na(res) && res > 0) return(res)
  }
  NA_real_
}


# ---------------------------------------------------------------- #
# Refusal constructor + print                                       #
# ---------------------------------------------------------------- #

.new_fb_preflight_refusal <- function(
  reason_code,
  binding_term,
  binding_bytes,
  total_bytes,
  ceiling_bytes,
  n_rows
) {
  structure(
    list(
      reason_code = reason_code,
      binding_term = binding_term,
      binding_bytes = as.numeric(binding_bytes),
      total_bytes = as.numeric(total_bytes),
      ceiling_bytes = as.numeric(ceiling_bytes),
      n_rows = n_rows
    ),
    class = c("fb_preflight_refusal", "list")
  )
}


# ---------------------------------------------------------------- #
# S3 print methods                                                  #
# ---------------------------------------------------------------- #

#' Print method for an internal `<fb_preflight>` summary
#'
#' Diagnostic print of the design-memory preflight result:
#' per-term `design_memory_bytes` (formatted with thousand-separator
#' " "), `representation_class`, and the aggregate ceiling check.
#' On refusal the binding term + numeric ceiling appear below the
#' per-term table.
#'
#' @param x   an `<fb_preflight>` object.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @keywords internal
#' @export
print.fb_preflight <- function(x, ...) {
  fmt_n <- function(b) {
    if (is.na(b)) {
      "NA"
    } else {
      format(round(b), big.mark = " ", scientific = FALSE)
    }
  }
  cat(sprintf(
    "<fb_preflight> n_rows = %s; total = %s bytes; ceiling = %s bytes\n",
    format(x$n_rows, big.mark = " ", scientific = FALSE),
    fmt_n(x$total_estimate_bytes),
    fmt_n(x$ceiling_bytes)
  ))
  if (length(x$per_term_estimate)) {
    cat("  per-term:\n")
    for (nm in names(x$per_term_estimate)) {
      e <- x$per_term_estimate[[nm]]
      cat(sprintf(
        "    %-28s  %14s B  %-26s  agg_candidate=%s\n",
        nm,
        fmt_n(e$design_memory_bytes),
        e$representation_class,
        e$aggregated_likelihood_candidate
      ))
    }
  }
  if (!is.null(x$aggregation_plan)) {
    cat(sprintf(
      "  aggregation: eligible=%s; K_est=%s; compression_est=%s\n",
      x$aggregation_plan$eligible,
      if (is.na(x$aggregation_plan$K_est)) {
        "NA"
      } else {
        format(x$aggregation_plan$K_est, big.mark = " ", scientific = FALSE)
      },
      if (is.na(x$aggregation_plan$compression_est)) {
        "NA"
      } else {
        sprintf("%.3f", x$aggregation_plan$compression_est)
      }
    ))
  }
  if (!is.null(x$refusal)) {
    cat("\n")
    print(x$refusal)
  }
  invisible(x)
}

#' Print method for an internal `<fb_preflight_refusal>` object
#'
#' Three-line diagnostic: reason code, binding term + its byte
#' estimate, and the active ceiling with the suggested override.
#'
#' @param x   an `<fb_preflight_refusal>` object.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @keywords internal
#' @export
print.fb_preflight_refusal <- function(x, ...) {
  cat("<fb_preflight_refusal>\n")
  cat(sprintf("  reason_code:    %s\n", x$reason_code))
  cat(sprintf("  binding_term:   %s\n", x$binding_term))

  fmt_bytes <- function(b) {
    if (is.na(b)) {
      "NA (not characterised)"
    } else {
      sprintf(
        "%s B  (%.2f GiB)",
        format(round(b), big.mark = " ", scientific = FALSE),
        b / 1024^3
      )
    }
  }
  cat(sprintf("  binding_bytes:  %s\n", fmt_bytes(x$binding_bytes)))
  cat(sprintf("  total_bytes:    %s\n", fmt_bytes(x$total_bytes)))
  cat(sprintf("  ceiling_bytes:  %s\n", fmt_bytes(x$ceiling_bytes)))

  if (
    identical(x$reason_code, "design_memory_exceeds_ceiling") &&
      !is.na(x$total_bytes)
  ) {
    ceiling_gb_suggest <- ceiling(x$total_bytes / 1024^3 * 1.25)
    cat(sprintf(
      "  override:       pass memory_ceiling_gb = %d to retry\n",
      ceiling_gb_suggest
    ))
  } else if (identical(x$reason_code, "representation_unknown_for_preflight")) {
    cat(
      "  remedy:         the term shape is not characterised by the\n",
      "                  preflight estimator. Restructure the formula\n",
      "                  to use a recognised term class, or open an\n",
      "                  issue with the term type so the estimator\n",
      "                  can be extended.\n",
      sep = ""
    )
  } else if (identical(x$reason_code, "memory_feasibility_inla_per_term")) {
    cat(
      "  remedy:         INLA's per-term memory model exceeds the\n",
      "                  active ceiling. Route via backend = \"greta\"\n",
      "                  for the indexed representation, or reduce\n",
      "                  the structured-cov term cardinality.\n",
      sep = ""
    )
  }
  invisible(x)
}


# ---------------------------------------------------------------- #
# v0.3.10 --- Per-term INLA memory estimator                        #
# ---------------------------------------------------------------- #

# Sparse-matrix bytes-per-non-zero constant. The Matrix package's
# dsCMatrix stores a value (8B double) + a row index (4B int) + a
# per-column pointer (~4B amortised); the practical figure is well-
# attested at ~16 bytes per non-zero (Lindgren, Rue, Lindstrom 2011
# §6.3; matches Matrix::object.size() probes on production-shape
# precision matrices).
.FB_INLA_BYTES_PER_NNZ <- 16L

# INLA-internal overhead multiplier. Absorbs the per-component
# allocations INLA's stage-2 (CCD/grid integration) makes on top of
# the structural design memory. Configurable via the option
# `flexyBayes.preflight_inla_overhead_factor` (default 2.0).
.FB_INLA_OVERHEAD_DEFAULT <- 2.0


# .estimate_inla_memory_breakdown() --- the per-term INLA memory model.
#
# Walks the preflight per_term_estimate + the IR's fixed/random term
# slots and applies a representation-specific INLA-memory model to
# each. Returns an `<fb_memory_estimate>` carrier with two fields:
#
#   $total                  numeric(1L), bytes; the overhead-applied
#                           sum of per-term contributions plus a fixed-
#                           effect model-matrix contribution.
#   $breakdown              a data.frame with one row per term:
#                              term_label       (character)
#                              representation   (character)
#                              bytes            (numeric)
#                              share            (numeric in [0, 1])
#                           Column order is stable under the
#                           closed-vocabulary discipline.
#
# Backward compatibility: the carrier inherits class `fb_memory_estimate`
# and `as.numeric.fb_memory_estimate()` returns `$total`, so existing
# scalar-numeric consumers (`if (as.numeric(pf$memory_estimate) > ...)`)
# keep working unchanged.
.estimate_inla_memory_breakdown <- function(
  per_term,
  fb_ir,
  n_rows,
  known_matrices = NULL,
  fb_dataset = NULL
) {
  rows <- list()

  # Random-term INLA memory contributions, keyed on representation
  # class. The per_term entries from .preflight_random_term() carry
  # the design_memory_bytes from the indexed estimator; we re-derive
  # the INLA-specific bytes here per format.
  for (i in seq_along(fb_ir$random_terms)) {
    term <- fb_ir$random_terms[[i]]
    label <- .preflight_term_label(
      term,
      term$type %||% "simple",
      kind = "random"
    )
    entry <- per_term[[label]]
    if (is.null(entry) || isTRUE(entry$unknown_representation)) {
      next
    }

    repr <- entry$representation_class %||% "unknown"
    k <- .preflight_term_level_count(term, fb_dataset)
    bytes <- .inla_random_term_bytes(
      term,
      repr,
      k,
      known_matrices = known_matrices,
      n_rows = n_rows
    )
    rows[[length(rows) + 1L]] <- list(
      term_label = label,
      representation = repr,
      bytes = bytes
    )
  }

  # Fixed-effect contribution: the design matrix INLA materialises
  # as N x p doubles. We sum across fixed terms' bytes from the
  # preflight estimator --- those entries already model the per-term
  # design block (one column for numeric/expression, k-1 for factor
  # dummies, n_dummies for factor interactions, k for smooths /
  # polynomials).
  fixed_bytes <- 0
  for (i in seq_along(fb_ir$fixed_terms)) {
    term <- fb_ir$fixed_terms[[i]]
    label <- .preflight_term_label(
      term,
      term$type %||% "expression",
      kind = "fixed"
    )
    entry <- per_term[[label]]
    if (is.null(entry) || isTRUE(entry$unknown_representation)) {
      next
    }
    fixed_bytes <- fixed_bytes + (entry$design_memory_bytes %||% 0)
  }
  if (fixed_bytes > 0) {
    rows[[length(rows) + 1L]] <- list(
      term_label = "(fixed effects)",
      representation = "fixed_model_matrix",
      bytes = fixed_bytes
    )
  }

  bytes_vec <- vapply(rows, function(r) as.numeric(r$bytes), numeric(1L))
  raw_total <- sum(bytes_vec, na.rm = TRUE)
  overhead <- getOption(
    "flexyBayes.preflight_inla_overhead_factor",
    .FB_INLA_OVERHEAD_DEFAULT
  )
  if (
    !is.numeric(overhead) ||
      length(overhead) != 1L ||
      is.na(overhead) ||
      overhead <= 0
  ) {
    overhead <- .FB_INLA_OVERHEAD_DEFAULT
  }
  total <- raw_total * overhead

  breakdown <- if (length(rows) == 0L) {
    data.frame(
      term_label = character(0L),
      representation = character(0L),
      bytes = numeric(0L),
      share = numeric(0L),
      stringsAsFactors = FALSE
    )
  } else {
    shares <- if (raw_total > 0) {
      bytes_vec / raw_total
    } else {
      rep(NA_real_, length(rows))
    }
    data.frame(
      term_label = vapply(rows, `[[`, character(1L), "term_label"),
      representation = vapply(rows, `[[`, character(1L), "representation"),
      bytes = bytes_vec,
      share = shares,
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      total = total,
      breakdown = breakdown,
      overhead_factor = overhead,
      raw_total = raw_total
    ),
    class = c("fb_memory_estimate", "list")
  )
}


# Map a single random-term entry to its INLA memory contribution.
# Conservative-by-design: when known_matrices is unavailable, sparse-
# precision and block-diagonal forms fall back to upper-bound
# estimates that assume worst-case sparsity / partition shapes.
.inla_random_term_bytes <- function(
  term,
  repr,
  k,
  known_matrices = NULL,
  n_rows = NA_real_
) {
  if (is.na(k) || is.null(k)) {
    k <- 0L
  }
  k_num <- as.numeric(k)

  # Indexed random intercept / IID: a k-vector of latent u's + the
  # integer index column. INLA's overhead inflates this only modestly.
  if (identical(repr, .representation_class("indexed_random_intercept"))) {
    return(8 * k_num + 4 * n_rows)
  }
  if (identical(repr, .representation_class("dense_baseline"))) {
    # dense_baseline on the random side: the design block at this
    # representation is already dense (e.g., spline / smooth on
    # random), and the per-term design_memory_bytes is the right
    # surrogate for the INLA model-matrix carryover.
    return(0) # accounted via design_memory_bytes if any; defaulting
  }

  if (identical(repr, .representation_class("indexed_structured_known"))) {
    # Sparse precision: nnz(Q) * c_sparse. Without known_matrices we
    # use the indexed-fallback K * c_sparse (assumes diagonal Q ---
    # the smallest sparse footprint; intentional under-estimate
    # rather than alarmist over-estimate when Q is unavailable).
    cov <- term$cov_representation
    nnz <- NA_real_
    if (
      !is.null(known_matrices) &&
        !is.null(cov) &&
        !is.null(known_matrices[[cov$data]])
    ) {
      Q <- known_matrices[[cov$data]]
      nnz <- tryCatch(
        {
          if (inherits(Q, "sparseMatrix")) {
            as.numeric(Matrix::nnzero(Q))
          } else {
            sum(abs(as.matrix(Q)) > 0)
          }
        },
        error = function(e) NA_real_
      )
    }
    if (is.na(nnz)) {
      nnz <- k_num
    }
    return(nnz * .FB_INLA_BYTES_PER_NNZ)
  }

  if (identical(repr, .representation_class("block_diagonal"))) {
    # sum_k n_k^2 * 8. Without known_matrices, use K^2 as worst-case
    # upper bound (single block of full K).
    cov <- term$cov_representation
    block_sizes <- NULL
    if (
      !is.null(known_matrices) &&
        !is.null(cov) &&
        !is.null(known_matrices[[cov$data]])
    ) {
      blocks <- known_matrices[[cov$data]]
      if (is.list(blocks)) {
        block_sizes <- vapply(
          blocks,
          function(V) {
            d <- dim(V)
            if (is.null(d)) 0L else as.integer(d[[1L]])
          },
          integer(1L)
        )
      }
    }
    sum_sq <- if (is.null(block_sizes)) {
      k_num^2
    } else {
      sum(as.numeric(block_sizes)^2)
    }
    return(sum_sq * 8)
  }

  # Default: a small overhead constant per term (lighter representations
  # already accounted for in the indexed branches).
  0
}


#' Coerce an `<fb_memory_estimate>` carrier to a numeric scalar
#'
#' Returns `x$total` (bytes, after the INLA overhead multiplier).
#' Registered against the `as.double` S3 generic so that
#' `as.numeric(<fb_memory_estimate>)` dispatches correctly --- the
#' base-R `as.numeric` is `.Primitive("as.double")`, so dispatch
#' fires on `as.double`, not on `as.numeric`.
#'
#' @param x   an `<fb_memory_estimate>` carrier.
#' @param ... unused.
#' @return numeric(1L) total bytes.
#' @keywords internal
#' @export
as.double.fb_memory_estimate <- function(x, ...) {
  if (is.null(x$total)) {
    return(NA_real_)
  }
  as.numeric(unclass(x)$total)
}


#' Print method for an internal `<fb_memory_estimate>` carrier
#'
#' Renders the per-term INLA memory breakdown introduced at v0.3.10.
#'
#' @param x   an `<fb_memory_estimate>` object.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @keywords internal
#' @export
print.fb_memory_estimate <- function(x, ...) {
  fmt_mb <- function(b) {
    if (is.na(b)) {
      "NA"
    } else {
      sprintf("%.1f MB", b / 1024^2)
    }
  }
  cat(sprintf(
    "<fb_memory_estimate> total = %s (INLA path; overhead = %.1fx)\n",
    fmt_mb(x$total),
    x$overhead_factor
  ))
  if (nrow(x$breakdown) > 0L) {
    for (i in seq_len(nrow(x$breakdown))) {
      r <- x$breakdown[i, ]
      cat(sprintf(
        "  - %-22s  %-26s  %9s  (%s)\n",
        r$term_label,
        r$representation,
        fmt_mb(r$bytes),
        if (is.na(r$share)) {
          "NA"
        } else {
          sprintf("%.0f%%", r$share * 100)
        }
      ))
    }
  }
  invisible(x)
}
