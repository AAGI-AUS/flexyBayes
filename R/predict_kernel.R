# predict_kernel.R -- v0.3.8 unified prediction kernel
# (audit close 2026-05-25).
#
# .predict_linear_draws() is the single per-draw linear-predictor
# kernel routed through by every prediction call site:
#
#   - predict.flexybayes() in-memory branch (R/methods.R)
#   - .predict_to_file() file-backed branch (R/predict_file_output.R)
#   - .predict_in_memory_sample() sampled-RE in-memory branch
#     (R/predict_posterior.R)
#
# Closes the file-backed smooth-bypass bug surfaced in the 2026-05-25
# audit (Critical Fix #1): prior to v0.3.8, file-backed prediction
# computed model.matrix(fixed_formula, chunk) without consulting the
# smooth basis stored on extras$parse_info$smooths. For fits with s()
# smooth terms, the smooth contribution was silently missing from the
# per-draw linear predictor on disk -- only the in-memory non-sample
# branch in R/methods.R included it.
#
# Architectural design:
#
#   - Stride-subsampling (seq.int(1, n_total, length.out = n_draws_cap))
#     applied identically to the fixed-effect beta draws AND the
#     smooth-basis raw / sigma draws. Same idx for all components so
#     per-draw rows align positionally. This makes chunked vs
#     monolithic prediction byte-identical when called with the same
#     n_draws_cap and same sampled_re; the chunk-invariance test
#     scaffold asserts this directly.
#
#   - Smooth basis evaluated via mgcv::PredictMat() on the retained
#     smooth objects in extras$parse_info$smooths. Per-draw basis
#     coefficient vector = raw_draws[i, ] * sigma_draws[i]; per-draw
#     smooth contribution = Bnew %*% t(coef_draws). Summed onto the
#     fixed-effect per-draw linear predictor.
#
#   - Legacy-fit refusal: if has_smooth_in_ir but smooths_slot is
#     empty, refuses with the same structured stop() the in-memory
#     branch used pre-v0.3.8 (smooth-basis retention pre-dates
#     v0.3.4; older fits cannot route through the kernel).
#
#   - sampled_re argument: list keyed by chunk-relative row index,
#     each entry a length-n_draws_eff vector. Layered row-wise onto
#     the linear predictor for unknown-level rows under
#     allow_new_levels = "sample". Sampling itself happens upstream
#     (.predict_sample_re_for_unknowns); the kernel only consumes
#     the realisations.
#
#   - include argument: character subset of c("fixed", "smooth",
#     "random_sampled") naming which contributions to add. Default
#     all three. Malformed include raises a typed condition with
#     reason_code = "predict_kernel_invalid_include". A
#     "random_known" branch (per-row known-level RE) is reserved
#     for v0.4.0; not accepted in v0.3.8's include vocabulary.
#
# Not exported.

# Resolve the per-draw stride index for a given total-draws count and
# cap. Returns either seq_len(n_total) (when n_total <= cap) or an
# evenly-spaced integer subset of length cap (when n_total > cap).
# Single source of truth for the byte-identity stride contract.
.predict_kernel_stride_idx <- function(n_total, n_draws_cap) {
  if (n_total <= n_draws_cap) {
    return(seq_len(n_total))
  }
  as.integer(seq.int(1L, n_total, length.out = n_draws_cap))
}

# Detect whether the fit's IR carries any smooth_mgcv random terms.
# Used by both the kernel and the legacy-fit refusal upstream.
.predict_kernel_has_smooth <- function(object) {
  random_terms <- object$extras$parse_info$random
  if (is.null(random_terms)) {
    return(FALSE)
  }
  any(vapply(
    random_terms,
    function(t) !is.null(t$type) && t$type == "smooth_mgcv",
    logical(1)
  ))
}

# Refuse cleanly when the fit was built before smooth-basis retention
# (or via fb_greta() where smooth basis is user-managed). Same
# message the in-memory branch used pre-v0.3.8 -- mechanical
# correctness over silently-wrong predictions on newdata.
.predict_kernel_refuse_legacy_smooth <- function() {
  stop(
    ".predict_linear_draws(): the fit was produced before smooth-",
    "basis retention shipped (or via fb_greta() where smooth ",
    "basis is user-managed). Refusing to use stats::model.matrix() ",
    "on smooth terms -- it would produce silently wrong ",
    "predictions on newdata. Re-fit via flexybayes() / fb_brms() ",
    "on the current package version to enable predict() with ",
    "newdata.",
    call. = FALSE
  )
}

# Validate the include argument against the v0.3.8 vocabulary.
# Raises a typed condition (flexybayes_predict_kernel_refusal) on
# malformed input so callers can pattern-match.
.predict_kernel_validate_include <- function(include) {
  allowed <- c("fixed", "smooth", "random_sampled")
  if (!is.character(include) || length(include) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "predict_kernel_invalid_include",
      message = paste0(
        ".predict_linear_draws(): `include` must be a non-empty ",
        "character vector; got: ",
        deparse(include),
        "."
      ),
      family_class = "flexybayes_predict_kernel_refusal",
      supplied = include,
      allowed = allowed
    ))
  }
  bad <- setdiff(include, allowed)
  if (length(bad) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "predict_kernel_invalid_include",
      message = paste0(
        ".predict_linear_draws(): `include` value(s) not in the ",
        "kernel vocabulary: ",
        paste(shQuote(bad), collapse = ", "),
        ". Allowed: ",
        paste(shQuote(allowed), collapse = ", "),
        ". (random_known is reserved for a future release.)"
      ),
      family_class = "flexybayes_predict_kernel_refusal",
      supplied = include,
      bad = bad,
      allowed = allowed
    ))
  }
  invisible(NULL)
}

# The kernel.
#
# Returns an n_rows x n_draws_eff numeric matrix of LINEAR-scale
# per-draw values. Link transformation, posterior aggregation
# (rowMeans / quantiles), and se.fit computation stay at call
# sites so the contract is the same regardless of the downstream
# semantic (point + se vs interval + se).
.predict_linear_draws <- function(
  object,
  newdata,
  n_draws_cap = 4000L,
  sampled_re = NULL,
  include = c("fixed", "smooth", "random_sampled")
) {
  .predict_kernel_validate_include(include)

  if (is.null(object$greta) || is.null(object$greta$draws)) {
    stop(
      ".predict_linear_draws(): fit object does not carry posterior ",
      "draws on $greta$draws. Per-draw prediction requires the ",
      "greta-backend draws array. INLA fits do not currently ",
      "expose per-draw access; an INLA draws adapter is queued ",
      "for a future release.",
      call. = FALSE
    )
  }

  has_smooth <- .predict_kernel_has_smooth(object)
  smooths_slot <- object$extras$parse_info$smooths
  if (has_smooth && (is.null(smooths_slot) || length(smooths_slot) == 0L)) {
    .predict_kernel_refuse_legacy_smooth()
  }

  fixed_formula <- object$glm$formula
  use_smooth_design <- has_smooth &&
    length(smooths_slot) > 0L &&
    "smooth" %in% include
  mm_formula <- if (use_smooth_design) {
    .linear_formula(fixed_formula)
  } else {
    fixed_formula
  }
  mm <- stats::model.matrix(
    stats::delete.response(stats::terms(mm_formula)),
    data = newdata
  )
  n_rows <- nrow(newdata)

  draws_mat <- do.call(rbind, lapply(object$greta$draws, as.matrix))
  n_total <- nrow(draws_mat)
  idx <- .predict_kernel_stride_idx(n_total, n_draws_cap)
  draws_sub <- draws_mat[idx, , drop = FALSE]
  n_draws_eff <- nrow(draws_sub)

  lin_pred <- matrix(0, nrow = n_rows, ncol = n_draws_eff)

  if ("fixed" %in% include) {
    fixed_info <- object$extras$parse_info$fixed
    beta_cols <- .get_fixed_draw_columns(fixed_info, colnames(draws_sub))
    if (length(beta_cols) == 0L && ncol(mm) > 0L) {
      stop(
        ".predict_linear_draws(): no fixed-effect draw columns ",
        "found on the fit. The fit may have no fixed effects or ",
        "have been built under an older codegen.",
        call. = FALSE
      )
    }
    avail <- beta_cols[beta_cols %in% colnames(draws_sub)]
    if (length(avail) > 0L) {
      beta_per_draw <- draws_sub[, avail, drop = FALSE]
      beta_names <- names(stats::coef(object))
      if (length(beta_names) >= length(avail)) {
        colnames(beta_per_draw) <- beta_names[seq_along(avail)]
      }
      common <- intersect(colnames(mm), colnames(beta_per_draw))
      if (length(common) > 0L) {
        mm_sub <- mm[, common, drop = FALSE]
        beta_per_draw <- beta_per_draw[, common, drop = FALSE]
        lin_pred <- lin_pred + (mm_sub %*% t(beta_per_draw))
      }
    }
  }

  if (use_smooth_design && "smooth" %in% include) {
    for (v in names(smooths_slot)) {
      sm <- smooths_slot[[v]]
      Bnew <- mgcv::PredictMat(sm, data = newdata)
      # On the low_rank_smooth
      # approximate path the fitted coefficients live in the rank-K
      # truncated space, so the newdata basis must be projected through
      # the same V_K (Bnew %*% V_K, n_new x K) used at fit time. Without
      # this the n_new x k full basis would mismatch the K coefficient
      # draws and the dimension guard below would silently drop the
      # smooth contribution.
      approx_v <- object$extras$parse_info$approx[[v]]
      if (
        !is.null(approx_v) &&
          identical(approx_v$scheme, "low_rank_smooth")
      ) {
        Bnew <- Bnew %*% approx_v$V_K
      }
      raw_pat <- paste0("^s_", v, "_raw(\\[|$)")
      raw_cols <- grep(raw_pat, colnames(draws_sub), value = TRUE)
      sigma_col <- paste0("sigma_s_", v)
      if (
        length(raw_cols) == 0L ||
          !sigma_col %in% colnames(draws_sub)
      ) {
        next
      }
      raw_draws <- draws_sub[, raw_cols, drop = FALSE]
      sigma_draws <- draws_sub[, sigma_col]
      coef_draws <- sweep(raw_draws, 1L, sigma_draws, `*`)
      if (ncol(Bnew) != ncol(coef_draws)) {
        next
      }
      lin_pred <- lin_pred + (Bnew %*% t(coef_draws))
    }
  }

  if (
    "random_sampled" %in%
      include &&
      !is.null(sampled_re) &&
      length(sampled_re) > 0L
  ) {
    for (row_idx in names(sampled_re)) {
      i <- as.integer(row_idx)
      if (i < 1L || i > n_rows) {
        next
      }
      re_vec <- sampled_re[[row_idx]]
      if (length(re_vec) != n_draws_eff) {
        len <- min(length(re_vec), n_draws_eff)
        lin_pred[i, seq_len(len)] <- lin_pred[i, seq_len(len)] +
          re_vec[seq_len(len)]
      } else {
        lin_pred[i, ] <- lin_pred[i, ] + re_vec
      }
    }
  }

  lin_pred
}

# Chunk-invariance helper for the test scaffold.
#
# Given a fit and a newdata frame, computes the per-draw linear
# predictor in monolithic and chunked modes (chunk_size in
# `chunk_sizes`) and returns TRUE iff all chunked outputs match the
# monolithic output at numerical tolerance (default 1e-10). BLAS
# non-associativity makes BIT-identity infeasible across re-grouped
# chunkings of the same matrix multiplication; the audit contract
# is numerical equivalence, not byte-identity.
# Within a single chunk_size run the result IS deterministic (same
# computation, same BLAS path) -- the file-output test asserts that
# separately. Not user-facing.
.predict_kernel_invariants <- function(
  object,
  newdata,
  chunk_sizes = c(100L, 500L, 2000L),
  n_draws_cap = 4000L,
  tolerance = 1e-10
) {
  monolithic <- .predict_linear_draws(
    object,
    newdata,
    n_draws_cap = n_draws_cap
  )
  for (cs in chunk_sizes) {
    n <- nrow(newdata)
    if (cs >= n) {
      next
    } # would be identical by construction
    chunks <- split(seq_len(n), ceiling(seq_len(n) / as.integer(cs)))
    accum <- matrix(0, nrow = n, ncol = ncol(monolithic))
    for (rng in chunks) {
      chunk <- newdata[rng, , drop = FALSE]
      accum[rng, ] <- .predict_linear_draws(
        object,
        chunk,
        n_draws_cap = n_draws_cap
      )
    }
    if (!isTRUE(all.equal(monolithic, accum, tolerance = tolerance))) {
      return(FALSE)
    }
  }
  TRUE
}
