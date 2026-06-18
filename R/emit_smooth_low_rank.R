# Low-rank smooth approximation --- v0.4.0
# (first registered scheme `low_rank_smooth`)
#
# The truncation engine behind the `low_rank_smooth` approximation
# scheme. An s() smooth contributes an n x k dense basis matrix B
# (built by mgcv::smoothCon() at parse time and bound into the greta
# model as `B_s_<var>` --- see R/codegen.R). For a large basis
# dimension k the exact dense path costs k coefficients and an n x k
# design block; the low-rank scheme replaces B with its rank-K
# principal-component truncation B_K = B V_K, where V_K holds the top
# K right singular vectors of B. The model then carries K (rather than
# k) basis coefficients, and prediction projects the newdata basis
# through the same V_K (see R/predict_kernel.R).
#
# This file ships three things:
#
#   1. .truncate_smooth_basis()   the SVD-based rank-K truncation,
#                                 returning B_K, the projection V_K,
#                                 the full singular values, and the
#                                 realised Frobenius capture.
#   2. .validate_low_rank_rank()  the rank refusal contract (positive
#                                 integer; rank <= min(basis rank, n)).
#   3. .validate_low_rank_smooth() the per-scheme validation procedure
#                                 registered into .approximation_registry
#                                 as the `low_rank_smooth` validation_fn;
#                                 validate_approximation() dispatches
#                                 here (see R/validate_approximation.R).
#
# Bias bound. For the rank-K
# truncation B_K of B, the analytical bias bound is the relative
# squared Frobenius residual
#
#     bias(B, K) = ||B - B_K||_F^2 / ||B||_F^2 = 1 - capture(B, K),
#
# and with the SVD B = U D V^T this is exact in the singular values:
#
#     capture(B, K) = sum_{i<=K} d_i^2 / sum_i d_i^2.
#
# The default pass threshold is capture >= 0.99. The
# validation procedure reports the realised capture, the threshold, the
# pass/fail verdict, and the fallback hint --- the user keeps the
# judgement, the contract surfaces the number.
#
# Backend scope. The mgcv dense basis B exists only on the greta
# emit path; the INLA backend represents smooths via its own rw2
# random-walk path and exposes no truncatable dense basis (see
# R/emit_inla.R). The `low_rank_smooth` scheme is therefore a
# greta-backend approximation; routing a low-rank smooth to INLA
# refuses upstream in dispatch rather than silently fitting an
# unrelated rw2 smooth.

# --- rank refusal contract ---------------------------------------- #

# .validate_low_rank_rank() --- the rank refusal contract for a
# low_rank_smooth approximation of a k-column basis built on n rows.
# Refuses non-positive-integer ranks and ranks exceeding the
# truncation ceiling min(k, n) (a rank at or above the basis
# dimension is not an approximation --- it is the exact basis, which
# the user should request via the default exact path). Returns the
# validated rank as an integer.
.validate_low_rank_rank <- function(rank, k, n, var = NULL) {
  where <- if (is.null(var)) "" else paste0(" for smooth s(", var, ")")

  if (
    length(rank) != 1L ||
      is.na(rank) ||
      !is.numeric(rank) ||
      rank != as.integer(rank) ||
      rank < 1L
  ) {
    stop(.fb_refusal_condition(
      reason_code = "low_rank_rank_invalid",
      message = paste0(
        "low_rank_smooth rank",
        where,
        " must be a single positive ",
        "integer; got ",
        .format_rank_for_message(rank),
        "."
      ),
      family_class = "flexybayes_low_rank_rank_refusal"
    ))
  }

  rank <- as.integer(rank)
  max_rank <- min(as.integer(k), as.integer(n))

  if (rank > max_rank) {
    stop(.fb_refusal_condition(
      reason_code = "low_rank_rank_exceeds_basis",
      message = paste0(
        "low_rank_smooth rank",
        where,
        " (",
        rank,
        ") exceeds the ",
        "truncation ceiling min(basis dimension k = ",
        k,
        ", n = ",
        n,
        ") = ",
        max_rank,
        ". A rank at or above the ",
        "basis dimension reproduces the exact basis; request the ",
        "exact smooth (drop the approximation) instead, or choose ",
        "rank <= ",
        max_rank,
        "."
      ),
      family_class = "flexybayes_low_rank_rank_refusal"
    ))
  }

  rank
}

# .format_rank_for_message() --- compact, NA-safe rendering of a
# rejected rank for the refusal message (avoids deparse noise on
# vectors / lists).
.format_rank_for_message <- function(rank) {
  if (length(rank) != 1L) {
    return(paste0("a length-", length(rank), " value"))
  }
  if (is.na(rank)) {
    return("NA")
  }
  if (is.numeric(rank)) {
    return(format(rank))
  }
  paste0("'", as.character(rank), "'")
}


# --- truncation engine -------------------------------------------- #

# .truncate_smooth_basis() --- rank-K principal-component truncation
# of the n x k smooth basis matrix `X`. Returns the truncated basis
# B_K = X V_K (n x K), the projection V_K (k x K, the top-K right
# singular vectors), the full singular-value spectrum, and the
# realised Frobenius capture. The economy SVD gives min(n, k)
# singular values and the k x min(n, k) right-singular-vector matrix;
# the capture denominator sums over the full spectrum so the bound is
# exact regardless of K.
#
# Arguments
#   X      the n x k dense smooth basis (term$X from mgcv::smoothCon()).
#   rank   the validated truncation rank K (use .validate_low_rank_rank()
#          on the caller side before calling this).
#   var    optional smooth variable name, threaded into refusal text.
#
# Returns a list with elements
#   B_K                n x K truncated basis (the design block the
#                      greta model carries).
#   V_K                k x K projection (newdata basis projects through
#                      this at predict time).
#   singular_values    the full singular-value spectrum of X (length
#                      min(n, k)).
#   frobenius_capture  sum(d[1:K]^2) / sum(d^2), in [0, 1].
#   rank               K (integer).
#   k                  the full basis dimension (integer).
.truncate_smooth_basis <- function(X, rank, var = NULL) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(
      ".truncate_smooth_basis(): `X` must be a numeric matrix.",
      call. = FALSE
    )
  }

  n <- nrow(X)
  k <- ncol(X)
  rank <- .validate_low_rank_rank(rank, k = k, n = n, var = var)

  sv <- svd(X)
  d <- sv$d
  v_k <- sv$v[, seq_len(rank), drop = FALSE]
  b_k <- X %*% v_k

  total <- sum(d^2)
  capture <- if (total > 0) sum(d[seq_len(rank)]^2) / total else 0
  capture <- min(max(capture, 0), 1)

  list(
    B_K = b_k,
    V_K = v_k,
    singular_values = d,
    frobenius_capture = capture,
    rank = rank,
    k = k
  )
}


# --- per-scheme validation procedure ------------------------------ #

# .validate_low_rank_smooth() --- the `low_rank_smooth` validation
# procedure, registered into .approximation_registry as the scheme's
# validation_fn at .onLoad() and dispatched to by
# validate_approximation() (R/validate_approximation.R). Reads the
# truncation metadata recorded on the fit at emit time
# (fit$extras$parse_info$approx, keyed by smooth variable) and reports
# the realised Frobenius capture against the pass threshold.
#
# A fit may carry more than one low-rank smooth; the procedure returns
# one per-smooth row plus an overall verdict (all smooths must clear
# the threshold for the fit to pass). The return value is the
# <fb_approximation_validation> classed object built by
# .new_approximation_validation() (R/validate_approximation.R).
#
# Arguments
#   fit        a fitted flexybayes object carrying at least one
#              low_rank_smooth approximation.
#   threshold  the Frobenius-capture pass threshold; defaults to
#              the value 0.99.
#   ...        ignored; present for the generic's signature.
.validate_low_rank_smooth <- function(fit, threshold = 0.99, ...) {
  approx <- fit$extras$parse_info$approx
  schemes <- vapply(
    approx,
    function(a) a$scheme %||% NA_character_,
    character(1)
  )
  low_rank <- approx[schemes == "low_rank_smooth"]

  if (length(low_rank) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "approximation_absent",
      message = paste0(
        "validate_approximation(): the fit carries no low_rank_smooth ",
        "smooth to validate. The scheme registered but no smooth term ",
        "was routed through the truncation path."
      ),
      family_class = "flexybayes_approximation_absent"
    ))
  }

  entry <- .lookup_approximation("low_rank_smooth")

  per_smooth <- lapply(names(low_rank), function(v) {
    a <- low_rank[[v]]
    capture <- a$frobenius_capture
    list(
      smooth = v,
      scheme = "low_rank_smooth",
      rank = a$rank,
      k = a$k,
      frobenius_capture = capture,
      bias_bound = 1 - capture,
      threshold = threshold,
      pass = isTRUE(capture >= threshold)
    )
  })
  names(per_smooth) <- names(low_rank)

  .new_approximation_validation(
    scheme = "low_rank_smooth",
    per_smooth = per_smooth,
    pass = all(vapply(per_smooth, `[[`, logical(1), "pass")),
    threshold = threshold,
    fallback_hint = entry$fallback_hint
  )
}
