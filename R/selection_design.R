# Selection design --- the protected reduction and the term algebra
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 5.2, 5.3, 6
# and 12. This file carries oracles 5 and 6 of the correctness ledger.
#
#   Oracle 5  the protected reduction. "After accounting for X0" is not an
#             implementation, so the contract names one: Frisch-Waugh-Lovell
#             on the whitened design. The provider is handed the residualised
#             candidate block and nothing else, and the protected
#             coefficients come back by GLS conditional on whatever the
#             provider returned.
#
#   Oracle 6  term resolution and contrast invariance. A factor's dummy
#             columns depend on its reference level and its contrast matrix,
#             so a coefficient column is not a stable thing to select. The
#             term is. Each term block is carried into an orthonormal basis
#             before it is penalised, because ordinary group lasso is
#             invariant to orthogonal transformations and to nothing wider,
#             while a contrast map is any invertible matrix at all.
#
# The two are one file because the block transforms act on the residualised
# design the reduction produces: the norm the group penalty sees is the norm
# of the fitted contribution AFTER the protected effects are accounted for,
# not before, and separating them would invite a caller to transform the raw
# block and lose that.
#
# Nothing here fits anything. There is no solver in this file and none is
# wanted: section 14 rules out reimplementing glmnet or grpreg, and the
# design side is complete when the provider can be handed a problem whose
# answer means the same thing in every parameterisation of the model.


# --- estimable bases -------------------------------------------------- #

# .fb_sel_estimable_basis() --- a rank-revealing thin factorisation.
#
# Returns an orthonormal Q spanning the column space of X together with the
# R that maps coefficients into that basis, and the NAMES of the columns
# that were kept and dropped. The names are the point. A silent rank drop
# is the defect this function exists to prevent: the caller can then refuse,
# or report, but not proceed in ignorance.
#
# The pivoting convention is R's own --- `qr()` with the dqrdc2 tolerance,
# which is what `lm()` uses to decide which coefficients come back NA. That
# makes "the declared estimable parameterisation" of section 5.3 something a
# user has already seen rather than a private convention of this module:
# earlier columns are kept, later collinear ones are declared aliased.
#
# Q and R satisfy  X[, basis_columns] = Q %*% R, so for a coefficient vector
# over the estimable columns, gamma = R beta and X beta = Q gamma. The back
# transform is a triangular solve, never a pseudo-inverse of a dense product.
.fb_sel_estimable_basis <- function(X, tol = 1e-7) {
  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  nms <- colnames(X)
  if (is.null(nms)) {
    nms <- if (ncol(X) > 0L) paste0("V", seq_len(ncol(X))) else character(0)
  }

  if (ncol(X) == 0L) {
    return(list(
      Q = matrix(0, nrow = nrow(X), ncol = 0L),
      R = matrix(0, nrow = 0L, ncol = 0L),
      basis_columns = integer(0),
      basis_names = character(0),
      aliased_columns = integer(0),
      aliased_names = character(0),
      rank = 0L,
      n_columns = 0L,
      condition = NA_real_,
      singular_values = numeric(0),
      tol = tol
    ))
  }

  qx <- qr(X, tol = tol)
  keep <- sort(qx$pivot[seq_len(qx$rank)])
  drop_cols <- setdiff(seq_len(ncol(X)), keep)

  # Re-factor the estimable columns on their own. Taking Q from the pivoted
  # factorisation of the full matrix would carry the arbitrary orthonormal
  # completion that sits beyond the rank, and that completion is not a basis
  # of anything the data determines.
  qe <- qr(X[, keep, drop = FALSE])
  ord <- qe$pivot
  basis_columns <- keep[ord]

  # The singular values are an independent read on the rank. `qr()`'s
  # pivoting is a heuristic that stops at the first negligible column it
  # reaches; the spectrum is not, so it is recorded and the caller can see
  # when the two would disagree.
  sv <- svd(X, nu = 0L, nv = 0L)$d
  cond <- if (length(sv) == 0L || min(sv) <= 0) Inf else max(sv) / min(sv)

  list(
    Q = qr.Q(qe),
    R = qr.R(qe),
    basis_columns = basis_columns,
    basis_names = nms[basis_columns],
    aliased_columns = drop_cols,
    aliased_names = nms[drop_cols],
    rank = qx$rank,
    n_columns = ncol(X),
    condition = cond,
    singular_values = sv,
    tol = tol
  )
}


# --- the protected reduction (oracle 5) ------------------------------- #

# .fb_sel_protected_reduction() --- Frisch-Waugh-Lovell on the whitened
# design.
#
#   y*  = W y        X0* = W X0        Xp* = W Xp
#   y_tilde  = y*  - Q0 (Q0' y*)
#   Xp_tilde = Xp* - Q0 (Q0' Xp*)
#
# with Q0 an orthonormal basis of the ESTIMABLE column space of X0*. The
# residual maker I - Q0 Q0' is never formed: it is n by n, and at the 73,341
# rows the deposit carries that is about 43 GB for a matrix whose only use is
# to be multiplied by something narrow.
#
# Why the provider must not see X0 at all, rather than seeing it with a zero
# penalty factor: `glmnet` renormalises `penalty.factor` so the vector sums
# to the number of variables, so appending unpenalised columns rescales every
# other entry and the same printed lambda is a different problem. A zero in
# that vector is not a hole in the penalty, it is a change of units.
#
# `offset` is subtracted here and nowhere else, BEFORE whitening. The order
# matters even though subtraction and whitening are both linear: subtracting
# the raw offset from an already-whitened response leaves W y - o, which is
# not W (y - o) unless W is the identity.
#
# Arguments
#   y                  numeric response, length n.
#   X0                 protected design, n by p0. May have zero columns.
#   Xp                 candidate design, n by pp. May have zero columns.
#   whitener           an <fb_sel_whitener>, or NULL when the caller has
#                      already whitened or V is the identity.
#   offset             numeric length n, or NULL.
#   on_rank_deficient  "refuse" (default) or "report". Section 5.2 asks for
#                      an explicit estimable-space result rather than a
#                      silent drop; both of those are explicit, and which
#                      one is wanted is the caller's decision, not this
#                      function's.
#   tol                rank tolerance passed to the estimable basis.
.fb_sel_protected_reduction <- function(
  y,
  X0,
  Xp,
  whitener = NULL,
  offset = NULL,
  on_rank_deficient = c("refuse", "report"),
  tol = 1e-7
) {
  on_rank_deficient <- match.arg(on_rank_deficient)

  y <- as.numeric(y)
  n <- length(y)
  X0 <- .fb_sel_as_design(X0, n, "X0")
  Xp <- .fb_sel_as_design(Xp, n, "Xp")

  if (!is.null(offset)) {
    offset <- as.numeric(offset)
    if (length(offset) != n) {
      stop(
        "`offset` must have one value per observation: ", n, " expected, ",
        length(offset), " supplied.",
        call. = FALSE
      )
    }
    if (anyNA(offset) || any(!is.finite(offset))) {
      stop("`offset` must be finite at every observation.", call. = FALSE)
    }
    y <- y - offset
  }

  shared <- intersect(colnames(X0), colnames(Xp))
  if (length(shared) > 0L) {
    stop(
      "a column may be protected or a candidate, not both: ",
      paste(shared, collapse = ", "),
      call. = FALSE
    )
  }

  # Whiten once, together, so the three objects cannot fall out of step.
  if (is.null(whitener)) {
    y_star <- y
    X0_star <- X0
    Xp_star <- Xp
  } else {
    if (!inherits(whitener, "fb_sel_whitener")) {
      stop("`whitener` must be an <fb_sel_whitener> or NULL.", call. = FALSE)
    }
    y_star <- drop(.fb_sel_whiten(whitener, matrix(y, ncol = 1L)))
    X0_star <- .fb_sel_whiten(whitener, X0)
    Xp_star <- .fb_sel_whiten(whitener, Xp)
    colnames(X0_star) <- colnames(X0)
    colnames(Xp_star) <- colnames(Xp)
  }

  basis <- .fb_sel_estimable_basis(X0_star, tol = tol)
  if (basis$rank < basis$n_columns) {
    if (identical(on_rank_deficient, "refuse")) {
      stop(.fb_refusal_condition(
        reason_code = "selection_rank_deficient_protected",
        message = paste0(
          "fb_select(): the protected design is rank deficient after ",
          "whitening -- rank ", basis$rank, " over ", basis$n_columns,
          " columns. The estimable space is spanned by: ",
          paste(basis$basis_names, collapse = ", "),
          ". These columns add nothing to it and carry no separately ",
          "identified coefficient: ",
          paste(basis$aliased_names, collapse = ", "),
          ". Drop them, or pass `on_rank_deficient = \"report\"` to proceed ",
          "on the estimable space with the drop recorded."
        ),
        rank = basis$rank,
        estimable = basis$basis_names,
        aliased = basis$aliased_names
      ))
    }
  }

  Q0 <- basis$Q
  # Held rather than recomputed: these two are what the protected
  # coefficients are recovered from, and they are r0 by 1 and r0 by pp, so
  # holding them costs the rank of the protected block rather than a copy of
  # the candidate design. y* and Xp* are exactly recoverable from them and
  # the residualised pair, which is why neither is stored.
  Q0ty_star <- crossprod(Q0, matrix(y_star, ncol = 1L))
  Q0tXp_star <- crossprod(Q0, Xp_star)

  y_tilde <- drop(y_star - Q0 %*% Q0ty_star)
  Xp_tilde <- Xp_star - Q0 %*% Q0tXp_star
  colnames(Xp_tilde) <- colnames(Xp)

  structure(
    list(
      n = n,
      y_tilde = y_tilde,
      Xp_tilde = Xp_tilde,
      Q0 = Q0,
      R0 = basis$R,
      Q0ty_star = drop(Q0ty_star),
      Q0tXp_star = Q0tXp_star,
      X0_star = X0_star,
      protected = colnames(X0),
      protected_estimable = basis$basis_names,
      protected_aliased = basis$aliased_names,
      protected_rank = basis$rank,
      protected_condition = basis$condition,
      candidates = colnames(Xp),
      rank_deficient = basis$rank < basis$n_columns,
      on_rank_deficient = on_rank_deficient,
      whitened = !is.null(whitener),
      whitener_builder = if (is.null(whitener)) NA_character_ else
        whitener$builder,
      offset = offset,
      offset_applied = !is.null(offset),
      tol = tol
    ),
    class = "fb_sel_reduction_fwl"
  )
}


# .fb_sel_as_design() --- coerce and check one design block.
.fb_sel_as_design <- function(X, n, what) {
  if (is.null(X)) {
    return(matrix(0, nrow = n, ncol = 0L))
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("`", what, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) != n) {
    stop(
      "`", what, "` must have one row per observation: ", n, " expected, ",
      nrow(X), " supplied.",
      call. = FALSE
    )
  }
  if (ncol(X) > 0L && is.null(colnames(X))) {
    colnames(X) <- paste0(what, "_", seq_len(ncol(X)))
  }
  X
}


# .fb_sel_provider_payload() --- everything the provider is allowed to see.
#
# The assertion inside is not decoration. The failure this module is most
# exposed to is a protected column reaching the solver by some later
# convenience --- a cbind for a diagnostic, an intercept re-added by a
# default argument --- and the symptom would be a plausible path with the
# wrong lambda scale rather than an error. So the payload is built in one
# place and checked against the protected names on the way out.
.fb_sel_provider_payload <- function(red) {
  stopifnot(inherits(red, "fb_sel_reduction_fwl"))
  leaked <- intersect(colnames(red$Xp_tilde), red$protected)
  if (length(leaked) > 0L) {
    stop(
      "internal contract violation: protected column(s) reached the ",
      "provider payload: ", paste(leaked, collapse = ", "),
      call. = FALSE
    )
  }
  list(
    x = red$Xp_tilde,
    y = red$y_tilde,
    intercept = FALSE,
    n = red$n,
    columns = colnames(red$Xp_tilde)
  )
}


# .fb_sel_reduction_response() --- y* rebuilt from what the object holds.
.fb_sel_reduction_response <- function(red) {
  drop(red$y_tilde + red$Q0 %*% matrix(red$Q0ty_star, ncol = 1L))
}


# .fb_sel_reduction_candidates() --- Xp* rebuilt from what the object holds.
.fb_sel_reduction_candidates <- function(red) {
  out <- red$Xp_tilde + red$Q0 %*% red$Q0tXp_star
  colnames(out) <- red$candidates
  out
}


# .fb_sel_recover_protected() --- beta0 by GLS given a candidate solution.
#
# Profiling the protected block out of the objective is exact and has a
# closed form: with beta_p held, the minimising beta0 solves the whitened
# normal equations of the protected block against y* - Xp* beta_p, which in
# the estimable basis is a single triangular solve.
#
# Aliased columns come back NA rather than 0, the way `lm()` reports them. A
# zero would read as an estimate of no effect; NA reads as what it is, a
# coefficient the design does not identify.
.fb_sel_recover_protected <- function(red, beta_candidate = NULL) {
  stopifnot(inherits(red, "fb_sel_reduction_fwl"))
  beta_candidate <- .fb_sel_check_beta(beta_candidate, red$candidates)

  rhs <- matrix(red$Q0ty_star, ncol = 1L)
  if (length(beta_candidate) > 0L) {
    rhs <- rhs - red$Q0tXp_star %*% matrix(beta_candidate, ncol = 1L)
  }
  out <- stats::setNames(
    rep(NA_real_, length(red$protected)),
    red$protected
  )
  if (red$protected_rank > 0L) {
    est <- drop(backsolve(red$R0, rhs))
    out[red$protected_estimable] <- est
  }
  out
}


# .fb_sel_reduction_fitted() --- whitened-scale fitted values.
#
# On the whitened scale by design. Putting these back on the response scale
# needs V's off-diagonal structure and the fold's training rows, which is the
# prediction estimand of section 5.4 and belongs to that unit, not this one.
# Returning something called "fitted" that quietly meant a marginal
# prediction would be the more useful-looking of the two and the wrong one.
.fb_sel_reduction_fitted <- function(red, beta_candidate = NULL) {
  stopifnot(inherits(red, "fb_sel_reduction_fwl"))
  beta_candidate <- .fb_sel_check_beta(beta_candidate, red$candidates)
  beta0 <- .fb_sel_recover_protected(red, beta_candidate)

  out <- rep(0, red$n)
  if (red$protected_rank > 0L) {
    keep <- red$protected_estimable
    out <- out + drop(
      red$X0_star[, keep, drop = FALSE] %*%
        matrix(beta0[keep], ncol = 1L)
    )
  }
  if (length(beta_candidate) > 0L) {
    out <- out + drop(
      .fb_sel_reduction_candidates(red) %*% matrix(beta_candidate, ncol = 1L)
    )
  }
  out
}


# .fb_sel_reduction_objective() --- the profiled objective.
#
# 1/2 ||y_tilde - Xp_tilde beta_p||^2 + penalty. This is not an
# approximation of the full objective with beta0 at its optimum, it IS that
# objective: (I - Q0 Q0') is idempotent and symmetric, so profiling the
# protected block out leaves the residualised least-squares term exactly.
# The penalty enters unchanged because it is zero on every protected
# coordinate.
.fb_sel_reduction_objective <- function(red, beta_candidate = NULL,
                                        penalty = 0) {
  stopifnot(inherits(red, "fb_sel_reduction_fwl"))
  beta_candidate <- .fb_sel_check_beta(beta_candidate, red$candidates)
  resid <- red$y_tilde
  if (length(beta_candidate) > 0L) {
    resid <- resid - drop(red$Xp_tilde %*% matrix(beta_candidate, ncol = 1L))
  }
  0.5 * sum(resid^2) + penalty
}


# .fb_sel_check_beta() --- a candidate solution in the module's own order.
#
# Matched by NAME when it carries names. A positional match between a
# provider's coefficient vector and the candidate block is one column
# reordering away from attaching every estimate to the wrong covariate, and
# nothing downstream would look wrong.
.fb_sel_check_beta <- function(beta, candidates) {
  if (is.null(beta)) {
    return(stats::setNames(rep(0, length(candidates)), candidates))
  }
  beta <- as.numeric(beta)
  if (length(beta) != length(candidates)) {
    stop(
      "a candidate solution must carry one coefficient per candidate ",
      "column: ", length(candidates), " expected, ", length(beta),
      " supplied.",
      call. = FALSE
    )
  }
  beta
}


# --- term resolution (oracle 6a) -------------------------------------- #

# .fb_sel_term_map() --- the typed resolver: targets are terms, not columns.
#
# The map comes from `terms()` and the `assign` attribute `model.matrix()`
# writes, which is the only place in R where the correspondence between a
# term and the columns it generated is recorded. Reconstructing it from
# column names is the tempting shortcut and it is wrong: a factor level
# called "x1" and a numeric called "x" give column names that no parser can
# separate reliably, and an interaction's name depends on the order the
# formula happened to list its variables.
#
# One numeric main effect is one term and one column. A factor is one term
# and ALL of its contrast columns. An interaction is one term and every
# column it generates. The intercept is assign 0 and is reported separately:
# it belongs in the protected block, never in a candidate group.
.fb_sel_term_map <- function(formula, data, contrasts = NULL) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  # na.pass, deliberately: the map's row order must stay aligned with `data`
  # so a caller can pair it with a response and a fold index. Dropping
  # incomplete rows here would silently renumber every row downstream, and
  # missing responses are section 5.5's business, not this function's.
  mf <- stats::model.frame(
    formula, data,
    na.action = stats::na.pass, drop.unused.levels = FALSE
  )
  mt <- stats::terms(mf)
  X <- stats::model.matrix(mt, mf, contrasts.arg = contrasts)

  assign <- attr(X, "assign")
  labels <- attr(mt, "term.labels")
  orders <- attr(mt, "order")
  factors <- attr(mt, "factors")
  variables <- vapply(
    seq_along(labels),
    function(k) paste(rownames(factors)[factors[, k] > 0], collapse = ", "),
    character(1L)
  )

  intercept_column <- which(assign == 0L)
  term_map <- lapply(
    seq_along(labels),
    function(k) which(assign == k)
  )
  names(term_map) <- labels

  # A term that generated no column at all --- an interaction fully aliased
  # by an emptied level combination, say --- is kept in the table with zero
  # columns rather than dropped, so the caller can see it went missing.
  n_columns <- lengths(term_map)

  involves_factor <- vapply(
    seq_along(labels),
    function(k) {
      vars <- rownames(factors)[factors[, k] > 0]
      any(vapply(
        vars,
        function(v) {
          col <- mf[[v]]
          is.factor(col) || is.character(col) || is.logical(col)
        },
        logical(1L)
      ))
    },
    logical(1L)
  )

  type <- ifelse(
    orders > 1L, "interaction",
    ifelse(involves_factor, "factor", "numeric")
  )

  columns <- data.frame(
    column = colnames(X),
    index = seq_len(ncol(X)),
    term = ifelse(assign == 0L, "(Intercept)", labels[pmax(assign, 1L)]),
    term_index = assign,
    stringsAsFactors = FALSE
  )

  terms_table <- data.frame(
    term = labels,
    term_index = seq_along(labels),
    order = orders,
    n_columns = as.integer(n_columns),
    involves_factor = involves_factor,
    type = type,
    variables = variables,
    stringsAsFactors = FALSE
  )

  list(
    formula = formula,
    terms = mt,
    X = X,
    assign = assign,
    has_intercept = length(intercept_column) > 0L,
    intercept_column = intercept_column,
    term_labels = labels,
    term_map = term_map,
    # The group map is the term map by default, and that default IS the
    # contract: the unit of selection is the term. It is a separate field
    # because a later grammar layer may merge terms with `as_group()`, and a
    # merge must be visible as a departure from the default rather than an
    # edit to it.
    group_map = term_map,
    columns = columns,
    terms_table = terms_table,
    contrasts = attr(X, "contrasts"),
    xlevels = stats::.getXlevels(mt, mf)
  )
}


# .fb_sel_partition_design() --- where the vertical line falls in the design.
#
# Section 2's partition principle in one function: the same model matrix
# becomes the base model, the extended base or the covariate study depending
# only on which terms are protected. The index arithmetic that splits the
# columns and renumbers the group map lives here and nowhere else, because a
# group map whose indices are one column out of step would attach the wrong
# penalty to the wrong term and nothing downstream would look wrong.
#
# The intercept is always protected. It is not a candidate under any
# partition: removing it is a change of model, not a selection.
.fb_sel_partition_design <- function(map, protected = character(0)) {
  protected <- as.character(protected)
  unknown <- setdiff(protected, map$term_labels)
  if (length(unknown) > 0L) {
    stop(
      "these protected terms are not in the model: ",
      paste(unknown, collapse = ", "),
      ". The model's terms are: ",
      paste(map$term_labels, collapse = ", "),
      call. = FALSE
    )
  }

  candidate_terms <- setdiff(map$term_labels, protected)
  protected_columns <- sort(c(
    map$intercept_column,
    unlist(map$term_map[protected], use.names = FALSE)
  ))
  candidate_columns <- sort(unlist(
    map$term_map[candidate_terms],
    use.names = FALSE
  ))
  if (is.null(candidate_columns)) {
    candidate_columns <- integer(0)
  }

  X0 <- map$X[, protected_columns, drop = FALSE]
  Xp <- map$X[, candidate_columns, drop = FALSE]

  # Renumber against Xp, by matching the ORIGINAL column positions rather
  # than by assuming the terms come out in blocks. They do here, but the
  # assumption is exactly the kind that survives every test until a caller
  # protects a term from the middle of the formula.
  group_map <- lapply(
    map$group_map[candidate_terms],
    function(idx) match(idx, candidate_columns)
  )
  names(group_map) <- candidate_terms

  list(
    X0 = X0,
    Xp = Xp,
    protected_terms = protected,
    candidate_terms = candidate_terms,
    protected_columns = protected_columns,
    candidate_columns = candidate_columns,
    group_map = group_map,
    has_intercept = map$has_intercept
  )
}


# --- coordinate semantics --------------------------------------------- #

# .fb_sel_penalty_component_traits() --- which components act coordinate by
# coordinate.
#
# The distinction is whether the penalty can tell one column of a block from
# another. An L2 group norm cannot: rotate the block and the norm is
# unchanged. An L1 sum can, and so its answer depends on which contrast
# matrix happened to be in `options("contrasts")` when the model matrix was
# built. That is not a tuning question, it is a different estimand.
.fb_sel_penalty_component_traits <- function() {
  data.frame(
    component = c(
      "ridge", "lasso", "elastic_net", "group_lasso", "exclusive_lasso",
      "sparse_group_lasso", "group_exclusive", "adaptive"
    ),
    coordinate_wise_l1 = c(
      FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )
}


# .fb_sel_targets() --- the typed targets, each carrying its coordinate
# semantics.
#
#   intrinsic        one column whose meaning survives any reparameterisation
#                    of the model: a numeric main effect, or an interaction
#                    of numerics. Scale is a separate matter, handled by the
#                    loadings.
#   term_subspace    a multi-column block addressed as a whole. Its meaning
#                    is the subspace it spans, which every full-rank contrast
#                    map preserves.
#   basis_dependent  a single contrast column of a multi-column term. Its
#                    meaning changes when the factor is releveled although
#                    the model has not.
#
# `resolution = "coefficient"` is available and it is the reason the third
# level exists. It is not forbidden --- a caller may have a canonical
# scientific basis, and `canonical_basis` is how that is declared --- but it
# is typed so that nothing downstream can treat it as if it were the same
# kind of object as a term.
.fb_sel_targets <- function(map, resolution = c("term", "coefficient"),
                            canonical_basis = FALSE) {
  resolution <- match.arg(resolution)
  tt <- map$terms_table

  if (identical(resolution, "term")) {
    semantics <- ifelse(
      tt$n_columns == 1L & !tt$involves_factor, "intrinsic", "term_subspace"
    )
    out <- data.frame(
      target = tt$term,
      term = tt$term,
      resolution = "term",
      n_columns = tt$n_columns,
      coordinate_semantics = semantics,
      stringsAsFactors = FALSE
    )
  } else {
    cols <- map$columns[map$columns$term_index > 0L, , drop = FALSE]
    hit <- match(cols$term, tt$term)
    single_numeric <- tt$n_columns[hit] == 1L & !tt$involves_factor[hit]
    out <- data.frame(
      target = cols$column,
      term = cols$term,
      resolution = "coefficient",
      n_columns = 1L,
      coordinate_semantics = ifelse(
        single_numeric, "intrinsic", "basis_dependent"
      ),
      stringsAsFactors = FALSE
    )
  }

  out$canonical_basis <- rep(as.logical(canonical_basis), nrow(out))
  out
}


# .fb_sel_assert_coordinate_semantics() --- the normative restriction of
# section 5.3.
#
# This release does not promise contrast-invariant selection WITHIN a factor
# block, and the refusal says so rather than returning an answer whose
# support would change if someone releveled the factor. The check is on the
# pair, not on the target alone: a sparse-group or exclusive component makes
# a multi-column term basis dependent even when the target was resolved as a
# whole term, because the L1 part of the penalty reaches inside the block.
#
# `canonical_basis = TRUE` on a target is the declared escape the contract
# allows. It is recorded on the returned table, so a certificate can say that
# the invariance guarantee was set aside on purpose and by whom.
.fb_sel_assert_coordinate_semantics <- function(targets, components) {
  components <- as.character(components)
  traits <- .fb_sel_penalty_component_traits()
  unknown <- setdiff(components, traits$component)
  if (length(unknown) > 0L) {
    stop(
      "internal contract violation: unknown penalty component(s) ",
      paste(unknown, collapse = ", "),
      ". The grammar layer validates the vocabulary before this point.",
      call. = FALSE
    )
  }

  has_l1 <- any(traits$coordinate_wise_l1[traits$component %in% components])
  effective <- targets$coordinate_semantics
  # A coordinate-wise component reaching inside a multi-column block re-types
  # the target: whatever it was called, what is being selected now is a
  # contrast column. The argument is sharper than "contrasts differ". A block
  # of rank r > 1 has an orthonormal basis Q_g that is itself determined only
  # up to an r by r rotation, so an L1 sum over the coordinates of gamma_g is
  # not merely parameterisation-sensitive, it is not a well-defined function
  # of the model at all. At rank 1 the basis is fixed up to sign, |gamma| is
  # well defined, and the restriction does not bite.
  #
  # The test is on the column count rather than the rank, because the target
  # table is built from the term map, before any block has been factorised.
  # A block whose columns are collinear therefore refuses on a column count
  # its rank does not justify. That is the conservative direction, and the
  # caller who wants it can resolve the block and declare a canonical basis.
  if (has_l1) {
    effective[targets$n_columns > 1L] <- "basis_dependent"
  }
  targets$effective_semantics <- effective

  offending <- targets[
    effective == "basis_dependent" & !targets$canonical_basis, ,
    drop = FALSE
  ]
  if (has_l1 && nrow(offending) > 0L) {
    l1 <- components[
      components %in% traits$component[traits$coordinate_wise_l1]
    ]
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_coordinate_basis_dependent",
      message = paste0(
        "fb_select(): an L1 component (", paste(l1, collapse = ", "),
        ") acts one contrast column at a time, and these targets are ",
        "basis dependent: ", paste(offending$target, collapse = ", "),
        ". Their dummy columns depend on the factor's reference level and ",
        "contrast matrix, so the support would change under a relevel that ",
        "changes no model. This release does not promise contrast-invariant ",
        "within-factor selection. Penalise the term as a group, or declare ",
        "a canonical scientific basis for the target."
      ),
      targets = offending$target,
      components = l1
    ))
  }

  invisible(targets)
}


# --- block transforms (oracle 6b) ------------------------------------- #

# .fb_sel_block_transforms() --- one orthonormal basis per residualised term
# block.
#
# Column-by-column scaling does not buy contrast invariance. Standardise the
# columns of a factor block and ||beta_g||_2 still moves when the contrast
# matrix changes, because the change of basis is not diagonal and is not
# orthogonal either. Ordinary group lasso is invariant to orthogonal
# transformations and to nothing wider, so the block is carried into an
# orthonormal basis before it is penalised:
#
#   X_tilde_g = Q_g R_g,   gamma_g = R_g beta_g,   ||gamma_g||_2 = ||fit_g||_2
#
# and the quantity penalised is the norm of the block's FITTED CONTRIBUTION
# in the GLS-residualised metric. That is a property of the model, not of the
# parameterisation, which is exactly why replacing X_g by X_g A for any
# invertible A leaves the fit, the risk and the support alone.
#
# The weight is sqrt(rank_g), not sqrt(p_g). A block with a collinear column
# has fewer directions than columns, and weighting by the dummy count would
# tax it for a degree of freedom it does not have.
.fb_sel_block_transforms <- function(X_tilde, groups, tol = 1e-7) {
  X_tilde <- as.matrix(X_tilde)
  if (!is.list(groups) || is.null(names(groups))) {
    stop("`groups` must be a named list of column indices.", call. = FALSE)
  }
  nms <- colnames(X_tilde)
  if (is.null(nms)) {
    nms <- paste0("V", seq_len(ncol(X_tilde)))
    colnames(X_tilde) <- nms
  }

  seen <- unlist(groups, use.names = FALSE)
  if (length(seen) > 0L) {
    if (any(seen < 1L | seen > ncol(X_tilde))) {
      stop(
        "`groups` indexes columns outside the design: ",
        paste(sort(unique(seen[seen < 1L | seen > ncol(X_tilde)])),
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    if (anyDuplicated(seen) > 0L) {
      stop(
        "a column may belong to one group only; repeated: ",
        paste(nms[unique(seen[duplicated(seen)])], collapse = ", "),
        call. = FALSE
      )
    }
  }

  blocks <- lapply(names(groups), function(g) {
    cols <- as.integer(groups[[g]])
    basis <- .fb_sel_estimable_basis(
      X_tilde[, cols, drop = FALSE],
      tol = tol
    )
    list(
      group = g,
      columns = cols,
      column_names = nms[cols],
      Q = basis$Q,
      R = basis$R,
      basis_columns = cols[basis$basis_columns],
      basis_names = nms[cols[basis$basis_columns]],
      aliased_columns = cols[basis$aliased_columns],
      aliased_names = nms[cols[basis$aliased_columns]],
      rank = basis$rank,
      n_columns = length(cols),
      condition = basis$condition,
      # sqrt(rank), never sqrt(p): section 5.1.
      weight = sqrt(basis$rank)
    )
  })
  names(blocks) <- names(groups)

  structure(
    list(
      blocks = blocks,
      n = nrow(X_tilde),
      group_names = names(groups),
      n_columns = ncol(X_tilde),
      column_names = nms,
      tol = tol
    ),
    class = "fb_sel_block_transforms"
  )
}


# .fb_sel_block_design() --- the group-orthonormalised design the provider
# actually sees, with the index that says which gamma belongs to which group.
.fb_sel_block_design <- function(bt) {
  stopifnot(inherits(bt, "fb_sel_block_transforms"))
  ranks <- vapply(bt$blocks, function(b) b$rank, integer(1L))
  Q <- do.call(cbind, c(
    lapply(bt$blocks, function(b) b$Q),
    list(matrix(0, nrow = bt$n, ncol = 0L))
  ))
  colnames(Q) <- unlist(
    lapply(bt$blocks, function(b) {
      if (b$rank == 0L) character(0) else paste0(b$group, ".g", seq_len(b$rank))
    }),
    use.names = FALSE
  )
  list(
    Q = Q,
    group = rep(bt$group_names, times = ranks),
    group_index = rep(seq_along(ranks), times = ranks),
    ranks = ranks,
    weights = vapply(bt$blocks, function(b) b$weight, numeric(1L))
  )
}


# .fb_sel_group_norms() --- ||gamma_g||_2 per group.
#
# The support of a group-penalised fit is read off these, and they are the
# invariant quantity: ||gamma_g||_2 is the norm of the fitted contribution,
# so it does not move when the block is reparameterised.
.fb_sel_group_norms <- function(bt, gamma) {
  split_gamma <- .fb_sel_split_gamma(bt, gamma)
  vapply(split_gamma, function(g) sqrt(sum(g^2)), numeric(1L))
}


# .fb_sel_back_transform() --- gamma back to beta in the declared estimable
# parameterisation.
#
# beta_g = R_g^-1 gamma_g on the estimable columns, NA on the aliased ones,
# which is `lm()`'s convention and the one the reduction already uses for the
# protected block. A triangular solve, not a pseudo-inverse: R_g is upper
# triangular and full rank by construction, so there is nothing to regularise
# and no tolerance to choose.
.fb_sel_back_transform <- function(bt, gamma) {
  stopifnot(inherits(bt, "fb_sel_block_transforms"))
  split_gamma <- .fb_sel_split_gamma(bt, gamma)
  out <- stats::setNames(
    rep(NA_real_, bt$n_columns), bt$column_names
  )
  for (g in bt$group_names) {
    b <- bt$blocks[[g]]
    if (b$rank > 0L) {
      out[b$basis_columns] <- drop(backsolve(b$R, split_gamma[[g]]))
    }
  }
  out
}


# .fb_sel_block_fitted() --- the fitted contribution, sum_g Q_g gamma_g.
.fb_sel_block_fitted <- function(bt, gamma) {
  des <- .fb_sel_block_design(bt)
  gamma <- .fb_sel_flatten_gamma(bt, gamma)
  if (length(gamma) == 0L) {
    return(rep(0, bt$n))
  }
  drop(des$Q %*% matrix(gamma, ncol = 1L))
}


# .fb_sel_split_gamma() --- gamma as one vector or one entry per group.
.fb_sel_split_gamma <- function(bt, gamma) {
  stopifnot(inherits(bt, "fb_sel_block_transforms"))
  ranks <- vapply(bt$blocks, function(b) b$rank, integer(1L))
  if (is.list(gamma)) {
    if (!setequal(names(gamma), bt$group_names)) {
      stop(
        "a gamma supplied as a list must carry one entry per group: ",
        paste(bt$group_names, collapse = ", "),
        call. = FALSE
      )
    }
    out <- gamma[bt$group_names]
    bad <- lengths(out) != ranks
    if (any(bad)) {
      stop(
        "gamma for group(s) ",
        paste(bt$group_names[bad], collapse = ", "),
        " must have one entry per basis direction (the group's rank).",
        call. = FALSE
      )
    }
    return(out)
  }
  gamma <- as.numeric(gamma)
  if (length(gamma) != sum(ranks)) {
    stop(
      "gamma must have one entry per basis direction: ", sum(ranks),
      " expected, ", length(gamma), " supplied.",
      call. = FALSE
    )
  }
  ends <- cumsum(ranks)
  starts <- ends - ranks + 1L
  out <- lapply(seq_along(ranks), function(i) {
    if (ranks[[i]] == 0L) numeric(0) else gamma[starts[[i]]:ends[[i]]]
  })
  names(out) <- bt$group_names
  out
}


# .fb_sel_flatten_gamma() --- the inverse of the split, in group order.
.fb_sel_flatten_gamma <- function(bt, gamma) {
  unlist(.fb_sel_split_gamma(bt, gamma), use.names = FALSE)
}
