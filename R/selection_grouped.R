# Grouped penalties on the penalised-optimiser provider.
#
# Contract: final specification row B9 (group lasso, sparse group lasso,
# exclusive lasso and the group-exclusive composite), on the canonical
# objective R/selection_bridges.R states and the whitened,
# protected-residualised substrate R/selection_glmnet.R builds.
#
# The canonical objective is unchanged in shape and wider by two numbers and
# a grouping:
#
#   0.5 ||z - A beta||^2
#     + l1 ||beta||_1 + (l2 / 2) ||beta||_2^2
#     + lg sum_g w_g ||beta_g||_2
#     + (le / 2) sum_g (||beta_g||_1)^2
#
# One proximal operator and one optimality residual cover all seven families,
# which is the same economy the two-number objective already had: ridge,
# laplace and elastic net are the corner where every group is a singleton and
# lg = le = 0, and the residual below reduces to
# `.fb_sel_kkt_residual_lambda()` there exactly rather than approximately.
#
# **Why the grouping is declared and never inferred.** A group is a statement
# about the model, not about the data, and two different scientific questions
# want two different groupings: a factor's contrast columns belong together
# because dropping one level is not a model anybody wrote, while a cluster of
# correlated measurements belongs together because the caller says those
# variables compete. Both are written in the grammar --- `group = "term"` and
# `group = "block"` --- both carry a family default, and both are printed on
# the certificate. Nothing here reads a correlation matrix and decides.
#
# **What the exclusive penalty does and does not promise.** It encourages
# competition inside a group; it does not guarantee exactly one non-zero
# member, and for a singleton group it is exactly a ridge term. The review of
# 2026-09-11 withdrew the "recovers one member per active group" claim and put
# calibrated recovery curves in its place, and the tests follow that: they
# measure the direction and the curve, and assert no guarantee the penalty
# does not carry.

# --- the vocabulary ------------------------------------------------------- #

# .fb_sel_grouped_families() --- the families whose penalty reads a grouping.
.fb_sel_grouped_families <- function() {
  c(
    "group_lasso",
    "sparse_group_lasso",
    "exclusive_lasso",
    "group_exclusive"
  )
}


# .fb_sel_group_rules() --- the two ways a group is declared, and what each
# one says, in the words the certificate prints.
.fb_sel_group_rules <- function() {
  c(
    term = paste0(
      "each model term is one group: a factor's contrast columns move ",
      "together, and a numeric main effect is a group of one"
    ),
    block = paste0(
      "the specification's own target set is one group: the terms it names ",
      "are declared competitors"
    )
  )
}


# .fb_sel_group_index() --- the declared grouping of a column set.
#
# `columns` are the design column names in their declared order and `terms`
# the model term each column came from, which is the map
# `.fb_sel_term_map()` writes and `R/selection_candidate_space.R` already
# carries per fold. Under `rule = "block"` every column named by the
# specification is one group, whatever term it came from; `block_label` names
# it on the certificate.
#
# The weight is sqrt(group size). On this substrate that is also sqrt(rank):
# the fold's candidate space has already dropped every column its training
# rows leave aliased (P2), so the columns reaching a group here are linearly
# independent there, and the sqrt(rank) rule `.fb_sel_block_transforms()`
# states is the same number. Where a fold drops part of a group the size is
# the surviving size, so the weight is recomputed rather than carried.
.fb_sel_group_index <- function(
  columns,
  terms,
  rule = c("term", "block"),
  block_label = "block"
) {
  rule <- match.arg(rule)
  columns <- as.character(columns)
  terms <- as.character(terms)
  if (length(columns) != length(terms)) {
    stop(
      ".fb_sel_group_index(): `columns` and `terms` must name the same ",
      "number of design columns: ",
      length(columns),
      " and ",
      length(terms),
      ".",
      call. = FALSE
    )
  }
  if (!is.character(block_label) || length(block_label) != 1L) {
    stop(
      ".fb_sel_group_index(): `block_label` must be a single string.",
      call. = FALSE
    )
  }

  labels <- if (identical(rule, "term")) {
    terms
  } else {
    rep(block_label, length(columns))
  }
  # The group order follows the declared column order, so a certificate reads
  # in the order the model was written rather than alphabetically.
  names_out <- unique(labels)
  index <- match(labels, names_out)
  sizes <- stats::setNames(
    as.integer(tabulate(index, nbins = length(names_out))),
    names_out
  )

  list(
    rule = rule,
    rule_description = unname(.fb_sel_group_rules()[[rule]]),
    index = as.integer(index),
    names = names_out,
    sizes = sizes,
    weights = stats::setNames(sqrt(as.numeric(sizes)), names_out),
    weight_rule = "sqrt(group size on the fold)",
    declared = data.frame(
      column = columns,
      term = terms,
      group = labels,
      stringsAsFactors = FALSE
    )
  )
}


# .fb_sel_group_parts() --- the (index, weights) pair every solver here reads.
#
# Accepts a full `.fb_sel_group_index()` result or the bare list a caller
# builds by hand, so the residual can be exercised on a grouping that was not
# derived from a design.
.fb_sel_group_parts <- function(groups, p) {
  if (!is.list(groups) || is.null(groups$index)) {
    stop(
      "`groups` must carry an integer `index`, one entry per penalised ",
      "column.",
      call. = FALSE
    )
  }
  index <- as.integer(groups$index)
  if (length(index) != p) {
    stop(
      "`groups$index` must name one group per penalised column: ",
      p,
      " expected, ",
      length(index),
      " supplied.",
      call. = FALSE
    )
  }
  n_groups <- if (length(index) == 0L) 0L else max(index)
  weights <- if (is.null(groups$weights)) {
    sqrt(as.numeric(tabulate(index, nbins = n_groups)))
  } else {
    as.numeric(groups$weights)
  }
  if (length(weights) != n_groups) {
    stop(
      "`groups$weights` must carry one weight per group: ",
      n_groups,
      " expected, ",
      length(weights),
      " supplied.",
      call. = FALSE
    )
  }
  list(index = index, weights = weights, n_groups = n_groups)
}


# .fb_sel_penalty_parts() --- the four numbers, from a canonical penalty or a
# bare list.
.fb_sel_penalty_parts <- function(penalty) {
  if (!is.list(penalty)) {
    stop("`penalty` must be a canonical penalty list.", call. = FALSE)
  }
  out <- list(
    l1 = as.numeric(penalty$lambda_l1 %||% 0),
    l2 = as.numeric(penalty$lambda_l2 %||% 0),
    lg = as.numeric(penalty$lambda_group %||% 0),
    le = as.numeric(penalty$lambda_exclusive %||% 0)
  )
  for (nm in names(out)) {
    value <- out[[nm]]
    if (length(value) != 1L || !is.finite(value) || value < 0) {
      stop(
        "the canonical penalty's `",
        nm,
        "` must be a single finite non-negative number.",
        call. = FALSE
      )
    }
  }
  out
}


# --- the proximal operator ------------------------------------------------ #

# .fb_sel_prox_group_block() --- the proximal operator of one group's penalty.
#
# Within a group the penalty is a function of the magnitudes alone, so every
# coordinate of the minimiser carries the sign of the matching coordinate of
# `v` and the problem reduces to the nonnegative orthant. Writing a = |v| and
# u = |x|, stationarity is
#
#   u_j (1 + t l2 + g / r) = a_j - theta,   theta = t l1 + t le S
#
# with S = sum(u), r = ||u||_2 and g = t lg w. Eliminating r gives the block
# factor in closed form: with m = 1 + t l2 and R = ||(a - theta)_+||_2,
#
#   u = (1 / m) max(0, 1 - g / R) (a - theta)_+
#
# so the whole operator is one coordinate threshold followed by one block
# shrinkage, exactly as a sparse group lasso is --- except that when le > 0
# the threshold depends on the answer, because the exclusive penalty charges
# a group for its own total.
#
# Three regimes, and the split is for accuracy and cost rather than taste:
#
#   le == 0              theta = t l1 is constant. Closed form.
#   le > 0, lg == 0      the block factor is one, so theta solves a piecewise
#                        linear equation exactly, by sorting.
#   le > 0, lg > 0       theta solves a scalar equation whose left side is
#                        strictly increasing, so bisection brackets it from
#                        the start and converges to the double's precision.
#
# The naive alternative --- threshold for the exclusive part first, then
# shrink the block --- computes the threshold at the wrong point and
# over-thresholds. It is the falsifier in
# tests/testthat/test-selection-grouped-penalties.R, and it is wrong by order
# one rather than by a tolerance.
.fb_sel_prox_group_block <- function(v, t, l1, l2, lg, le, weight) {
  v <- as.numeric(v)
  if (length(v) == 0L) {
    return(v)
  }
  a <- abs(v)
  m <- 1 + t * l2
  g <- t * lg * weight
  c_exclusive <- t * le

  # kappa and the shifted magnitudes at a given coordinate threshold.
  shifted_at <- function(theta) pmax(a - theta, 0)
  kappa_at <- function(shifted) {
    radius <- sqrt(sum(shifted^2))
    if (radius > 0) max(0, 1 - g / radius) else 0
  }

  if (c_exclusive == 0) {
    theta <- t * l1
  } else if (g == 0) {
    # kappa is one, so theta - t l1 = (t le / m) sum((a - theta)_+), which is
    # piecewise linear in theta and solved exactly on the sorted magnitudes.
    coefficient <- c_exclusive / m
    b <- sort(a - t * l1, decreasing = TRUE)
    totals <- cumsum(b)
    phi <- 0
    for (k in seq_along(b)) {
      candidate <- coefficient * totals[[k]] / (1 + k * coefficient)
      if (candidate < 0) {
        candidate <- 0
      }
      above <- b[[k]] > candidate
      below <- if (k < length(b)) b[[k + 1L]] <= candidate else TRUE
      if (above && below) {
        phi <- candidate
        break
      }
    }
    theta <- t * l1 + phi
  } else {
    # theta - t l1 - (t le / m) kappa(theta) sum((a - theta)_+) is continuous
    # and strictly increasing, negative at t l1 and positive at t l1 + max(a).
    root_at <- function(theta) {
      shifted <- shifted_at(theta)
      theta - t * l1 - (c_exclusive / m) * kappa_at(shifted) * sum(shifted)
    }
    lower <- t * l1
    upper <- t * l1 + max(a)
    for (iteration in seq_len(80L)) {
      middle <- 0.5 * (lower + upper)
      if (root_at(middle) > 0) upper <- middle else lower <- middle
    }
    theta <- 0.5 * (lower + upper)
  }

  shifted <- shifted_at(theta)
  sign(v) * (kappa_at(shifted) / m) * shifted
}


# .fb_sel_prox_grouped() --- the proximal operator over every group.
.fb_sel_prox_grouped <- function(v, t, parts, group_parts) {
  out <- numeric(length(v))
  for (k in seq_len(group_parts$n_groups)) {
    idx <- which(group_parts$index == k)
    out[idx] <- .fb_sel_prox_group_block(
      v[idx],
      t,
      parts$l1,
      parts$l2,
      parts$lg,
      parts$le,
      group_parts$weights[[k]]
    )
  }
  out
}


# --- the optimality residual ---------------------------------------------- #

# .fb_sel_penalty_value_grouped() --- pen(beta) at a point.
.fb_sel_penalty_value_grouped <- function(beta, parts, group_parts) {
  total <- parts$l1 * sum(abs(beta)) + 0.5 * parts$l2 * sum(beta^2)
  for (k in seq_len(group_parts$n_groups)) {
    block <- beta[group_parts$index == k]
    total <- total +
      parts$lg * group_parts$weights[[k]] * sqrt(sum(block^2)) +
      0.5 * parts$le * sum(abs(block))^2
  }
  total
}


# .fb_sel_objective_grouped() --- the canonical objective at a point.
.fb_sel_objective_grouped <- function(A, z, beta, parts, group_parts) {
  residual <- as.numeric(z) - drop(A %*% beta)
  0.5 *
    sum(residual^2) +
    .fb_sel_penalty_value_grouped(beta, parts, group_parts)
}


# .fb_sel_kkt_residual_grouped() --- how far a point is from the canonical
# optimum, in the max norm of the subdifferential.
#
# The smooth part has gradient A'(A beta - z) + l2 beta. The penalty's
# subdifferential at beta splits by group, and the exclusive term contributes
# a threshold the group sets for itself: with S_g = ||beta_g||_1 the
# coordinate threshold is tau_g = l1 + le S_g, because d/d|beta_j| of
# (le / 2) S_g^2 is le S_g.
#
#   group non-zero   the L2 norm is differentiable there, so an active
#                    coordinate reads |grad_j + tau_g sign(beta_j) +
#                    lg w_g beta_j / ||beta_g||_2| and a zero coordinate
#                    inside an active group reads max(|grad_j| - tau_g, 0)
#   group all zero   S_g = 0, so tau_g = l1, and the condition is that
#                    -grad_g lies in the sum of an L-infinity ball of radius
#                    l1 and an L2 ball of radius lg w_g, whose distance is
#                    max(0, ||soft(grad_g, l1)||_2 - lg w_g)
#
# With lg = le = 0 both branches are the two lines
# `.fb_sel_kkt_residual_lambda()` already computes, which is why the widened
# residual can replace it rather than sit beside it. The test asserts that
# identity with `expect_equal()` on the value, not on the formula.
#
# The measurement is of the point, not of the procedure that produced it.
.fb_sel_kkt_residual_grouped <- function(A, z, beta, penalty, groups) {
  A <- as.matrix(A)
  z <- as.numeric(z)
  beta <- as.numeric(beta)
  if (ncol(A) == 0L || length(beta) == 0L) {
    return(0)
  }
  if (nrow(A) != length(z)) {
    stop(
      "`A` and `z` must agree on the number of rows: ",
      nrow(A),
      " and ",
      length(z),
      ".",
      call. = FALSE
    )
  }
  if (ncol(A) != length(beta)) {
    stop(
      "`beta` must carry one coefficient per column of `A`: ",
      ncol(A),
      " expected, ",
      length(beta),
      " supplied.",
      call. = FALSE
    )
  }
  parts <- .fb_sel_penalty_parts(penalty)
  group_parts <- .fb_sel_group_parts(groups, ncol(A))

  gradient <- drop(crossprod(A, drop(A %*% beta) - z)) + parts$l2 * beta
  worst <- 0
  for (k in seq_len(group_parts$n_groups)) {
    idx <- which(group_parts$index == k)
    block <- beta[idx]
    block_gradient <- gradient[idx]
    weight <- group_parts$weights[[k]]
    radius <- sqrt(sum(block^2))
    if (radius > 0) {
      tau <- parts$l1 + parts$le * sum(abs(block))
      full <- block_gradient + parts$lg * weight * block / radius
      active <- abs(block) > 0
      residual <- numeric(length(block))
      residual[active] <- abs(full[active] + tau * sign(block[active]))
      residual[!active] <- pmax(abs(full[!active]) - tau, 0)
      worst <- max(worst, max(residual))
      next
    }
    soft <- .fb_sel_soft_threshold(block_gradient, parts$l1)
    worst <- max(
      worst,
      max(0, sqrt(sum(soft^2)) - parts$lg * weight)
    )
  }
  worst
}


# --- the certified mode --------------------------------------------------- #

# .fb_sel_prox_mode_grouped() --- the canonical mode by FISTA, with a KKT
# certificate and a restricted polish.
#
# The same three ingredients as `.fb_sel_prox_mode_lambda()`, for the same
# reasons: accelerated proximal gradient at step 1 / L, adaptive restart
# because FISTA's momentum is not monotone, and a polish because first-order
# methods identify a support long before they resolve its values.
#
# The polish differs. With a group norm and an exclusive term the restricted
# problem is not linear in the active coefficients --- the threshold and the
# block factor both depend on the answer --- so the polish is a short fixed
# point on the restricted stationarity condition rather than one Cholesky
# solve. It is accepted only when the recomputed KKT residual improves, which
# is a measurement and not an assertion about the iteration.
#
# A point whose residual is still above `tol` when the budget runs out is
# refused with `selection_prox_mode_not_converged`, exactly as the ungrouped
# route refuses: a mode that fails its own optimality check is not a mode.
.fb_sel_prox_mode_grouped <- function(
  A,
  z,
  penalty,
  groups,
  tol = 1e-9,
  max_iter = 2e5
) {
  A <- as.matrix(A)
  z <- as.numeric(z)
  if (nrow(A) != length(z)) {
    stop(
      "`A` and `z` must agree on the number of rows: ",
      nrow(A),
      " and ",
      length(z),
      ".",
      call. = FALSE
    )
  }
  max_iter <- as.integer(max_iter)
  if (is.na(max_iter) || max_iter < 1L) {
    stop("`max_iter` must be at least one iteration.", call. = FALSE)
  }
  parts <- .fb_sel_penalty_parts(penalty)
  p <- ncol(A)
  group_parts <- .fb_sel_group_parts(groups, p)
  family <- penalty$family %||% "canonical"

  if (p == 0L) {
    return(list(
      beta = numeric(0),
      kkt = 0,
      iterations = 0L,
      converged = TRUE,
      family = family,
      objective = 0.5 * sum(z^2),
      polished = FALSE
    ))
  }

  gram <- crossprod(A)
  aty <- drop(crossprod(A, z))
  lipschitz <- max(eigen(gram, symmetric = TRUE, only.values = TRUE)$values) +
    parts$l2
  if (!is.finite(lipschitz) || lipschitz <= 0) {
    lipschitz <- 1
  }
  step <- 1 / lipschitz

  beta <- numeric(p)
  momentum <- beta
  weight <- 1
  objective <- .fb_sel_objective_grouped(A, z, beta, parts, group_parts)
  best <- .fb_sel_kkt_residual_grouped(A, z, beta, penalty, groups)
  polished <- FALSE
  iterations <- 0L
  check_every <- 20L

  for (iteration in seq_len(max_iter)) {
    iterations <- iteration
    gradient <- drop(gram %*% momentum) - aty
    proposal <- .fb_sel_prox_grouped(
      momentum - step * gradient,
      step,
      parts,
      group_parts
    )
    proposal_objective <- .fb_sel_objective_grouped(
      A,
      z,
      proposal,
      parts,
      group_parts
    )
    if (proposal_objective > objective) {
      # Adaptive restart: drop the momentum and step from the incumbent.
      weight <- 1
      gradient <- drop(gram %*% beta) - aty
      proposal <- .fb_sel_prox_grouped(
        beta - step * gradient,
        step,
        parts,
        group_parts
      )
      proposal_objective <- .fb_sel_objective_grouped(
        A,
        z,
        proposal,
        parts,
        group_parts
      )
    }
    next_weight <- (1 + sqrt(1 + 4 * weight^2)) / 2
    momentum <- proposal + ((weight - 1) / next_weight) * (proposal - beta)
    beta <- proposal
    weight <- next_weight
    objective <- proposal_objective

    if (iteration %% check_every == 0L || iteration == max_iter) {
      best <- .fb_sel_kkt_residual_grouped(A, z, beta, penalty, groups)
      if (best <= tol) {
        break
      }
      candidate <- .fb_sel_polish_grouped(
        gram,
        aty,
        beta,
        parts,
        group_parts
      )
      if (!is.null(candidate)) {
        candidate_kkt <- .fb_sel_kkt_residual_grouped(
          A,
          z,
          candidate,
          penalty,
          groups
        )
        if (candidate_kkt < best) {
          beta <- candidate
          momentum <- candidate
          weight <- 1
          objective <- .fb_sel_objective_grouped(
            A,
            z,
            beta,
            parts,
            group_parts
          )
          best <- candidate_kkt
          polished <- TRUE
        }
        if (best <= tol) {
          break
        }
      }
    }
  }

  if (best > tol) {
    stop(.fb_refusal_condition(
      reason_code = "selection_prox_mode_not_converged",
      message = paste0(
        "fb_select(): the grouped proximal-gradient mode for family `",
        family,
        "` stopped after ",
        iterations,
        " iteration(s) with a KKT residual of ",
        format(best, digits = 3),
        ", above the tolerance ",
        format(tol, digits = 3),
        ". No point is returned: a mode that fails its own optimality check ",
        "is not a mode."
      ),
      kkt = best,
      tol = tol,
      iterations = iterations
    ))
  }

  list(
    beta = stats::setNames(beta, colnames(A)),
    kkt = best,
    iterations = iterations,
    converged = TRUE,
    family = family,
    objective = .fb_sel_objective_grouped(A, z, beta, parts, group_parts),
    polished = polished
  )
}


# .fb_sel_polish_grouped() --- a short fixed point on the restricted
# stationarity condition.
#
# On the current support, and holding the signs and the group scalars from the
# incumbent, the condition is linear:
#
#   (A_S'A_S + l2 I + diag(lg w_g / r_g)) beta_S = A_S'z - tau_g sign(beta_S)
#
# Re-deriving r_g and tau_g from each solve and repeating a few times is the
# fixed point. It returns NULL when there is nothing to polish, when the
# restricted system has no Cholesky factor, or when a sign flips --- a flipped
# sign means the support assumption the solve rested on was wrong, so the
# point is not a candidate at all. The caller recomputes the residual and
# keeps whichever point wins.
.fb_sel_polish_grouped <- function(gram, aty, beta, parts, group_parts) {
  support <- which(abs(beta) > 0)
  if (length(support) == 0L) {
    return(NULL)
  }
  signs <- sign(beta[support])
  index <- group_parts$index[support]
  current <- beta[support]
  system_base <- gram[support, support, drop = FALSE]

  for (round in seq_len(8L)) {
    scale_vector <- numeric(length(support))
    shift <- numeric(length(support))
    for (k in unique(index)) {
      in_group <- which(index == k)
      block <- current[in_group]
      radius <- sqrt(sum(block^2))
      if (radius <= 0) {
        return(NULL)
      }
      scale_vector[in_group] <- parts$lg * group_parts$weights[[k]] / radius
      shift[in_group] <- parts$l1 + parts$le * sum(abs(block))
    }
    system <- system_base
    diag(system) <- diag(system) + parts$l2 + scale_vector
    factor_system <- tryCatch(chol(system), error = function(e) NULL)
    if (is.null(factor_system)) {
      return(NULL)
    }
    rhs <- aty[support] - shift * signs
    solved <- drop(backsolve(
      factor_system,
      forwardsolve(t(factor_system), rhs)
    ))
    if (any(sign(solved) != signs)) {
      return(NULL)
    }
    if (max(abs(solved - current)) <= 1e-14 * max(1, max(abs(solved)))) {
      current <- solved
      break
    }
    current <- solved
  }

  out <- numeric(length(beta))
  out[support] <- current
  out
}


# --- the path anchor ------------------------------------------------------ #

# .fb_sel_group_fold_index() --- one fold's grouping, from the declared one.
#
# A fold fits the candidate columns its own training rows identify
# (R/selection_candidate_space.R). A group whose members are partly dropped
# there is therefore a smaller group on that fold, and this is the one place
# that is decided: the group is fitted on its surviving members and its weight
# is recomputed from the surviving size, because sqrt(group size) is a
# statement about how many directions the penalty is charging for and the fold
# is charging for fewer. A group with no surviving member is absent from the
# fold rather than present at weight zero.
#
# The alternative -- dropping the whole group when any member is dropped --
# was not taken. P2 already chose column-level drops for this route, and
# discarding an identified factor level because a different one is aliased
# removes information the fold has. What the choice must not be is invisible,
# so `complete` travels to the certificate and the test asserts it.
.fb_sel_group_fold_index <- function(grouping, kept_names, fold_label = NA) {
  declared <- grouping$declared
  kept <- declared$column %in% kept_names
  labels <- declared$group[kept]
  names_out <- unique(labels)
  index <- match(labels, names_out)
  sizes <- stats::setNames(
    as.integer(tabulate(index, nbins = length(names_out))),
    names_out
  )
  declared_sizes <- table(declared$group)

  membership <- data.frame(
    fold = rep(as.character(fold_label), length(names_out)),
    group = names_out,
    declared_size = as.integer(declared_sizes[names_out]),
    fitted_size = as.integer(sizes),
    weight = sqrt(as.numeric(sizes)),
    stringsAsFactors = FALSE
  )
  membership$complete <- membership$fitted_size == membership$declared_size
  rownames(membership) <- NULL

  list(
    rule = grouping$rule,
    index = as.integer(index),
    names = names_out,
    sizes = sizes,
    weights = stats::setNames(sqrt(as.numeric(sizes)), names_out),
    columns = declared$column[kept],
    membership = membership
  )
}


# .fb_sel_group_support_history() --- the support over the tuning grid, read
# at the group's resolution.
#
# The column history `.fb_sel_glmnet_support_history()` writes is the finer
# statement and it is kept; this is the same membership summed to the group,
# because a group lasso's answer is about groups. The two statements
# U-SEL-PATH separated are separated here too: `available` is whether the
# fold could estimate any member of the group at all, and `in_support` is
# whether the group carried a non-zero coefficient at that point. A group no
# fold could estimate is reported unavailable rather than never selected.
.fb_sel_group_support_history <- function(points, grouping) {
  empty <- data.frame(
    fold = character(0),
    lambda_index = integer(0),
    s = numeric(0),
    group = character(0),
    available = logical(0),
    in_support = logical(0),
    n_available = integer(0),
    n_in_support = integer(0),
    stringsAsFactors = FALSE
  )
  if (is.null(grouping) || is.null(points) || nrow(points) == 0L) {
    return(empty)
  }
  group_of <- grouping$declared$group[
    match(points$column, grouping$declared$column)
  ]
  points$group <- group_of
  keyed <- split(
    points,
    list(points$fold, points$lambda_index, points$group),
    drop = TRUE
  )
  rows <- lapply(keyed, function(part) {
    data.frame(
      fold = part$fold[[1L]],
      lambda_index = part$lambda_index[[1L]],
      s = part$s[[1L]],
      group = part$group[[1L]],
      available = any(part$available),
      in_support = any(part$in_support),
      n_available = sum(part$available),
      n_in_support = sum(part$in_support),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$fold, out$lambda_index, out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}


# .fb_sel_grouped_lambda_max() --- the smallest multiple of the declared
# penalty at which every coefficient is zero, or Inf where there is none.
#
# A family with a finite null anchor is one whose penalty is non-smooth at the
# origin, so the subdifferential there has interior: the group norm and the L1
# sum both do, and the anchor is the largest scaling at which -A'z still sits
# inside. The exclusive penalty does NOT: (le / 2) S_g^2 has gradient zero at
# the origin, so it behaves like a ridge in the limit and reaches zero only as
# the scaling grows without bound. Reporting a finite anchor for it would put
# the grid where the sequence happened to stop rather than where the family
# has a null, which is the defect `.fb_sel_support_semantics()` already
# refuses to repeat for ridge.
.fb_sel_grouped_lambda_max <- function(A, z, penalty, groups) {
  A <- as.matrix(A)
  if (ncol(A) == 0L) {
    return(0)
  }
  parts <- .fb_sel_penalty_parts(penalty)
  group_parts <- .fb_sel_group_parts(groups, ncol(A))
  if (parts$l1 == 0 && parts$lg == 0) {
    return(Inf)
  }
  correlation <- drop(crossprod(A, as.numeric(z)))
  anchor <- 0
  for (k in seq_len(group_parts$n_groups)) {
    idx <- which(group_parts$index == k)
    block <- correlation[idx]
    weight <- group_parts$weights[[k]]
    # The smallest s with ||soft(block, s l1)||_2 <= s lg w. With lg = 0 this
    # is max|block| / l1; with l1 = 0 it is ||block||_2 / (lg w); in between
    # the left side is convex and decreasing in s, so a bracketed bisection
    # is exact to the double's precision.
    if (parts$lg == 0) {
      anchor <- max(anchor, max(abs(block)) / parts$l1)
      next
    }
    if (parts$l1 == 0) {
      anchor <- max(anchor, sqrt(sum(block^2)) / (parts$lg * weight))
      next
    }
    gap_at <- function(s) {
      soft <- .fb_sel_soft_threshold(block, s * parts$l1)
      sqrt(sum(soft^2)) - s * parts$lg * weight
    }
    upper <- max(abs(block)) / parts$l1
    lower <- 0
    for (iteration in seq_len(80L)) {
      middle <- 0.5 * (lower + upper)
      if (gap_at(middle) > 0) lower <- middle else upper <- middle
    }
    anchor <- max(anchor, 0.5 * (lower + upper))
  }
  anchor
}
