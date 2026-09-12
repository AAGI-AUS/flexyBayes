# =============================================================================
# Selection design: oracles 5 and 6 of the correctness ledger.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 5.2, 5.3, 6
# and 12.
#
#   Oracle 5  the protected reduction. The Frisch-Waugh-Lovell residualised
#             solve must agree with a direct penalised GLS on the whole
#             design -- on the coefficients, on the fitted values and on the
#             objective -- and the provider must never see a protected
#             column.
#
#   Oracle 6  term resolution and contrast invariance. Replacing a term's
#             columns by any full-rank reparameterisation of them, including
#             a different contrast matrix, a different reference level, or a
#             random invertible A, must leave the fit, the risk and the term
#             support exactly where they were.
#
# Nothing here fits a model. Both oracles are deterministic linear algebra on
# small well-conditioned doubles, so the reference is an independent dense
# computation rather than a second engine, and the bar is the contract's:
# relative 1e-8. The measured agreement is four to seven orders tighter than
# that on every check in this file, so the tolerance is a contract floor and
# not a fitted one.
#
# Every gate below ships with the defect that makes it red. Two gates the
# consolidated audit examined had passed for every possible input because
# nobody had ever made them fail on purpose; the negative controls in this
# file exist so that class does not recur here.
# =============================================================================


# --- fixtures ---------------------------------------------------------- #

# An <fb_sel_whitener> built by hand from a known V. Building it here rather
# than fitting a model keeps these oracles independent of the covariance
# builder: if the block-Cholesky route were wrong, that is oracle 1 and 2's
# business, and it must not be able to turn oracle 5 red as well.
fbsd_whitener <- function(V, block) {
  rows_list <- split(seq_along(block), block)
  structure(
    list(
      builder = "block_cholesky",
      estimand = "test_fixture",
      n = nrow(V),
      block_factor = "block",
      n_blocks = length(rows_list),
      factors = lapply(rows_list, function(rows) {
        list(rows = rows, L = t(chol(V[rows, rows, drop = FALSE])))
      }),
      components = list(residual = NA_real_, terms = numeric(0))
    ),
    class = "fb_sel_whitener"
  )
}


# A correlated Gaussian case with a protected block, a candidate block and a
# block-diagonal V: exactly the shape section 5.2 is written for.
fbsd_case <- function(seed = 8121L, n = 90L, block_size = 9L) {
  set.seed(seed)
  block <- rep(seq_len(n / block_size), each = block_size)
  V <- diag(0.6, n) + 1.4 * outer(block, block, "==")
  X0 <- cbind(
    `(Intercept)` = rep(1, n), z1 = rnorm(n), z2 = rnorm(n)
  )
  Xp <- cbind(c1 = rnorm(n), c2 = rnorm(n), c3 = rnorm(n))
  beta <- c(2, -1.2, 0.6, 1.1, 0, -0.8)
  y <- as.numeric(
    cbind(X0, Xp) %*% beta + crossprod(chol(V), rnorm(n))
  )
  list(
    n = n, block = block, V = V, X0 = X0, Xp = Xp, y = y,
    whitener = fbsd_whitener(V, block)
  )
}


# The independent oracle for oracle 5: penalised GLS on the whole design, in
# the ORIGINAL basis, with V inverted densely. Deliberately the computation
# the module refuses to perform -- forming V^-1 and a joint normal-equation
# solve -- so that agreement is evidence about the reduction rather than
# about a shared implementation.
fbsd_gls <- function(y, X0, Xp, V, lambda = 0) {
  X <- cbind(X0, Xp)
  Vi <- solve(V)
  pen <- diag(
    c(rep(0, ncol(X0)), rep(lambda, ncol(Xp))),
    nrow = ncol(X), ncol = ncol(X)
  )
  beta <- drop(solve(t(X) %*% Vi %*% X + pen, t(X) %*% Vi %*% y))
  names(beta) <- colnames(X)
  resid <- y - drop(X %*% beta)
  list(
    beta = beta,
    beta0 = beta[seq_len(ncol(X0))],
    beta_p = beta[ncol(X0) + seq_len(ncol(Xp))],
    fitted = drop(X %*% beta),
    objective = 0.5 * drop(t(resid) %*% Vi %*% resid) +
      0.5 * lambda * sum(beta[ncol(X0) + seq_len(ncol(Xp))]^2)
  )
}


# Block-coordinate descent for the group lasso, used ONLY as an oracle for
# oracle 6. With an orthonormal Q_g the block update is an exact group soft
# threshold, so this converges to the global optimum of a convex problem and
# needs no line search, no step size and no random start. It is not a
# provider and nothing in R/ calls it: section 14 rules out shipping a
# solver, and this one exists to make "identical support" a statement about
# something rather than a statement about nothing.
fbsd_group_lasso <- function(y, Qs, weights, lambda,
                             max_iter = 10000L, tol = 1e-13) {
  gamma <- lapply(Qs, function(Q) rep(0, ncol(Q)))
  fitted <- rep(0, length(y))
  for (iter in seq_len(max_iter)) {
    delta <- 0
    for (i in seq_along(Qs)) {
      Q <- Qs[[i]]
      if (ncol(Q) == 0L) {
        next
      }
      contrib <- drop(Q %*% gamma[[i]])
      z <- drop(crossprod(Q, y - fitted + contrib))
      nz <- sqrt(sum(z^2))
      thr <- lambda * weights[[i]]
      new <- if (nz <= thr) z * 0 else z * (1 - thr / nz)
      delta <- max(delta, max(abs(new - gamma[[i]])))
      fitted <- fitted - contrib + drop(Q %*% new)
      gamma[[i]] <- new
    }
    if (delta < tol) {
      break
    }
  }
  list(gamma = gamma, fitted = fitted, iterations = iter)
}


# The oracle-6 data: two numerics, a four-level factor, a three-level factor,
# a factor-by-numeric interaction, and a block-correlated V so the invariance
# claim is tested after whitening rather than in the iid corner where a
# contrast map is easier to get right by accident.
fbsd_terms_case <- function(seed = 4021L, n = 150L, block_size = 10L) {
  set.seed(seed)
  d <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    f = factor(sample(letters[1:4], n, TRUE)),
    h = factor(sample(c("lo", "mid", "hi"), n, TRUE))
  )
  d$block <- rep(seq_len(n / block_size), each = block_size)
  d$y <- 1 + 0.9 * d$x1 - 1.6 * d$x2 +
    2.4 * (d$f == "c") - 1.8 * (d$f == "d") +
    0.7 * d$x1 * (d$f == "b") + rnorm(n, sd = 0.8)
  V <- diag(0.5, n) + 0.9 * outer(d$block, d$block, "==")
  list(data = d, V = V, whitener = fbsd_whitener(V, d$block))
}


# One pass of the whole design-side pipeline: term map, partition,
# residualisation, block transforms, group lasso. `transform` optionally
# replaces one term's columns by X_g %*% A, which is the metamorphic move
# section 5.3 asks for.
fbsd_pipeline <- function(case, data = NULL, contrasts_opt = NULL,
                          transform = NULL, term = NULL, lambda = 2) {
  if (!is.null(contrasts_opt)) {
    # on.exit, not a bare restore afterwards: a failing expectation in the
    # caller would otherwise leave the session's contrasts changed and take
    # every later test file down with it.
    old <- options(contrasts = contrasts_opt)
    on.exit(options(old), add = TRUE)
  }
  d <- if (is.null(data)) case$data else data
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, d)
  part <- .fb_sel_partition_design(map, protected = character(0))
  Xp <- part$Xp
  if (!is.null(transform)) {
    cols <- part$group_map[[term]]
    Xp[, cols] <- Xp[, cols, drop = FALSE] %*% transform
  }
  red <- .fb_sel_protected_reduction(
    d$y, part$X0, Xp,
    whitener = case$whitener
  )
  bt <- .fb_sel_block_transforms(red$Xp_tilde, part$group_map)
  des <- .fb_sel_block_design(bt)
  gl <- fbsd_group_lasso(
    red$y_tilde,
    lapply(bt$blocks, function(b) b$Q),
    des$weights,
    lambda = lambda
  )
  norms <- .fb_sel_group_norms(bt, gl$gamma)
  list(
    map = map, part = part, red = red, bt = bt, design = des,
    gamma = gl$gamma,
    fitted = gl$fitted,
    risk = sqrt(mean((red$y_tilde - gl$fitted)^2)),
    norms = norms,
    support = names(norms)[norms > 1e-9],
    beta = .fb_sel_back_transform(bt, gl$gamma),
    iterations = gl$iterations
  )
}


# =============================================================================
# ORACLE 5 --- the protected reduction (contract section 5.2)
# =============================================================================

# --- the fixture is sound before anything is asked of it --------------- #

test_that("the hand-built whitener satisfies W V W' = I", {
  # A precondition, not an oracle. If the fixture's whitener were wrong,
  # every disagreement below would be blamed on the reduction.
  case <- fbsd_case()
  WV <- .fb_sel_whiten(case$whitener, case$V)
  WVWt <- t(.fb_sel_whiten(case$whitener, t(WV)))
  expect_equal(WVWt, diag(1, case$n), tolerance = 1e-10)
})


# --- ORACLE 5 (i): the residualised solve equals direct GLS ------------ #

test_that("oracle 5: the residualised solve equals direct GLS at lambda 0", {
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  direct <- fbsd_gls(case$y, case$X0, case$Xp, case$V, lambda = 0)

  # The reduced problem is an ordinary least-squares solve on the
  # residualised pair, which needs no solver at all.
  beta_p <- drop(solve(
    crossprod(red$Xp_tilde),
    crossprod(red$Xp_tilde, red$y_tilde)
  ))
  names(beta_p) <- red$candidates

  # Measured: candidate coefficients agree to 3.7e-16 relative, protected to
  # 1.2e-16, fitted values to 2.9e-16, objective to 1.4e-16. The bar is 1e-8.
  expect_equal(beta_p, direct$beta_p, tolerance = 1e-8)
  expect_equal(
    .fb_sel_recover_protected(red, beta_p), direct$beta0,
    tolerance = 1e-8
  )
  expect_equal(
    .fb_sel_reduction_fitted(red, beta_p),
    drop(.fb_sel_whiten(case$whitener, matrix(direct$fitted, ncol = 1L))),
    tolerance = 1e-8
  )
  expect_equal(
    .fb_sel_reduction_objective(red, beta_p), direct$objective,
    tolerance = 1e-8
  )
})


test_that("oracle 5: the agreement holds under a live penalty", {
  # Ridge on the candidates only has a closed form on both sides, so a
  # genuinely penalised comparison costs nothing and closes the acceptance
  # clause properly: at lambda 0 a reduction that mishandled the penalty
  # partition would still pass.
  case <- fbsd_case()
  lambda <- 3.5
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  direct <- fbsd_gls(case$y, case$X0, case$Xp, case$V, lambda = lambda)

  p <- ncol(case$Xp)
  beta_p <- drop(solve(
    crossprod(red$Xp_tilde) + lambda * diag(p),
    crossprod(red$Xp_tilde, red$y_tilde)
  ))
  names(beta_p) <- red$candidates

  expect_equal(beta_p, direct$beta_p, tolerance = 1e-8)
  expect_equal(
    .fb_sel_recover_protected(red, beta_p), direct$beta0,
    tolerance = 1e-8
  )
  expect_equal(
    .fb_sel_reduction_objective(
      red, beta_p,
      penalty = 0.5 * lambda * sum(beta_p^2)
    ),
    direct$objective,
    tolerance = 1e-8
  )
})


test_that("the profiled objective is exact away from the optimum too", {
  # Profiling the protected block out is an identity, not a first-order
  # approximation, so it must hold at an arbitrary candidate vector. A
  # reduction that only agreed at the optimum would be a coincidence of two
  # correct solvers rather than a correct reduction.
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  arbitrary <- c(c1 = 0.3, c2 = -1.1, c3 = 2)

  y_star <- .fb_sel_reduction_response(red)
  fitted <- .fb_sel_reduction_fitted(red, arbitrary)
  expect_equal(
    .fb_sel_reduction_objective(red, arbitrary),
    0.5 * sum((y_star - fitted)^2),
    tolerance = 1e-8
  )
})


test_that("oracle 5 falsifier: skipping the candidate residualisation is red", {
  # The defect that survives casual review: residualise the response, forget
  # the candidate block. The coefficients are then wrong by the protected
  # block's projection, and nothing about the shape of the answer changes.
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  direct <- fbsd_gls(case$y, case$X0, case$Xp, case$V, lambda = 0)

  Xp_unresidualised <- .fb_sel_reduction_candidates(red)
  beta_wrong <- drop(solve(
    crossprod(Xp_unresidualised),
    crossprod(Xp_unresidualised, red$y_tilde)
  ))
  expect_false(isTRUE(all.equal(
    unname(beta_wrong), unname(direct$beta_p),
    tolerance = 1e-6
  )))
})


test_that("oracle 5 falsifier: ignoring the whitener is red", {
  # Running the same reduction on the unwhitened design is the other silent
  # failure: it is ordinary least squares wearing the name of GLS.
  case <- fbsd_case()
  red_ok <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  red_iid <- .fb_sel_protected_reduction(case$y, case$X0, case$Xp)
  direct <- fbsd_gls(case$y, case$X0, case$Xp, case$V, lambda = 0)

  solve_for <- function(red) {
    drop(solve(crossprod(red$Xp_tilde), crossprod(red$Xp_tilde, red$y_tilde)))
  }
  expect_equal(unname(solve_for(red_ok)), unname(direct$beta_p),
    tolerance = 1e-8
  )
  expect_false(isTRUE(all.equal(
    unname(solve_for(red_iid)), unname(direct$beta_p),
    tolerance = 1e-4
  )))
})


# --- ORACLE 5 (ii): rank deficiency is explicit ------------------------ #

test_that("a rank-deficient protected block refuses, naming the columns", {
  case <- fbsd_case()
  X0 <- cbind(case$X0, z1_again = case$X0[, "z1"])

  expect_error(
    .fb_sel_protected_reduction(
      case$y, X0, case$Xp,
      whitener = case$whitener
    ),
    "rank deficient after whitening"
  )
  # The reason code, not just any error: an unregistered code surfaces as a
  # contract violation, on which a bare expect_error() would also pass.
  expect_error(
    .fb_sel_protected_reduction(
      case$y, X0, case$Xp,
      whitener = case$whitener
    ),
    class = "flexybayes_refusal_selection_rank_deficient_protected"
  )
  # And it names both halves: what is estimable, and what is not.
  err <- tryCatch(
    .fb_sel_protected_reduction(
      case$y, X0, case$Xp,
      whitener = case$whitener
    ),
    error = function(e) e
  )
  expect_match(conditionMessage(err), "z1_again")
  expect_match(conditionMessage(err), "\\(Intercept\\), z1, z2")
  expect_identical(err$aliased, "z1_again")
  expect_identical(err$rank, 3L)
})


test_that("the estimable space is reported when the caller asks for it", {
  case <- fbsd_case()
  X0 <- cbind(case$X0, z1_again = case$X0[, "z1"])
  red <- .fb_sel_protected_reduction(
    case$y, X0, case$Xp,
    whitener = case$whitener, on_rank_deficient = "report"
  )

  expect_true(red$rank_deficient)
  expect_identical(red$protected_rank, 3L)
  expect_identical(red$protected_estimable, c("(Intercept)", "z1", "z2"))
  expect_identical(red$protected_aliased, "z1_again")
  # Not a silent drop: the aliased coefficient comes back NA, the way lm()
  # reports one, rather than a zero that would read as an estimate.
  beta0 <- .fb_sel_recover_protected(red, c(c1 = 1, c2 = 0, c3 = -1))
  expect_true(is.na(beta0[["z1_again"]]))
  expect_false(anyNA(beta0[c("(Intercept)", "z1", "z2")]))
})


test_that("rank-deficiency falsifier: a full-rank block must not refuse", {
  # The gate has to be able to say no. A rank check wired to a constant
  # would pass the refusal test above and fail here.
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  expect_false(red$rank_deficient)
  expect_identical(red$protected_aliased, character(0))
})


# --- ORACLE 5 (iv): a duplicate protected column changes nothing ------- #

test_that("a duplicate protected column does not change predictions", {
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  red_dup <- .fb_sel_protected_reduction(
    case$y, cbind(case$X0, z1_again = case$X0[, "z1"]), case$Xp,
    whitener = case$whitener, on_rank_deficient = "report"
  )

  beta_p <- c(c1 = 0.7, c2 = -0.2, c3 = 1.3)
  expect_equal(red_dup$y_tilde, red$y_tilde, tolerance = 1e-10)
  expect_equal(
    unname(red_dup$Xp_tilde), unname(red$Xp_tilde),
    tolerance = 1e-10
  )
  expect_equal(
    .fb_sel_reduction_fitted(red_dup, beta_p),
    .fb_sel_reduction_fitted(red, beta_p),
    tolerance = 1e-10
  )
  expect_equal(
    .fb_sel_reduction_objective(red_dup, beta_p),
    .fb_sel_reduction_objective(red, beta_p),
    tolerance = 1e-10
  )
})


# --- ORACLE 5 (iii): the provider sees no protected column ------------- #

test_that("the provider payload carries the candidate block and nothing else", {
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  payload <- .fb_sel_provider_payload(red)

  expect_identical(colnames(payload$x), colnames(case$Xp))
  expect_identical(ncol(payload$x), ncol(case$Xp))
  expect_length(intersect(colnames(payload$x), colnames(case$X0)), 0L)
  # The intercept is disabled at the provider: it is a protected column and
  # the residualisation has already removed its direction.
  expect_false(payload$intercept)
  # Asserted on the numbers, not only on the names. Every candidate column
  # handed over is orthogonal to the protected column space, so no protected
  # direction can be hiding in the payload under another name.
  expect_lt(max(abs(crossprod(red$Q0, payload$x))), 1e-10)
  # And nothing else of the protected block travels with it.
  expect_setequal(names(payload), c("x", "y", "intercept", "n", "columns"))
})


test_that("a protected column smuggled in under a new name is annihilated", {
  # The name check alone would let this through: the same numbers, a
  # different label. The reduction removes it because the reduction works on
  # the column space, not on the names.
  case <- fbsd_case()
  Xp <- cbind(case$Xp, z1_in_disguise = case$X0[, "z1"])
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, Xp,
    whitener = case$whitener
  )
  payload <- .fb_sel_provider_payload(red)
  expect_lt(max(abs(payload$x[, "z1_in_disguise"])), 1e-10)
})


test_that("payload falsifier: a leaked protected column is caught on the way out", {
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  # Simulate the defect the guard exists for: a later convenience widens the
  # candidate matrix with a protected column.
  red$Xp_tilde <- cbind(red$Xp_tilde, z1 = red$X0_star[, "z1"])
  expect_error(
    .fb_sel_provider_payload(red),
    "protected column\\(s\\) reached the provider payload"
  )
})


test_that("a column cannot be protected and a candidate at once", {
  case <- fbsd_case()
  expect_error(
    .fb_sel_protected_reduction(
      case$y, case$X0, cbind(case$Xp, z1 = case$X0[, "z1"]),
      whitener = case$whitener
    ),
    "a column may be protected or a candidate, not both: z1"
  )
})


# --- ORACLE 5 (v): the offset is subtracted exactly once --------------- #

test_that("an offset is subtracted once, before the reduction", {
  case <- fbsd_case()
  offset <- 0.5 + 0.25 * seq_len(case$n)

  with_offset <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener, offset = offset
  )
  pre_subtracted <- .fb_sel_protected_reduction(
    case$y - offset, case$X0, case$Xp,
    whitener = case$whitener
  )

  expect_equal(with_offset$y_tilde, pre_subtracted$y_tilde, tolerance = 1e-12)
  expect_equal(with_offset$Xp_tilde, pre_subtracted$Xp_tilde, tolerance = 1e-12)
  expect_true(with_offset$offset_applied)
  expect_identical(with_offset$offset, offset)
  expect_false(pre_subtracted$offset_applied)
  expect_null(pre_subtracted$offset)

  # A zero offset is the identity, so the argument cannot be doing something
  # on its own account.
  zero <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener, offset = rep(0, case$n)
  )
  no_offset <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  expect_equal(zero$y_tilde, no_offset$y_tilde, tolerance = 1e-12)
})


test_that("offset falsifier: subtracting twice, or after whitening, is red", {
  case <- fbsd_case()
  offset <- 0.5 + 0.25 * seq_len(case$n)
  once <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener, offset = offset
  )

  # Twice. The negative control section 5.5 asks for, in the offset's own
  # terms.
  twice <- .fb_sel_protected_reduction(
    case$y - offset, case$X0, case$Xp,
    whitener = case$whitener, offset = offset
  )
  expect_false(isTRUE(all.equal(
    twice$y_tilde, once$y_tilde,
    tolerance = 1e-6
  )))

  # After whitening rather than before. Both operations are linear, which is
  # exactly why this defect looks harmless: the error is that the offset is
  # not whitened with the response, so W y - o is compared with W(y - o).
  y_star <- drop(.fb_sel_whiten(case$whitener, matrix(case$y, ncol = 1L)))
  late <- .fb_sel_protected_reduction(y_star - offset, case$X0, case$Xp)
  expect_false(isTRUE(all.equal(
    late$y_tilde, once$y_tilde,
    tolerance = 1e-6
  )))
})


test_that("a malformed offset refuses rather than recycling", {
  case <- fbsd_case()
  expect_error(
    .fb_sel_protected_reduction(
      case$y, case$X0, case$Xp,
      whitener = case$whitener, offset = c(1, 2)
    ),
    "must have one value per observation"
  )
  expect_error(
    .fb_sel_protected_reduction(
      case$y, case$X0, case$Xp,
      whitener = case$whitener, offset = c(NA, rep(0, case$n - 1L))
    ),
    "must be finite at every observation"
  )
})


# --- the reduction's own arithmetic ------------------------------------ #

test_that("y* and Xp* are recoverable from what the object holds", {
  # The object stores the residualised pair plus two small projections
  # instead of a second copy of the whitened design, so the reconstruction is
  # part of the contract rather than a convenience.
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  expect_equal(
    .fb_sel_reduction_response(red),
    drop(.fb_sel_whiten(case$whitener, matrix(case$y, ncol = 1L))),
    tolerance = 1e-10
  )
  expect_equal(
    .fb_sel_reduction_candidates(red),
    .fb_sel_whiten(case$whitener, case$Xp),
    tolerance = 1e-10
  )
})


test_that("a candidate solution of the wrong length refuses", {
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, case$X0, case$Xp,
    whitener = case$whitener
  )
  expect_error(
    .fb_sel_recover_protected(red, c(1, 2)),
    "one coefficient per candidate column"
  )
})


test_that("an empty protected block reduces to the identity", {
  case <- fbsd_case()
  red <- .fb_sel_protected_reduction(
    case$y, NULL, case$Xp,
    whitener = case$whitener
  )
  expect_identical(red$protected_rank, 0L)
  expect_equal(
    red$y_tilde,
    drop(.fb_sel_whiten(case$whitener, matrix(case$y, ncol = 1L))),
    tolerance = 1e-12
  )
  expect_length(.fb_sel_recover_protected(red, c(1, 0, 0)), 0L)
})


# =============================================================================
# ORACLE 6a --- term resolution (contract sections 5.3 and 6)
# =============================================================================

test_that("the resolver groups columns by term, not by column", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)

  expect_identical(map$term_labels, c("x2", "f", "h", "f:x1"))
  # One numeric main effect: one term, one column.
  expect_length(map$term_map[["x2"]], 1L)
  # A four-level factor: one term, all three of its contrast columns.
  expect_length(map$term_map[["f"]], 3L)
  expect_identical(
    colnames(map$X)[map$term_map[["f"]]], c("fb", "fc", "fd")
  )
  # A three-level factor: one term, both of its columns.
  expect_length(map$term_map[["h"]], 2L)
  # The interaction: every column it generates, which is four here because
  # the factor's main effect is already in the model.
  expect_length(map$term_map[["f:x1"]], 4L)
  # The intercept is assign 0 and belongs to no candidate term.
  expect_true(map$has_intercept)
  expect_identical(map$intercept_column, 1L)
  expect_false(1L %in% unlist(map$term_map, use.names = FALSE))
  # Every non-intercept column is claimed exactly once.
  claimed <- sort(unlist(map$term_map, use.names = FALSE))
  expect_identical(claimed, setdiff(seq_len(ncol(map$X)), 1L))

  expect_identical(
    map$terms_table$type, c("numeric", "factor", "factor", "interaction")
  )
  expect_identical(
    map$terms_table$n_columns, c(1L, 3L, 2L, 4L)
  )
  # The group map is the term map: the unit of selection is the term.
  expect_identical(map$group_map, map$term_map)
})


test_that("resolver falsifier: column names cannot do the assignment map's job", {
  # A numeric called `fb` and a factor `f` with a level `b` both produce a
  # column called "fb". A resolver that parsed names would put two different
  # terms in one group, or one term in two. The assign attribute does not
  # have the problem, and this is the case that proves the shortcut fails.
  set.seed(517L)
  d <- data.frame(
    fb = rnorm(40),
    f = factor(sample(c("a", "b"), 40, TRUE))
  )
  map <- .fb_sel_term_map(~ fb + f, d)

  expect_true(anyDuplicated(colnames(map$X)) > 0L)
  expect_identical(colnames(map$X), c("(Intercept)", "fb", "fb"))
  # Two distinct terms, one column each, correctly separated.
  expect_identical(map$term_map[["fb"]], 2L)
  expect_identical(map$term_map[["f"]], 3L)
  # The naive alternative, and the damage it does.
  naive <- which(colnames(map$X) == "fb")
  expect_length(naive, 2L)
  expect_false(identical(naive, map$term_map[["fb"]]))
})


test_that("the resolver keeps row alignment with the data it was given", {
  # na.pass, deliberately. Dropping incomplete rows here would renumber every
  # row downstream and silently misalign the response and the folds.
  case <- fbsd_terms_case()
  d <- case$data
  d$x2[c(3L, 17L)] <- NA_real_
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, d)
  expect_identical(nrow(map$X), nrow(d))
  expect_true(all(is.na(map$X[c(3L, 17L), "x2"])))
})


test_that("the partition puts the line where the caller asks, intercept aside", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)

  none <- .fb_sel_partition_design(map)
  expect_identical(colnames(none$X0), "(Intercept)")
  expect_identical(none$candidate_terms, c("x2", "f", "h", "f:x1"))
  expect_identical(ncol(none$Xp), ncol(map$X) - 1L)

  # Protecting a term from the middle of the formula is the case the index
  # arithmetic has to survive.
  some <- .fb_sel_partition_design(map, protected = "f")
  expect_identical(ncol(some$X0), 4L)
  expect_identical(some$candidate_terms, c("x2", "h", "f:x1"))
  expect_identical(names(some$group_map), c("x2", "h", "f:x1"))
  # The renumbered group map must still point at the same columns.
  for (g in names(some$group_map)) {
    expect_identical(
      colnames(some$Xp)[some$group_map[[g]]],
      colnames(map$X)[map$term_map[[g]]]
    )
  }
  expect_error(
    .fb_sel_partition_design(map, protected = "not_a_term"),
    "these protected terms are not in the model: not_a_term"
  )
})


# =============================================================================
# ORACLE 6b --- block transforms and contrast invariance (section 5.3)
# =============================================================================

test_that("the block transform reproduces the block it came from", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  part <- .fb_sel_partition_design(map)
  red <- .fb_sel_protected_reduction(
    case$data$y, part$X0, part$Xp,
    whitener = case$whitener
  )
  bt <- .fb_sel_block_transforms(red$Xp_tilde, part$group_map)

  for (g in bt$group_names) {
    b <- bt$blocks[[g]]
    expect_equal(
      unname(red$Xp_tilde[, b$basis_columns, drop = FALSE]),
      unname(b$Q %*% b$R),
      tolerance = 1e-10
    )
    expect_equal(crossprod(b$Q), diag(1, b$rank), tolerance = 1e-10)
  }
})


test_that("group weights are sqrt(rank), never sqrt(column count)", {
  # Section 5.1: a block with a collinear column has fewer directions than
  # columns, and taxing it for the missing degree of freedom would shrink it
  # harder than a block that really carries that direction.
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  part <- .fb_sel_partition_design(map)
  red <- .fb_sel_protected_reduction(
    case$data$y, part$X0, part$Xp,
    whitener = case$whitener
  )

  bt <- .fb_sel_block_transforms(red$Xp_tilde, part$group_map)
  expect_equal(
    unname(vapply(bt$blocks, function(b) b$rank, integer(1L))),
    c(1L, 3L, 2L, 4L)
  )
  expect_equal(
    unname(.fb_sel_block_design(bt)$weights),
    sqrt(c(1, 3, 2, 4)),
    tolerance = 1e-12
  )

  # Now make one column of `f` a copy of another. The block has four columns
  # and rank three.
  Xp_dup <- cbind(part$Xp, f_extra = part$Xp[, part$group_map[["f"]][[1L]]])
  groups <- part$group_map
  groups[["f"]] <- c(groups[["f"]], ncol(Xp_dup))
  red_dup <- .fb_sel_protected_reduction(
    case$data$y, part$X0, Xp_dup,
    whitener = case$whitener
  )
  bt_dup <- .fb_sel_block_transforms(red_dup$Xp_tilde, groups)

  expect_identical(bt_dup$blocks[["f"]]$n_columns, 4L)
  expect_identical(bt_dup$blocks[["f"]]$rank, 3L)
  expect_equal(bt_dup$blocks[["f"]]$weight, sqrt(3), tolerance = 1e-12)
  # The falsifier for this gate: sqrt(p) would be 2 here, sqrt(rank) is
  # sqrt(3). They must not coincide.
  expect_false(isTRUE(all.equal(
    bt_dup$blocks[["f"]]$weight, sqrt(4),
    tolerance = 1e-6
  )))
  expect_identical(bt_dup$blocks[["f"]]$aliased_names, "f_extra")
})


test_that("the back transform returns beta in the estimable parameterisation", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  part <- .fb_sel_partition_design(map)
  red <- .fb_sel_protected_reduction(
    case$data$y, part$X0, part$Xp,
    whitener = case$whitener
  )
  bt <- .fb_sel_block_transforms(red$Xp_tilde, part$group_map)

  set.seed(31L)
  gamma <- rnorm(sum(vapply(bt$blocks, function(b) b$rank, integer(1L))))
  beta <- .fb_sel_back_transform(bt, gamma)

  expect_identical(names(beta), colnames(red$Xp_tilde))
  expect_false(anyNA(beta))
  # X beta and Q gamma are the same fitted contribution.
  expect_equal(
    drop(red$Xp_tilde %*% matrix(beta, ncol = 1L)),
    .fb_sel_block_fitted(bt, gamma),
    tolerance = 1e-9
  )
})


test_that("an aliased column comes back NA, not zero", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  part <- .fb_sel_partition_design(map)
  Xp_dup <- cbind(part$Xp, f_extra = part$Xp[, part$group_map[["f"]][[1L]]])
  groups <- part$group_map
  groups[["f"]] <- c(groups[["f"]], ncol(Xp_dup))
  red <- .fb_sel_protected_reduction(
    case$data$y, part$X0, Xp_dup,
    whitener = case$whitener
  )
  bt <- .fb_sel_block_transforms(red$Xp_tilde, groups)
  gamma <- rep(1, sum(vapply(bt$blocks, function(b) b$rank, integer(1L))))
  beta <- .fb_sel_back_transform(bt, gamma)

  expect_true(is.na(beta[["f_extra"]]))
  expect_false(anyNA(beta[setdiff(names(beta), "f_extra")]))
})


test_that("a column in two groups refuses", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  part <- .fb_sel_partition_design(map)
  groups <- part$group_map
  groups[["h"]] <- c(groups[["h"]], groups[["x2"]])
  expect_error(
    .fb_sel_block_transforms(part$Xp, groups),
    "a column may belong to one group only"
  )
})


# --- ORACLE 6: the invariance gate ------------------------------------- #

test_that("oracle 6: fit, risk and support survive every reparameterisation", {
  case <- fbsd_terms_case()
  base <- fbsd_pipeline(case)

  # The gate is only worth anything if the support is non-trivial: at this
  # lambda three of the four terms are in and one is out.
  expect_identical(base$support, c("x2", "f", "f:x1"))
  expect_false(base$norms[["h"]] > 1e-9)

  variants <- list(
    sum = c("contr.sum", "contr.poly"),
    helmert = c("contr.helmert", "contr.poly"),
    sas = c("contr.SAS", "contr.poly")
  )
  for (nm in names(variants)) {
    v <- fbsd_pipeline(case, contrasts_opt = variants[[nm]])
    expect_equal(v$fitted, base$fitted, tolerance = 1e-8, info = nm)
    expect_equal(v$risk, base$risk, tolerance = 1e-8, info = nm)
    expect_identical(v$support, base$support, info = nm)
    expect_equal(unname(v$norms), unname(base$norms),
      tolerance = 1e-8, info = nm
    )
  }

  # A different reference level. The model has not changed; only the labels
  # of the columns have.
  relevelled <- case$data
  relevelled$f <- factor(relevelled$f, levels = c("c", "a", "d", "b"))
  relevelled$h <- factor(relevelled$h, levels = c("mid", "hi", "lo"))
  v <- fbsd_pipeline(case, data = relevelled)
  expect_equal(v$fitted, base$fitted, tolerance = 1e-8)
  expect_equal(v$risk, base$risk, tolerance = 1e-8)
  expect_identical(v$support, base$support)

  # A random well-conditioned full-rank A on the factor block, which no
  # contrast argument could produce. This is the clause that separates
  # invariance from "the four contrast schemes we happened to try".
  set.seed(9109L)
  A <- matrix(rnorm(9), 3L, 3L) + 3 * diag(3L)
  expect_lt(kappa(A, exact = TRUE), 10)
  v <- fbsd_pipeline(case, transform = A, term = "f")
  expect_equal(v$fitted, base$fitted, tolerance = 1e-8)
  expect_equal(v$risk, base$risk, tolerance = 1e-8)
  expect_identical(v$support, base$support)

  # After back-transformation the coefficients are related by A, exactly:
  # the block was reparameterised, so beta was too.
  cols <- base$part$group_map[["f"]]
  expect_equal(
    drop(A %*% matrix(v$beta[cols], ncol = 1L)),
    unname(base$beta[cols]),
    tolerance = 1e-8
  )

  # And on the interaction block as well, which has four columns and is the
  # one a term resolver is most likely to mis-group.
  set.seed(2207L)
  A4 <- matrix(rnorm(16), 4L, 4L) + 4 * diag(4L)
  expect_lt(kappa(A4, exact = TRUE), 10)
  v4 <- fbsd_pipeline(case, transform = A4, term = "f:x1")
  expect_equal(v4$fitted, base$fitted, tolerance = 1e-8)
  expect_equal(v4$risk, base$risk, tolerance = 1e-8)
  expect_identical(v4$support, base$support)
})


test_that("oracle 6 falsifier: column standardisation loses the invariance", {
  # The plausible wrong answer, and the one the contract names: scale each
  # column to unit norm and penalise the resulting coefficient vector. It is
  # invariant to a diagonal change of basis and to nothing else, so a
  # contrast map moves the fit, the risk and, at the wrong lambda, the
  # support.
  case <- fbsd_terms_case()
  bad_pipeline <- function(contrasts_opt = NULL, transform = NULL,
                           term = NULL, lambda = 2) {
    if (!is.null(contrasts_opt)) {
      old <- options(contrasts = contrasts_opt)
      on.exit(options(old))
    }
    d <- case$data
    map <- .fb_sel_term_map(~ x2 + f + h + x1:f, d)
    part <- .fb_sel_partition_design(map)
    Xp <- part$Xp
    if (!is.null(transform)) {
      cols <- part$group_map[[term]]
      Xp[, cols] <- Xp[, cols, drop = FALSE] %*% transform
    }
    red <- .fb_sel_protected_reduction(
      d$y, part$X0, Xp,
      whitener = case$whitener
    )
    Qs <- lapply(part$group_map, function(cols) {
      block <- red$Xp_tilde[, cols, drop = FALSE]
      sweep(block, 2L, sqrt(colSums(block^2)), "/")
    })
    gl <- fbsd_group_lasso(
      red$y_tilde, Qs,
      sqrt(lengths(part$group_map)),
      lambda = lambda
    )
    norms <- vapply(gl$gamma, function(g) sqrt(sum(g^2)), numeric(1L))
    list(
      fitted = gl$fitted,
      risk = sqrt(mean((red$y_tilde - gl$fitted)^2)),
      support = names(norms)[norms > 1e-9]
    )
  }

  base_bad <- bad_pipeline()
  set.seed(9109L)
  A <- matrix(rnorm(9), 3L, 3L) + 3 * diag(3L)
  moved <- bad_pipeline(transform = A, term = "f")
  helmert <- bad_pipeline(contrasts_opt = c("contr.helmert", "contr.poly"))

  expect_false(isTRUE(all.equal(
    moved$fitted, base_bad$fitted,
    tolerance = 1e-6
  )))
  expect_false(isTRUE(all.equal(
    moved$risk, base_bad$risk,
    tolerance = 1e-6
  )))
  expect_false(isTRUE(all.equal(
    helmert$fitted, base_bad$fitted,
    tolerance = 1e-6
  )))

  # And the same comparison through the block transform is invariant, so the
  # difference is the transform and not the fixture.
  ok_base <- fbsd_pipeline(case)
  ok_moved <- fbsd_pipeline(case, transform = A, term = "f")
  expect_equal(ok_moved$fitted, ok_base$fitted, tolerance = 1e-8)
})


test_that("oracle 6 falsifier: a mis-grouped block loses the invariance", {
  # Splitting a factor's contrast columns into separate groups is the
  # column-not-term defect the contract's section 6 is about. It gives a
  # different answer under a relevel that changes no model.
  case <- fbsd_terms_case()
  split_pipeline <- function(data) {
    map <- .fb_sel_term_map(~ x2 + f + h + x1:f, data)
    part <- .fb_sel_partition_design(map)
    groups <- part$group_map
    f_cols <- groups[["f"]]
    groups[["f"]] <- NULL
    for (i in seq_along(f_cols)) {
      groups[[paste0("f_col", i)]] <- f_cols[[i]]
    }
    red <- .fb_sel_protected_reduction(
      data$y, part$X0, part$Xp,
      whitener = case$whitener
    )
    bt <- .fb_sel_block_transforms(red$Xp_tilde, groups)
    des <- .fb_sel_block_design(bt)
    gl <- fbsd_group_lasso(
      red$y_tilde, lapply(bt$blocks, function(b) b$Q), des$weights,
      lambda = 2
    )
    gl$fitted
  }
  relevelled <- case$data
  relevelled$f <- factor(relevelled$f, levels = c("c", "a", "d", "b"))

  expect_false(isTRUE(all.equal(
    split_pipeline(relevelled), split_pipeline(case$data),
    tolerance = 1e-6
  )))
  # The term-grouped pipeline is invariant to the same relevel.
  expect_equal(
    fbsd_pipeline(case, data = relevelled)$fitted,
    fbsd_pipeline(case)$fitted,
    tolerance = 1e-8
  )
})


# =============================================================================
# Coordinate semantics --- the normative restriction of section 5.3
# =============================================================================

test_that("every target carries a coordinate semantics", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)

  by_term <- .fb_sel_targets(map, "term")
  expect_identical(by_term$target, c("x2", "f", "h", "f:x1"))
  expect_identical(
    by_term$coordinate_semantics,
    c("intrinsic", "term_subspace", "term_subspace", "term_subspace")
  )
  expect_true(all(
    by_term$coordinate_semantics %in%
      c("intrinsic", "basis_dependent", "term_subspace")
  ))

  by_coef <- .fb_sel_targets(map, "coefficient")
  expect_identical(nrow(by_coef), ncol(map$X) - 1L)
  expect_identical(
    by_coef$coordinate_semantics[by_coef$target == "x2"], "intrinsic"
  )
  expect_true(all(
    by_coef$coordinate_semantics[by_coef$term == "f"] == "basis_dependent"
  ))
})


test_that("an L1 component over a basis-dependent target refuses by name", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  by_term <- .fb_sel_targets(map, "term")

  expect_error(
    .fb_sel_assert_coordinate_semantics(by_term, "lasso"),
    "basis dependent"
  )
  expect_error(
    .fb_sel_assert_coordinate_semantics(by_term, "lasso"),
    class = "flexybayes_refusal_selection_penalty_coordinate_basis_dependent"
  )
  err <- tryCatch(
    .fb_sel_assert_coordinate_semantics(by_term, c("group_lasso", "lasso")),
    error = function(e) e
  )
  # It names every offending target -- all three multi-column terms, `h`
  # included, because a three-level factor contributes two contrast columns
  # and an L1 sum can tell them apart -- and the component that offended.
  expect_setequal(err$targets, c("f", "h", "f:x1"))
  expect_identical(err$components, "lasso")
  expect_match(conditionMessage(err), "f, h, f:x1")
  # The single-column numeric main effect is not among them.
  expect_false("x2" %in% err$targets)
  expect_match(conditionMessage(err), "within-factor selection")

  # Sparse group and exclusive reach inside the block too, so they refuse on
  # the same grounds even though the target was resolved as a whole term.
  for (component in c("sparse_group_lasso", "exclusive_lasso",
                      "group_exclusive", "elastic_net", "adaptive")) {
    expect_error(
      .fb_sel_assert_coordinate_semantics(by_term, component),
      "basis dependent",
      info = component
    )
  }
})


test_that("a group-only penalty over the same targets is accepted", {
  # The falsifier for the refusal above: a gate that refused everything
  # would pass the test above and be useless. Group L2 and ridge are
  # invariant to the rotation the block transform leaves free, so they must
  # pass.
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  by_term <- .fb_sel_targets(map, "term")

  for (component in c("ridge", "group_lasso")) {
    out <- .fb_sel_assert_coordinate_semantics(by_term, component)
    expect_s3_class(out, "data.frame")
    expect_true(all(out$effective_semantics != "basis_dependent"))
  }

  # An L1 over single-column intrinsic targets is fine: there is no contrast
  # matrix to depend on.
  numeric_only <- .fb_sel_targets(
    .fb_sel_term_map(~ x1 + x2, case$data), "term"
  )
  expect_s3_class(
    .fb_sel_assert_coordinate_semantics(numeric_only, "lasso"),
    "data.frame"
  )
})


test_that("a declared canonical basis is the documented escape", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  declared <- .fb_sel_targets(map, "term", canonical_basis = TRUE)

  out <- .fb_sel_assert_coordinate_semantics(declared, "lasso")
  expect_s3_class(out, "data.frame")
  # The escape does not rewrite the typing: the targets are still recorded as
  # basis dependent, so a certificate can say the guarantee was set aside.
  expect_true(any(out$effective_semantics == "basis_dependent"))
  expect_true(all(out$canonical_basis))
})


test_that("an unknown penalty component is an internal contract violation", {
  case <- fbsd_terms_case()
  map <- .fb_sel_term_map(~ x2 + f + h + x1:f, case$data)
  expect_error(
    .fb_sel_assert_coordinate_semantics(
      .fb_sel_targets(map, "term"), "bridge_regression"
    ),
    "unknown penalty component\\(s\\) bridge_regression"
  )
})


# =============================================================================
# Registry
# =============================================================================

test_that("every refusal this module can throw is registered", {
  # An unregistered reason code raises a contract violation rather than the
  # refusal, so a bare expect_error() would pass on the wrong error. This is
  # the direct check that the codes above are real.
  codes <- c(
    "selection_rank_deficient_protected",
    "selection_penalty_coordinate_basis_dependent"
  )
  expect_true(all(codes %in% fb_refusals()$reason_code))
})
