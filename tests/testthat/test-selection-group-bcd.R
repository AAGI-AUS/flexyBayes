# Group lasso and group elastic net by block coordinate descent, and the
# relative KKT threshold of a path point.
#
# Oracle: grpreg 3.6.0 (Breheny and Huang), called one lambda at a time. On a
# design whose groups are centred and orthonormal (crossprod(X_g) / n = I) its
# objective, measured, is
#
#   (1 / 2n) ||y - X b||^2 + lambda sum_g m_g [alpha ||b_g|| + (1 - alpha) / 2 ||b_g||^2]
#
# with m_g its `group.multiplier`, so in this package's units lambda_group =
# n lambda alpha and lambda_l2 = n lambda (1 - alpha) m_g, which is one number
# only when every m_g is equal. With alpha < 1 grpreg leaves out some groups
# whose score lies between lambda alpha m_g and lambda m_g, so it is compared
# only at points where its own answer meets the stationarity condition to
# 1e-8; the blocks assert those points are the majority.
#
# Tolerance: over 20 seeds and 1,549 such points the largest coefficient
# difference measured was 5.9e-12; the blocks use 1e-10.

fbgb_orthonormal <- function(groups, n = 80L, seed = 1L) {
  .fb_sel_with_seed(seed, {
    p <- length(groups)
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    for (pass in seq_len(2L)) {
      x <- scale(x, center = TRUE, scale = FALSE)
      for (k in unique(groups)) {
        idx <- which(groups == k)
        x[, idx] <- qr.Q(qr(x[, idx, drop = FALSE])) * sqrt(n)
      }
    }
    x <- matrix(as.numeric(x), nrow = n, ncol = p)
    colnames(x) <- paste0("v", seq_len(p))
    truth <- stats::rnorm(p) * (groups %in% c(1, 3, 4))
    z <- drop(x %*% truth) + stats::rnorm(n)
    list(A = x, z = z - mean(z), n = n, groups = groups)
  })
}


fbgb_grouping <- function(design, weights = NULL) {
  out <- .fb_sel_group_index(
    colnames(design$A), paste0("t", design$groups), rule = "term"
  )
  if (!is.null(weights)) {
    out$weights[] <- weights
  }
  out
}


# grpreg at one lambda, and its answer's stationarity residual under the
# measured convention, per observation.
fbgb_grpreg_at <- function(design, grouping, lambda, alpha, multiplier) {
  fit <- grpreg::grpreg(
    design$A, design$z,
    group = design$groups, penalty = "grLasso", alpha = alpha,
    lambda = lambda, group.multiplier = multiplier,
    eps = 1e-12, max.iter = 1e7
  )
  beta <- fit$beta[-1L, 1L]
  penalty <- list(
    lambda_l1 = 0,
    lambda_l2 = design$n * lambda * (1 - alpha) * multiplier[[1L]],
    lambda_group = design$n * lambda * alpha,
    lambda_exclusive = 0
  )
  list(
    beta = unname(beta),
    penalty = penalty,
    kkt = .fb_sel_kkt_residual_grouped(
      design$A, design$z, beta, penalty, grouping
    ) / design$n
  )
}


fbgb_compare <- function(design, grouping, alpha, multiplier) {
  lambdas <- grpreg::grpreg(
    design$A, design$z,
    group = design$groups, penalty = "grLasso", alpha = alpha,
    group.multiplier = multiplier
  )$lambda[seq(2L, 100L, by = 7L)]
  rows <- lapply(lambdas, function(lambda) {
    theirs <- fbgb_grpreg_at(design, grouping, lambda, alpha, multiplier)
    mine <- .fb_sel_bcd_grouped(
      design$A, design$z, theirs$penalty, grouping,
      tol = 1e-12 * max(abs(crossprod(design$A, design$z)))
    )
    data.frame(
      lambda = lambda,
      grpreg_kkt = theirs$kkt,
      difference = max(abs(unname(mine$beta) - theirs$beta)),
      n_zero_groups = sum(tapply(theirs$beta == 0, design$groups, all))
    )
  })
  do.call(rbind, rows)
}


test_that("the path KKT threshold is 1e-6 of max|A'z| and scales with the payload", {
  design <- fbgb_orthonormal(c(1, 1, 2, 2, 2))
  expected <- 1e-6 * max(abs(crossprod(design$A, design$z)))
  expect_identical(.fb_sel_path_kkt_tolerance(design$A, design$z), expected)
  expect_equal(
    .fb_sel_path_kkt_tolerance(design$A, 1000 * design$z),
    1000 * expected
  )
  expect_identical(
    .fb_sel_path_kkt_tolerance(design$A[, 0L, drop = FALSE], design$z),
    0
  )
  # Ridge, laplace and elastic net keep the absolute threshold at any scale.
  expect_identical(
    .fb_sel_path_kkt_tolerance(design$A, 1000 * design$z, grouped = FALSE),
    1e-6
  )
})


test_that("the block update solves its subproblem exactly", {
  .fb_sel_with_seed(5L, {
    for (trial in seq_len(50L)) {
      k <- sample(2:8, 1L)
      root <- matrix(stats::rnorm(k * k), k)
      hessian <- crossprod(root) + diag(stats::runif(1L, 1e-3, 1), k)
      r <- stats::rnorm(k, sd = 3)
      threshold <- stats::runif(1L, 0, 1.2) * sqrt(sum(r^2))
      decomposition <- eigen(hessian, symmetric = TRUE)
      b <- .fb_sel_group_block_update(
        r, decomposition$values, decomposition$vectors, threshold
      )
      radius <- sqrt(sum(b^2))
      if (sqrt(sum(r^2)) <= threshold) {
        expect_identical(b, numeric(k))
        next
      }
      # Stationarity: H b - r + threshold b / ||b|| = 0.
      stationarity <- drop(hessian %*% b) - r + threshold * b / radius
      expect_lt(max(abs(stationarity)), 1e-10 * max(1, sqrt(sum(r^2))))
    }
  })
})


test_that("block coordinate descent reproduces grpreg's group lasso and group elastic net", {
  skip_if_not_installed("grpreg")
  # Equal group sizes, so grpreg's ridge multiplier is one number.
  design <- fbgb_orthonormal(rep(1:5, each = 4L), seed = 3L)
  grouping <- fbgb_grouping(design)
  for (alpha in c(1, 0.6, 0.3)) {
    table <- fbgb_compare(design, grouping, alpha, rep(2, 5L))
    certified <- table$grpreg_kkt <= 1e-8
    expect_gt(mean(certified), 0.5, label = paste("alpha", alpha))
    expect_lt(max(table$difference[certified]), 1e-10)
    # The comparison covers points with groups both in and out.
    expect_true(any(table$n_zero_groups[certified] > 0))
    if (alpha == 1) {
      expect_true(all(certified))
    }
  }
})


test_that("unequal groups at unit weights reproduce grpreg too", {
  skip_if_not_installed("grpreg")
  groups <- rep(1:5, times = c(2L, 3L, 4L, 5L, 8L))
  design <- fbgb_orthonormal(groups, seed = 4L)
  grouping <- fbgb_grouping(design, weights = 1)
  for (alpha in c(1, 0.5)) {
    table <- fbgb_compare(design, grouping, alpha, rep(1, 5L))
    certified <- table$grpreg_kkt <= 1e-8
    expect_gt(mean(certified), 0.5, label = paste("alpha", alpha))
    expect_lt(max(table$difference[certified]), 1e-10)
  }
})


test_that("the ridge part equals appended rows, checked against grpreg's group lasso", {
  skip_if_not_installed("grpreg")
  # The package's own group elastic net: unequal groups at sqrt(size)
  # weights and one unweighted ridge. Appending rows sqrt(l2) M, with
  # M'M = I and columns of M summing to zero, turns it into a group lasso
  # that grpreg (alpha = 1) fits on a design still centred and orthonormal
  # within groups, at (n + l2) instead of n.
  groups <- rep(1:5, times = c(2L, 3L, 4L, 5L, 8L))
  design <- fbgb_orthonormal(groups, seed = 6L)
  grouping <- fbgb_grouping(design)
  n <- design$n
  p <- ncol(design$A)
  basis <- qr.Q(qr(cbind(1, diag(p + 1L)[, seq_len(p)])))[, -1L]
  for (l2 in c(5, 40)) {
    rows <- p + 1L
    augmented <- rbind(design$A, sqrt(l2) * basis)
    z_augmented <- c(design$z, numeric(rows))
    total <- n + rows
    # grpreg divides each group by sqrt((n + l2) / total) before penalising.
    shrink <- sqrt((n + l2) / total)
    for (lambda in c(0.4, 0.15, 0.05)) {
      fit <- grpreg::grpreg(
        augmented, z_augmented,
        group = groups, penalty = "grLasso", lambda = lambda,
        eps = 1e-12, max.iter = 1e7
      )
      theirs <- unname(fit$beta[-1L, 1L])
      penalty <- list(
        lambda_l1 = 0,
        lambda_l2 = l2,
        lambda_group = total * lambda * shrink,
        lambda_exclusive = 0
      )
      mine <- .fb_sel_bcd_grouped(
        design$A, design$z, penalty, grouping,
        tol = 1e-12 * max(abs(crossprod(design$A, design$z)))
      )
      expect_lt(max(abs(unname(mine$beta) - theirs)), 1e-10)
      # And the package's two routes to the same point: the ridge inside the
      # block Hessian, and a group lasso on the appended rows.
      appended <- .fb_sel_bcd_grouped(
        augmented, z_augmented,
        list(lambda_l1 = 0, lambda_l2 = 0,
             lambda_group = penalty$lambda_group, lambda_exclusive = 0),
        grouping,
        tol = 1e-12 * max(abs(crossprod(augmented, z_augmented)))
      )
      expect_lt(max(abs(appended$beta - mine$beta)), 1e-10)
    }
  }
})


test_that("block coordinate descent and the proximal route agree on correlated blocks", {
  design <- .fb_sel_with_seed(7L, {
    n <- 120L
    latent <- stats::rnorm(n)
    a <- cbind(
      sapply(seq_len(4L), function(j) latent + stats::rnorm(n, sd = 0.1)),
      sapply(seq_len(3L), function(j) stats::rnorm(n))
    )
    colnames(a) <- paste0("a", seq_len(7L))
    z <- drop(a %*% c(0.5, 0, 0, 0, 0.3, 0, 0)) + stats::rnorm(n)
    list(A = a, z = z)
  })
  grouping <- .fb_sel_group_index(
    colnames(design$A), c(rep("g1", 4L), rep("g2", 3L))
  )
  for (l2 in c(0, 3)) {
    penalty <- list(lambda_l1 = 0, lambda_l2 = l2, lambda_group = 20,
                    lambda_exclusive = 0)
    bcd <- .fb_sel_bcd_grouped(design$A, design$z, penalty, grouping, tol = 1e-11)
    prox <- .fb_sel_prox_mode_grouped(
      design$A, design$z, penalty, grouping, tol = 1e-11
    )
    expect_lt(max(abs(bcd$beta - prox$beta)), 1e-10)
    expect_lt(
      .fb_sel_kkt_residual_grouped(design$A, design$z, bcd$beta, penalty, grouping),
      1e-11
    )
  }
  # Sparse group and exclusive penalties are not block problems.
  expect_error(
    .fb_sel_bcd_grouped(
      design$A, design$z,
      list(lambda_l1 = 1, lambda_group = 20), grouping
    ),
    "internal contract violation"
  )
  expect_error(
    .fb_sel_bcd_grouped(design$A, design$z, list(lambda_group = 20), grouping,
      tol = 1e-14, max_sweeps = 1L),
    class = "flexybayes_refusal_selection_prox_mode_not_converged"
  )
})


test_that("group_lasso() takes alpha, and alpha below one adds the ridge part", {
  default <- fb_penalty(term(everything()) ~ group_lasso(scale = 0.5))
  expect_identical(default$specs[[1L]]$args$alpha, 1)
  mixed <- fb_penalty(term(everything()) ~ group_lasso(scale = 0.5, alpha = 0.4))
  expect_identical(mixed$specs[[1L]]$args$alpha, 0.4)
  positional <- fb_penalty(term(everything()) ~ group_lasso(0.5, "block"))
  expect_identical(positional$specs[[1L]]$args$group, "block")
  expect_identical(positional$specs[[1L]]$args$alpha, 1)
  for (bad in list(0, 1.2, -0.1)) {
    expect_error(
      fb_penalty(term(everything()) ~ group_lasso(scale = 0.5, alpha = bad)),
      class = "flexybayes_refusal_penalty_hyperparameter_invalid"
    )
  }

  lasso <- .fb_sel_canonical_penalty(
    "group_lasso", list(scale = 0.5, group = "term", alpha = 1), n = 40
  )
  expect_identical(c(lasso$lambda_group, lasso$lambda_l2), c(2, 0))
  net <- .fb_sel_canonical_penalty(
    "group_lasso", list(scale = 0.5, group = "term", alpha = 0.4), n = 40
  )
  expect_equal(net$lambda_group, 0.4 / 0.5)
  expect_equal(net$lambda_l2, 0.6 / 0.5)
  expect_identical(net$lambda_l1, 0)
  expect_true(net$finite_anchor)
})


# --- the public verb -------------------------------------------------------- #

fbgb_small_study <- function() {
  fbgs_study(
    n_trial = 12L,
    sizes = c(2L, 3L, 4L),
    signal = c(g1_1 = 0.6, g3_2 = -0.4),
    seed = 11L
  )
}


fbgb_fit <- function(d, family_call, sizes = c(2L, 3L, 4L)) {
  fb_select(
    fbgs_base(d),
    candidates = fbgs_candidates(d),
    penalty = fbgs_penalty(sizes, family_call),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~trial, v = 3L),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 5L
  )
}


test_that("fb_select() fits a group elastic net whose groups enter and leave whole", {
  d <- fbgb_small_study()
  ratios <- numeric(0)
  original <- .fb_sel_bcd_grouped
  local_mocked_bindings(
    .fb_sel_bcd_grouped = function(A, z, penalty, groups, ...) {
      ratios <<- c(ratios, penalty$lambda_l2 / penalty$lambda_group)
      original(A, z, penalty, groups, ...)
    }
  )
  fit <- fbgb_fit(d, quote(group_lasso(scale = 1, group = "block", alpha = 0.5)))
  # alpha = 0.5 puts equal weight on the group norm and the ridge at every
  # point the solver receives.
  expect_gt(length(ratios), 0L)
  expect_equal(ratios, rep(1, length(ratios)))
  arm <- fit$benchmark$glmnet
  expect_identical(arm$family, "group_lasso")
  expect_gt(nrow(arm$kkt), 0L)
  expect_true(all(arm$kkt$kkt <= arm$kkt$kkt_tolerance))
  expect_true(all(arm$kkt$route == "internal_bcd"))
  history <- arm$support$group_history
  expect_gt(nrow(history), 0L)
  sizes <- c(block_1 = 2L, block_2 = 3L, block_3 = 4L)
  # Every group is held with all of its columns or with none.
  expect_true(all(
    history$n_in_support == 0L |
      history$n_in_support == sizes[history$group]
  ))
  expect_true(any(history$n_in_support > 0L))
  expect_true(any(history$n_in_support == 0L))
})


test_that("group lasso path points use block coordinate descent unless the option asks for prox", {
  d <- fbgb_small_study()
  call_count <- 0L
  original <- .fb_sel_bcd_grouped
  local_mocked_bindings(
    .fb_sel_bcd_grouped = function(...) {
      call_count <<- call_count + 1L
      original(...)
    }
  )
  default <- fbgb_fit(d, quote(group_lasso(scale = 1, group = "block")))
  expect_gt(call_count, 0L)
  expect_gt(nrow(default$benchmark$glmnet$kkt), 0L)
  expect_true(all(default$benchmark$glmnet$kkt$route == "internal_bcd"))

  call_count <- 0L
  withr::local_options(flexyBayes.selection_grouped_solver = "prox")
  prox <- fbgb_fit(d, quote(group_lasso(scale = 1, group = "block")))
  expect_identical(call_count, 0L)
  expect_gt(nrow(prox$benchmark$glmnet$kkt), 0L)
  expect_true(all(prox$benchmark$glmnet$kkt$route == "internal_prox"))
  # Both solvers meet the same threshold, and the supports they report agree.
  expect_true(all(prox$benchmark$glmnet$kkt$kkt <= prox$benchmark$glmnet$kkt$kkt_tolerance))
  expect_identical(
    default$benchmark$glmnet$support$group_history$in_support,
    prox$benchmark$glmnet$support$group_history$in_support
  )

  # A sparse group lasso has an L1 part and stays on the proximal route.
  withr::local_options(flexyBayes.selection_grouped_solver = NULL)
  sparse <- fbgb_fit(d, quote(sparse_group_lasso(scale = 1, alpha = 0.3, group = "block")))
  expect_identical(call_count, 0L)
  expect_gt(nrow(sparse$benchmark$glmnet$kkt), 0L)
  expect_true(all(sparse$benchmark$glmnet$kkt$route == "internal_prox"))
})


test_that("an elastic-net path keeps the absolute 1e-6 at every point", {
  skip_if_not_installed("glmnet")
  d <- fbgb_small_study()
  d$y <- 1000 * d$y
  fit <- fb_select(
    fbgs_base(d),
    candidates = fbgs_candidates(d),
    penalty = fb_penalty(
      term(everything()) ~ elastic_net(lambda = 0.1, alpha = 0.5)
    ),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~trial, v = 3L),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 5L
  )
  kkt <- fit$benchmark$glmnet$kkt
  expect_gt(nrow(kkt), 0L)
  expect_true(all(kkt$kkt_tolerance == 1e-6))
  expect_false(any(kkt$route == "internal_bcd"))
})
