# =============================================================================
# The exclusive lasso on blocks of strongly correlated candidates, such as one
# soil property measured at several depths.
#
# Oracle: the exact optimum of
#
#   0.5 ||z - A b||^2 + (le / 2) sum_g (||b_g||_1)^2
#
# from an active-set quadratic programme (`quadprog`) in one sign orthant,
# refined until the subgradient residual written in this file is below
# 1e-10. On the orthant with signs s, b = s u with u >= 0 and ||b_g||_1 is the
# sum of u over the group, so the problem is a positive-definite quadratic
# programme in u.
#
# Bounds below are from review/units/U-EXCL_2026-09-24. A polish that carries
# the exclusive term from its previous solve diverges on these blocks and
# fails every block of this file.
# =============================================================================

# --- fixtures ------------------------------------------------------------- #

# Blocks of `sizes` columns whose members correlate at `rho`; the first
# member of each block carries the signal.
fbxb_matrix <- function(seed, rho, n = 200L, sizes = c(5L, 4L, 3L)) {
  .fb_sel_with_seed(seed, {
    columns <- list()
    for (g in seq_along(sizes)) {
      latent <- stats::rnorm(n)
      for (j in seq_len(sizes[[g]])) {
        columns[[length(columns) + 1L]] <- sqrt(rho) * latent +
          sqrt(1 - rho) * stats::rnorm(n)
      }
    }
    a <- scale(do.call(cbind, columns), scale = FALSE)
    colnames(a) <- paste0("v", seq_len(ncol(a)))
    first <- cumsum(sizes) - sizes + 1L
    z <- drop(a[, first] %*% rep(0.5, length(sizes))) + stats::rnorm(n)
    list(
      A = a,
      z = z - mean(z),
      groups = list(index = rep(seq_along(sizes), sizes))
    )
  })
}

fbxb_penalty <- function(le) {
  list(lambda_l1 = 0, lambda_l2 = 0, lambda_group = 0, lambda_exclusive = le,
    family = "exclusive_lasso")
}

# The subgradient residual of the exclusive lasso, from its definition.
fbxb_residual <- function(gram, atz, b, le, index) {
  g <- drop(gram %*% b) - atz
  threshold <- le * tapply(abs(b), index, sum)[index]
  max(ifelse(b != 0, abs(g + threshold * sign(b)), pmax(abs(g) - threshold, 0)))
}

fbxb_exact <- function(A, z, le, index) {
  gram <- crossprod(A)
  atz <- drop(crossprod(A, z))
  p <- ncol(A)
  same_group <- outer(index, index, "==") * 1
  signs <- ifelse(atz >= 0, 1, -1)
  for (pass in seq_len(20L)) {
    D <- (signs %o% signs) * gram + le * same_group
    sol <- quadprog::solve.QP(D, signs * atz, diag(p), rep(0, p))
    u <- pmax(sol$solution, 0)
    u[sol$iact[sol$iact > 0]] <- 0
    b <- signs * u
    if (fbxb_residual(gram, atz, b, le, index) < 1e-10) {
      return(b)
    }
    g <- drop(gram %*% b) - atz
    threshold <- le * tapply(abs(b), index, sum)[index]
    flip <- b == 0 & abs(g) > threshold
    signs[flip] <- ifelse(g[flip] > 0, -1, 1)
  }
  stop("the quadratic programme did not certify an optimum")
}

# A multi-trial design shaped like site-level soil covariates: localities
# hold two trials each, every candidate takes one value per trial, and each
# block's members correlate at `rho`. Only `clay_1` carries signal.
fbxb_study <- function(seed = 11L, n_loc = 12L, per = 25L,
                       sizes = c(clay = 5L, ph = 4L, soc = 3L),
                       rho = c(clay = 0.97, ph = 0.9, soc = 0.8)) {
  .fb_sel_with_seed(seed, {
    n_trial <- 2L * n_loc
    loc_of_trial <- rep(seq_len(n_loc), each = 2L)
    trial_index <- rep(seq_len(n_trial), each = per)
    d <- data.frame(
      loc = factor(sprintf("L%02d", loc_of_trial[trial_index])),
      trial = factor(sprintf("t%02d", trial_index)),
      x = stats::rnorm(length(trial_index))
    )
    for (b in names(sizes)) {
      latent <- sqrt(0.7) * stats::rnorm(n_loc)[loc_of_trial] +
        sqrt(0.3) * stats::rnorm(n_trial)
      for (j in seq_len(sizes[[b]])) {
        member <- sqrt(rho[[b]]) * latent +
          sqrt(1 - rho[[b]]) * stats::rnorm(n_trial)
        d[[paste0(b, "_", j)]] <- member[trial_index]
      }
    }
    d$y <- 2 + 0.5 * d$x + 0.6 * d$clay_1 +
      stats::rnorm(n_trial, sd = 0.3)[trial_index] +
      stats::rnorm(nrow(d), sd = 1)
    d
  })
}

fbxb_stub_base <- function(d) {
  fb <- fb_from_asreml(y ~ x, random = ~trial, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = y ~ x, random = ~trial, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}


# --- the restricted polish ------------------------------------------------ #

test_that("the polish solves a correlated exclusive block exactly", {
  skip_if_not_installed("quadprog")
  for (seed in 1:4) {
    f <- fbxb_matrix(seed, rho = 0.99)
    exact <- fbxb_exact(f$A, f$z, le = 2, f$groups$index)
    start <- exact
    held <- exact != 0
    start[held] <- exact[held] * (1 + 1e-4)
    polished <- .fb_sel_polish_grouped(
      crossprod(f$A),
      drop(crossprod(f$A, f$z)),
      start,
      .fb_sel_penalty_parts(fbxb_penalty(2)),
      .fb_sel_group_parts(f$groups, ncol(f$A))
    )
    expect_false(is.null(polished))
    expect_lt(max(abs(polished - exact)), 1e-10)
  }
})


test_that("the grouped mode certifies correlated exclusive blocks by polish", {
  skip_if_not_installed("quadprog")
  for (seed in 1:4) {
    for (le in c(1, 2, 4)) {
      f <- fbxb_matrix(seed, rho = 0.99)
      solved <- .fb_sel_prox_mode_grouped(
        f$A, f$z, fbxb_penalty(le), f$groups, tol = 1e-12, max_iter = 2e4
      )
      exact <- fbxb_exact(f$A, f$z, le, f$groups$index)
      expect_true(solved$polished)
      # At most 160 measured; 8,820 with the exclusive term carried over.
      expect_lte(solved$iterations, 400L)
      expect_lt(max(abs(solved$beta - exact)), 1e-9)
      expect_identical(unname(solved$beta != 0), unname(exact != 0))
    }
  }
})


# --- through fb_select() -------------------------------------------------- #

test_that("fb_select holds every exclusive path point at the exact optimum", {
  skip_if_not_installed("quadprog")
  withr::local_options(
    flexyBayes.selection_workers = 1L,
    flexyBayes.silence_default_prior_note = TRUE
  )
  d <- fbxb_study()
  columns <- grep("^(clay|ph|soc)_", names(d), value = TRUE)
  calls <- list()
  solver <- .fb_sel_prox_mode_grouped
  testthat::local_mocked_bindings(
    .fb_sel_prox_mode_grouped = function(A, z, penalty, groups, ...) {
      out <- solver(A, z, penalty, groups, ...)
      calls[[length(calls) + 1L]] <<- list(A = A, z = z, penalty = penalty,
        groups = groups, beta = out$beta, polished = out$polished)
      out
    }
  )
  fit <- fb_select(
    fbxb_stub_base(d),
    candidates = stats::reformulate(columns),
    penalty = fb_penalty(
      term(matches("^clay_")) ~ exclusive_lasso(scale = 0.7, group = "block"),
      term(matches("^ph_")) ~ exclusive_lasso(scale = 0.7, group = "block"),
      term(matches("^soc_")) ~ exclusive_lasso(scale = 0.7, group = "block")
    ),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~loc, v = 3L),
    decision = fb_decision(method = "full", metric = "rmse", rule = "one_se"),
    seed = 5L
  )
  arm <- fit$benchmark$glmnet
  expect_true(all(arm$kkt$kkt <= arm$kkt$kkt_tolerance))
  expect_identical(length(calls), nrow(arm$kkt))
  polished <- vapply(calls, function(call) call$polished, logical(1L))
  # 65 of 66 measured; 30 of 66 with the exclusive term carried over.
  expect_gt(mean(polished), 0.9)
  for (call in calls) {
    exact <- fbxb_exact(call$A, call$z, call$penalty$lambda_exclusive,
      call$groups$index)
    expect_identical(unname(call$beta != 0), unname(exact != 0))
    if (call$polished) {
      # 3e-14 measured; 4.4e-5 with the exclusive term carried over.
      expect_lt(max(abs(call$beta - exact)), 1e-10)
    }
  }

  # Competition inside a block: at the chosen point every block holds one or
  # two of its members on every fold, and more members enter as the penalty
  # falls.
  coef <- as.data.frame(arm$coefficients)
  coef <- coef[coef$arm == "benchmark" & coef$role == "candidate", ]
  coef$block <- sub("_[0-9]+$", "", coef$term)
  chosen <- coef[coef$lambda_index == arm$selection$selected_lambda_index, ]
  held <- tapply(chosen$estimate != 0, list(chosen$fold, chosen$block), sum)
  expect_true(all(held >= 1L & held <= 2L))
  smallest <- coef[coef$lambda_index == max(coef$lambda_index), ]
  expect_gt(sum(smallest$estimate != 0), sum(chosen$estimate != 0))
})
