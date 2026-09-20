# =============================================================================
# The support rule of the penalised optimiser: which penalised coefficients a
# path point holds, read through fb_select(), tidy() and the stability record.
#
# At s = 1 the canonical lambda is the fold's lambda_max, where the exact
# solution is the zero vector. A solver leaves values of order 1e-16 there on
# the column attaining the anchor, and read as held they are a selection the
# optimum does not contain. Two oracles, neither the package's:
#
#   The closed form at the anchor. At s = 1 every penalised coefficient is
#   zero, so nothing is held, retained or entered there.
#
#   The exact optimum elsewhere. Each fold's whitened, residualised problem is
#   rebuilt from base R and solved by an active-set quadratic programme
#   (`quadprog`) in a sign orthant, refined until the subgradient residual
#   written in this file is below 1e-8. The path's support at every other
#   grid point must equal that optimum's.
#
# What each block holds:
#
#   The chosen support. On the two data sets below, the one-standard-error
#   rule chooses s = 1; the chosen support, `tidy()$retained`, the fold
#   counts, the stored coefficients and the history at that point are empty.
#
#   The exact optimum. Lasso and elastic-net supports at every point of every
#   fold equal the quadratic programme's.
#
#   One rule for every read-out. The history, the entrance columns of
#   `tidy()`, the group profile and the stability frequencies agree with the
#   stored coefficient tables read as `estimate != 0`, and none of them holds
#   a column at s = 1.
#
#   The proximal provider. A grouped penalty holds nothing at s = 1.
#
#   The rule itself. At the anchor it returns zeros; elsewhere it removes a
#   coefficient only when the point still meets the path's KKT threshold
#   without it.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

# Twenty-four sites of five rows, one protected covariate, ten candidates of
# which the first three carry a weak signal, and a site random effect. At
# seeds 7 and 3 the one-standard-error rule chooses s = 1.
fbsz_data <- function(seed, n_site = 24L, per = 5L, p = 10L, signal = 0.3) {
  .fb_sel_with_seed(seed, {
    n <- n_site * per
    site <- factor(rep(sprintf("s%02d", seq_len(n_site)), each = per))
    d <- data.frame(site = site, x = stats::rnorm(n))
    d$region <- factor(paste0("r", (as.integer(site) - 1L) %% 3L + 1L))
    for (j in seq_len(p)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    u <- stats::rnorm(n_site, sd = 1)
    d$y <- 2 + 0.5 * d$x +
      signal * (1.0 * d$z1 - 0.9 * d$z2 + 0.8 * d$z3) +
      u[as.integer(site)] + stats::rnorm(n, sd = 1)
    d
  })
}

fbsz_candidates <- stats::as.formula(
  paste("~", paste0("z", 1:10, collapse = " + "))
)

# A base carrying the parsed model, the data and the call record the fold
# routes read, with no engine behind it.
fbsz_stub_base <- function(d) {
  fb <- fb_from_asreml(y ~ x, random = ~site, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = y ~ x, random = ~site, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

fbsz_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

fbsz_select <- function(d,
                        penalty = fb_penalty(
                          term(everything()) ~ laplace(scale = 0.15)
                        ),
                        decision = fb_decision(method = "full",
                                               metric = "rmse")) {
  fb_select(
    fbsz_stub_base(d), fbsz_candidates,
    penalty = penalty,
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 6L),
    decision = decision,
    seed = 1L
  )
}

# The stored coefficient table of the benchmark grid, candidate rows only.
fbsz_grid_coefficients <- function(sel) {
  coef <- as.data.frame(sel$benchmark$glmnet$coefficients)
  coef[coef$arm == "benchmark" & coef$role == "candidate", ]
}

# --- the exact-optimum oracle ------------------------------------------- #

# The fold's whitened, protected-residualised problem, from base R alone:
# forwardsolve against the fold's Cholesky factor and `qr.resid()` for the
# residualisation, sharing no code with the reduction under test.
fbsz_hand_payload <- function(substrate) {
  train_vec <- substrate$train_rows
  lower <- t(substrate$blocks$factor_tt)
  whiten <- function(m) forwardsolve(lower, m)
  kept_vec <- colnames(substrate$x_candidate)
  kept_vec <- kept_vec[kept_vec %in% substrate$candidate_kept]
  x0 <- whiten(substrate$x_protected[train_vec, , drop = FALSE])
  xp <- whiten(substrate$x_candidate[train_vec, kept_vec, drop = FALSE])
  yw <- drop(whiten(matrix(substrate$response[train_vec], ncol = 1L)))
  qr0 <- qr(x0)
  xt <- qr.resid(qr0, xp)
  colnames(xt) <- kept_vec
  list(x = xt, y = qr.resid(qr0, yw))
}

# The subgradient residual of 0.5 ||z - A b||^2 + l1 ||b||_1 + 0.5 l2 ||b||^2,
# written from its definition.
fbsz_residual <- function(gram, atz, b, l1, l2) {
  g <- drop(gram %*% b) - atz + l2 * b
  max(ifelse(b != 0, abs(g + l1 * sign(b)), pmax(abs(g) - l1, 0)))
}

# The exact optimum by an active-set quadratic programme in one sign orthant,
# the orthant refined from the gradient until the residual certifies it.
fbsz_exact <- function(A, z, l1, l2) {
  gram <- crossprod(A)
  atz <- drop(crossprod(A, z))
  p <- ncol(A)
  signs <- ifelse(atz >= 0, 1, -1)
  for (pass in seq_len(20L)) {
    D <- (signs %o% signs) * gram + diag(l2, p)
    sol <- quadprog::solve.QP(D, signs * atz - l1, diag(p), rep(0, p))
    u <- pmax(sol$solution, 0)
    u[sol$iact[sol$iact > 0]] <- 0
    b <- stats::setNames(signs * u, colnames(A))
    if (fbsz_residual(gram, atz, b, l1, l2) < 1e-8) {
      return(b)
    }
    g <- drop(gram %*% b) - atz + l2 * b
    flip <- b == 0 & abs(g) > l1
    signs[flip] <- ifelse(g[flip] > 0, -1, 1)
  }
  stop("the quadratic programme did not certify an optimum")
}


# --- (i) the chosen support at s = 1 -------------------------------------- #

test_that("a one-standard-error choice of s = 1 reports the empty support the exact solution has", {
  skip_if_not_installed("glmnet")
  fbsz_quiet()
  for (seed in c(7L, 3L)) {
    sel <- fbsz_select(fbsz_data(seed))
    arm <- sel$benchmark$glmnet
    chosen <- arm$selection
    # The case under test: the rule chose the null anchor.
    expect_identical(chosen$selected_s, 1)
    expect_identical(chosen$selected_lambda_index, 1L)

    expect_identical(chosen$selected_terms, character(0))
    expect_identical(sel$selection$chosen_terms, character(0))
    expect_identical(sel$selection$chosen_size, 0L)
    td <- tidy(sel)
    expect_false(any(td$retained[td$role == "candidate"]))
    expect_true(all(chosen$candidate_availability$n_selected == 0L))

    # The stored coefficients, the support rows and the complexity the rule
    # compared at that point.
    coef <- fbsz_grid_coefficients(sel)
    at_anchor <- coef$estimate[coef$lambda_index == 1L]
    expect_identical(
      at_anchor[!is.na(at_anchor)],
      rep(0, sum(!is.na(at_anchor)))
    )
    rows <- as.data.frame(arm$support$by_lambda)
    rows <- rows[rows$arm == "benchmark" & rows$lambda_index == 1L, ]
    expect_identical(nrow(rows), 6L)
    expect_true(all(rows$n_nonzero == 0L))
    expect_true(all(rows$terms == ""))
    scores <- as.data.frame(chosen$scores)
    expect_identical(scores$complexity[scores$lambda_index == 1L], 0)

    # Nothing enters at the anchor.
    history <- arm$support$history
    expect_false(any(history$in_support[history$s == 1]))
    candidates <- td[td$role == "candidate", ]
    expect_true(all(is.na(candidates$path_entry_s) |
                      candidates$path_entry_s < 1))
  }
})


# --- (ii) the exact optimum at every point ------------------------------- #

test_that("the path's support equals the exact optimum at every point of every fold", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("quadprog")
  d <- fbsz_data(7L)
  cv <- .fb_sel_cv(
    d, "y", ~x, fbsz_candidates, random = "site",
    resamples = fb_resamples(unit = ~site, v = 6L), seed = 1L
  )
  blocks <- list(
    laplace = list(list(
      family = "laplace", columns = 1:10,
      hyperparameters = list(scale = 0.15)
    )),
    elastic_net = list(list(
      family = "elastic_net", columns = 1:10,
      hyperparameters = list(lambda = 0.2, alpha = 0.6)
    ))
  )
  for (family in names(blocks)) {
    path <- .fb_sel_glmnet_path(
      cv$substrates, blocks[[family]], folds = cv$folds, seed = 1L
    )
    alpha <- path$penalty$alpha_composite
    history <- path$support_by_resolution$history
    compared <- 0L
    anchors <- 0L
    for (label in names(cv$substrates)) {
      payload <- fbsz_hand_payload(cv$substrates[[label]])
      # The oracle's problem is the path's: the same null anchor.
      expect_equal(
        max(abs(drop(crossprod(payload$x, payload$y)))) / alpha,
        path$selection$lambda_max[[label]],
        tolerance = 1e-10
      )
      rows <- history[history$fold == label, ]
      for (index in sort(unique(rows$lambda_index))) {
        at_point <- rows[rows$lambda_index == index, ]
        held <- sort(at_point$column[at_point$in_support])
        if (at_point$s[[1L]] >= 1) {
          # The closed form: lambda_max is where zero becomes the solution.
          exact <- character(0)
          anchors <- anchors + 1L
        } else {
          lambda <- at_point$lambda_canonical[[1L]]
          b <- fbsz_exact(payload$x, payload$y, lambda * alpha,
                          lambda * (1 - alpha))
          exact <- sort(names(b)[b != 0])
        }
        expect_identical(held, exact,
                         label = paste(family, label, "point", index))
        compared <- compared + 1L
      }
    }
    expect_identical(compared, 20L * length(cv$substrates))
    expect_identical(anchors, length(cv$substrates))
  }
})


# --- (iii) one rule for every read-out ------------------------------------ #

test_that("the history, the entrance columns and stability read the same rule and hold nothing at s = 1", {
  skip_if_not_installed("glmnet")
  fbsz_quiet()
  sel <- fbsz_select(
    fbsz_data(7L),
    decision = fb_decision(method = "full", metric = "rmse",
                           stability_subsamples = 8L)
  )
  arm <- sel$benchmark$glmnet

  # The cross-validated path: membership is the stored coefficient read as
  # `estimate != 0`, at every point, fold and column.
  coef <- fbsz_grid_coefficients(sel)
  history <- arm$support$history
  key_history <- paste(history$fold, history$lambda_index, history$column)
  key_coef <- paste(coef$fold, coef$lambda_index, coef$term)
  estimate <- coef$estimate[match(key_history, key_coef)]
  expect_false(anyNA(match(key_history, key_coef)))
  expect_identical(
    history$in_support,
    !is.na(estimate) & estimate != 0
  )
  expect_false(any(history$in_support[history$s == 1]))

  # The entrance columns and the group profile start below the anchor.
  td <- tidy(sel)
  candidates <- td[td$role == "candidate", ]
  for (field in c("path_entry_s", "entry_s_median", "entry_s_q25",
                  "entry_s_q75")) {
    expect_true(all(is.na(candidates[[field]]) | candidates[[field]] < 1),
                label = field)
  }
  tg <- tidy(sel, level = "group")
  expect_true(all(is.na(tg$path_entry_s) | tg$path_entry_s < 1))

  # Stability: every frequency is the stored per-subsample coefficient read
  # as `estimate != 0` over the B subsamples, and none is positive at s = 1.
  stability <- arm$stability
  expect_true(isTRUE(stability$ran))
  sub_coef <- stability$coefficients
  expect_identical(
    sub_coef$estimate[sub_coef$s == 1 & !is.na(sub_coef$estimate)],
    rep(0, sum(sub_coef$s == 1 & !is.na(sub_coef$estimate)))
  )
  frequency <- stability$path[stability$path$level == "column", ]
  hand <- vapply(seq_len(nrow(frequency)), function(i1) {
    e <- sub_coef$estimate[sub_coef$term == frequency$name[[i1]] &
                             sub_coef$lambda_index ==
                               frequency$lambda_index[[i1]]]
    sum(!is.na(e) & e != 0) / 8
  }, numeric(1L))
  expect_identical(frequency$frequency, hand)
  expect_true(all(stability$path$n_selected[stability$path$s == 1] == 0L))
  expect_identical(stability$stable, character(0))
})


# --- (iv) the proximal provider ----------------------------------------- #

test_that("a grouped penalty solved by the proximal provider holds nothing at s = 1", {
  fbsz_quiet()
  sel <- fbsz_select(
    fbsz_data(7L),
    penalty = fb_penalty(
      term(z1, z2, z3, z4, z5) ~ group_lasso(scale = 0.3, group = "block"),
      term(z6, z7, z8, z9, z10) ~ group_lasso(scale = 0.3, group = "block")
    )
  )
  arm <- sel$benchmark$glmnet
  expect_identical(sel$provider, "prox")
  coef <- fbsz_grid_coefficients(sel)
  at_anchor <- coef$estimate[coef$lambda_index == 1L]
  expect_identical(
    at_anchor[!is.na(at_anchor)],
    rep(0, sum(!is.na(at_anchor)))
  )
  history <- arm$support$history
  expect_false(any(history$in_support[history$s == 1]))
  groups <- arm$support$group_history
  expect_false(any(groups$in_support[groups$s == 1]))
})


# --- (v) the rule itself ------------------------------------------------- #

test_that("the rule zeroes the anchor and removes a coefficient only while the point stays certified", {
  tol <- .fb_sel_path_kkt_tolerance()
  basis <- .fb_sel_with_seed(3L, qr.Q(qr(matrix(stats::rnorm(100L), 50L))))
  kkt <- function(A, z, l1) {
    function(b) .fb_sel_kkt_residual_lambda(A, z, b, l1, 0)
  }

  # At the null anchor every coefficient is zero, and the names are kept.
  beta <- c(a = -7.6e-17, b = 0.3)
  expect_identical(
    .fb_sel_path_support(beta, basis, 0, kkt(basis, numeric(50L), 0.1),
                         at_null_anchor = TRUE),
    c(a = 0, b = 0)
  )

  # Orthonormal columns and a second column on its threshold (A_2'z = l1),
  # so the exact optimum is (0.5, 0). A solver answer of (0.5, 3e-7) meets the
  # threshold; the rule returns the exact optimum's support.
  l1 <- 0.1
  z <- drop(basis %*% c(0.5 + l1, l1))
  beta <- c(a = 0.5, b = 3e-7)
  expect_lte(kkt(basis, z, l1)(beta), tol)
  out <- .fb_sel_path_support(beta, basis, 0, kkt(basis, z, l1))
  expect_identical(out[["b"]], 0)
  expect_identical(out[["a"]], beta[["a"]])
  expect_lte(kkt(basis, z, l1)(out), tol)
  # Above the threshold the coefficient is kept.
  beta_large <- c(a = 0.5, b = 2e-6)
  expect_identical(
    .fb_sel_path_support(beta_large, basis, 0, kkt(basis, z, l1)),
    beta_large
  )
  # The threshold is on the gradient's scale: on a column of squared norm
  # 1e-4 on its threshold, a coefficient of 3e-3 moves its own condition by
  # 3e-7 and is removed, although it is far above 1e-6 as a number.
  narrow <- cbind(basis[, 1L], 0.01 * basis[, 2L])
  z_narrow <- drop(basis %*% c(0.5 + l1, l1 / 0.01))
  beta_narrow <- c(a = 0.5, b = 3e-3)
  expect_lte(kkt(narrow, z_narrow, l1)(beta_narrow), tol)
  expect_identical(
    .fb_sel_path_support(beta_narrow, narrow, 0, kkt(narrow, z_narrow, l1)),
    c(a = 0.5, b = 0)
  )

  # A column of squared norm 1e6 correlated 0.9 with a unit column, and an
  # exact optimum whose second coefficient is 5e-7: its own condition moves by
  # 5e-7 when it is removed, but the first column's moves by 4.5e-4, so the
  # point without it is not certified and the solved point is kept.
  A <- cbind(1000 * basis[, 1L], 0.9 * basis[, 1L] + sqrt(0.19) * basis[, 2L])
  beta <- c(a = 0.2, b = 5e-7)
  gram <- crossprod(A)
  z <- drop(A %*% (beta + l1 * solve(gram, sign(beta))))
  expect_lte(kkt(A, z, l1)(beta), tol)
  removed <- beta
  removed[["b"]] <- 0
  expect_gt(kkt(A, z, l1)(removed), tol)
  expect_identical(.fb_sel_path_support(beta, A, 0, kkt(A, z, l1)), beta)
})


# --- (vi) where the rule is applied -------------------------------------- #

test_that("every point of a selecting path passes through the rule, flagged as the anchor exactly where s >= 1", {
  skip_if_not_installed("glmnet")
  cv <- .fb_sel_cv(
    fbsz_data(7L), "y", ~x, fbsz_candidates, random = "site",
    resamples = fb_resamples(unit = ~site, v = 3L), seed = 1L
  )
  block <- function(family, hyper, columns = 1:10) {
    list(family = family, columns = columns, hyperparameters = hyper)
  }
  specs <- list(
    laplace = list(block("laplace", list(scale = 0.15))),
    group_lasso = list(
      block("group_lasso", list(scale = 0.3, group = "block"), 1:5),
      block("group_lasso", list(scale = 0.3, group = "block"), 6:10)
    ),
    exclusive_lasso = list(
      block("exclusive_lasso", list(scale = 0.7, group = "block"), 1:5),
      block("exclusive_lasso", list(scale = 0.7, group = "block"), 6:10)
    ),
    ridge = list(block("ridge", list(scale = 0.6)))
  )
  flags <- logical(0)
  local_mocked_bindings(
    .fb_sel_path_support = function(beta, A, lambda_l2, residual,
                                    at_null_anchor = FALSE) {
      flags <<- c(flags, at_null_anchor)
      beta
    }
  )
  for (family in names(specs)) {
    flags <- logical(0)
    path <- .fb_sel_glmnet_path(
      cv$substrates, specs[[family]], folds = cv$folds,
      decision = list(n_lambda = 6L), seed = 1L
    )
    if (identical(family, "ridge")) {
      expect_length(flags, 0L)
      next
    }
    # The coordinates table is written once per point in the order the
    # points are solved, so it lines up with the recorded calls.
    coords <- path$path_coordinates
    expect_identical(flags, !is.na(coords$s) & coords$s >= 1, label = family)
    first <- coords$arm == "benchmark" & coords$lambda_index == 1L
    if (identical(family, "exclusive_lasso")) {
      expect_false(any(flags), label = family)
    } else {
      expect_identical(flags[first], rep(TRUE, 3L), label = family)
    }
  }
})
