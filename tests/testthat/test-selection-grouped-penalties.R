# =============================================================================
# Grouped penalties on the penalised-optimiser provider (unit P8).
#
# Contract: final specification row B9 --- group lasso, sparse group lasso,
# exclusive lasso and the group-exclusive composite, with three acceptance
# items: the exclusive proximal operator against a brute-force convex solve to
# 1e-8, a KKT certificate for every fitted solution, and known-truth recovery
# curves on designs whose group structure is known.
#
# The canonical objective is the one R/selection_bridges.R already states,
# widened from two numbers to four plus a grouping:
#
#   0.5 ||z - A beta||^2
#     + l1 ||beta||_1 + (l2 / 2) ||beta||_2^2
#     + lg sum_g w_g ||beta_g||_2
#     + (le / 2) sum_g (||beta_g||_1)^2
#
# Nothing about a grouped solution is taken from a solver's convergence flag.
# Three oracles grade the work and none of them shares a source with the code
# it grades:
#
#   quadprog   the exclusive prox restricted to the nonnegative orthant is a
#              positive-definite quadratic programme, so an active-set dual
#              method solves it exactly. This is the 1e-8 acceptance.
#   grpreg     an independent group-lasso implementation, at the objective
#              convention MEASURED in this file rather than read from its
#              documentation.
#   subgradient  the optimality residual computed from the definition of each
#              component's subdifferential, at whatever point it is handed.
#
# Injected defects that turn these gates red, all standing in this file:
#   G1  the naive composite prox -- exclusive thresholding followed by block
#       soft thresholding -- which uses a coordinate threshold computed
#       without the group shrinkage and so thresholds too hard
#   G2  the sqrt(|g|) group weight replaced by 1, which grpreg's measured
#       convention rejects by two orders of magnitude
#   G3  a lasso path on the same folds as a group-lasso path, so the
#       "a term enters as a whole" assertion is not vacuous
# =============================================================================

# --- fixtures ------------------------------------------------------------- #

# Ten sites of six rows. Two factors of three levels each (two contrast
# columns apiece) and four numeric candidates, of which z1 and z2 carry
# signal. The factor is the group the group lasso acts on; the numeric cluster
# is the group the exclusive lasso acts on.
fbgp_data <- function(n_cluster = 10L, per = 6L, seed = 51L) {
  .fb_sel_with_seed(seed, {
    n <- n_cluster * per
    site <- factor(rep(seq_len(n_cluster), each = per))
    u <- stats::rnorm(n_cluster, sd = 1.2)
    d <- data.frame(
      site = site,
      x = stats::rnorm(n),
      f = factor(rep(c("a", "b", "c"), length.out = n)),
      g = factor(rep(c("p", "q", "r"), each = 2L, length.out = n))
    )
    for (j in seq_len(4L)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    d$y <- 3 +
      0.8 * d$x +
      1.3 * d$z1 -
      1.0 * d$z2 +
      1.1 * as.integer(d$f == "b") +
      u[as.integer(site)] +
      stats::rnorm(n, sd = 0.5)
    d
  })
}


fbgp_substrates <- function(data = fbgp_data(), v = 3L, seed = 77L) {
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = v), data, seed = seed)
  substrates <- lapply(folds, function(split) {
    .fb_sel_fold_substrate(
      data,
      split,
      response = "y",
      protected = ~x,
      candidates = ~ f + g + z1 + z2 + z3 + z4,
      random = "site"
    )
  })
  names(substrates) <- vapply(folds, function(f) f$label, character(1L))
  list(data = data, folds = folds, substrates = substrates)
}


# The candidate design of `fbgp_substrates()` is f (2 contrast columns), g (2)
# and z1..z4, so eight columns in declared order.
fbgp_block <- function(family, hyper, columns = 1:8) {
  list(list(family = family, columns = columns, hyperparameters = hyper))
}


# Accumulate a worst case under a named seed.
#
# `.fb_sel_with_seed()` forces its argument as a promise, so the loop below
# evaluates in THIS frame and a `<<-` inside it would begin its search in the
# frame's parent -- skipping the accumulator entirely and writing to the
# global environment, while the local one stayed at its initial value and
# every `expect_lt()` against it passed vacuously. The seeded block therefore
# returns its value rather than assigning one.
fbgp_worst <- function(seed, trials, one) {
  .fb_sel_with_seed(seed, {
    worst <- 0
    for (trial in seq_len(trials)) {
      worst <- max(worst, one())
    }
    worst
  })
}


# A dense matrix problem with a declared grouping, for the solver-level
# blocks. Four groups of sizes 3, 2, 2 and 1.
fbgp_matrix <- function(n = 70L, seed = 8L) {
  .fb_sel_with_seed(seed, {
    groups <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L)
    p <- length(groups)
    a <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    colnames(a) <- paste0("v", seq_len(p))
    truth <- c(1.4, -0.9, 0.6, 0, 0, 1.1, 0, 0.5)
    z <- drop(a %*% truth) + stats::rnorm(n, sd = 0.5)
    list(A = a, z = z, n = n, p = p, groups = groups, truth = truth)
  })
}


# The independent convex oracle for the exclusive proximal operator.
#
# prox_t(v) = argmin_x 0.5 ||x - v||^2 + (t le / 2) (||x||_1)^2. Every
# coordinate of the minimiser carries the sign of the matching coordinate of
# `v`, so with u = |x| >= 0 and a = |v| the problem is
#
#   min_{u >= 0} 0.5 u' (I + t le 11') u - a'u
#
# a positive-definite quadratic programme. `quadprog::solve.QP` is Goldfarb
# and Idnani's dual active-set method: a different algorithm by a different
# author, and exact rather than iterative, which is what lets the comparison
# be stated at 1e-8 rather than at a solver tolerance.
fbgp_prox_exclusive_qp <- function(v, t, le) {
  a <- abs(v)
  p <- length(a)
  sol <- quadprog::solve.QP(
    Dmat = diag(p) + (t * le) * matrix(1, p, p),
    dvec = a,
    Amat = diag(p),
    bvec = rep(0, p)
  )
  sign(v) * pmax(sol$solution, 0)
}


# G1. The naive composite: threshold for the exclusive part without the group
# shrinkage, then shrink the block. It is the obvious wrong answer and it is
# here so the subgradient block below cannot pass vacuously.
fbgp_naive_composite <- function(v, t, lg, le, weight) {
  inner <- .fb_sel_prox_group_block(
    v, t, l1 = 0, l2 = 0, lg = 0, le = le, weight = weight
  )
  radius <- sqrt(sum(inner^2))
  if (radius == 0) {
    return(inner)
  }
  inner * max(0, 1 - t * lg * weight / radius)
}


# The subgradient residual of the prox subproblem, written from the definition
# rather than from the formula under test.
fbgp_prox_residual <- function(x, v, t, l1, l2, lg, le, weight) {
  base <- (x - v) + t * l2 * x
  tau <- t * (l1 + le * sum(abs(x)))
  radius <- sqrt(sum(x^2))
  if (radius > 0) {
    gradient <- base + t * lg * weight * x / radius
    active <- abs(x) > 0
    out <- numeric(length(x))
    out[active] <- abs(gradient[active] + tau * sign(x[active]))
    out[!active] <- pmax(abs(gradient[!active]) - tau, 0)
    return(max(out))
  }
  soft <- sign(base) * pmax(abs(base) - t * l1, 0)
  max(0, sqrt(sum(soft^2)) - t * lg * weight)
}


# A design whose groups are already centred and orthonormal in grpreg's own
# metric (crossprod(X_g) / n = I), so grpreg's internal standardisation and
# within-group orthogonalisation are the identity and its returned
# coefficients are on the scale this package states its objective on.
fbgp_orthonormal_design <- function(n = 60L, groups = c(1, 1, 1, 2, 2, 3, 3, 4),
                                    seed = 99L) {
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
    truth <- c(1.2, -0.7, 0.3, 0, 0, 0.9, 0.4, 0)
    z <- drop(x %*% truth) + stats::rnorm(n, sd = 0.4)
    list(A = x, z = z - mean(z), n = n, groups = groups)
  })
}


# =============================================================================
# The grammar
# =============================================================================

test_that("the four grouped families are admitted and declare their grouping", {
  for (spec in list(
    list(call = quote(group_lasso(scale = 1)), family = "group_lasso",
         rule = "term"),
    list(call = quote(sparse_group_lasso(scale = 1, alpha = 0.5)),
         family = "sparse_group_lasso", rule = "term"),
    list(call = quote(exclusive_lasso(scale = 1)), family = "exclusive_lasso",
         rule = "block"),
    list(call = quote(group_exclusive(group_scale = 1, exclusive_scale = 2)),
         family = "group_exclusive", rule = "term")
  )) {
    formula <- stats::as.formula(call("~", quote(term(everything())), spec$call))
    pen <- fb_penalty(formula)
    expect_s3_class(pen, "fb_penalty")
    expect_identical(pen$specs[[1L]]$family, spec$family)
    # The grouping rule is never guessed: the family table states a default,
    # the specification carries it, and the certificate prints it.
    expect_identical(pen$specs[[1L]]$args$group, spec$rule)
  }

  # The default is a default, not a fixture: naming the other rule takes.
  named <- fb_penalty(term(everything()) ~ group_lasso(scale = 1, group = "block"))
  expect_identical(named$specs[[1L]]$args$group, "block")
})


test_that("a grouped hyperparameter outside its domain refuses by name", {
  expect_error(
    fb_penalty(term(everything()) ~ group_lasso(scale = 0)),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
  expect_error(
    fb_penalty(term(everything()) ~ sparse_group_lasso(scale = 1, alpha = 1)),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
  expect_error(
    fb_penalty(term(everything()) ~ exclusive_lasso(scale = -1)),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
  expect_error(
    fb_penalty(term(everything()) ~ group_lasso(scale = 1, group = "factor")),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
  # `adaptive()` has no provider and stays withdrawn: admitting four families
  # is not admitting five.
  expect_error(
    fb_penalty(term(everything()) ~ adaptive(scale = 1)),
    class = "flexybayes_refusal_selection_penalty_family_withdrawn"
  )
})


test_that("each grouped family maps to the canonical four-number penalty", {
  gl <- .fb_sel_canonical_penalty(
    "group_lasso", list(scale = 0.5, group = "term"), n = 40
  )
  expect_equal(gl$lambda_group, 2)
  expect_equal(gl$lambda_l1, 0)
  expect_equal(gl$lambda_exclusive, 0)
  expect_true(gl$grouped)
  expect_identical(gl$group_rule, "term")

  sgl <- .fb_sel_canonical_penalty(
    "sparse_group_lasso", list(scale = 0.5, alpha = 0.25, group = "term"),
    n = 40
  )
  expect_equal(sgl$lambda_l1, 0.5)
  expect_equal(sgl$lambda_group, 1.5)

  # The exclusive penalty is quadratic in the coefficients, so its scale is
  # stated the way ridge states one: le = 1 / scale^2. That is what makes a
  # singleton group exactly a ridge, which the solver block below checks.
  ex <- .fb_sel_canonical_penalty(
    "exclusive_lasso", list(scale = 0.5, group = "block"), n = 40
  )
  expect_equal(ex$lambda_exclusive, 4)
  expect_equal(ex$lambda_l1, 0)
  expect_identical(ex$group_rule, "block")

  ge <- .fb_sel_canonical_penalty(
    "group_exclusive",
    list(group_scale = 0.5, exclusive_scale = 0.25, group = "term"), n = 40
  )
  expect_equal(ge$lambda_group, 2)
  expect_equal(ge$lambda_exclusive, 16)

  # The three families that were already here keep their maps untouched.
  expect_equal(
    .fb_sel_canonical_penalty("ridge", list(scale = 0.8), n = 10)$lambda_l2,
    1 / 0.8^2
  )
  expect_false(
    .fb_sel_canonical_penalty("laplace", list(scale = 0.4), n = 10)$grouped
  )
})


# =============================================================================
# Acceptance 1: the exclusive prox against a brute-force convex solve
# =============================================================================

test_that("the exclusive prox matches an independent convex solve to 1e-8", {
  skip_if_not_installed("quadprog")
  worst <- fbgp_worst(2026L, 1500L, function() {
    p <- sample.int(8L, 1L)
    v <- stats::rnorm(p, sd = stats::runif(1L, 0.1, 4))
    if (stats::runif(1L) < 0.25) {
      v[sample.int(p, 1L)] <- 0
    }
    t <- 10^stats::runif(1L, -2, 1)
    le <- 10^stats::runif(1L, -3, 2)
    mine <- .fb_sel_prox_group_block(
      v, t, l1 = 0, l2 = 0, lg = 0, le = le, weight = 1
    )
    max(abs(mine - fbgp_prox_exclusive_qp(v, t, le)))
  })
  # The measurement has to be a measurement: a worst case of exactly zero
  # would mean the loop never ran, which is how this block passed before the
  # accumulator was fixed.
  expect_gt(worst, 0)
  # The acceptance is 1e-8. The closed form is exact up to floating point, so
  # the margin is reported rather than the threshold being met narrowly.
  expect_lt(worst, 1e-8)
  expect_lt(worst, 1e-11)
})


test_that("a singleton exclusive group is exactly a ridge", {
  # (le / 2) * (|beta|)^2 is (le / 2) beta^2 when the group has one member, so
  # the prox must be the ridge prox v / (1 + t le). This is the exact-limit
  # identity that fixes the family's scale convention.
  worst <- fbgp_worst(31L, 500L, function() {
    v <- stats::rnorm(1L, sd = 3)
    t <- 10^stats::runif(1L, -2, 1)
    le <- 10^stats::runif(1L, -3, 2)
    mine <- .fb_sel_prox_group_block(
      v, t, l1 = 0, l2 = 0, lg = 0, le = le, weight = 1
    )
    abs(mine - v / (1 + t * le))
  })
  expect_gt(worst, 0)
  expect_lt(worst, 1e-12)
})


test_that("every family's prox satisfies its own subgradient condition", {
  families <- list(
    ridge = c(l1 = 0, l2 = 1.3, lg = 0, le = 0),
    laplace = c(l1 = 0.7, l2 = 0, lg = 0, le = 0),
    elastic_net = c(l1 = 0.5, l2 = 0.9, lg = 0, le = 0),
    group_lasso = c(l1 = 0, l2 = 0, lg = 1.1, le = 0),
    sparse_group_lasso = c(l1 = 0.4, l2 = 0, lg = 0.9, le = 0),
    exclusive_lasso = c(l1 = 0, l2 = 0, lg = 0, le = 1.7),
    group_exclusive = c(l1 = 0, l2 = 0, lg = 0.8, le = 1.4)
  )
  for (nm in names(families)) {
    f <- families[[nm]]
    worst <- fbgp_worst(7L, 400L, function() {
      p <- sample.int(6L, 1L)
      v <- stats::rnorm(p, sd = stats::runif(1L, 0.05, 4))
      t <- 10^stats::runif(1L, -2, 1)
      w <- sqrt(p)
      x <- .fb_sel_prox_group_block(
        v, t, f[["l1"]], f[["l2"]], f[["lg"]], f[["le"]], w
      )
      fbgp_prox_residual(
        x, v, t, f[["l1"]], f[["l2"]], f[["lg"]], f[["le"]], w
      )
    })
    expect_lt(worst, 1e-10, label = paste0("prox residual for ", nm))
  }
})


test_that("the naive composite prox is wrong, so the check above is not vacuous", {
  # G1. Thresholding for the exclusive part before the group shrinkage is
  # applied uses a threshold computed at the wrong point, and over-thresholds.
  # The draw keeps the group alive: a penalty strong enough to zero the whole
  # block makes both operators return zero and agree for a reason that has
  # nothing to do with the composition being right.
  measured <- .fb_sel_with_seed(13L, {
    worst_deviation <- 0
    worst_residual <- 0
    for (trial in seq_len(200L)) {
      p <- sample(2:5, 1L)
      v <- stats::rnorm(p, sd = 2)
      t <- 10^stats::runif(1L, -1.5, -0.5)
      lg <- 10^stats::runif(1L, -1, 0)
      le <- 10^stats::runif(1L, -0.5, 0.5)
      w <- sqrt(p)
      mine <- .fb_sel_prox_group_block(v, t, 0, 0, lg, le, w)
      if (all(mine == 0)) {
        next
      }
      naive <- fbgp_naive_composite(v, t, lg, le, w)
      worst_deviation <- max(worst_deviation, max(abs(naive - mine)))
      worst_residual <- max(
        worst_residual, fbgp_prox_residual(naive, v, t, 0, 0, lg, le, w)
      )
    }
    list(deviation = worst_deviation, residual = worst_residual)
  })
  expect_gt(measured$deviation, 1e-3)
  expect_gt(measured$residual, 1e-3)
})


test_that("the grouped KKT residual collapses to the ungrouped one", {
  # With no group and no exclusive term the widened residual must BE the
  # residual R/selection_bridges.R already computed, not merely close to it.
  fixture <- fbgp_matrix()
  singletons <- seq_len(fixture$p)
  .fb_sel_with_seed(4L, {
    for (trial in seq_len(25L)) {
      beta <- stats::rnorm(fixture$p) * stats::rbinom(fixture$p, 1L, 0.6)
      for (pair in list(c(0.8, 0), c(0, 1.2), c(0.5, 0.7))) {
        penalty <- list(
          lambda_l1 = pair[[1L]], lambda_l2 = pair[[2L]],
          lambda_group = 0, lambda_exclusive = 0
        )
        expect_equal(
          .fb_sel_kkt_residual_grouped(
            fixture$A, fixture$z, beta, penalty,
            list(index = singletons, weights = rep(1, fixture$p))
          ),
          .fb_sel_kkt_residual_lambda(
            fixture$A, fixture$z, beta, pair[[1L]], pair[[2L]]
          )
        )
      }
    }
  })
})


# =============================================================================
# Acceptance 2: a certificate for every fitted solution
# =============================================================================

test_that("the grouped mode returns a certified point or refuses", {
  fixture <- fbgp_matrix()
  groups <- .fb_sel_group_index(
    colnames(fixture$A), paste0("t", fixture$groups), rule = "term"
  )
  expect_identical(groups$sizes, c(t1 = 3L, t2 = 2L, t3 = 2L, t4 = 1L))
  expect_equal(unname(groups$weights), sqrt(c(3, 2, 2, 1)))

  for (family in c("group_lasso", "sparse_group_lasso", "exclusive_lasso",
                   "group_exclusive")) {
    hyper <- switch(
      family,
      group_lasso = list(scale = 0.4, group = "term"),
      sparse_group_lasso = list(scale = 0.4, alpha = 0.4, group = "term"),
      exclusive_lasso = list(scale = 0.6, group = "term"),
      group_exclusive = list(group_scale = 0.5, exclusive_scale = 0.7,
        group = "term")
    )
    penalty <- .fb_sel_canonical_penalty(family, hyper, n = fixture$n)
    solved <- .fb_sel_prox_mode_grouped(
      fixture$A, fixture$z, penalty, groups, tol = 1e-9
    )
    expect_true(solved$converged, label = family)
    expect_lt(solved$kkt, 1e-9, label = paste0("kkt for ", family))
    # The residual is recomputed here, from the definition, at the point the
    # solver returned. A solver that reported its own flag would pass on its
    # own word; this does not.
    expect_lt(
      .fb_sel_kkt_residual_grouped(
        fixture$A, fixture$z, solved$beta, penalty, groups
      ),
      1e-9,
      label = paste0("independent kkt for ", family)
    )
  }
})


test_that("optimiser paths use one declared KKT acceptance threshold", {
  # The provider bridge has always been accepted against 1e-6. Grouped path
  # points must use the same independently recomputed certificate, while the
  # lower-level solver remains available at its stricter 1e-9 default.
  expect_identical(.fb_sel_path_kkt_tolerance(), 1e-6)
  expect_identical(formals(.fb_sel_prox_mode_grouped)$tol, 1e-9)
  expect_identical(formals(.fb_sel_prox_mode_lambda)$tol, 1e-10)
})


test_that("a mode that cannot be certified refuses rather than returning", {
  fixture <- fbgp_matrix()
  groups <- .fb_sel_group_index(
    colnames(fixture$A), paste0("t", fixture$groups), rule = "term"
  )
  penalty <- .fb_sel_canonical_penalty(
    "group_lasso", list(scale = 0.4, group = "term"), n = fixture$n
  )
  expect_error(
    .fb_sel_prox_mode_grouped(
      fixture$A, fixture$z, penalty, groups, tol = 1e-14, max_iter = 3L
    ),
    class = "flexybayes_refusal_selection_prox_mode_not_converged"
  )
})


# =============================================================================
# The independent package oracle
# =============================================================================

test_that("the group lasso reproduces grpreg at its measured convention", {
  skip_if_not_installed("grpreg")
  design <- fbgp_orthonormal_design()
  groups <- .fb_sel_group_index(
    colnames(design$A), paste0("t", design$groups), rule = "term"
  )

  # grpreg's objective convention is MEASURED here, not read: the weight that
  # drives its own KKT residual to its convergence floor is sqrt(|g|), and the
  # two competing weights are asserted to lose by orders of magnitude (G2).
  lambda_grpreg <- 0.30
  fit <- grpreg::grpreg(
    design$A, design$z, group = design$groups, penalty = "grLasso",
    lambda = lambda_grpreg, family = "gaussian"
  )
  theirs <- as.numeric(stats::coef(fit))[-1L]

  residual_at <- function(weight_of) {
    gradient <- drop(crossprod(
      design$A, drop(design$A %*% theirs) - design$z
    )) / design$n
    out <- 0
    for (k in unique(design$groups)) {
      idx <- which(design$groups == k)
      block <- theirs[idx]
      grad <- gradient[idx]
      w <- weight_of(length(idx))
      norm <- sqrt(sum(block^2))
      out <- max(out, if (norm > 1e-12) {
        max(abs(grad + lambda_grpreg * w * block / norm))
      } else {
        max(sqrt(sum(grad^2)) - lambda_grpreg * w, 0)
      })
    }
    out
  }
  expect_lt(residual_at(sqrt), 1e-5)
  expect_gt(residual_at(function(m) 1), 1e-2)
  expect_gt(residual_at(function(m) m), 1e-2)

  # The canonical objective carries n in front of grpreg's lambda, because
  # grpreg states its loss per observation and this package states it as a sum.
  penalty <- list(
    lambda_l1 = 0, lambda_l2 = 0,
    lambda_group = design$n * lambda_grpreg, lambda_exclusive = 0,
    grouped = TRUE
  )
  mine <- .fb_sel_prox_mode_grouped(
    design$A, design$z, penalty, groups, tol = 1e-10
  )
  expect_equal(unname(mine$beta), theirs, tolerance = 1e-5)
  # The supports agree exactly, which is the statement a selection makes.
  expect_identical(unname(abs(mine$beta) > 1e-8), unname(abs(theirs) > 1e-8))
  # Not every coefficient is non-zero, so the agreement above is about a
  # selection and not about two dense vectors.
  expect_true(any(abs(theirs) <= 1e-8))
})


# =============================================================================
# What the user asked for
# =============================================================================

test_that("a group lasso moves a factor's contrast columns together", {
  skip_if_not_installed("glmnet")
  built <- fbgp_substrates()
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fbgp_block("group_lasso", list(scale = 0.3, group = "term")),
    folds = built$folds,
    decision = list(n_lambda = 8L),
    seed = 5L
  )
  expect_identical(path$support_by_resolution$semantics, "group")
  expect_true(path$support_by_resolution$selects)
  grouped_columns <- path$grouping$declared$column[
    path$grouping$declared$term == "f"
  ]
  expect_identical(length(grouped_columns), 2L)
  # `f` has three treatment-coded levels. The reference level is represented
  # by the intercept and has no candidate column; only the two non-reference
  # contrasts belong to the penalised group.
  expect_false(any(grepl("fa$", grouped_columns)))

  # The property: at every tuning point of every fold, the two contrast
  # columns of a factor are both held or both out. A coefficient-wise penalty
  # cannot promise this, which is what G3 below shows.
  coefs <- path$coefficients_original_scale
  candidates <- coefs[coefs$role == "candidate", ]
  for (arm_name in unique(candidates$arm)) {
    for (index in unique(candidates$lambda_index)) {
      for (fold in unique(candidates$fold)) {
        rows <- candidates[
          candidates$arm == arm_name & candidates$lambda_index == index &
            candidates$fold == fold & candidates$term %in% grouped_columns,
        ]
        if (nrow(rows) != 2L || anyNA(rows$estimate)) {
          next
        }
        expect_true(
          all(abs(rows$estimate) > 0) || all(abs(rows$estimate) == 0),
          label = "a factor's contrast columns move together"
        )
      }
    }
  }

  # G3. A coefficient-wise lasso can split two coordinates. This tiny
  # orthogonal problem is a deterministic negative control for the predicate
  # above; using another path on the random fixture made the control depend on
  # whether a coarse lambda grid happened to visit the split interval.
  x_lasso <- diag(2)
  colnames(x_lasso) <- c("contrast_1", "contrast_2")
  split <- .fb_sel_glmnet_fit(
    x_lasso,
    y = c(2, 0.1),
    alpha = 1,
    lambda = 0.2,
    thresh = 1e-14,
    maxit = 1e5
  )
  expect_identical(sum(abs(split) > 1e-8), 1L)
})


test_that("an exclusive lasso thins a correlated cluster from inside", {
  # The user-facing need is that where a cluster of highly associated
  # variables exists, one member enters and the rest stay out. The exclusive
  # penalty encourages that competition; it does NOT guarantee exactly one
  # non-zero member, and the package must not claim it does. What is asserted
  # here is the direction: within an active cluster the exclusive penalty
  # holds strictly fewer members than a ridge of matched strength.
  .fb_sel_with_seed(17L, {
    n <- 90L
    base <- stats::rnorm(n)
    a <- cbind(
      c1 = base + stats::rnorm(n, sd = 0.12),
      c2 = base + stats::rnorm(n, sd = 0.12),
      c3 = base + stats::rnorm(n, sd = 0.12),
      c4 = base + stats::rnorm(n, sd = 0.12)
    )
    a <- scale(a, center = TRUE, scale = FALSE)
    a <- matrix(as.numeric(a), nrow = n, dimnames = list(NULL, colnames(a)))
    z <- 1.5 * a[, 1L] + stats::rnorm(n, sd = 0.3)

    groups <- .fb_sel_group_index(
      colnames(a), rep("cluster", 4L), rule = "block"
    )
    expect_identical(unname(groups$sizes), 4L)

    exclusive <- .fb_sel_prox_mode_grouped(
      a, z,
      .fb_sel_canonical_penalty(
        "exclusive_lasso", list(scale = 0.25, group = "block"), n = n
      ),
      groups, tol = 1e-9
    )
    held <- sum(abs(exclusive$beta) > 1e-8)
    expect_lt(held, 4L)
    expect_gte(held, 1L)
    # A ridge at the same canonical strength holds every member: the thinning
    # is the exclusive penalty's doing, not the data's.
    ridge <- .fb_sel_prox_mode_grouped(
      a, z,
      list(lambda_l1 = 0, lambda_l2 = 1 / 0.25^2, lambda_group = 0,
        lambda_exclusive = 0, grouped = FALSE),
      groups, tol = 1e-9
    )
    expect_identical(sum(abs(ridge$beta) > 1e-8), 4L)
  })
})


# =============================================================================
# Acceptance 3: known-truth recovery curves
# =============================================================================

test_that("recovery curves behave on truths with 0, 1 and several per group", {
  # Six groups of four correlated members. Three truth regimes, as the review
  # of 2026-09-11 required in place of the withdrawn "one member per active
  # group" claim: a group with no active member, one with exactly one, and one
  # with three. The curve is read over a canonical lambda sweep.
  curve <- function(active_per_group, family, hyper_at) {
    .fb_sel_with_seed(23L, {
      n <- 120L
      n_group <- 6L
      per <- 4L
      truth <- numeric(n_group * per)
      cols <- list()
      for (k in seq_len(n_group)) {
        shared <- stats::rnorm(n)
        for (j in seq_len(per)) {
          cols[[length(cols) + 1L]] <- shared + stats::rnorm(n, sd = 0.25)
        }
        if (active_per_group > 0L) {
          idx <- (k - 1L) * per + seq_len(active_per_group)
          truth[idx] <- 1.2
        }
      }
      a <- do.call(cbind, cols)
      colnames(a) <- paste0("g", rep(seq_len(n_group), each = per), "_",
        rep(seq_len(per), times = n_group))
      a <- matrix(as.numeric(scale(a, center = TRUE, scale = FALSE)),
        nrow = n, dimnames = list(NULL, colnames(a)))
      z <- drop(a %*% truth) + stats::rnorm(n, sd = 0.6)
      groups <- .fb_sel_group_index(
        colnames(a), paste0("g", rep(seq_len(n_group), each = per)),
        rule = "term"
      )
      out <- NULL
      # Strongest penalty first: these families state a SCALE, so a small
      # scale is a large lambda. The sweep therefore runs from the strongest
      # penalty to the weakest and the held set grows along it.
      for (scale_value in c(0.004, 0.015, 0.05, 0.2, 1)) {
        solved <- .fb_sel_prox_mode_grouped(
          a, z, .fb_sel_canonical_penalty(family, hyper_at(scale_value), n = n),
          groups, tol = 1e-8
        )
        held <- abs(solved$beta) > 1e-8
        out <- rbind(out, data.frame(
          scale = scale_value,
          true_positive = sum(held & truth != 0) / max(1L, sum(truth != 0)),
          false_positive = sum(held & truth == 0) / max(1L, sum(truth == 0)),
          n_held = sum(held)
        ))
      }
      out
    })
  }

  gl <- curve(1L, "group_lasso", function(s) list(scale = s, group = "term"))
  ex <- curve(1L, "exclusive_lasso",
    function(s) list(scale = s, group = "term"))

  # A curve is a curve: the held set grows as the penalty weakens.
  expect_gt(gl$n_held[nrow(gl)], gl$n_held[1L])
  expect_gt(ex$n_held[nrow(ex)], ex$n_held[1L])
  # At its weakest penalty every truly active coefficient is recovered.
  expect_equal(gl$true_positive[nrow(gl)], 1)
  expect_equal(ex$true_positive[nrow(ex)], 1)
  # The two families differ in a way the support semantics already claim. A
  # group norm is non-smooth at the origin, so a strong enough group lasso
  # holds nothing; the exclusive penalty is smooth there and has no finite
  # null anchor, so however strong it gets, a group keeps its dominant
  # member. This is the `null_anchor` row of .fb_sel_support_semantics(),
  # measured rather than asserted.
  expect_identical(gl$n_held[1L], 0L)
  expect_gt(ex$n_held[1L], 0L)

  # With one active member per group the exclusive penalty, which competes
  # inside a group, admits fewer false positives at matched strength than the
  # group lasso, which moves a whole group at once.
  matched <- merge(gl, ex, by = "scale", suffixes = c("_group", "_excl"))
  expect_true(any(
    matched$false_positive_excl < matched$false_positive_group
  ))

  # A truth with no active member in any group: the strongest penalty holds
  # nothing, so the curve starts where it should.
  empty <- curve(0L, "group_lasso", function(s) list(scale = s, group = "term"))
  expect_identical(empty$n_held[1L], 0L)

  # A truth with three active members per group: the group lasso recovers them
  # all once the group is admitted, which is the behaviour that makes it the
  # right penalty for a factor.
  several <- curve(3L, "group_lasso",
    function(s) list(scale = s, group = "term"))
  expect_equal(several$true_positive[nrow(several)], 1)
})


# =============================================================================
# The fold interaction (P2) and the support history (U-SEL-PATH)
# =============================================================================

test_that("a group partly dropped on a fold is fitted on its survivors", {
  skip_if_not_installed("glmnet")
  # Trial-constant candidates: nine of them across twelve trials, so a fold
  # whose training rows cover eight trials cannot identify all nine and P2's
  # candidate space drops some for that fold.
  #
  # The grouping here is `block` rather than `term`, and that is forced by the
  # tree rather than chosen. A partly dropped FACTOR group is not reachable:
  # a contrast column is dropped only when its level is missing from the
  # training rows, and `.fb_sel_check_fold_levels()` refuses that case with
  # `selection_fold_level_absent` before the candidate space is computed. The
  # reachable partial drop is a numeric candidate aliased on a fold, which is
  # what a declared block of competitors is made of.
  d <- .fb_sel_with_seed(11L, {
    n_loc <- 6L
    trials_per <- 2L
    plots <- 8L
    n_trial <- n_loc * trials_per
    dd <- data.frame(
      trial = factor(rep(sprintf("T%02d", seq_len(n_trial)), each = plots)),
      location = factor(rep(
        rep(sprintf("L%d", seq_len(n_loc)), each = trials_per),
        each = plots
      ))
    )
    n <- nrow(dd)
    tr <- as.integer(dd$trial)
    for (j in seq_len(9L)) {
      dd[[paste0("c", j)]] <- stats::rnorm(n_trial)[tr]
    }
    dd$x <- stats::rnorm(n)
    dd$y <- 1 + 0.5 * dd$x + 0.4 * dd$c1 +
      stats::rnorm(n_trial, sd = 0.5)[tr] + stats::rnorm(n, sd = 0.4)
    dd
  })
  folds <- .fb_sel_folds(fb_resamples(unit = ~location, v = 3L), d, seed = 3L)
  substrates <- lapply(folds, function(split) {
    .fb_sel_fold_substrate(
      d, split, response = "y", protected = ~x,
      candidates = ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9,
      random = "trial"
    )
  })
  names(substrates) <- vapply(folds, function(f) f$label, character(1L))

  dropped_somewhere <- any(vapply(
    substrates,
    function(s) !all(s$candidate_space$declared$column %in% s$candidate_kept),
    logical(1L)
  ))
  expect_true(dropped_somewhere)

  path <- .fb_sel_glmnet_path(
    substrates,
    fbgp_block(
      "group_lasso", list(scale = 0.4, group = "block"), columns = 1:9
    ),
    folds = folds,
    decision = list(benchmark = FALSE)
  )
  membership <- path$support_by_resolution$group_membership
  expect_true(is.data.frame(membership))
  expect_true(all(c(
    "fold", "group", "declared_size", "fitted_size", "weight", "complete"
  ) %in% names(membership)))

  # The defined behaviour: the group is fitted on the members its fold
  # identifies, its weight is recomputed from that size, and the fold says so
  # rather than the drop being invisible.
  partial <- membership[membership$fitted_size < membership$declared_size, ]
  expect_gt(nrow(partial), 0L)
  expect_false(any(partial$complete))
  expect_equal(partial$weight, sqrt(partial$fitted_size))
  complete_rows <- membership[membership$complete, ]
  expect_equal(complete_rows$weight, sqrt(complete_rows$declared_size))
  # A group with no surviving member is absent from the fold, not weight zero.
  expect_true(all(membership$fitted_size > 0L))
})


test_that("grouped support reaches the per-fold support history", {
  skip_if_not_installed("glmnet")
  built <- fbgp_substrates()
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fbgp_block("group_lasso", list(scale = 0.3, group = "term")),
    folds = built$folds,
    decision = list(n_lambda = 6L),
    seed = 9L
  )
  history <- path$support_by_resolution$history
  by_fold <- path$support_by_resolution$history_by_fold
  expect_true(is.data.frame(history))
  expect_gt(nrow(history), 0L)
  expect_true(all(c(
    "fold", "lambda_index", "s", "column", "available", "in_support"
  ) %in% names(history)))

  # The group history is the column history read at the group's resolution,
  # and it carries the same availability statement U-SEL-PATH established:
  # a column a fold could not estimate was never on its path.
  group_history <- path$support_by_resolution$group_history
  expect_true(is.data.frame(group_history))
  expect_true(all(c(
    "fold", "lambda_index", "s", "group", "available", "in_support"
  ) %in% names(group_history)))
  expect_setequal(
    unique(group_history$group),
    unique(built$substrates[[1L]]$candidate_space$declared$term)
  )
  expect_true(all(c("n_points", "n_in_support", "first_entry_s") %in%
    names(by_fold)))
})


# =============================================================================
# The refusals
# =============================================================================

test_that("a grouping that names no group refuses by name", {
  skip_if_not_installed("glmnet")
  built <- fbgp_substrates()
  # Every candidate here is a single-column numeric, so `group = "term"` makes
  # every group a singleton: the group lasso would silently be a lasso and the
  # exclusive lasso silently a ridge. That is the "never guessed" case.
  numeric_only <- lapply(built$substrates, function(s) {
    keep <- s$candidate_space$declared$term %in% paste0("z", 1:4)
    s$x_candidate <- s$x_candidate[, keep, drop = FALSE]
    s$candidate_space$declared <- s$candidate_space$declared[keep, ]
    s$candidate_kept <- intersect(
      s$candidate_kept, s$candidate_space$declared$column
    )
    s
  })
  expect_error(
    .fb_sel_glmnet_path(
      numeric_only,
      fbgp_block(
        "group_lasso", list(scale = 0.4, group = "term"), columns = 1:4
      ),
      folds = built$folds,
      decision = list(benchmark = FALSE)
    ),
    class = "flexybayes_refusal_selection_penalty_group_structure_absent"
  )

  # The same columns under `group = "block"` are one declared group, so the
  # grouping does something and the path runs.
  expect_s3_class(
    .fb_sel_glmnet_path(
      numeric_only,
      fbgp_block(
        "exclusive_lasso", list(scale = 0.4, group = "block"), columns = 1:4
      ),
      folds = built$folds,
      decision = list(benchmark = FALSE)
    ),
    "fb_sel_canonical_path"
  )
})


test_that("a grouped family on the Bayesian route refuses by name", {
  skip_if_not_installed("brms")
  # The grouped families have an optimiser provider and no brms one. Asking
  # for the Bayesian workflow must say that, rather than reaching the compiler
  # and stopping on an unknown route.
  map <- .fb_sel_term_map(~ f + z1 + z2, fbgp_data())
  blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ group_lasso(scale = 1)),
    map,
    protected = character(0)
  )
  expect_identical(blocks[[1L]]$route, "optimizer_only")
  expect_false(blocks[[1L]]$provider_capabilities$brms)
  expect_error(
    .fb_pen_compile_brms(blocks, c("fb", "fc", "z1", "z2"), n_train = 40),
    class = "flexybayes_refusal_selection_penalty_family_optimizer_only"
  )
})


test_that("a grouped path is a prox path and needs no glmnet", {
  # The grouped families are not expressible as a glmnet problem at all, so
  # the certificate names the provider that did the work rather than the one
  # the module is named after.
  built <- fbgp_substrates()
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fbgp_block("exclusive_lasso", list(scale = 0.4, group = "block")),
    folds = built$folds,
    decision = list(benchmark = FALSE)
  )
  expect_identical(path$provider, "prox")
  expect_false(path$fully_bayesian)
  expect_identical(
    path$equivalence_level$level, "map_exact_fixed_covariance"
  )
  # The route is recorded on the objective tables, not on the coordinates.
  # Reading it off `path_coordinates` compares against a column that is not
  # there, which yields `logical(0)`, and `all(logical(0))` is TRUE -- so that
  # spelling of this assertion passed whatever the path had done. Each table
  # is asserted non-empty before its contents are read.
  expect_gt(nrow(path$canonical_objective), 0L)
  expect_true(all(path$canonical_objective$route == "internal_prox"))
  expect_true(all(path$provider_objective$route == "internal_prox"))
  expect_gt(nrow(path$path_coordinates), 0L)
  expect_true(all(is.na(path$path_coordinates$lambda_glmnet)))
  # `glmnet` solved nothing here, so it has no objective to report.
  expect_true(all(is.na(path$provider_objective$provider_objective)))

  # The provider gate is the one binding that asks, so absence can be staged
  # without uninstalling anything.
  testthat::local_mocked_bindings(
    .fb_sel_provider_installed = function(package) FALSE
  )
  expect_s3_class(
    .fb_sel_glmnet_path(
      built$substrates,
      fbgp_block("exclusive_lasso", list(scale = 0.4, group = "block")),
      folds = built$folds,
      decision = list(benchmark = FALSE)
    ),
    "fb_sel_canonical_path"
  )
})


test_that("the certificate records the grouping rule and every group", {
  skip_if_not_installed("glmnet")
  built <- fbgp_substrates()
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fbgp_block("group_exclusive",
      list(group_scale = 0.5, exclusive_scale = 0.8, group = "term")),
    folds = built$folds,
    decision = list(benchmark = FALSE)
  )
  grouping <- path$grouping
  expect_true(is.list(grouping))
  expect_identical(grouping$rule, "term")
  expect_true(nzchar(grouping$rule_description))
  expect_true(is.data.frame(grouping$declared))
  expect_true(all(c("column", "term", "group") %in% names(grouping$declared)))
  # The factor's two contrast columns are one group; each numeric is its own.
  expect_identical(
    grouping$declared$group[grouping$declared$column %in% c("fb", "fc")],
    c("f", "f")
  )
  expect_identical(grouping$weight_rule, "sqrt(group size on the fold)")
  expect_identical(path$family, "group_exclusive")
  expect_equal(path$penalty$lambda_group, 1 / 0.5)
  expect_equal(path$penalty$lambda_exclusive, 1 / 0.8^2)
})


# =============================================================================
# The public route (release integration P8-INT).
#
# The solver-level contract above is necessary but not sufficient: 0.13.0 is
# the first release in which these families are promised to fb_select() users.
# This block therefore starts at the exported verb and is deliberately not a
# second set of proximal-operator tests.
# =============================================================================

fbgp_stub_base <- function(d) {
  fb <- fb_from_asreml(y ~ x, random = ~site, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(
          fixed = y ~ x,
          random = ~site,
          family = "gaussian"
        )
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}


fbgp_public_penalties <- function() {
  list(
    group_lasso = fb_penalty(
      term(everything()) ~ group_lasso(scale = 0.7)
    ),
    sparse_group_lasso = fb_penalty(
      term(everything()) ~ sparse_group_lasso(scale = 0.7, alpha = 0.35)
    ),
    exclusive_lasso = fb_penalty(
      term(everything()) ~ exclusive_lasso(scale = 0.7)
    ),
    group_exclusive = fb_penalty(
      term(everything()) ~ group_exclusive(
        group_scale = 0.7, exclusive_scale = 0.9
      )
    )
  )
}


test_that("all grouped families resolve through the public optimiser plan", {
  d <- fbgp_data(n_cluster = 8L, per = 6L)
  base <- fbgp_stub_base(d)
  for (family in names(fbgp_public_penalties())) {
    plan <- fb_select(
      base,
      candidates = ~ f + g + z1 + z2 + z3 + z4,
      penalty = fbgp_public_penalties()[[family]],
      workflow = "penalized_optimizer",
      resamples = fb_resamples(unit = ~site, v = 2L),
      seed = 91L,
      plan = TRUE
    )
    expect_s3_class(plan, "fb_select_plan")
    expect_identical(plan$selection$provider, "prox")
    expect_identical(plan$selection$prior_route, "optimizer_only")
    expect_false(plan$selection$stanvars)
  }
})


test_that("all grouped families execute through fb_select and certify prox", {
  d <- fbgp_data(n_cluster = 8L, per = 6L)
  base <- fbgp_stub_base(d)
  for (family in names(fbgp_public_penalties())) {
    fit <- fb_select(
      base,
      candidates = ~ f + g + z1 + z2 + z3 + z4,
      penalty = fbgp_public_penalties()[[family]],
      workflow = "penalized_optimizer",
      resamples = fb_resamples(unit = ~site, v = 2L),
      decision = fb_decision(method = "full", metric = "rmse"),
      seed = 91L
    )
    expect_s3_class(fit, "fb_select")
    expect_identical(fit$workflow, "penalized_optimizer")
    expect_identical(fit$provider, "prox")
    expect_identical(fit$benchmark$glmnet$provider, "prox")
    expect_identical(fit$prior, family)
    expect_true(is.data.frame(tidy(fit)))
    expect_identical(nrow(glance(fit)), 1L)
    expect_gt(nrow(fit$benchmark$glmnet$kkt), 0L)
    expect_true(all(is.finite(fit$benchmark$glmnet$kkt$kkt)))
    expect_true(is.data.frame(
      fit$benchmark$glmnet$support$group_history
    ))
  }
})


test_that("a grouped family still refuses on the public Bayesian route", {
  skip_if_not_installed("brms")
  d <- fbgp_data(n_cluster = 8L, per = 6L)
  expect_error(
    fb_select(
      fbgp_stub_base(d),
      candidates = ~ f + g + z1 + z2 + z3 + z4,
      penalty = fbgp_public_penalties()$group_lasso,
      workflow = "bayes",
      resamples = fb_resamples(unit = ~site, v = 2L),
      decision = fb_decision(method = "full", metric = "rmse"),
      seed = 92L,
      plan = TRUE
    ),
    class = "flexybayes_refusal_selection_penalty_family_optimizer_only"
  )
})


test_that("the public grouped route needs no glmnet but remains non-Bayesian", {
  d <- fbgp_data(n_cluster = 8L, per = 6L)
  base <- fbgp_stub_base(d)
  penalty <- fbgp_public_penalties()$group_lasso

  testthat::local_mocked_bindings(
    .fb_sel_provider_installed = function(package) FALSE
  )
  fit <- fb_select(
    base,
    candidates = ~ f + g + z1 + z2 + z3 + z4,
    penalty = penalty,
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 2L),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 92L
  )
  expect_identical(fit$provider, "prox")
  expect_false(fit$benchmark$glmnet$fully_bayesian)

  expect_error(
    fb_select(
      base,
      candidates = ~ f + g + z1 + z2 + z3 + z4,
      penalty = penalty,
      workflow = "bayes",
      resamples = fb_resamples(unit = ~site, v = 2L),
      decision = fb_decision(method = "full", metric = "rmse"),
      seed = 92L,
      plan = TRUE
    ),
    class = "flexybayes_refusal_selection_fallback_not_authorised"
  )
})
