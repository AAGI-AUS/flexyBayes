# =============================================================================
# Exact-limit bridges: the mode oracles of the benchmark contract.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 4.3, 8.1 and
# 8.2, on the whitened substrate of SPEC_penalised_selection_v3_2026-09-11.md
# section 5.1.
#
# Three oracles live in this file and none of them shares a source with the
# code it grades.
#
#   The ridge closed form is checked against two independent dense routes: an
#   LU solve of the same normal equations, and a QR solve of the augmented
#   least-squares system, which is a different derivation rather than a
#   different spelling of the same one.
#
#   The KKT residual is checked at a point known to be optimal and at the
#   same point moved by hand. A residual that is small at the optimum proves
#   nothing on its own -- a function returning zero would pass it -- so the
#   perturbed point is a standing part of the gate.
#
#   The proximal mode is checked by the residual it certifies itself with, on
#   problems chosen to be hard for a first-order method: correlated columns,
#   and more columns than rows.
#
# Injected defects that turn these gates red, all standing in the file:
#   D1  the optimum moved by 0.05 in one coordinate           (KKT gate)
#   D2  the ridge bridge asked for lambda = 1/s instead of    (ridge gate)
#       1/s^2, graded against the independent dense solve
#   D3  the proximal mode given one iteration                 (certificate)
#   D5  a scale-mixture family put to an optimiser            (equivalence)
# =============================================================================

# --- fixtures ---------------------------------------------------------- #

fb_bridge_design <- function(n = 60L, p = 6L, seed = 4L, mean_shift = 2) {
  .fb_sel_with_seed(seed, {
    A <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    colnames(A) <- paste0("a", seq_len(p))
    z <- mean_shift +
      drop(A %*% c(2, -1.2, 0, 0, 0.6, 0)[seq_len(p)]) +
      stats::rnorm(n, sd = 0.7)
    list(A = A, z = z, n = n, p = p)
  })
}


fb_bridge_correlated <- function(n = 70L, seed = 12L) {
  .fb_sel_with_seed(seed, {
    core <- matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L)
    mixed <- core %*%
      matrix(stats::rnorm(9L), nrow = 3L, ncol = 3L) +
      0.03 * matrix(stats::rnorm(n * 3L), nrow = n, ncol = 3L)
    A <- cbind(core, mixed, matrix(stats::rnorm(n * 6L), nrow = n, ncol = 6L))
    colnames(A) <- paste0("c", seq_len(ncol(A)))
    z <- drop(A %*% c(1.5, 0, -1, rep(0, ncol(A) - 3L))) +
      stats::rnorm(n, sd = 0.5)
    list(A = A, z = z, n = n, p = ncol(A))
  })
}


fb_bridge_wide <- function(n = 30L, p = 60L, seed = 33L) {
  .fb_sel_with_seed(seed, {
    A <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    colnames(A) <- paste0("w", seq_len(p))
    z <- drop(A %*% c(2, -1.5, rep(0, p - 2L))) + stats::rnorm(n, sd = 0.4)
    list(A = A, z = z, n = n, p = p)
  })
}


# --- the ridge closed form --------------------------------------------- #

test_that("the ridge bridge matches two independent dense solves", {
  fixture <- fb_bridge_design()
  scale <- 0.75
  lambda_l2 <- 1 / scale^2

  bridged <- .fb_sel_bridge_ridge(fixture$A, fixture$z, scale)

  # Route one: an LU solve of the normal equations. Same equations, different
  # factorisation, and it forms the system matrix explicitly.
  by_lu <- unname(drop(solve(
    crossprod(fixture$A) + diag(lambda_l2, fixture$p),
    crossprod(fixture$A, fixture$z)
  )))

  # Route two: the augmented least-squares system. Ridge is ordinary least
  # squares on [A ; sqrt(lambda) I] against [z ; 0], which is a different
  # derivation rather than a rearrangement of the first.
  augmented <- rbind(fixture$A, diag(sqrt(lambda_l2), fixture$p))
  by_qr <- unname(drop(qr.solve(augmented, c(fixture$z, rep(0, fixture$p)))))

  expect_equal(unname(bridged), by_lu, tolerance = 1e-10)
  expect_equal(unname(bridged), by_qr, tolerance = 1e-10)
  expect_identical(names(bridged), colnames(fixture$A))
})


test_that("the ridge gate is red on the wrong lambda map (defect D2)", {
  # D2. `ridge(scale = s)` is Normal(0, s^2), so the canonical penalty
  # coefficient is 1 / s^2. Reading the prior scale as the penalty -- an easy
  # transcription -- gives a different estimator, and the independent dense
  # solve at the correct lambda has to see it.
  fixture <- fb_bridge_design()
  scale <- 0.75
  correct <- unname(drop(solve(
    crossprod(fixture$A) + diag(1 / scale^2, fixture$p),
    crossprod(fixture$A, fixture$z)
  )))
  defective <- unname(drop(solve(
    crossprod(fixture$A) + diag(1 / scale, fixture$p),
    crossprod(fixture$A, fixture$z)
  )))

  expect_equal(
    unname(.fb_sel_bridge_ridge(fixture$A, fixture$z, scale)),
    correct,
    tolerance = 1e-10
  )
  expect_gt(max(abs(correct - defective)), 1e-3)
})


test_that("the ridge bridge refuses a scale that is not a scale", {
  fixture <- fb_bridge_design()
  expect_error(
    .fb_sel_bridge_ridge(fixture$A, fixture$z, 0),
    "finite positive number",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_bridge_ridge(fixture$A, fixture$z[-1L], 1),
    "must agree on the number of rows",
    fixed = TRUE
  )
})


# --- the KKT residual --------------------------------------------------- #

test_that("the KKT residual is zero at the ridge optimum and large beside it", {
  fixture <- fb_bridge_design()
  scale <- 0.75
  optimum <- .fb_sel_bridge_ridge(fixture$A, fixture$z, scale)

  at_optimum <- .fb_sel_kkt_residual(
    fixture$A,
    fixture$z,
    optimum,
    "ridge",
    list(scale = scale)
  )
  expect_lt(at_optimum, 1e-8)

  # D1. The same point moved by 0.05 in one coordinate. A gate nobody has
  # seen fail is not evidence, so the perturbation is a standing assertion
  # rather than a comment.
  perturbed <- optimum
  perturbed[[2L]] <- perturbed[[2L]] + 0.05
  moved <- .fb_sel_kkt_residual(
    fixture$A,
    fixture$z,
    perturbed,
    "ridge",
    list(scale = scale)
  )
  expect_gt(moved, 1)
  expect_gt(moved / max(at_optimum, .Machine$double.eps), 1e6)
})


test_that("the KKT residual reads the subgradient at a zero coordinate", {
  fixture <- fb_bridge_design()
  # A penalty large enough that the all-zero point is the lasso optimum: it
  # is optimal exactly when max |A'z| is inside lambda_1.
  lambda_l1 <- max(abs(drop(crossprod(fixture$A, fixture$z)))) * 1.5
  scale <- 1 / lambda_l1
  zeroes <- rep(0, fixture$p)

  expect_lt(
    .fb_sel_kkt_residual(
      fixture$A,
      fixture$z,
      zeroes,
      "laplace",
      list(scale = scale)
    ),
    1e-10
  )

  # Below the threshold the same all-zero point is no longer optimal, and the
  # residual is the amount by which the gradient leaves the subdifferential.
  loose <- 1 / (lambda_l1 / 3)
  expect_gt(
    .fb_sel_kkt_residual(
      fixture$A,
      fixture$z,
      zeroes,
      "laplace",
      list(scale = loose)
    ),
    1
  )
})


test_that("the canonical penalty maps each family as section 4.3 fixes it", {
  ridge <- .fb_sel_canonical_penalty("ridge", list(scale = 0.5), n = 40)
  expect_equal(ridge$lambda_l1, 0)
  expect_equal(ridge$lambda_l2, 1 / 0.5^2)

  laplace <- .fb_sel_canonical_penalty("laplace", list(scale = 0.25), n = 40)
  expect_equal(laplace$lambda_l1, 1 / 0.25)
  expect_equal(laplace$lambda_l2, 0)

  elastic <- .fb_sel_canonical_penalty(
    "elastic_net",
    list(lambda = 0.1, alpha = 0.4),
    n = 40
  )
  expect_equal(elastic$lambda_l1, 40 * 0.1 * 0.4)
  expect_equal(elastic$lambda_l2, 40 * 0.1 * 0.6)
  expect_equal(elastic$alpha_composite, 0.4)

  expect_error(
    .fb_sel_canonical_penalty("ridge", list(scale = -1), 10),
    "finite positive number",
    fixed = TRUE
  )
})


# --- P6: the lasso endpoint of elastic_net ------------------------------- #
#
# `alpha = 1` is the lasso corner of the elastic-net path (final
# specification B6): lambda_l2 is exactly zero there, not merely small, and
# the resulting (l1, l2) pair is the same canonical coordinate laplace()
# already reaches at alpha_composite = 1. `alpha = 0` stays refused -- that
# corner already has its own family, ridge().
test_that("elastic_net admits alpha = 1 as an exact pure-lasso point (P6)", {
  lasso_corner <- .fb_sel_canonical_penalty(
    "elastic_net",
    list(lambda = 1, alpha = 1),
    n = 10
  )
  expect_equal(lasso_corner$lambda_l1, 10)
  expect_equal(lasso_corner$lambda_l2, 0)
  expect_equal(lasso_corner$alpha_composite, 1)

  # The same (l1, l2) point a matching laplace() scale reaches: alpha = 1
  # is not a new numeric corner, only a newly admitted name for one already
  # exercised.
  matching_laplace <- .fb_sel_canonical_penalty(
    "laplace",
    list(scale = 1 / 10),
    n = 10
  )
  expect_equal(lasso_corner$lambda_l1, matching_laplace$lambda_l1)
  expect_equal(lasso_corner$lambda_l2, matching_laplace$lambda_l2)
})


test_that("elastic_net's alpha stays refused outside (0, 1] (P6)", {
  expect_error(
    .fb_sel_canonical_penalty("elastic_net", list(lambda = 1, alpha = 0), 10),
    "must be a single number inside (0, 1].",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_canonical_penalty(
      "elastic_net", list(lambda = 1, alpha = 1.01), 10
    ),
    "must be a single number inside (0, 1].",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_canonical_penalty(
      "elastic_net", list(lambda = 1, alpha = -0.1), 10
    ),
    "must be a single number inside (0, 1].",
    fixed = TRUE
  )
})


test_that("sparse_group_lasso's alpha is untouched by the elastic_net widening (P6)", {
  # A falsifier for the fix: if the widening had touched the shared
  # `unit_open` domain instead of adding a family-specific one, this would
  # now admit alpha = 1 too, and sparse_group_lasso(alpha = 1) is a
  # different, already-named family (group_lasso), not a corner of this one.
  expect_error(
    fb_penalty(term(everything()) ~ sparse_group_lasso(scale = 1, alpha = 1)),
    "strictly inside (0, 1)",
    fixed = TRUE
  )
})


# --- the proximal mode --------------------------------------------------- #

test_that("the proximal mode certifies a lasso optimum on a plain design", {
  fixture <- fb_bridge_design()
  # Half of the all-zero threshold, so some coordinates are certain to be at
  # zero and the certificate is exercised on both branches of the
  # subdifferential rather than only on the smooth one.
  null_anchor <- max(abs(drop(crossprod(fixture$A, fixture$z))))
  solved <- .fb_sel_prox_mode(
    fixture$A,
    fixture$z,
    "laplace",
    list(scale = 1 / (0.5 * null_anchor)),
    tol = 1e-10
  )
  expect_true(solved$converged)
  expect_lt(solved$kkt, 1e-8)
  expect_true(any(solved$beta == 0))
  expect_identical(names(solved$beta), colnames(fixture$A))

  # A looser penalty keeps more coordinates, which is the monotone direction
  # a lasso path has to move in.
  loose <- .fb_sel_prox_mode(
    fixture$A,
    fixture$z,
    "laplace",
    list(scale = 1 / (0.02 * null_anchor)),
    tol = 1e-10
  )
  expect_gte(sum(loose$beta != 0), sum(solved$beta != 0))
  expect_lt(loose$kkt, 1e-8)
})


test_that("the proximal mode certifies a lasso optimum on correlated columns", {
  fixture <- fb_bridge_correlated()
  for (scale in c(0.05, 0.2, 1)) {
    solved <- .fb_sel_prox_mode(
      fixture$A,
      fixture$z,
      "laplace",
      list(scale = scale),
      tol = 1e-10
    )
    expect_lt(solved$kkt, 1e-8)
  }
})


test_that("the proximal mode certifies an elastic-net optimum with p > n", {
  fixture <- fb_bridge_wide()
  for (alpha in c(0.3, 0.7)) {
    solved <- .fb_sel_prox_mode(
      fixture$A,
      fixture$z,
      "elastic_net",
      list(lambda = 0.2, alpha = alpha),
      tol = 1e-10
    )
    expect_lt(solved$kkt, 1e-8)
    expect_lt(sum(solved$beta != 0), fixture$p)
  }
  lasso <- .fb_sel_prox_mode(
    fixture$A,
    fixture$z,
    "laplace",
    list(scale = 0.15),
    tol = 1e-10
  )
  expect_lt(lasso$kkt, 1e-8)
})


test_that("the proximal mode refuses rather than returning a near-miss (defect D3)", {
  # D3. One iteration is not enough to reach 1e-12 on this problem, and the
  # point the iteration is standing on is not a mode. Probe B4 found an
  # optimiser reporting convergence 0.16 away from a lasso optimum; the
  # refusal here is what stops this module doing the same.
  fixture <- fb_bridge_correlated()
  expect_error(
    .fb_sel_prox_mode(
      fixture$A,
      fixture$z,
      "laplace",
      list(scale = 0.2),
      tol = 1e-12,
      max_iter = 1L
    ),
    class = "flexybayes_refusal_selection_prox_mode_not_converged"
  )
  condition <- tryCatch(
    .fb_sel_prox_mode(
      fixture$A,
      fixture$z,
      "laplace",
      list(scale = 0.2),
      tol = 1e-12,
      max_iter = 1L
    ),
    flexybayes_refusal_selection_prox_mode_not_converged = function(e) e
  )
  expect_gt(condition$kkt, 1e-12)
  expect_identical(condition$iterations, 1L)
})


test_that("the no-penalty limit is generalised least squares", {
  # Section 8.2's Tier-1 no-penalty limit. With the prior scale taken to
  # infinity the canonical penalty vanishes and the mode is the ordinary
  # least-squares solution on the whitened design -- which is generalised
  # least squares in the original basis, and is what `.fb_sel_gls_fit()`
  # computes by a rank-revealing QR.
  fixture <- fb_bridge_design()
  gls <- unname(drop(qr.solve(fixture$A, fixture$z)))

  by_prox <- .fb_sel_prox_mode(
    fixture$A,
    fixture$z,
    "laplace",
    list(scale = 1e12),
    tol = 1e-9
  )
  expect_equal(unname(by_prox$beta), gls, tolerance = 1e-7)

  by_ridge <- .fb_sel_bridge_ridge(fixture$A, fixture$z, 1e8)
  expect_equal(unname(by_ridge), gls, tolerance = 1e-7)

  # And the module's own generalised-least-squares route, on the same design
  # with an identity covariance, is the same answer by a third algorithm.
  by_module <- .fb_sel_gls_fit(
    chol(diag(1, fixture$n)),
    fixture$A,
    fixture$z
  )
  expect_equal(unname(by_module$beta), gls, tolerance = 1e-8)
})


test_that("the proximal mode handles the corners a provider will not", {
  fixture <- fb_bridge_design(p = 1L)
  one_column <- .fb_sel_prox_mode(
    fixture$A,
    fixture$z,
    "laplace",
    list(scale = 0.3),
    tol = 1e-10
  )
  expect_lt(one_column$kkt, 1e-8)

  empty <- .fb_sel_prox_mode(
    matrix(0, nrow = 10L, ncol = 0L),
    stats::rnorm(10L),
    "laplace",
    list(scale = 1)
  )
  expect_length(empty$beta, 0L)
  expect_true(empty$converged)
})


# --- the equivalence level ----------------------------------------------- #

test_that("the equivalence level table reads as section 8.1 writes it", {
  exact <- .fb_sel_equivalence_level("ridge", "glmnet", TRUE, TRUE)
  expect_identical(exact$level, "map_exact_fixed_covariance")

  for (family in c("laplace", "elastic_net", "none")) {
    expect_identical(
      .fb_sel_equivalence_level(family, "glmnet", TRUE, TRUE)$level,
      "map_exact_fixed_covariance"
    )
  }

  # Integrating the covariance, or leaving the residual scale free, changes
  # the estimand: predictive risk on common folds still compares, coefficients
  # do not -- and "on common folds" is a condition, not a figure of speech,
  # so the middle level needs the shared table stated.
  expect_identical(
    .fb_sel_equivalence_level(
      "ridge", "glmnet", FALSE, TRUE, shared_predictive_table = TRUE
    )$level,
    "posterior_predictive_comparable"
  )
  expect_identical(
    .fb_sel_equivalence_level(
      "laplace", "glmnet", TRUE, FALSE, shared_predictive_table = TRUE
    )$level,
    "posterior_predictive_comparable"
  )

  # D5. A scale-mixture prior put to a penalised optimiser. `glmnet` has no
  # representation of it at all, so no identity is claimed at any setting of
  # the two flags.
  for (fixed in list(c(TRUE, TRUE), c(FALSE, TRUE), c(FALSE, FALSE))) {
    verdict <- .fb_sel_equivalence_level(
      "regularized_horseshoe",
      "glmnet",
      fixed[[1L]],
      fixed[[2L]]
    )
    expect_identical(verdict$level, "no_identity_expected")
    expect_true(grepl("no representation", verdict$reason, fixed = TRUE))
  }

  # The unpenalised benchmark expresses the flat corner and nothing else.
  expect_identical(
    .fb_sel_equivalence_level("none", "gls_reml_plugin", TRUE, TRUE)$level,
    "map_exact_fixed_covariance"
  )
  expect_identical(
    .fb_sel_equivalence_level("laplace", "gls_reml_plugin", TRUE, TRUE)$level,
    "no_identity_expected"
  )
  expect_identical(
    .fb_sel_equivalence_level("laplace", "closed_form", TRUE, TRUE)$level,
    "no_identity_expected"
  )

  expect_error(
    .fb_sel_equivalence_level("ridge", "sparsegl", TRUE, TRUE),
    "should be one of",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_equivalence_level("ridge", "glmnet", NA, TRUE),
    "must be TRUE or FALSE",
    fixed = TRUE
  )
})


test_that("the middle level is not reached without a shared out-of-fold score", {
  # The defect this closes. `posterior_predictive_comparable` says the arms
  # were scored "on the same folds, target, loss, unit and candidate
  # design". The Bayesian penalised route wrote no out-of-fold score at all,
  # so the label was printed about a comparison nobody had made. The flag
  # defaults to FALSE: a caller who has not built the shared table cannot
  # reach the level by forgetting.
  no_table <- .fb_sel_equivalence_level("ridge", "glmnet", FALSE, FALSE)
  expect_identical(no_table$level, "no_identity_expected")
  expect_match(
    no_table$reason, "no common out-of-fold score", fixed = TRUE
  )
  expect_match(no_table$reason, "benchmark = TRUE", fixed = TRUE)
  expect_false(no_table$shared_predictive_table)

  with_table <- .fb_sel_equivalence_level(
    "ridge", "glmnet", FALSE, FALSE, shared_predictive_table = TRUE
  )
  expect_identical(with_table$level, "posterior_predictive_comparable")

  # The flag cannot raise a comparison past what the family supports: a
  # scale mixture glmnet has no representation of stays where it was.
  expect_identical(
    .fb_sel_equivalence_level(
      "regularized_horseshoe", "glmnet", FALSE, FALSE,
      shared_predictive_table = TRUE
    )$level,
    "no_identity_expected"
  )
  # Nor past the exact limit: with both conditions held the level is the
  # Tier-1 one whatever the flag says.
  expect_identical(
    .fb_sel_equivalence_level(
      "ridge", "glmnet", TRUE, TRUE, shared_predictive_table = FALSE
    )$level,
    "map_exact_fixed_covariance"
  )
  expect_error(
    .fb_sel_equivalence_level(
      "ridge", "glmnet", TRUE, TRUE, shared_predictive_table = NA
    ),
    "must be TRUE or FALSE",
    fixed = TRUE
  )
})


# --- the provider lambda map --------------------------------------------- #

test_that("the canonical-to-provider map inverts exactly", {
  convention <- list(
    response_scale_rule = "rms",
    l1_scale_rule = "none",
    response_scale = 2.3
  )
  for (alpha in c(0, 0.25, 0.5, 0.9, 1)) {
    for (lambda in c(0.01, 1, 250)) {
      forward <- .fb_sel_canonical_to_glmnet(lambda, alpha, 80, convention)
      back <- .fb_sel_glmnet_to_canonical(
        forward$lambda,
        forward$alpha,
        80,
        convention
      )
      expect_equal(back$lambda_canonical, lambda, tolerance = 1e-12)
      expect_equal(back$alpha, alpha, tolerance = 1e-12)
    }
  }

  # The L1 corner carries no response scale and the L2 corner carries it
  # once. That asymmetry is the measured convention, and the arithmetic has
  # to show it.
  pure_l1 <- .fb_sel_canonical_to_glmnet(10, 1, 50, convention)
  expect_equal(pure_l1$lambda, 10 / 50)
  expect_equal(pure_l1$alpha, 1)
  pure_l2 <- .fb_sel_canonical_to_glmnet(10, 0, 50, convention)
  expect_equal(pure_l2$lambda, 10 * 2.3 / 50)
  expect_equal(pure_l2$alpha, 0)
})


test_that("the lambda map refuses a convention with no response scale", {
  expect_error(
    .fb_sel_canonical_to_glmnet(
      1,
      0.5,
      10,
      list(response_scale_rule = "rms", l1_scale_rule = "none")
    ),
    "response_scale",
    fixed = TRUE
  )
})


test_that("the response scale rules are what they are named", {
  z <- c(1, 2, 3, 4, 10)
  expect_equal(.fb_sel_response_scale(z, "none"), 1)
  expect_equal(.fb_sel_response_scale(z, "rms"), sqrt(mean(z^2)))
  expect_equal(
    .fb_sel_response_scale(z, "sd_population"),
    sqrt(mean((z - mean(z))^2))
  )
  expect_equal(.fb_sel_response_scale(z, "sd_sample"), stats::sd(z))
  # The rules separate only on a response with a non-zero mean, which is why
  # the convention measurement gives its internal design one.
  expect_gt(
    abs(
      .fb_sel_response_scale(z, "rms") -
        .fb_sel_response_scale(z, "sd_population")
    ),
    1e-3
  )
})
