# test-known-matrix-alignment.R -- v0.3.8 audit Critical Fix #2
# (known-matrix dim/level alignment). Covers the new validators in
# R/structured_cov.R: .check_known_matrix_dim() +
# .check_known_matrix_dimnames(), and the typed refusal classes they
# raise via .stop_structured_cov_refusal():
#
#   reason_code = "known_matrix_dim_mismatch"
#     matrix dim does not match the grouping factor's level count.
#   reason_code = "known_matrix_dimnames_mismatch"
#     row dimnames != column dimnames (matrix is not square in name
#     space even if square in dim space).
#   reason_code = "known_matrix_level_mismatch"
#     dimnames carry levels not in (or missing from) the factor, OR
#     dimnames carry the correct levels but in a different order
#     than levels(<group>) -- INLA's generic0 + greta's t(chol(V))
#     both require positional alignment.
#
# Plus the three happy paths:
#
#   dimnames-present-and-equal-and-ordered -> validator returns NULL
#   dimnames-absent                         -> validator returns NULL
#                                              (dim check still runs)
#   expected_n = NULL                       -> validator no-ops
#                                              (dispatch could not
#                                               supply the level count)
#
# Tests are at the validator level (no MCMC fit needed) -- the
# refusals fire structurally before any fit attempt.

suppressPackageStartupMessages({
  library(testthat)
})


# ---------------------------------------------------------------- #
# Helpers                                                            #
# ---------------------------------------------------------------- #

# Build a small PD symmetric matrix for the precision tests.
mk_pd_matrix <- function(n, names_ = NULL, seed = 1L) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), nrow = n)
  Q <- crossprod(A) + diag(n) * 0.5
  if (!is.null(names_)) {
    dimnames(Q) <- list(names_, names_)
  }
  Q
}

# Build a small lower-triangular Cholesky factor.
mk_chol_matrix <- function(n, names_ = NULL, seed = 1L) {
  Q <- mk_pd_matrix(n, seed = seed)
  L <- t(chol(Q))
  if (!is.null(names_)) {
    dimnames(L) <- list(names_, names_)
  }
  L
}


# ---------------------------------------------------------------- #
# Dim alignment: refuses when matrix dim != expected_n               #
# ---------------------------------------------------------------- #

test_that("known_matrix_dim_mismatch: precision matrix dim != factor level count", {
  Q <- mk_pd_matrix(5L)
  err <- expect_error(
    flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "g",
      expected_n = 4L,
      fit_levels = NULL
    ),
    "carries 4 levels"
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_identical(err$reason_code, "known_matrix_dim_mismatch")
  expect_identical(err$expected_n, 4L)
  expect_identical(err$actual_dim, c(5L, 5L))
})

test_that("known_matrix_dim_mismatch: chol factor dim != factor level count", {
  L <- mk_chol_matrix(5L)
  err <- expect_error(
    flexyBayes:::.validate_chol_input(
      L,
      name = "L",
      group_var = "g",
      expected_n = 4L,
      fit_levels = NULL
    ),
    "carries 4 levels"
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_identical(err$reason_code, "known_matrix_dim_mismatch")
})


# ---------------------------------------------------------------- #
# Dimnames refusal: row != col names                                 #
# ---------------------------------------------------------------- #

test_that("known_matrix_dimnames_mismatch: rownames != colnames", {
  Q <- mk_pd_matrix(4L)
  dimnames(Q) <- list(c("a", "b", "c", "d"), c("a", "b", "c", "X"))
  err <- expect_error(
    flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "g",
      expected_n = 4L,
      fit_levels = c("a", "b", "c", "d")
    ),
    "different rownames and colnames"
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_identical(err$reason_code, "known_matrix_dimnames_mismatch")
})


# ---------------------------------------------------------------- #
# Level refusal: set-different levels                                #
# ---------------------------------------------------------------- #

test_that("known_matrix_level_mismatch: matrix carries unknown levels", {
  Q <- mk_pd_matrix(4L, names_ = c("a", "b", "c", "Z"))
  err <- expect_error(
    flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "g",
      expected_n = 4L,
      fit_levels = c("a", "b", "c", "d")
    ),
    "do not match levels"
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_identical(err$reason_code, "known_matrix_level_mismatch")
})


# ---------------------------------------------------------------- #
# Level refusal: permuted dimnames                                   #
# ---------------------------------------------------------------- #

test_that("known_matrix_level_mismatch: correct level set but permuted order", {
  # Matrix carries levels in alphabetical order, factor in the
  # natural insertion order. INLA's generic0 requires positional
  # alignment, so this refuses.
  Q <- mk_pd_matrix(4L, names_ = c("a", "b", "c", "d"))
  err <- expect_error(
    flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "g",
      expected_n = 4L,
      fit_levels = c("d", "c", "b", "a")
    ),
    "ordered differently"
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_identical(err$reason_code, "known_matrix_level_mismatch")
  # Refusal message names the perm fix verbatim so the user can
  # copy-paste.
  expect_match(conditionMessage(err), "perm <- match\\(levels")
})


# ---------------------------------------------------------------- #
# Happy paths                                                        #
# ---------------------------------------------------------------- #

test_that("dimnames-present-and-equal-and-ordered: validator returns NULL", {
  Q <- mk_pd_matrix(4L, names_ = c("a", "b", "c", "d"))
  expect_silent(
    invisible(flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "g",
      expected_n = 4L,
      fit_levels = c("a", "b", "c", "d")
    ))
  )
})

test_that("dimnames-absent: dim check passes, level check is no-op (caution surfaces upstream)", {
  Q <- mk_pd_matrix(4L)
  # Validator passes even though we cannot enforce level alignment.
  # The future <fb_plan> object will surface this case as a caution
  # on backend_decision; here we just confirm the validator does
  # not refuse for the dimnames-absent case.
  expect_silent(
    invisible(flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "g",
      expected_n = 4L,
      fit_levels = c("a", "b", "c", "d")
    ))
  )
})

test_that("expected_n = NULL: dim check no-ops (parser-bypass test fixture path)", {
  Q <- mk_pd_matrix(5L)
  # Even with dim != typical level count, no expected_n means no
  # dim refusal.
  expect_silent(
    invisible(flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "g",
      expected_n = NULL,
      fit_levels = NULL
    ))
  )
})


# ---------------------------------------------------------------- #
# Sparse-precision PD path uses Matrix::Cholesky (no dense coerce)   #
# ---------------------------------------------------------------- #

test_that("sparse precision PD probe routes through Matrix::Cholesky (no dense coerce)", {
  testthat::skip_if_not_installed("Matrix")
  # PD sparse via a small symmetric SPD construction. The validator
  # picks the sparse-native branch (inherits 'Matrix') and uses
  # Matrix::Cholesky(forceSymmetric(Q), perm = FALSE) -- no
  # as.matrix() coercion of Q to the dense representation.
  Q <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 4L, 5L),
    j = c(1L, 2L, 3L, 4L, 5L),
    x = c(2, 2, 2, 2, 2),
    dims = c(5L, 5L),
    symmetric = TRUE
  )
  expect_silent(
    invisible(flexyBayes:::.validate_precision_input_pd(
      Q,
      name = "Q",
      group_var = "g"
    ))
  )

  # Non-PD sparse: a symmetric matrix with a negative diagonal entry.
  # Built as a dsCMatrix so the validator's sparseMatrix branch runs;
  # Matrix::Cholesky(perm = FALSE) on the symmetric form refuses
  # because a non-pivoted simplicial Cholesky encounters
  # sqrt(non-positive) at the first negative pivot.
  Q_neg <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 4L, 5L),
    j = c(1L, 2L, 3L, 4L, 5L),
    x = c(-1, 2, 2, 2, 2),
    dims = c(5L, 5L),
    symmetric = TRUE
  )
  err <- expect_error(
    flexyBayes:::.validate_precision_input_pd(
      Q_neg,
      name = "Q_neg",
      group_var = "g"
    ),
    "positive-definite"
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_identical(err$reason_code, "precision_not_positive_definite")
})


# ---------------------------------------------------------------- #
# flexyBayes.trust_pd option skips the PD probe                      #
# ---------------------------------------------------------------- #

test_that("flexyBayes.trust_pd = TRUE skips the PD probe", {
  # A symmetric-but-indefinite dense matrix: eigenvalues 3 and -1.
  # Force double storage so Matrix::isSymmetric does not get tripped
  # up by integer / numeric type mixing.
  Q <- matrix(c(2, 1, 1, 2), nrow = 2L) - 2 * diag(2L)
  storage.mode(Q) <- "double"
  stopifnot(isSymmetric(Q)) # sanity: symmetric

  withr::with_options(
    list(flexyBayes.trust_pd = TRUE),
    expect_silent(
      invisible(flexyBayes:::.validate_precision_input(
        Q,
        name = "Q",
        group_var = "g",
        expected_n = 2L,
        fit_levels = NULL
      ))
    )
  )
  # Without the option, the same matrix refuses.
  withr::with_options(
    list(flexyBayes.trust_pd = FALSE),
    {
      err <- expect_error(
        flexyBayes:::.validate_precision_input(
          Q,
          name = "Q",
          group_var = "g",
          expected_n = 2L,
          fit_levels = NULL
        ),
        "positive-definite"
      )
      expect_identical(err$reason_code, "precision_not_positive_definite")
    }
  )
})
