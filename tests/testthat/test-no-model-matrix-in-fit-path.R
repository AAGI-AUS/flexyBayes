# Tests for the the design spec representation acceptance:
# the greta fit-time codegen path must NOT call `model.matrix()`
# for supported indexed term classes. The indexed-emit path binds
# either a per-row integer level-index vector or a parse-time basis
# matrix into the model environment and references it from the
# emitted code; either way, no `model.matrix()` evaluation happens
# inside the fit-time code, so generated-code size stays nearly
# constant in n.
#
# Strategy: trace `flexybayes(..., return_code = TRUE)` output --
# the literal string of greta code that the fit-time pipeline would
# evaluate -- and assert `grepl("model\\.matrix\\(", code)` returns
# FALSE. Equivalent reasoning applies for any other variant of the
# call (`model.matrix(` with parens, `stats::model.matrix(`,
# `base::model.matrix(`).
#
# Coverage cells:
#   (a) simple fixed factor (no interaction)         y ~ g
#   (b) simple random intercept (1|g)                y ~ x + (1|g)
#   (c) smooth s(x)                                  y ~ 1, ~ s(x)
#   (d) factor_numeric_interaction (Stage 1A)        y ~ g:x
#   (e) uncorrelated random slopes (x || g) (1B)     y ~ x + (x || g)
#
# Cells (d) and (e) are skipped on the stage-1c worktree because the
# R/ machinery for those term classes lives on the parallel Stage 1A
# (`R/codegen.R`, `R/parse_formula.R`) and Stage 1B branches; once
# those merge into main, the skips become active subtests.


skip_if_no_mgcv <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mgcv")
}

mk_nomm_data <- function(n = 60L, n_group = 6L, seed = 20260523L) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_group), length.out = n))
  x <- stats::rnorm(n)
  y <- 1 +
    0.4 * x +
    stats::rnorm(n_group, 0, 0.3)[as.integer(g)] +
    stats::rnorm(n, 0, 0.5)
  data.frame(y = y, x = x, g = g)
}

# Centralised assertion. Treats the captured code as a single
# character string; failures print a short excerpt around the
# offending match for diagnostic clarity.
expect_no_model_matrix_in_code <- function(code, label) {
  testthat::expect_type(code, "character")
  testthat::expect_true(
    nzchar(code),
    info = paste0(label, ": codegen returned empty string")
  )
  # Match any of: model.matrix(, stats::model.matrix(, base::model.matrix(.
  pat <- "(stats::|base::)?model\\.matrix\\("
  hit <- regmatches(code, regexpr(pat, code))
  testthat::expect_false(
    grepl(pat, code),
    info = paste0(
      label,
      ": fit-time code contains a model.matrix() ",
      "call (match: '",
      paste(hit, collapse = ", "),
      "')"
    )
  )
}

# ---------------------------------------------------------------- #
# (a) Simple fixed factor                                          #
# ---------------------------------------------------------------- #

test_that("greta fit-time code path emits no model.matrix() call for simple fixed factor", {
  skip_if_no_greta()
  d <- mk_nomm_data()
  code <- suppressMessages(flexybayes(
    fixed = y ~ g,
    data = d,
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_no_model_matrix_in_code(code, "fixed_factor_simple")
})

# ---------------------------------------------------------------- #
# (b) Simple random intercept                                      #
# ---------------------------------------------------------------- #

test_that("greta fit-time code path emits no model.matrix() call for random intercept (1|g)", {
  skip_if_no_greta()
  d <- mk_nomm_data()
  code <- suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_no_model_matrix_in_code(code, "random_intercept")
})

# ---------------------------------------------------------------- #
# (c) Smooth s(x)                                                  #
# ---------------------------------------------------------------- #

test_that("greta fit-time code path emits no model.matrix() call for smooth s(x)", {
  skip_if_no_greta()
  skip_if_no_mgcv()
  d <- mk_nomm_data()
  code <- suppressMessages(flexybayes(
    fixed = y ~ 1,
    random = ~ s(x, k = 6L),
    data = d,
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_no_model_matrix_in_code(code, "smooth")
})

# ---------------------------------------------------------------- #
# (c.1) n-scaling regression: code size stays nearly constant       #
# ---------------------------------------------------------------- #
#
# Companion check on the indexed-emit invariant: code size should
# not blow up linearly in n. We compare code length at n = 100 and
# n = 5000 for the same model and require the larger one to stay
# within a 2x band of the smaller -- the indexed path should give
# near-identical sizes, but we allow some slack because seed-text
# embedding may shift by a few bytes.

test_that("indexed emit keeps return_code size bounded across n for random intercept", {
  skip_if_no_greta()
  d_small <- mk_nomm_data(n = 100L)
  d_big <- mk_nomm_data(n = 5000L)
  code_small <- suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d_small,
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  code_big <- suppressMessages(flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d_big,
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  n_small <- nchar(code_small)
  n_big <- nchar(code_big)
  # Allow 2x because of incidental seed-text or repr drift; the
  # important property is that code size is not linear in n. A
  # dense fall-back at n = 5000 would inline a 5000-row matrix
  # literal that would be orders of magnitude larger.
  expect_lt(n_big, 2L * n_small + 200L)
})

# ---------------------------------------------------------------- #
# (d) factor_numeric_interaction (Stage 1A)                         #
# ---------------------------------------------------------------- #

test_that("greta fit-time code path emits no model.matrix() call for factor:numeric interaction", {
  skip_if_no_greta()
  d <- mk_nomm_data()
  code <- suppressMessages(flexybayes(
    fixed = y ~ g:x,
    data = d,
    return_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_no_model_matrix_in_code(code, "factor_numeric_interaction")
})

# ---------------------------------------------------------------- #
# (e) Uncorrelated random slopes (x || g) (Stage 1B)                #
# ---------------------------------------------------------------- #

test_that("greta fit-time code path emits no model.matrix() call for (x || g)", {
  skip_if_no_greta()
  d <- mk_nomm_data()
  # fb_brms() is the brms-format ingest; the (x || g) syntax lands
  # on that entry point. The expectation is identical: indexed emit,
  # no model.matrix() call in the greta-side codegen.
  code <- suppressMessages(fb(
    y ~ x + (x || g),
    data = d,
    backend = "greta",
    return_code = TRUE,
    verbose = FALSE
  ))
  expect_no_model_matrix_in_code(code, "uncorrelated_random_slopes")
})
