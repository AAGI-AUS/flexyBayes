# Tests for the low_rank_smooth approximation path --- ADR 0030 C5 +
# ADR 0027 (v0.4.0 Wave 1 Phase 1B). Covers the rank-K PCA truncation
# engine, the s(x, representation = ...) parse-time interception, the
# codegen substitution of the truncated basis, the greta-only routing
# guard, and (stress-gated) the end-to-end greta fit + predict
# projection. The exact dense smooth path is unchanged --- regression
# checks confirm a plain s(x) carries no approximation.

skip_if_no_mgcv <- function() skip_if_not_installed("mgcv")
mk_lr_data <- function() {
  set.seed(2026L)
  n <- 60L
  x <- sort(runif(n, 0, 10))
  data.frame(x = x, y = sin(x) + rnorm(n, sd = 0.2))
}

# ---------------------------------------------------------------- #
# (a) truncation engine numerics                                    #
# ---------------------------------------------------------------- #

test_that(".truncate_smooth_basis() is an exact rank-K PCA truncation", {
  skip_if_no_mgcv()
  d <- mk_lr_data()
  sm <- mgcv::smoothCon(
    mgcv::s(x, k = 12L),
    data = d,
    absorb.cons = TRUE,
    scale.penalty = TRUE
  )[[1L]]
  B <- sm$X
  tr <- flexyBayes:::.truncate_smooth_basis(B, rank = 5L, var = "x")

  expect_identical(dim(tr$B_K), c(nrow(B), 5L))
  expect_identical(dim(tr$V_K), c(ncol(B), 5L))
  expect_identical(tr$rank, 5L)
  expect_identical(tr$k, ncol(B))
  # B_K = B V_K
  expect_equal(unname(tr$B_K), unname(B %*% tr$V_K))
  # capture == sum top-K singular^2 / sum all singular^2
  d_sv <- svd(B)$d
  expect_equal(tr$frobenius_capture, sum(d_sv[seq_len(5L)]^2) / sum(d_sv^2))
  expect_true(tr$frobenius_capture > 0 && tr$frobenius_capture <= 1)
})

# ---------------------------------------------------------------- #
# (b) rank refusal contract                                         #
# ---------------------------------------------------------------- #

test_that(".validate_low_rank_rank() enforces the rank contract", {
  # positive integer
  expect_error(
    flexyBayes:::.validate_low_rank_rank(0L, k = 10L, n = 50L),
    class = "flexybayes_low_rank_rank_refusal"
  )
  expect_error(
    flexyBayes:::.validate_low_rank_rank(2.5, k = 10L, n = 50L),
    class = "flexybayes_low_rank_rank_refusal"
  )
  expect_error(
    flexyBayes:::.validate_low_rank_rank(c(1L, 2L), k = 10L, n = 50L),
    class = "flexybayes_low_rank_rank_refusal"
  )
  # ceiling min(k, n)
  err <- tryCatch(
    flexyBayes:::.validate_low_rank_rank(11L, k = 10L, n = 50L),
    flexybayes_low_rank_rank_refusal = function(e) e
  )
  expect_s3_class(err, "flexybayes_low_rank_rank_refusal")
  expect_match(conditionMessage(err), "exceeds the truncation ceiling")
  expect_identical(err$reason_code, "low_rank_rank_exceeds_basis")
  # valid rank returns an integer
  expect_identical(
    flexyBayes:::.validate_low_rank_rank(4, k = 10L, n = 50L),
    4L
  )
})

# ---------------------------------------------------------------- #
# (c) parse-time interception: representation -> approx_spec         #
# ---------------------------------------------------------------- #

test_that("s(x, representation = ...) truncates the basis at parse time", {
  skip_if_no_mgcv()
  d <- mk_lr_data()
  fb <- flexyBayes:::fb_from_asreml(
    fixed = y ~ 1,
    random = ~ s(
      x,
      k = 10L,
      representation = list(scheme = "low_rank_smooth", rank = 4L)
    ),
    data = d
  )
  rt <- fb$random_terms[[1L]]
  expect_identical(rt$type, "smooth_mgcv")
  expect_false(is.null(rt$approx_spec))
  expect_identical(rt$approx_spec$scheme, "low_rank_smooth")
  expect_identical(rt$approx_spec$rank, 4L)
  expect_identical(dim(rt$X_K), c(nrow(d), 4L))
  expect_identical(dim(rt$approx_spec$V_K), c(rt$k, 4L))
  # the full basis is retained for validation
  expect_identical(ncol(rt$X), rt$k)
  expect_true(rt$approx_spec$frobenius_capture > 0)
})

test_that("plain s(x) carries no approximation (regression)", {
  skip_if_no_mgcv()
  d <- mk_lr_data()
  fb <- flexyBayes:::fb_from_asreml(fixed = y ~ 1, random = ~ s(x), data = d)
  rt <- fb$random_terms[[1L]]
  expect_identical(rt$type, "smooth_mgcv")
  expect_null(rt$approx_spec)
  expect_null(rt$X_K)
})

# ---------------------------------------------------------------- #
# (d) parse-time refusals                                           #
# ---------------------------------------------------------------- #

test_that("parse-time refuses bad rank / unknown scheme / bad spec", {
  skip_if_no_mgcv()
  d <- mk_lr_data()
  # rank exceeds ceiling
  expect_error(
    flexyBayes:::fb_from_asreml(
      fixed = y ~ 1,
      random = ~ s(
        x,
        k = 10L,
        representation = list(scheme = "low_rank_smooth", rank = 99L)
      ),
      data = d
    ),
    class = "flexybayes_low_rank_rank_refusal"
  )
  # unknown scheme
  expect_error(
    flexyBayes:::fb_from_asreml(
      fixed = y ~ 1,
      random = ~ s(x, representation = list(scheme = "nope", rank = 3L)),
      data = d
    ),
    regexp = "not a registered approximation scheme"
  )
  # spec without a scheme
  expect_error(
    flexyBayes:::fb_from_asreml(
      fixed = y ~ 1,
      random = ~ s(x, representation = list(rank = 3L)),
      data = d
    ),
    class = "flexybayes_approximation_spec_invalid"
  )
})

# ---------------------------------------------------------------- #
# (e) codegen substitutes the truncated basis (dim = K)             #
# ---------------------------------------------------------------- #

test_that("codegen emits dim = K on the low-rank path", {
  skip_if_no_greta()
  skip_if_no_mgcv()
  d <- mk_lr_data()
  code <- flexybayes(
    fixed = y ~ 1,
    random = ~ s(
      x,
      k = 10L,
      representation = list(scheme = "low_rank_smooth", rank = 4L)
    ),
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_match(code, "s_x_raw <- normal\\(0, 1, dim = 4\\)")
  expect_match(code, "B_g_s_x <- as_data\\(B_s_x\\)")
})

# ---------------------------------------------------------------- #
# (f) greta-only routing guard                                      #
# ---------------------------------------------------------------- #

test_that("low_rank smooth refuses on explicit inla / brms backends", {
  skip_if_no_mgcv()
  d <- mk_lr_data()
  expect_error(
    flexybayes(
      fixed = y ~ 1,
      random = ~ s(
        x,
        representation = list(scheme = "low_rank_smooth", rank = 4L)
      ),
      data = d,
      backend = "inla",
      verbose = FALSE
    ),
    class = "flexybayes_low_rank_requires_greta"
  )
})

# ---------------------------------------------------------------- #
# (g) end-to-end greta fit + predict projection (stress-gated)      #
# ---------------------------------------------------------------- #
#
# The full MCMC fit is gated behind FLEXYBAYES_RUN_STRESS to keep the
# routine tally fast; it exercises the exactness label, the V_K slot,
# validate_approximation() on a real fit, and the predict-side V_K
# projection (the smooth must contribute, not collapse to flat zero).

test_that("end-to-end low-rank greta fit, validation, and prediction", {
  skip_if_no_greta()
  skip_if_no_mgcv()
  if (!identical(Sys.getenv("FLEXYBAYES_RUN_STRESS"), "true")) {
    skip("stress-gated end-to-end greta fit (set FLEXYBAYES_RUN_STRESS=true)")
  }

  d <- mk_lr_data()
  fit <- flexybayes(
    fixed = y ~ 1,
    random = ~ s(
      x,
      k = 10L,
      representation = list(scheme = "low_rank_smooth", rank = 4L)
    ),
    data = d,
    backend = "greta",
    n_samples = 200L,
    warmup = 200L,
    chains = 2L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )

  expect_identical(fit$exactness, "approximate_low_rank_smooth")
  ap <- fit$extras$parse_info$approx
  expect_false(is.null(ap$x))
  expect_identical(dim(ap$x$V_K), c(ap$x$k, 4L))

  v <- validate_approximation(fit)
  expect_s3_class(v, "fb_approximation_validation")
  expect_identical(v$scheme, "low_rank_smooth")
  expect_true(
    v$per_smooth$x$frobenius_capture > 0 &&
      v$per_smooth$x$frobenius_capture <= 1
  )

  nd <- data.frame(x = seq(0.5, 9.5, length.out = 25L))
  pr <- predict(fit, newdata = nd)
  pe <- as.numeric(
    if (is.list(pr)) {
      (pr$fit %||% pr$prediction %||% pr[[1L]])
    } else {
      pr
    }
  )
  expect_length(pe, nrow(nd))
  expect_true(all(is.finite(pe)))
  # the projected smooth contributes -- not the silent flat-zero path
  expect_gt(stats::sd(pe), 1e-6)
})
