# Tests for the low_rank_smooth approximation path --- ADR 0030 C5 +
# ADR 0027 (v0.4.0 Wave 1 Phase 1B). Covers the rank-K PCA truncation
# engine and the s(x, representation = ...) parse-time interception.
# The scheme was built for a native engine withdrawn entirely in 0.9.3
# (see NEWS.md); neither active engine consumes the truncated basis, so
# the scheme now refuses unconditionally at dispatch
# (low_rank_smooth_unsupported) rather than routing anywhere -- the
# former codegen-substitution and end-to-end fit + predict sections are
# removed, since there is no active engine left to emit or fit them.
# The exact dense smooth path is unchanged --- regression checks
# confirm a plain s(x) carries no approximation.

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
# (e) low_rank_smooth has no active-engine consumer                #
# ---------------------------------------------------------------- #
#
# 0.9.3: the scheme was built for a native engine withdrawn entirely
# (see NEWS.md). Sections (e) "codegen substitutes the truncated
# basis" and (g) "end-to-end fit + predict projection" from the
# pre-0.9.3 file are removed here -- both asserted on that engine's own
# generated-code text or fit shape, and neither active engine has an
# emit path for the truncated basis at all (dispatch refuses before
# either would run). The truncation itself still runs during formula
# parsing (see R/emit_smooth_low_rank.R); section (c) above already
# covers that it is intercepted at parse time, and (a)/(b) cover the
# truncation engine's own numerics directly.

test_that("low_rank smooth refuses on every backend -- no active engine consumes it", {
  skip_if_no_mgcv()
  d <- mk_lr_data()
  for (be in c("inla", "brms", "auto")) {
    err <- tryCatch(
      flexybayes(
        fixed = y ~ 1,
        random = ~ s(
          x,
          representation = list(scheme = "low_rank_smooth", rank = 4L)
        ),
        data = d,
        backend = be,
        verbose = FALSE
      ),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_low_rank_smooth_unsupported_refusal")
    expect_identical(err$reason_code, "low_rank_smooth_unsupported", label = be)
  }
})
