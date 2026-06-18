# Tests for the Stage 3A foundation: .fb_aggregate_gaussian()
# sufficient-statistics extractor + algebraic log-likelihood
# equivalence (ADR 0022 / v0.3.2).
#
# These tests gate the core algebraic invariant of Stage 3A: the
# aggregated form must produce a bit-exact log-likelihood (to within
# floating-point reassociation tolerance) for every in-scope IR. If
# this property holds, every downstream emit (greta, INLA) can be
# verified against the per-row form via the same identity.
#
# Out of scope for this test file (deferred to the Stage 3A backend
# wiring commit): the greta + INLA emit paths, the dispatch branch
# that consumes the agg_candidate flag, the <flexybayes>$exactness
# slot. This commit lands the math.

test_that(".fb_aggregate_gaussian(): basic shape + counts", {
  set.seed(1L)
  N <- 200L
  J <- 10L
  df <- data.frame(
    y = rnorm(N),
    f = factor(sample(c("a", "b", "c"), N, replace = TRUE)),
    g = factor(sample.int(J, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = df)
  ds <- .fb_dataset(df)
  agg <- .fb_aggregate_gaussian(fb, ds)

  expect_s3_class(agg, "fb_aggregated")
  expect_identical(agg$N, N)
  expect_true(agg$K >= 1L)
  expect_true(agg$K <= N)
  expect_identical(sum(agg$sufficient_stats$n_k), N)
  # S1_k accumulates to the raw sum of y exactly
  expect_lt(abs(sum(agg$sufficient_stats$S1_k) - sum(df$y)), 1e-10)
  # S2_k accumulates to the raw sum of squares exactly
  expect_lt(abs(sum(agg$sufficient_stats$S2_k) - sum(df$y * df$y)), 1e-10)
})

test_that("intercept-only random-intercept model: K equals nlevels(g)", {
  set.seed(2L)
  N <- 100L
  J <- 5L
  df <- data.frame(
    y = rnorm(N),
    g = factor(sample.int(J, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ (1 | g), data = df)
  agg <- .fb_aggregate_gaussian(fb, .fb_dataset(df))

  # With only `(Intercept)` as fixed and `g` as cell-key, K = J
  expect_identical(agg$K, J)
  expect_identical(sum(agg$sufficient_stats$n_k), N)
})

test_that("bit-exact log-lik equivalence: balanced y ~ f + (1|g)", {
  set.seed(3L)
  N <- 500L
  J <- 20L
  df <- data.frame(
    y = rnorm(N),
    f = factor(sample(c("a", "b", "c"), N, replace = TRUE)),
    g = factor(sample.int(J, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = df)
  agg <- .fb_aggregate_gaussian(fb, .fb_dataset(df))

  # Per-row design matrix rebuilt locally for the equivalence check
  # (audit P1.4 -- the aggregator no longer carries the N x p form).
  X_row <- stats::model.matrix(~f, data = df)
  beta <- rnorm(ncol(X_row))
  u <- setNames(rnorm(J, sd = 0.5), levels(df$g))
  sigma <- 0.8

  mu_raw <- as.vector(X_row %*% beta) + u[as.character(df$g)]
  ll_raw <- .gaussian_loglik_raw(df$y, mu_raw, sigma)

  X_cells <- agg$cell_design
  ss <- agg$sufficient_stats
  mu_cells <- as.vector(X_cells %*% beta) +
    u[as.character(ss$g)]
  ll_agg <- .gaussian_loglik_aggregated(
    ss$n_k,
    ss$S1_k,
    ss$S2_k,
    mu_cells,
    sigma
  )

  expect_lt(abs(ll_raw - ll_agg) / abs(ll_raw), 1e-12)
})

test_that("property-based: 20 random in-scope IRs all hit bit-exact equivalence", {
  # Stage 3A's correctness rests on the algebraic identity holding
  # for EVERY in-scope IR + every (beta, u, sigma) triple. Test
  # 20 randomly-generated cases at three N scales.
  set.seed(4L)
  fail_count <- 0L
  max_rel_err <- 0
  for (rep in seq_len(20L)) {
    N <- sample(c(100L, 500L, 2000L), 1L)
    J <- sample(c(5L, 15L, 30L), 1L)
    n_f_levels <- sample(2L:4L, 1L)
    df <- data.frame(
      y = rnorm(N, sd = runif(1L, 0.5, 2)),
      f = factor(sample(letters[1:n_f_levels], N, replace = TRUE)),
      g = factor(sample.int(J, N, replace = TRUE))
    )
    fb <- fb_from_brms(y ~ f + (1 | g), data = df)
    agg <- .fb_aggregate_gaussian(fb, .fb_dataset(df))

    X_row <- stats::model.matrix(~f, data = df)
    beta <- rnorm(ncol(X_row))
    u <- setNames(rnorm(J, sd = 0.5), levels(df$g))
    sigma <- runif(1L, 0.3, 1.5)

    mu_raw <- as.vector(X_row %*% beta) + u[as.character(df$g)]
    ll_raw <- .gaussian_loglik_raw(df$y, mu_raw, sigma)

    X_cells <- agg$cell_design
    ss <- agg$sufficient_stats
    mu_cells <- as.vector(X_cells %*% beta) +
      u[as.character(ss$g)]
    ll_agg <- .gaussian_loglik_aggregated(
      ss$n_k,
      ss$S1_k,
      ss$S2_k,
      mu_cells,
      sigma
    )

    rel_err <- abs(ll_raw - ll_agg) / abs(ll_raw)
    max_rel_err <- max(max_rel_err, rel_err)
    if (rel_err > 1e-12) fail_count <- fail_count + 1L
  }
  expect_identical(fail_count, 0L)
  # Across 20 random IRs at three N scales: max rel error stays
  # comfortably under the 1e-12 acceptance threshold
  expect_lt(max_rel_err, 1e-12)
})

test_that("cell-mean shortcut produces a DIFFERENT log-lik", {
  # The cell-mean shortcut drops the sigma-dependent within-cell
  # sum-of-squares term and biases the residual-variance posterior.
  # Confirm the property-based gate above would catch a cell-mean
  # regression.
  set.seed(5L)
  N <- 300L
  J <- 10L
  df <- data.frame(
    y = rnorm(N),
    f = factor(sample(c("a", "b", "c"), N, replace = TRUE)),
    g = factor(sample.int(J, N, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = df)
  agg <- .fb_aggregate_gaussian(fb, .fb_dataset(df))

  X_row <- stats::model.matrix(~f, data = df)
  beta <- rnorm(ncol(X_row))
  u <- setNames(rnorm(J, sd = 0.5), levels(df$g))
  sigma <- 0.7

  mu_raw <- as.vector(X_row %*% beta) + u[as.character(df$g)]
  ll_raw <- .gaussian_loglik_raw(df$y, mu_raw, sigma)

  X_cells <- agg$cell_design
  ss <- agg$sufficient_stats
  mu_cells <- as.vector(X_cells %*% beta) +
    u[as.character(ss$g)]
  ll_cm <- .gaussian_loglik_cellmean(
    ss$n_k,
    ss$S1_k,
    mu_cells,
    sigma
  )

  # The cell-mean log-lik should differ materially -- well above
  # the property-based gate's 1e-12 threshold
  expect_gt(abs(ll_raw - ll_cm) / abs(ll_raw), 1e-6)
})

test_that("compression ratio reporting: balanced factor design", {
  # 5 factor levels x 10 reps each = 50 rows = 5 cells; compression
  # 10:1
  set.seed(6L)
  df <- data.frame(
    y = rnorm(50L),
    f = factor(rep(c("a", "b", "c", "d", "e"), 10L))
  )
  fb <- fb_from_brms(y ~ f, data = df)
  agg <- .fb_aggregate_gaussian(fb, .fb_dataset(df))

  expect_identical(agg$K, 5L)
  expect_identical(agg$N, 50L)
  expect_lt(abs(agg$compression - 5 / 50), 1e-12)

  printed <- capture.output(print(agg))
  expect_true(any(grepl("compression 10", printed, fixed = TRUE)))
  expect_true(any(grepl("N = 50", printed, fixed = TRUE)))
  expect_true(any(grepl("K = 5", printed, fixed = TRUE)))
})


# ---------------------------------------------------------------- #
# Out-of-scope refusals                                             #
# ---------------------------------------------------------------- #

# Helper for the out-of-scope tests: hand-build a minimal IR with a
# specific term shape (avoids the brms walker's pre-ingest refusals).
.test_make_aggregate_ir <- function(
  family = "gaussian",
  link = "identity",
  fixed_terms = list(),
  random_terms = list()
) {
  structure(
    list(
      response = "y",
      family = family,
      link = link,
      intercept = TRUE,
      fixed_terms = fixed_terms,
      random_terms = random_terms,
      rcov_terms = list(),
      addition_terms = list(),
      priors = list(),
      data_summary = list(n = 100L),
      capabilities = character(),
      source = "test"
    ),
    class = c("fb_terms", "list")
  )
}

test_that("refusal: non-gaussian family", {
  df <- data.frame(y = rpois(50L, 2), x = rnorm(50L))
  fb <- fb_from_brms(y ~ x, data = df, family = "poisson")
  ds <- .fb_dataset(df)

  err <- tryCatch(
    .fb_aggregate_gaussian(fb, ds),
    flexybayes_aggregate_out_of_scope = function(c) c
  )
  expect_s3_class(err, "flexybayes_aggregate_out_of_scope")
  expect_identical(err$reason_code, "non_gaussian_family")
})

test_that("refusal: smooth fixed term", {
  fb <- .test_make_aggregate_ir(
    fixed_terms = list(list(type = "smooth", var = "x", k = 10L))
  )
  df <- data.frame(y = rnorm(50L), x = rnorm(50L))
  ds <- .fb_dataset(df)

  err <- tryCatch(
    .fb_aggregate_gaussian(fb, ds),
    flexybayes_aggregate_out_of_scope = function(c) c
  )
  expect_s3_class(err, "flexybayes_aggregate_out_of_scope")
  expect_identical(err$reason_code, "smooth_term_not_aggregatable")
})

test_that("refusal: uncorrelated random slope (x || g)", {
  set.seed(7L)
  df <- data.frame(
    y = rnorm(50L),
    x = rnorm(50L),
    g = factor(rep(letters[1:5], 10L))
  )
  fb <- fb_from_brms(y ~ (x || g), data = df)
  ds <- .fb_dataset(df)

  err <- tryCatch(
    .fb_aggregate_gaussian(fb, ds),
    flexybayes_aggregate_out_of_scope = function(c) c
  )
  expect_s3_class(err, "flexybayes_aggregate_out_of_scope")
  expect_identical(err$reason_code, "random_slope_not_aggregatable")
})

test_that("refusal: structured-covariance random term (audit P2.8)", {
  fb <- .test_make_aggregate_ir(
    random_terms = list(list(type = "vm", var = "geno", var_n = 50L))
  )
  df <- data.frame(y = rnorm(50L), geno = factor(rep(letters[1:5], 10L)))
  ds <- .fb_dataset(df)

  err <- tryCatch(
    .fb_aggregate_gaussian(fb, ds),
    flexybayes_aggregate_out_of_scope = function(c) c
  )
  expect_s3_class(err, "flexybayes_aggregate_out_of_scope")
  expect_identical(err$reason_code, "structured_random_not_aggregatable")
})

test_that("refusal: metadata-only dataset (no y to sum)", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e7,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(100L)))
  )
  fb <- fb_from_brms(y ~ (1 | g), data = NULL, carry_n_rows = 1e7)

  expect_error(.fb_aggregate_gaussian(fb, ds), regexp = "data-backed")
})


# ---------------------------------------------------------------- #
# Internal-only contract                                            #
# ---------------------------------------------------------------- #

test_that(".fb_aggregate_gaussian + helpers are internal", {
  ns <- asNamespace("flexyBayes")
  exp <- getNamespaceExports(ns)
  expect_false(".fb_aggregate_gaussian" %in% exp)
  expect_false("fb_aggregate_gaussian" %in% exp)
  expect_false(".gaussian_loglik_raw" %in% exp)
  expect_false(".gaussian_loglik_aggregated" %in% exp)
  expect_false(".gaussian_loglik_cellmean" %in% exp)
  # All present internally
  expect_true(exists(".fb_aggregate_gaussian", envir = ns, inherits = FALSE))
  expect_true(exists(".gaussian_loglik_raw", envir = ns, inherits = FALSE))
  expect_true(exists(
    ".gaussian_loglik_aggregated",
    envir = ns,
    inherits = FALSE
  ))
})
