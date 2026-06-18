# Streaming out-of-core aggregation: the accumulated sufficient
# statistics must equal a single-pass aggregation exactly (up to
# floating-point summation order), be invariant to the chunk boundary,
# and refuse the same out-of-scope models the in-memory path refuses.

test_that("streamed gaussian sufficient stats match the in-memory path", {
  set.seed(101L)
  n <- 5e4L
  df <- data.frame(
    env = factor(sample(letters[1:5], n, replace = TRUE)),
    geno = factor(sample(1:30, n, replace = TRUE)),
    y = stats::rnorm(n, 1, 2)
  )

  ds <- .fb_dataset(df)
  fb <- fb_from_asreml(y ~ env, random = ~geno, data = df, family = "gaussian")
  ref <- .fb_aggregate_gaussian(fb, ds)

  agg <- flexybayes_stream(
    y ~ env,
    random = ~geno,
    source = df,
    family = "gaussian",
    chunk_rows = 7000,
    fit = FALSE,
    verbose = FALSE
  )

  expect_s3_class(agg, "fb_aggregated")
  expect_true(isTRUE(agg$streamed))
  expect_identical(agg$K, ref$K)
  expect_identical(as.numeric(agg$N), as.numeric(ref$N))

  key <- c(ref$fixed_cols, ref$random_cols)
  rs <- data.table::as.data.table(ref$sufficient_stats)
  as_ <- data.table::as.data.table(agg$sufficient_stats)
  data.table::setorderv(rs, key)
  data.table::setorderv(as_, key)

  expect_equal(as.numeric(rs$n_k), as_$n_k)
  expect_lt(max(abs(rs$S1_k - as_$S1_k) / pmax(abs(rs$S1_k), 1)), 1e-8)
  expect_lt(max(abs(rs$S2_k - as_$S2_k) / pmax(abs(rs$S2_k), 1)), 1e-8)
  expect_equal(sum(as_$n_k), as.numeric(n))
})

test_that("the streamed aggregation is invariant to the chunk boundary", {
  set.seed(202L)
  n <- 3e4L
  df <- data.frame(
    a = factor(sample(1:4, n, replace = TRUE)),
    b = factor(sample(1:10, n, replace = TRUE)),
    y = stats::rnorm(n)
  )
  one <- flexybayes_stream(
    y ~ a,
    random = ~b,
    source = df,
    chunk_rows = n,
    fit = FALSE,
    verbose = FALSE
  )
  many <- flexybayes_stream(
    y ~ a,
    random = ~b,
    source = df,
    chunk_rows = 4321,
    fit = FALSE,
    verbose = FALSE
  )
  key <- c(one$fixed_cols, one$random_cols)
  o <- data.table::as.data.table(one$sufficient_stats)
  m <- data.table::as.data.table(many$sufficient_stats)
  data.table::setorderv(o, key)
  data.table::setorderv(m, key)
  expect_identical(o$n_k, m$n_k)
  expect_equal(o$S1_k, m$S1_k, tolerance = 1e-8)
  expect_equal(o$S2_k, m$S2_k, tolerance = 1e-8)
})

test_that("count-family sufficient statistics are exact additive sums", {
  # `g` is a RANDOM grouping factor so its raw level is retained as a
  # cell-key column on the sufficient-statistics table (a fixed factor
  # is stored as model-matrix contrast columns instead). The reference
  # per-cell sums use base tapply / table to stay backend-independent.
  set.seed(303L)
  n <- 2e4L
  df <- data.frame(
    g = factor(sample(1:8, n, replace = TRUE)),
    y = stats::rpois(n, 3),
    m = sample(1:5, n, replace = TRUE)
  )
  ap <- flexybayes_stream(
    y ~ 1,
    random = ~g,
    source = df,
    family = "poisson",
    exposure = "m",
    chunk_rows = 2500,
    fit = FALSE,
    verbose = FALSE
  )
  ss_p <- as.data.frame(ap$sufficient_stats)
  ss_p <- ss_p[order(ss_p$g), ]
  ref_cnt <- tapply(df$y, df$g, sum)
  ref_e <- tapply(df$m, df$g, sum)
  expect_equal(ss_p$count_k, as.numeric(ref_cnt[as.character(ss_p$g)]))
  expect_equal(ss_p$expo_k, as.numeric(ref_e[as.character(ss_p$g)]))

  set.seed(404L)
  db <- data.frame(
    g = factor(sample(1:6, n, replace = TRUE)),
    y = stats::rbinom(n, 1, 0.4)
  )
  ab <- flexybayes_stream(
    y ~ 1,
    random = ~g,
    source = db,
    family = "binomial",
    chunk_rows = 3000,
    fit = FALSE,
    verbose = FALSE
  )
  ss_b <- as.data.frame(ab$sufficient_stats)
  ss_b <- ss_b[order(ss_b$g), ]
  ref_s <- tapply(db$y, db$g, sum)
  ref_tr <- as.numeric(table(db$g))
  names(ref_tr) <- names(table(db$g))
  expect_equal(ss_b$succ_k, as.numeric(ref_s[as.character(ss_b$g)]))
  expect_equal(ss_b$trials_k, as.numeric(ref_tr[as.character(ss_b$g)]))
})

test_that("an .fst source streams and reports the same K and N", {
  skip_if_not_installed("fst")
  set.seed(505L)
  n <- 4e4L
  df <- data.frame(
    env = factor(sample(letters[1:6], n, replace = TRUE)),
    geno = factor(sample(1:20, n, replace = TRUE)),
    y = stats::rnorm(n)
  )
  path <- withr::local_tempfile(fileext = ".fst")
  fst::write_fst(df, path)

  from_df <- flexybayes_stream(
    y ~ env,
    random = ~geno,
    source = df,
    chunk_rows = 9999,
    fit = FALSE,
    verbose = FALSE
  )
  from_fst <- flexybayes_stream(
    y ~ env,
    random = ~geno,
    source = path,
    chunk_rows = 9999,
    fit = FALSE,
    verbose = FALSE
  )
  expect_identical(from_fst$K, from_df$K)
  expect_identical(as.numeric(from_fst$N), as.numeric(from_df$N))

  key <- c(from_df$fixed_cols, from_df$random_cols)
  a <- data.table::as.data.table(from_df$sufficient_stats)
  b <- data.table::as.data.table(from_fst$sufficient_stats)
  data.table::setorderv(a, key)
  data.table::setorderv(b, key)
  expect_equal(a$S1_k, b$S1_k, tolerance = 1e-8)
})

test_that("the streamed row total is overflow-safe (carried as a double)", {
  # A partitioned dataset can exceed the 2^31 integer ceiling. The total
  # row count and the per-cell counts must therefore be doubles (exact to
  # 2^53), not integers -- this test pins the type contract without
  # materialising billions of rows.
  set.seed(909L)
  chunks <- list(
    data.frame(g = factor(c(1, 2, 1, 2)), y = stats::rnorm(4)),
    data.frame(g = factor(c(2, 1, 2, 1)), y = stats::rnorm(4))
  )
  gen <- function(i) if (i <= length(chunks)) chunks[[i]] else NULL

  agg <- flexybayes_stream(
    y ~ 1,
    random = ~g,
    source = gen,
    family = "gaussian",
    fit = FALSE,
    verbose = FALSE
  )
  expect_type(agg$N, "double")
  expect_type(agg$sufficient_stats$n_k, "double")
  expect_equal(agg$N, 8)
})

test_that("out-of-scope models are refused before any fit", {
  set.seed(606L)
  n <- 1e3L
  df <- data.frame(
    x = stats::rnorm(n),
    geno = factor(sample(1:5, n, replace = TRUE)),
    y = stats::rnorm(n)
  )
  # Continuous fixed effect -> one cell per row, refused.
  expect_error(
    flexybayes_stream(
      y ~ x,
      random = ~geno,
      source = df,
      fit = FALSE,
      verbose = FALSE
    ),
    "continuous|out of scope|cell"
  )
  # Intercept-only -> no factor cell key, refused.
  expect_error(
    flexybayes_stream(y ~ 1, source = df, fit = FALSE, verbose = FALSE),
    "factor cell-key|does not compress"
  )
})
