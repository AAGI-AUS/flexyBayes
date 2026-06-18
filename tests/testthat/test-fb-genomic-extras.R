# Genomic output contract (G0b) --- the `.fb_genomic_summary()`
# constructor that every genomic / MET fit carries on
# `fit$extras$genomic`. Tested with synthetic draws of known structure
# so the heritability, GEBV, reliability, and marker-effect summaries
# are checked against truth the test authored.

test_that(".fb_genomic_summary() recovers heritability from known-mean draws", {
  set.seed(1L)
  # sigma_g^2 ~ around 3, sigma_e^2 ~ around 1 -> h2 around 0.75.
  var_g <- stats::rgamma(2000L, shape = 100, rate = 100 / 3)
  var_e <- stats::rgamma(2000L, shape = 100, rate = 100 / 1)
  s <- flexyBayes:::.fb_genomic_summary(var_g, var_e)
  expect_s3_class(s, "fb_genomic_summary")
  expect_equal(s$heritability[["mean"]], 0.75, tolerance = 0.03)
  expect_equal(s$genetic_variance[["mean"]], 3, tolerance = 0.2)
  expect_equal(s$residual_variance[["mean"]], 1, tolerance = 0.1)
  expect_true(s$heritability[["q2.5"]] < s$heritability[["q97.5"]])
  expect_equal(s$n_genotypes, 0L)
  expect_null(s$gebv)
})

test_that(".fb_genomic_summary() computes GEBVs and reliability from draws", {
  set.seed(2L)
  n_draws <- 1500L
  n_geno <- 8L
  true_bv <- seq(-2, 2, length.out = n_geno)
  # A precise genotype (small posterior SD) and a noisy one to separate
  # reliability.
  sds <- c(0.1, rep(0.7, n_geno - 2L), 1.4)
  gebv_draws <- vapply(seq_len(n_geno), function(k) {
    stats::rnorm(n_draws, true_bv[k], sds[k])
  }, numeric(n_draws))
  var_g <- stats::rgamma(n_draws, shape = 200, rate = 200 / 2)
  var_e <- stats::rgamma(n_draws, shape = 200, rate = 200 / 2)

  s <- flexyBayes:::.fb_genomic_summary(
    var_g, var_e, gebv_draws = gebv_draws,
    labels = paste0("line", seq_len(n_geno))
  )
  expect_equal(s$n_genotypes, n_geno)
  expect_equal(s$gebv$mean, true_bv, tolerance = 0.15)
  expect_equal(s$gebv$genotype, paste0("line", seq_len(n_geno)))
  # Reliability is in [0, 1] and the precise genotype is more reliable
  # than the noisy one.
  expect_true(all(s$gebv$reliability >= 0 & s$gebv$reliability <= 1))
  expect_gt(s$gebv$reliability[1L], s$gebv$reliability[n_geno])
})

test_that(".fb_genomic_summary() summarises marker effects with retention probability", {
  set.seed(3L)
  n_draws <- 1200L
  n_markers <- 6L
  # Three real effects, three near-zero, with spike-and-slab-like draws.
  centres <- c(1.2, -0.9, 0.6, 0, 0, 0)
  marker_draws <- vapply(seq_len(n_markers), function(k) {
    if (centres[k] == 0) {
      ifelse(stats::runif(n_draws) < 0.9, 0, stats::rnorm(n_draws, 0, 0.3))
    } else {
      stats::rnorm(n_draws, centres[k], 0.2)
    }
  }, numeric(n_draws))
  var_g <- rep(1, n_draws)
  var_e <- rep(1, n_draws)
  s <- flexyBayes:::.fb_genomic_summary(
    var_g, var_e, marker_draws = marker_draws, inclusion_eps = 1e-6
  )
  expect_equal(s$n_markers, n_markers)
  expect_equal(s$marker_effects$mean, centres, tolerance = 0.1)
  # The real effects are retained in (almost) every draw; the spike
  # markers far less often.
  expect_true(all(s$marker_effects$prob_retained[1:3] > 0.95))
  expect_true(all(s$marker_effects$prob_retained[4:6] < 0.5))
})

test_that(".fb_genomic_summary() degenerate draws give NA heritability not NaN", {
  s <- flexyBayes:::.fb_genomic_summary(
    var_g_draws = c(0, 1, 2),
    var_e_draws = c(0, 1, 2)
  )
  # The first draw is 0/0; it must not poison the summary as NaN.
  expect_false(is.nan(s$heritability[["mean"]]))
  expect_true(is.finite(s$heritability[["mean"]]))
})

test_that(".fb_genomic_summary() rejects malformed inputs", {
  expect_error(
    flexyBayes:::.fb_genomic_summary(c(1, 2), c(1, 2, 3)),
    "same length"
  )
  expect_error(
    flexyBayes:::.fb_genomic_summary(c(-1, 2), c(1, 2)),
    "non-negative"
  )
  expect_error(
    flexyBayes:::.fb_genomic_summary(c(1, NA), c(1, 2)),
    "non-finite"
  )
  expect_error(
    flexyBayes:::.fb_genomic_summary(
      c(1, 2), c(1, 2),
      gebv_draws = matrix(0, 3L, 2L)
    ),
    "one row per posterior draw"
  )
})

test_that("print.fb_genomic_summary() renders the headline quantities", {
  s <- flexyBayes:::.fb_genomic_summary(
    var_g_draws = stats::rgamma(500L, 100, 50),
    var_e_draws = stats::rgamma(500L, 100, 100),
    gebv_draws = matrix(stats::rnorm(500L * 3L), 500L, 3L)
  )
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("<fb_genomic_summary>", out, fixed = TRUE)))
  expect_true(any(grepl("heritability", out)))
  expect_true(any(grepl("GEBVs", out)))
})
