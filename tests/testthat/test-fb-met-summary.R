# Breeder MET summary (G2). The overall-performance / stability / G x E
# BLUP arithmetic is pinned exactly on a constructed fit with known
# realised-effect draws; the environment correlation and the end-to-end
# pipeline are checked on a real factor-analytic fit where overall
# performance must track the genotype means.

# Build a minimal greta-shaped fit carrying known g_mat (and optional
# Lambda / psi) posterior draws, so fb_met_summary() can be exercised
# without a real MCMC run.
.mk_fake_fa_fit <- function(g_mat_draws, term, lambda = NULL, psi = NULL) {
  # g_mat_draws: n_draws x (n_inner * n_outer) with column names
  # g_mat_<tag>[i,j]. Optional Lambda / psi columns appended.
  mat <- g_mat_draws
  if (!is.null(lambda)) mat <- cbind(mat, lambda)
  if (!is.null(psi)) mat <- cbind(mat, psi)
  draws <- coda::as.mcmc.list(list(coda::as.mcmc(mat)))
  structure(
    list(
      greta = list(draws = draws),
      extras = list(parse_info = list(random = list(term)))
    ),
    class = c("flexybayes", "list")
  )
}

test_that("fb_met_summary() computes OP, stability, and BLUPs exactly", {
  ni <- 4L
  no <- 3L
  tag <- "gen_env_fa2"
  term <- list(
    type = "fa_gxe", inner = "gen", outer = "env", k = 2L,
    n_inner = ni, n_outer = no
  )

  # Known realised effects: genotype i has base level b_i plus an
  # environment pattern, so OP_i is predictable and stability orders the
  # genotypes by how much they swing across environments.
  set.seed(1L)
  n_draws <- 400L
  base <- c(2, 1, 0, -1)
  swing <- c(0.2, 0.2, 2.0, 0.2) # genotype 3 is highly unstable
  cols <- character(0)
  mat <- matrix(NA_real_, n_draws, ni * no)
  col_k <- 1L
  for (j in seq_len(no)) {
    for (i in seq_len(ni)) {
      cols[col_k] <- sprintf("g_mat_%s[%d,%d]", tag, i, j)
      env_shift <- c(-1, 0, 1)[j] * swing[i]
      mat[, col_k] <- base[i] + env_shift + stats::rnorm(n_draws, 0, 0.05)
      col_k <- col_k + 1L
    }
  }
  colnames(mat) <- cols

  fit <- .mk_fake_fa_fit(mat, term)
  ms <- fb_met_summary(fit)

  expect_s3_class(ms, "fb_met_summary")
  # OP is the across-environment mean -> base (the env shifts are
  # symmetric and cancel).
  expect_equal(ms$op$mean, base, tolerance = 0.02)
  # Stability orders genotypes by across-environment swing; genotype 3 is
  # the least stable (largest spread).
  expect_equal(which.max(ms$stability$mean), 3L)
  # G x E BLUP matrix recovers the realised effects.
  expect_equal(dim(ms$gxe_blup), c(ni, no))
  expect_equal(ms$gxe_blup[3L, 3L], base[3L] + swing[3L], tolerance = 0.02)
})

test_that("fb_met_summary() refuses fits without a factor-analytic term", {
  expect_error(fb_met_summary(list()), "flexybayes")
  bad <- structure(
    list(extras = list(parse_info = list(random = list()))),
    class = c("flexybayes", "list")
  )
  expect_error(fb_met_summary(bad), "no factor-analytic")
})

test_that("fb_met_summary() on a real factor-analytic fit tracks genotype means", {
  skip_on_cran()
  skip_if_not_installed("agridat")
  skip_if_not_installed("greta")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  local_tf_seed(2026L)
  data(yan.winterwheat, package = "agridat")
  d <- yan.winterwheat

  fit <- flexybayes(
    yield ~ env, random = ~ fa(env, 2):gen, data = d,
    backend = "greta", n_samples = 500L, warmup = 500L, chains = 2L,
    verbose = FALSE, mcmc_verbose = FALSE
  )
  gen_levels <- levels(d$gen)
  env_levels <- levels(d$env)
  ms <- fb_met_summary(fit,
    genotype_levels = gen_levels, environment_levels = env_levels
  )

  expect_equal(nrow(ms$op), length(gen_levels))
  expect_equal(dim(ms$gxe_blup), c(length(gen_levels), length(env_levels)))

  # Overall performance should rank genotypes like their empirical mean
  # yield (the genotype main effect the realised effects carry).
  gen_mean <- tapply(d$yield, d$gen, mean)[gen_levels]
  rho <- stats::cor(ms$op$mean, as.numeric(gen_mean), method = "spearman")
  expect_gt(rho, 0.5)

  # The environment correlation matrix is a valid correlation matrix.
  expect_equal(dim(ms$env_cor), c(length(env_levels), length(env_levels)))
  expect_equal(unname(diag(ms$env_cor)), rep(1, length(env_levels)),
    tolerance = 1e-6)
  expect_true(all(ms$env_cor >= -1.001 & ms$env_cor <= 1.001))
})

test_that("print.fb_met_summary() renders the headline breeder quantities", {
  ni <- 3L
  no <- 2L
  tag <- "gen_env_fa1"
  term <- list(type = "fa_gxe", inner = "gen", outer = "env", k = 1L,
    n_inner = ni, n_outer = no)
  set.seed(2L)
  mat <- matrix(stats::rnorm(100L * ni * no), 100L, ni * no)
  cols <- character(0)
  ck <- 1L
  for (j in seq_len(no)) for (i in seq_len(ni)) {
    cols[ck] <- sprintf("g_mat_%s[%d,%d]", tag, i, j)
    ck <- ck + 1L
  }
  colnames(mat) <- cols
  ms <- fb_met_summary(.mk_fake_fa_fit(mat, term))
  out <- utils::capture.output(print(ms))
  expect_true(any(grepl("<fb_met_summary>", out, fixed = TRUE)))
  expect_true(any(grepl("overall perform", out)))
})
