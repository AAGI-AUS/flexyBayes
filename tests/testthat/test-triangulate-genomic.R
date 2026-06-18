# Genomic triangulation: triangulate_genomic() (heritability / variance
# components / breeding values across engines or against a field lens) and
# triangulate_gwas() (hit-set / effect agreement). The lens-list path is
# pinned exactly; the cross-backend path is checked on real fits.

# ---------------------------------------------------------------- #
# (a) triangulate_genomic() on generic lenses.                     #
# ---------------------------------------------------------------- #

test_that("triangulate_genomic() compares a Bayesian posterior against a REML point lens", {
  set.seed(1L)
  # Lens A: a Bayesian posterior (draws). Lens B: a REML point + SE.
  lens_bayes <- list(
    h2 = stats::rnorm(2000L, 0.55, 0.06),
    var_g = stats::rgamma(2000L, 100, 50),
    var_e = stats::rgamma(2000L, 100, 100),
    gebv = stats::setNames(c(2, 1, 0, -1, -2), paste0("g", 1:5)),
    label = "bayes"
  )
  lens_reml <- list(
    h2 = list(estimate = 0.52, se = 0.05),
    var_g = list(estimate = 2.0, se = 0.3),
    var_e = list(estimate = 1.0, se = 0.15),
    gebv = stats::setNames(c(1.9, 1.1, 0.1, -0.9, -2.1), paste0("g", 1:5)),
    label = "sommer"
  )

  tg <- triangulate_genomic(lens_bayes, lens_reml)
  expect_s3_class(tg, "triangulate_genomic_result")
  expect_equal(tg$label_a, "bayes")
  expect_equal(tg$label_b, "sommer")

  h2_row <- tg$components[tg$components$quantity == "heritability", ]
  expect_equal(h2_row$value_a, 0.55, tolerance = 0.02)
  expect_equal(h2_row$value_b, 0.52, tolerance = 1e-9)
  expect_true(h2_row$intervals_overlap) # 0.55 and 0.52 agree

  # GEBVs match by genotype label and are highly correlated.
  expect_equal(tg$gebv$n_common, 5L)
  expect_gt(tg$gebv$pearson, 0.98)
})

test_that("triangulate_genomic() flags disjoint heritability intervals", {
  lens_a <- list(h2 = list(estimate = 0.20, se = 0.03),
    var_g = list(estimate = 0.5, se = 0.1),
    var_e = list(estimate = 2, se = 0.2), gebv = NULL, label = "a")
  lens_b <- list(h2 = list(estimate = 0.70, se = 0.03),
    var_g = list(estimate = 3, se = 0.2),
    var_e = list(estimate = 1, se = 0.1), gebv = NULL, label = "b")
  tg <- triangulate_genomic(lens_a, lens_b)
  h2_row <- tg$components[tg$components$quantity == "heritability", ]
  expect_false(h2_row$intervals_overlap)
  expect_equal(h2_row$difference, -0.5, tolerance = 1e-9)
})

test_that("triangulate_genomic() carries the shared-upstream caveat unless declared independent", {
  lens_a <- list(h2 = list(estimate = 0.5, se = 0.05), gebv = NULL, label = "a")
  lens_b <- list(h2 = list(estimate = 0.5, se = 0.05), gebv = NULL, label = "b")
  expect_false(is.na(triangulate_genomic(lens_a, lens_b)$shared_upstream_caveat))
  expect_true(is.na(
    triangulate_genomic(lens_a, lens_b, data_independence = TRUE)$shared_upstream_caveat
  ))
})

test_that("triangulate_genomic() requires named GEBVs in a lens", {
  lens_a <- list(h2 = list(estimate = 0.5, se = 0.05),
    gebv = c(1, 2, 3), label = "a")
  lens_b <- list(h2 = list(estimate = 0.5, se = 0.05), gebv = NULL, label = "b")
  expect_error(triangulate_genomic(lens_a, lens_b), "named")
})

# ---------------------------------------------------------------- #
# (b) triangulate_gwas() on GWAS lenses.                           #
# ---------------------------------------------------------------- #

test_that("triangulate_gwas() reports hit-set Jaccard, top-K overlap, and effect correlation", {
  mk <- function(sig_idx, effects) {
    p <- rep(0.5, 20L)
    p[sig_idx] <- 1e-10
    data.frame(
      marker = paste0("snp", 1:20),
      p_value = p,
      p_bonferroni = pmin(1, p * 20L),
      effect = effects,
      stringsAsFactors = FALSE
    )
  }
  set.seed(2L)
  eff <- stats::rnorm(20L)
  a <- list(results = mk(c(3L, 7L, 12L), eff), lambda_gc = 1.02)
  # b shares 2 of 3 hits and has correlated effects.
  b <- list(results = mk(c(3L, 7L, 18L), eff + stats::rnorm(20L, 0, 0.1)),
    lambda_gc = 0.98)

  tg <- triangulate_gwas(a, b, top_k = 5L)
  expect_s3_class(tg, "triangulate_gwas_result")
  expect_equal(tg$n_sig_a, 3L)
  expect_equal(tg$n_sig_b, 3L)
  expect_equal(tg$n_sig_common, 2L)
  expect_equal(tg$jaccard, 2 / 4, tolerance = 1e-9) # |{3,7}| / |{3,7,12,18}|
  expect_gt(tg$effect_correlation, 0.95)
  expect_equal(tg$lambda_gc_a, 1.02)
})

# ---------------------------------------------------------------- #
# (c) Real cross-backend genomic triangulation.                    #
# ---------------------------------------------------------------- #

test_that("triangulate_genomic() agrees across greta and brms GBLUP fits", {
  skip_on_cran()
  skip_if_not_installed("greta")
  skip_if_not_installed("brms")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  local_tf_seed(2026L)
  G <- sim_kinship(n_geno = 30L, n_markers = 300L, seed = 31L)
  sim <- sim_gblup_pheno(G, var_g = 1, var_e = 1, n_rep = 4L, seed = 32L)
  dat <- sim$data

  fit_greta <- flexybayes(
    y ~ 1, random = ~ vm(geno, Gmat), data = dat,
    known_matrices = list(Gmat = G), backend = "greta",
    n_samples = 600L, warmup = 600L, chains = 2L,
    verbose = FALSE, mcmc_verbose = FALSE
  )
  fit_brms <- flexybayes(
    y ~ 1, random = ~ vm(geno, Gmat), data = dat,
    known_matrices = list(Gmat = G), backend = "brms",
    n_samples = 600L, warmup = 600L, chains = 2L,
    verbose = FALSE, mcmc_verbose = FALSE
  )

  tg <- triangulate_genomic(fit_greta, fit_brms)
  h2_row <- tg$components[tg$components$quantity == "heritability", ]
  expect_true(h2_row$intervals_overlap)
  expect_gt(tg$gebv$n_common, 25L)
  expect_gt(tg$gebv$pearson, 0.8)
})
