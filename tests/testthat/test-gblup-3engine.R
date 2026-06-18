# GBLUP across three engines (G1). The genomic relationship random
# effect vm(geno, G) now reaches all three backends: greta and brms via
# the dense covariance (brms's native gr(cov = K) = the K = L L'
# decorrelation), INLA via the precision carrier (generic0). GBLUP is
# therefore three-engine triangulatable. The fast tests pin the emit /
# gate behaviour; the gated recovery test fits all three engines on a
# simulated known-heritability dataset and checks recovery + agreement.

# ---------------------------------------------------------------- #
# (a) Fast: brms vm() emit + the relaxed capability gate.           #
# ---------------------------------------------------------------- #

.gblup_term <- function(var = "geno", fmt = "dense", sym = "Gmat") {
  list(
    type = "vm", var = var, mat = sym,
    cov_representation = list(format = fmt, data = sym, scheme = NULL)
  )
}

test_that(".fb_brms_data2() materialises each carrier to the relationship covariance", {
  set.seed(1L)
  M <- matrix(stats::rnorm(36L), 6L, 6L)
  G <- crossprod(M) + diag(6L)
  dat <- data.frame(geno = factor(paste0("g", rep(1:6, 2L))))

  fb_dense <- list(random_terms = list(.gblup_term(fmt = "dense")))
  d2 <- flexyBayes:::.fb_brms_data2(fb_dense, list(Gmat = G), dat)
  expect_equal(unname(d2$Gmat), G, ignore_attr = TRUE)

  L <- t(chol(G))
  fb_chol <- list(random_terms = list(.gblup_term(fmt = "chol", sym = "Lc")))
  d2c <- flexyBayes:::.fb_brms_data2(fb_chol, list(Lc = L), dat)
  expect_equal(unname(d2c$Lc), G, tolerance = 1e-8, ignore_attr = TRUE)

  Q <- solve(G)
  fb_prec <- list(random_terms = list(.gblup_term(fmt = "precision", sym = "Qp")))
  d2p <- flexyBayes:::.fb_brms_data2(fb_prec, list(Qp = Q), dat)
  expect_equal(unname(d2p$Qp), G, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that(".fb_brms_covname() refuses the greta/INLA-only carriers on brms", {
  expect_error(
    flexyBayes:::.fb_brms_covname(.gblup_term(fmt = "blocks")),
    "greta / INLA-only"
  )
  expect_error(
    flexyBayes:::.fb_brms_covname(.gblup_term(fmt = "low_rank")),
    "greta / INLA-only"
  )
  expect_equal(flexyBayes:::.fb_brms_covname(.gblup_term(fmt = "dense")), "Gmat")
})

test_that("brms backend now generates GBLUP Stan code via the known-covariance route", {
  skip_if_not_installed("brms")
  set.seed(2L)
  G <- {
    Z <- matrix(stats::rbinom(8L * 40L, 2L, 0.3), 8L, 40L)
    Zc <- scale(Z, scale = FALSE)
    tcrossprod(Zc) / 40 + diag(8L) * 1e-2
  }
  dimnames(G) <- list(paste0("g", 1:8), paste0("g", 1:8))
  dat <- data.frame(
    geno = factor(rep(paste0("g", 1:8), 3L)),
    y = stats::rnorm(24L)
  )
  code <- flexybayes(
    y ~ 1, random = ~ vm(geno, Gmat), data = dat,
    known_matrices = list(Gmat = G), backend = "brms",
    return_code = TRUE, verbose = FALSE
  )
  # The Cholesky-of-known-covariance group effect is the GBLUP structure.
  expect_true(grepl("Lcov", code, fixed = TRUE))
})

test_that("the capability gate allows dense vm/ped on brms but still refuses fa and block carriers", {
  expect_true(flexyBayes:::.capability_brms(
    list(random_terms = list(.gblup_term(fmt = "dense")))
  ))
  expect_true(flexyBayes:::.capability_brms(
    list(random_terms = list(.gblup_term(fmt = "precision", sym = "Qp")))
  ))
  expect_identical(
    flexyBayes:::.capability_brms(
      list(random_terms = list(.gblup_term(fmt = "blocks")))
    ),
    "stan_cannot_represent_structured_cov"
  )
  expect_identical(
    flexyBayes:::.capability_brms(
      list(random_terms = list(list(type = "fa_gxe", var = "env")))
    ),
    TRUE # fa_gxe is gated at emit, not here (only fa/us/ar1 are listed)
  )
})

test_that("genomic_summary() refuses a non-flexybayes object", {
  expect_error(genomic_summary(list()), "flexybayes object")
})

# ---------------------------------------------------------------- #
# (b) Gated: three-engine recovery on a simulated GBLUP.            #
# ---------------------------------------------------------------- #

test_that("GBLUP recovers heritability and breeding values across three engines", {
  skip_on_cran()
  skip_if_not_installed("greta")
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  local_tf_seed(2026L)

  G <- sim_kinship(n_geno = 30L, n_markers = 300L, seed = 21L)
  sim <- sim_gblup_pheno(G, var_g = 1, var_e = 1, n_rep = 4L, seed = 22L)
  true_h2 <- sim$h2 # 0.5
  dat <- sim$data
  Q <- solve(G)
  dimnames(Q) <- dimnames(G)

  fit_greta <- flexybayes(
    y ~ 1, random = ~ vm(geno, Gmat), data = dat,
    known_matrices = list(Gmat = G), backend = "greta",
    n_samples = 700L, warmup = 700L, chains = 2L,
    verbose = FALSE, mcmc_verbose = FALSE
  )
  fit_brms <- flexybayes(
    y ~ 1, random = ~ vm(geno, Gmat), data = dat,
    known_matrices = list(Gmat = G), backend = "brms",
    n_samples = 700L, warmup = 700L, chains = 2L,
    verbose = FALSE, mcmc_verbose = FALSE
  )
  fit_inla <- flexybayes(
    y ~ 1, random = ~ vm(geno, cov = fb_cov(Qprec, type = "precision")),
    data = dat, known_matrices = list(Qprec = Q), backend = "inla",
    verbose = FALSE
  )

  gs_greta <- genomic_summary(fit_greta)
  gs_brms <- genomic_summary(fit_brms)
  gs_inla <- genomic_summary(fit_inla)

  # Every engine's heritability credible interval contains the truth.
  for (gs in list(gs_greta, gs_brms, gs_inla)) {
    expect_s3_class(gs, "fb_genomic_summary")
    expect_gte(true_h2, gs$heritability[["q2.5"]])
    expect_lte(true_h2, gs$heritability[["q97.5"]])
  }

  # Cross-engine agreement (the triangulation thesis): the three
  # posterior-mean heritabilities are mutually close.
  h2_means <- c(
    gs_greta$heritability[["mean"]],
    gs_brms$heritability[["mean"]],
    gs_inla$heritability[["mean"]]
  )
  expect_lt(max(h2_means) - min(h2_means), 0.2)

  # GEBVs are present on all three engines and track the true breeding
  # values (positive rank correlation well clear of zero).
  for (gs in list(gs_greta, gs_brms, gs_inla)) {
    expect_false(is.null(gs$gebv))
    expect_equal(nrow(gs$gebv), 30L)
    rho <- stats::cor(gs$gebv$mean, sim$u_true, method = "spearman")
    expect_gt(rho, 0.5)
    expect_true(all(gs$gebv$reliability >= 0 & gs$gebv$reliability <= 1))
  }
})
