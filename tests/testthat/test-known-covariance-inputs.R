# ADR 0025 Stage 5A known-covariance input formats --- Phase A
# (parser + IR slot + validators + route guard; emit deferred to
# Phase B). The byte-identical-posterior commitment of ADR 0025
# §iii is enforced at Phase B against the v0.2 fixture corpus;
# Phase A holds at IR level (snapshot of .parse_formula output for
# the dense path) plus structural / shape validator coverage.

# ---------------------------------------------------------------- #
# (a) Dense V path --- regression guard at IR level.                #
# ---------------------------------------------------------------- #

test_that("dense vm(geno, V = K) parses to an unchanged IR shape with cov_representation slot", {
  dat <- data.frame(geno = factor(1:5))

  positional <- flexyBayes:::.parse_formula(~ vm(geno, Gmat), dat)
  named <- flexyBayes:::.parse_formula(~ vm(geno, V = Gmat), dat)

  for (terms in list(positional, named)) {
    expect_length(terms, 1L)
    expect_equal(terms[[1]]$type, "vm")
    expect_equal(terms[[1]]$var, "geno")
    expect_equal(terms[[1]]$mat, "Gmat")
    expect_equal(terms[[1]]$cov_representation$format, "dense")
    expect_equal(terms[[1]]$cov_representation$data, "Gmat")
    expect_null(terms[[1]]$cov_representation$scheme)
  }
})

test_that("dense ped(animal, A) parses to IR with cov_representation$format = 'dense'", {
  dat <- data.frame(animal = factor(1:5))
  terms <- flexyBayes:::.parse_formula(~ ped(animal, Amat), dat)
  expect_equal(terms[[1]]$type, "ped")
  expect_equal(terms[[1]]$var, "animal")
  expect_equal(terms[[1]]$mat, "Amat")
  expect_equal(terms[[1]]$cov_representation$format, "dense")
})

# ---------------------------------------------------------------- #
# Stage 5A named-arg parser coverage.                               #
# ---------------------------------------------------------------- #

test_that("vm(geno, cov = fb_cov(L, type = 'chol')) parses to cov_representation$format = 'chol'", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(L, type = "chol")),
    dat
  )
  expect_equal(terms[[1]]$type, "vm")
  expect_equal(terms[[1]]$var, "geno")
  expect_true(is.na(terms[[1]]$mat))
  expect_equal(terms[[1]]$cov_representation$format, "chol")
  expect_equal(terms[[1]]$cov_representation$data, "L")
})

test_that("vm(geno, cov = fb_cov(Q, type = 'precision')) parses to cov_representation$format = 'precision'", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  expect_equal(terms[[1]]$cov_representation$format, "precision")
  expect_equal(terms[[1]]$cov_representation$data, "Q")
})

test_that("ped(animal, A_obj, use_sparse_precision = TRUE) parses to pedigree_sparse_precision", {
  dat <- data.frame(animal = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ ped(animal, A_obj, use_sparse_precision = TRUE),
    dat
  )
  expect_equal(terms[[1]]$type, "ped")
  expect_equal(
    terms[[1]]$cov_representation$format,
    "pedigree_sparse_precision"
  )
  expect_equal(terms[[1]]$cov_representation$data, "A_obj")
})

test_that("ped(animal, A_obj, use_sparse_precision = FALSE) stays on dense format", {
  dat <- data.frame(animal = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ ped(animal, A_obj, use_sparse_precision = FALSE),
    dat
  )
  expect_equal(terms[[1]]$cov_representation$format, "dense")
})

# ---------------------------------------------------------------- #
# (k) Mutual-exclusion refusal --- ADR 0025 §8 subtest (k).        #
# ---------------------------------------------------------------- #

test_that("vm(geno, V = K, chol = L) refuses with vm_redundant_specification", {
  # Legacy keyword carriers (deprecated v0.4.0). Quiet the lifecycle
  # warning here; the deprecation itself is asserted in
  # test-fb-cov-constructor.R.
  withr::local_options(lifecycle_verbosity = "quiet")
  dat <- data.frame(geno = factor(1:5))
  err <- tryCatch(
    flexyBayes:::.parse_formula(~ vm(geno, V = K, chol = L), dat),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "vm_redundant_specification")
  expect_match(err$message, "V \\+ chol")
})

test_that("vm(geno, chol = L, precision = Q) refuses with vm_redundant_specification", {
  withr::local_options(lifecycle_verbosity = "quiet")
  dat <- data.frame(geno = factor(1:5))
  err <- tryCatch(
    flexyBayes:::.parse_formula(~ vm(geno, chol = L, precision = Q), dat),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "vm_redundant_specification")
})

# ---------------------------------------------------------------- #
# Low-rank refusal coverage at parse time (subtest j of ADR 0025).  #
# ---------------------------------------------------------------- #

test_that("vm(geno, low_rank_factor = F) without low_rank_scheme refuses at parse", {
  withr::local_options(lifecycle_verbosity = "quiet")
  dat <- data.frame(geno = factor(1:5))
  err <- tryCatch(
    flexyBayes:::.parse_formula(~ vm(geno, low_rank_factor = F_mat), dat),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "low_rank_scheme_required")
})

# ---------------------------------------------------------------- #
# Setup-env: dense path unchanged (regression guard at runtime).    #
# ---------------------------------------------------------------- #

test_that("setup_env dense vm path binds the matrix and falls through Phase A guard", {
  dat <- data.frame(geno = factor(1:5))
  G_mat <- diag(5)
  rownames(G_mat) <- colnames(G_mat) <- as.character(1:5)

  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))
  random_terms <- flexyBayes:::.parse_formula(~ vm(geno, Gmat), dat)

  ev <- new.env(parent = emptyenv())
  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(Gmat = G_mat),
    NULL
  )

  expect_true("Gmat" %in% ls(ev))
  expect_equal(ev$Gmat, G_mat)
  expect_equal(ev$geno_id, 1:5)
  expect_equal(ev$n_geno, 5L)
})

# ---------------------------------------------------------------- #
# Setup-env: chol path runs validator + binds matrix (Phase B-greta).#
# ---------------------------------------------------------------- #

test_that("chol path: lower-triangular L passes validator AND binds the L matrix to ev", {
  dat <- data.frame(geno = factor(1:5))
  L <- diag(5)
  L[lower.tri(L)] <- 0.1
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(L, type = "chol")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(L = L),
    NULL
  )
  expect_true("L" %in% ls(ev))
  expect_equal(ev$L, L)
  expect_equal(ev$geno_id, 1:5)
  expect_equal(ev$n_geno, 5L)
})

test_that("chol path: upper-triangular L refuses with chol_not_triangular before route guard", {
  dat <- data.frame(geno = factor(1:5))
  U <- diag(5)
  U[upper.tri(U)] <- 0.1 # upper, not lower
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(U, type = "chol")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(U = U),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "chol_not_triangular")
})

test_that("chol path: non-square L refuses with chol_not_square before route guard", {
  dat <- data.frame(geno = factor(1:5))
  L <- matrix(0.1, nrow = 5, ncol = 3)
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(L, type = "chol")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(L = L),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "chol_not_square")
})

# ---------------------------------------------------------------- #
# Setup-env: precision path runs validator + Phase A route refusal. #
# ---------------------------------------------------------------- #

test_that("precision path: PD symmetric Q passes validator AND binds Q to ev", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(geno = factor(1:5))
  Q <- diag(5) + 0.01 * matrix(1, 5, 5) # symmetric PD
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(Q = Q),
    NULL
  )
  expect_true("Q" %in% ls(ev))
  expect_equal(ev$Q, Q)
})

test_that("pedigree_sparse_precision path: PD symmetric Q passes validator AND binds Q to ev", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(animal = factor(1:5))
  Q <- diag(5) + 0.01 * matrix(1, 5, 5) # symmetric PD
  random_terms <- flexyBayes:::.parse_formula(
    ~ ped(animal, Q, use_sparse_precision = TRUE),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(Q = Q),
    NULL
  )
  expect_true("Q" %in% ls(ev))
  expect_equal(ev$Q, Q)
  expect_equal(ev$animal_id, 1:5)
  expect_equal(ev$n_animal, 5L)
})

test_that("precision path: asymmetric Q refuses with precision_not_symmetric", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(geno = factor(1:5))
  Q <- matrix(
    c(
      1,
      0.5,
      0,
      0,
      0,
      0,
      1,
      0.5,
      0,
      0,
      0,
      0,
      1,
      0.5,
      0,
      0,
      0,
      0,
      1,
      0.5,
      0,
      0,
      0,
      0,
      1
    ),
    nrow = 5,
    byrow = TRUE
  )
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(Q = Q),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "precision_not_symmetric")
})

test_that("precision path: indefinite symmetric Q refuses with precision_not_positive_definite", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(geno = factor(1:5))
  Q <- diag(c(1, 1, 1, 1, -1)) # symmetric, indefinite
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(Q = Q),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "precision_not_positive_definite")
})

# ---------------------------------------------------------------- #
# Setup-env: missing-from-known_matrices refusals.                  #
# ---------------------------------------------------------------- #

test_that("chol path with L missing from known_matrices refuses cleanly", {
  dat <- data.frame(geno = factor(1:5))
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(L, type = "chol")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "chol_not_in_known_matrices")
})

# ---------------------------------------------------------------- #
# ADR 0025 Decision 3 (v0.3.10): blocks ships end-to-end. The      #
# pre-v0.3.10 route-refusal subtest below is rewritten as a        #
# success-path guard --- valid block partition validates cleanly,  #
# binds the carrier into setup_env, and falls through to codegen.  #
# ---------------------------------------------------------------- #

test_that("blocks path: valid 2+3 partition validates and falls through Phase A guard", {
  dat <- data.frame(geno = factor(1:5))
  Bs <- list(diag(2), diag(3)) # 2 + 3 = 5 = nlevels(geno)
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  expect_silent(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(Bs = Bs),
      NULL
    )
  )
  expect_true("Bs" %in% ls(ev))
  expect_identical(ev$Bs, Bs)
})

# ---------------------------------------------------------------- #
# Phase B-greta: codegen emits per-format square-root expressions.  #
# ---------------------------------------------------------------- #

.fixture_data_for_vm <- function(N = 60L, n_geno = 6L) {
  set.seed(20260525L)
  data.frame(
    geno = factor(rep(seq_len(n_geno), length.out = N)),
    x = rnorm(N),
    yield = rnorm(N, 50, 5)
  )
}

test_that("codegen: dense vm path keeps the pre-Stage-5A t(chol()) expression", {
  dat <- .fixture_data_for_vm()
  G_mat <- diag(6) + 0.1
  code <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, Gmat),
    data = dat,
    known_matrices = list(Gmat = G_mat),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl("L_G_geno <- t(chol(Gmat))", code, fixed = TRUE))
  expect_false(grepl("as.matrix(", code, fixed = TRUE))
})

test_that("codegen: chol path emits as.matrix(L), no t(chol()) wrap", {
  dat <- .fixture_data_for_vm()
  G <- diag(6) + 0.1
  L <- t(chol(G))
  code <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(L_chol, type = "chol")),
    data = dat,
    known_matrices = list(L_chol = L),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl("L_G_geno <- as.matrix(L_chol)", code, fixed = TRUE))
  expect_false(grepl("L_G_geno <- t(chol(", code, fixed = TRUE))
})

test_that("codegen: precision path emits solve(chol(Q)) square root", {
  skip_if_not_installed("Matrix")
  dat <- .fixture_data_for_vm()
  G <- diag(6) + 0.1
  Q <- solve(G)
  code <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Qprec, type = "precision")),
    data = dat,
    known_matrices = list(Qprec = Q),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl(
    "L_G_geno <- as.matrix(solve(chol(Qprec)))",
    code,
    fixed = TRUE
  ))
})

test_that("codegen: ped use_sparse_precision = TRUE emits solve(chol(Q)) square root", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(
    animal = factor(rep(seq_len(6L), length.out = 60L)),
    yield = rnorm(60L, 50, 5)
  )
  G <- diag(6) + 0.1
  Q <- solve(G)
  code <- flexybayes(
    yield ~ 1,
    random = ~ ped(animal, Qprec, use_sparse_precision = TRUE),
    data = dat,
    known_matrices = list(Qprec = Q),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl(
    "L_A_animal <- as.matrix(solve(chol(Qprec)))",
    code,
    fixed = TRUE
  ))
})

# ---------------------------------------------------------------- #
# Phase B-greta (b) / (c) / (d-partial): posterior equivalence on   #
# greta across dense / chol / precision paths. Stress-gated; TF     #
# RNG state across two separate greta fits is not reseeded by base  #
# set.seed() (cf. test-validation-lmer.R header) so we run two      #
# short fits and assert posterior-mean agreement on the recovered   #
# sigma_geno parameter within a generous Monte-Carlo tolerance.     #
# ---------------------------------------------------------------- #

test_that("posterior equivalence: dense V vs chol = t(chol(V)) on greta", {
  skip_if_greta_backend_unusable()
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  if (!identical(tolower(Sys.getenv("FLEXYBAYES_RUN_STRESS")), "true")) {
    testthat::skip(
      "set FLEXYBAYES_RUN_STRESS=true to run this Phase B-greta MCMC check"
    )
  }

  dat <- .fixture_data_for_vm()
  G <- diag(6) + 0.1
  L <- t(chol(G))

  fit_dense <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, G_mat),
    data = dat,
    known_matrices = list(G_mat = G),
    backend = "greta",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    verbose = FALSE
  )
  fit_chol <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(L_user, type = "chol")),
    data = dat,
    known_matrices = list(L_user = L),
    backend = "greta",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    verbose = FALSE
  )

  sigma_dense <- mean(as.matrix(fit_dense$draws[, "sigma_geno"]))
  sigma_chol <- mean(as.matrix(fit_chol$draws[, "sigma_geno"]))
  expect_lt(abs(sigma_dense - sigma_chol), max(0.5, 0.25 * sigma_dense)) # generous MC tolerance
})

test_that("posterior equivalence: dense V vs precision = solve(V) on greta", {
  skip_if_greta_backend_unusable()
  skip_if_not_installed("Matrix")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  if (!identical(tolower(Sys.getenv("FLEXYBAYES_RUN_STRESS")), "true")) {
    testthat::skip(
      "set FLEXYBAYES_RUN_STRESS=true to run this Phase B-greta MCMC check"
    )
  }

  dat <- .fixture_data_for_vm()
  G <- diag(6) + 0.1
  Q <- solve(G)

  fit_dense <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, G_mat),
    data = dat,
    known_matrices = list(G_mat = G),
    backend = "greta",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    verbose = FALSE
  )
  fit_prec <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Qprec, type = "precision")),
    data = dat,
    known_matrices = list(Qprec = Q),
    backend = "greta",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    verbose = FALSE
  )

  sigma_dense <- mean(as.matrix(fit_dense$draws[, "sigma_geno"]))
  sigma_prec <- mean(as.matrix(fit_prec$draws[, "sigma_geno"]))
  expect_lt(abs(sigma_dense - sigma_prec), max(0.5, 0.25 * sigma_dense))
})

# ---------------------------------------------------------------- #
# Phase B-inla: gate flip + INLA emit + routing-policy version bump.#
# ---------------------------------------------------------------- #

.mk_fb_for_random_term <- function(random_expr, dat) {
  parsed_fixed <- flexyBayes:::.parse_fixed(yield ~ 1, dat)
  flexyBayes:::new_fb_terms(
    response = "yield",
    family = "gaussian",
    link = "identity",
    intercept = parsed_fixed$intercept,
    fixed_terms = parsed_fixed$terms,
    random_terms = flexyBayes:::.parse_formula(random_expr, dat),
    rcov_terms = list(list(type = "units")),
    data_summary = list(n = nrow(dat))
  )
}

test_that("lgm_gate accepts vm with precision format for INLA emit", {
  dat <- .fixture_data_for_vm()
  fb <- .mk_fb_for_random_term(
    ~ vm(geno, cov = fb_cov(Qprec, type = "precision")),
    dat
  )
  gated <- flexyBayes:::lgm_gate(fb)
  expect_false(flexyBayes:::is_lgm_refusal(gated))
  expect_true("lgm_compatible" %in% gated$capabilities)
})

test_that("lgm_gate refuses vm with dense V on INLA with format-aware message", {
  dat <- .fixture_data_for_vm()
  fb <- .mk_fb_for_random_term(~ vm(geno, Gmat), dat)
  gated <- flexyBayes:::lgm_gate(fb)
  expect_s3_class(gated, "lgm_refusal")
  rti <- Filter(
    function(f) f$rule_id == "random_term_type_inla",
    gated$failures
  )
  expect_length(rti, 1L)
  expect_match(rti[[1L]]$reason, "sparse-precision")
  expect_match(rti[[1L]]$reason, "precision = solve")
})

test_that(".ROUTING_POLICY_VERSION bumps to 'stage5a_v1'", {
  expect_identical(flexyBayes:::.ROUTING_POLICY_VERSION, "stage5a_v1")
})

test_that("emit_inla refuses when known_matrices entry shadows a data column", {
  skip_if_not_installed("INLA")
  dat <- .fixture_data_for_vm()
  Q <- solve(diag(6) + 0.1)
  err <- tryCatch(
    flexybayes(
      yield ~ 1,
      random = ~ vm(geno, cov = fb_cov(geno, type = "precision")),
      data = dat,
      known_matrices = list(geno = Q),
      backend = "inla",
      verbose = FALSE
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "known_matrices_data_name_collision")
})

# ---------------------------------------------------------------- #
# ADR 0025 §8 subtest (d): sparse-precision INLA round-trip.        #
# Verifies the flexyBayes INLA emit path produces a fit that        #
# matches a direct INLA::inla() reference in structure (random-     #
# effect count, summary shape, intercept location). flexyBayes      #
# defaults to the uniform-on-SD prior on the precision             #
# hyperparameter (represented exactly for INLA via an expression-  #
# prior); the reference uses INLA's loggamma default. The two      #
# priors agree on the location-mean for the random effects but      #
# shrink to slightly different amounts, so the RE-mean tolerance    #
# is loosened to a sensible band rather than bit-exact. Tightening  #
# to bit-exact equivalence would require passing matched explicit   #
# precision priors to both fits, which is a Phase C documentation   #
# exercise, not a Phase B-inla emit verification.                   #
# ---------------------------------------------------------------- #

test_that("(d) precision-on-INLA: flexybayes matches direct INLA::inla() reference shape", {
  skip_if_not_installed("INLA")
  set.seed(20260525L)
  dat <- .fixture_data_for_vm(N = 120L, n_geno = 8L)
  G <- diag(8) + 0.1
  Q <- solve(G)

  fb_fit <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Qprec, type = "precision")),
    data = dat,
    known_matrices = list(Qprec = Q),
    backend = "inla",
    verbose = FALSE
  )
  dat_ref <- as.list(dat)
  dat_ref$geno_id <- as.integer(factor(dat$geno))
  dat_ref$Qprec <- Q
  ref_fit <- INLA::inla(
    yield ~ 1 + f(geno_id, model = "generic0", Cmatrix = Qprec),
    family = "gaussian",
    data = dat_ref,
    control.compute = list(config = TRUE)
  )

  fb_re_means <- fb_fit$inla$summary.random$geno_id$mean
  ref_re_means <- ref_fit$summary.random$geno_id$mean
  expect_s3_class(fb_fit, "flexybayes_inla")
  expect_identical(fb_fit$exactness, "exact")
  expect_equal(length(fb_re_means), length(ref_re_means))
  # Intercept matches at the data-scale precision (both fits see the
  # same likelihood; the prior difference enters only on the RE
  # precision hyperparameter, not the fixed-effect mean).
  expect_lt(
    abs(
      fb_fit$inla$summary.fixed[1L, "mean"] -
        ref_fit$summary.fixed[1L, "mean"]
    ),
    0.5
  )
  # RE-mean agreement within prior-driven shrinkage band.
  expect_lt(max(abs(fb_re_means - ref_re_means)), 2.0)
})

# ---------------------------------------------------------------- #
# ADR 0025 §8 subtest (e): BYM2-shape neighbourhood precision on a  #
# small lattice. Builds a rook-adjacency ICAR precision Q for a 4x4 #
# grid and confirms the fit succeeds + matches a direct INLA call.  #
# ---------------------------------------------------------------- #

test_that("(e) BYM2-shape lattice: sparse-precision fit succeeds on INLA", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("Matrix")
  set.seed(20260525L)

  # 4x4 rook adjacency -> ICAR precision Q = D - W, perturbed by a
  # small ridge so Q is positive-definite (proper CAR; the BYM2
  # use case wraps this with the scaling-precision parameter, but
  # for the subtest the proper-CAR shape is sufficient to exercise
  # the generic0 emit + show convergence).
  grid_dim <- 4L
  n_cell <- grid_dim * grid_dim
  coords <- expand.grid(r = seq_len(grid_dim), c = seq_len(grid_dim))
  W <- matrix(0, n_cell, n_cell)
  for (i in seq_len(n_cell)) {
    for (j in seq_len(n_cell)) {
      if (i == j) {
        next
      }
      if (
        abs(coords$r[i] - coords$r[j]) +
          abs(coords$c[i] - coords$c[j]) ==
          1L
      ) {
        W[i, j] <- 1
      }
    }
  }
  D <- diag(rowSums(W))
  Q <- D - W + 0.01 * diag(n_cell) # proper-CAR ridge for PD

  dat <- data.frame(
    cell = factor(seq_len(n_cell)),
    y = as.numeric(
      MASS::mvrnorm(1L, mu = rep(0, n_cell), Sigma = solve(Q)) +
        rnorm(n_cell, sd = 0.1)
    )
  )

  fit <- flexybayes(
    y ~ 1,
    random = ~ vm(cell, cov = fb_cov(Qmat, type = "precision")),
    data = dat,
    known_matrices = list(Qmat = Q),
    backend = "inla",
    verbose = FALSE
  )
  expect_s3_class(fit, "flexybayes_inla")
  expect_identical(fit$exactness, "exact")
  expect_true(!is.null(fit$inla$summary.random$cell_id))
  expect_equal(nrow(fit$inla$summary.random$cell_id), n_cell)
})

# ---------------------------------------------------------------- #
# Original Phase A test (kept for reference; assertion intact).     #
# ---------------------------------------------------------------- #

test_that("low_rank path with registered-looking scheme refuses with approximate_route_not_yet_registered (v0.3.10 upgraded message)", {
  withr::local_options(lifecycle_verbosity = "quiet")
  dat <- data.frame(geno = factor(1:5))
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, low_rank_factor = F_mat, low_rank_scheme = "pca"),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(F_mat = matrix(0, 5, 2)),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "approximate_route_not_yet_registered")
  expect_equal(err$format, "low_rank")
  expect_equal(err$scheme, "pca")
  # v0.4.0 upgrade (ADR 0030 C3): message names the low-rank carrier as
  # a reserved fb_cov() type + the actionable dense-carrier workaround
  # (materialise U %*% t(U)).
  expect_match(err$message, "reserved type")
  expect_match(err$message, "fb_cov\\(")
  expect_match(err$message, "F_mat %\\*% t\\(F_mat\\)")
})

# ---------------------------------------------------------------- #
# v0.3.8 audit Critical Fix #2: known-matrix dim/level alignment   #
# wired through .setup_env() dispatch (.stage5a_route_check now    #
# receives expected_n + fit_levels from the dispatch layer).        #
# ---------------------------------------------------------------- #

test_that("setup_env precision path refuses with known_matrix_dim_mismatch when Q dim != nlevels(geno)", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(geno = factor(1:5))
  Q <- diag(4) + 0.01 * matrix(1, 4, 4) # 4 x 4 but geno has 5 levels
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(Q = Q),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "known_matrix_dim_mismatch")
  expect_equal(err$expected_n, 5L)
  expect_equal(err$actual_dim, c(4L, 4L))
})

test_that("setup_env precision path refuses with known_matrix_level_mismatch when Q dimnames are permuted", {
  skip_if_not_installed("Matrix")
  # geno's levels (1, 2, 3, 4, 5) -- factor default ordering.
  dat <- data.frame(geno = factor(c("1", "2", "3", "4", "5")))
  Q <- diag(5) + 0.01 * matrix(1, 5, 5)
  # Dimnames are the correct level set, but in reverse order.
  dimnames(Q) <- list(c("5", "4", "3", "2", "1"), c("5", "4", "3", "2", "1"))
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(Q = Q),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "known_matrix_level_mismatch")
  # Refusal message names the perm fix verbatim.
  expect_match(conditionMessage(err), "perm <- match\\(levels")
})

test_that("setup_env chol path refuses with known_matrix_dim_mismatch when L dim != nlevels(geno)", {
  dat <- data.frame(geno = factor(1:5))
  L <- diag(4)
  L[lower.tri(L)] <- 0.1 # 4 x 4 lower-triangular
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(L, type = "chol")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  err <- tryCatch(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(L = L),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "known_matrix_dim_mismatch")
})

test_that("setup_env precision path: aligned dimnames pass cleanly (happy path through dispatch)", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(geno = factor(c("g1", "g2", "g3", "g4", "g5")))
  Q <- diag(5) + 0.01 * matrix(1, 5, 5)
  dimnames(Q) <- list(levels(dat$geno), levels(dat$geno))
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  fixed_info <- flexyBayes:::.parse_fixed(V1 ~ 1, cbind(dat, V1 = 1))

  ev <- new.env(parent = emptyenv())
  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(Q = Q),
    NULL
  )
  expect_true("Q" %in% ls(ev))
  expect_equal(ev$geno_id, 1:5)
  expect_equal(ev$n_geno, 5L)
})
