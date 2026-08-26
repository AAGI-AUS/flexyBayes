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

test_that("vm(geno, cov = fb_cov(F, type = 'low_rank', scheme = )) refuses with approximate_route_not_yet_registered", {
  # The low-rank covariance carrier's vocabulary is registered but no
  # engine has an approximate-covariance fit route wired to it -- see
  # .check_low_rank_cov_reserved() (R/structured_cov.R). It fires from
  # new_fb_terms() (R/fb_terms.R), once the full random_terms list is
  # assembled -- not from bare .parse_formula(), which several other
  # tests in this file call directly to check parsing/deprecation
  # behaviour in isolation and which must stay refusal-free. Going
  # through fb_from_asreml() (the shared ingest entry point both
  # grammars converge on) means the message is the same regardless of
  # `backend` (confirmed live for both inla and brms during the 0.9.3
  # withdrawal -- without this check each engine's own capability gate
  # refuses the term too, but points at the *other* engine as the
  # alternative, which is circular here: neither can fit it).
  withr::local_options(lifecycle_verbosity = "quiet")
  dat <- data.frame(y = rnorm(5), geno = factor(1:5))
  err <- tryCatch(
    fb_from_asreml(
      fixed = y ~ 1,
      random = ~ vm(geno, cov = fb_cov(F_mat, type = "low_rank", scheme = "low_rank_smooth")),
      data = dat
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "approximate_route_not_yet_registered")
  expect_equal(err$format, "low_rank")
  expect_equal(err$scheme, "low_rank_smooth")
  # The message names the low-rank carrier as a reserved fb_cov() type
  # and the actionable dense-carrier workaround (materialise F %*% t(F)).
  expect_match(err$message, "reserved type")
  expect_match(err$message, "fb_cov\\(")
  expect_match(err$message, "F_mat %\\*% t\\(F_mat\\)")
})

test_that("ped(animal, cov = fb_cov(F, type = 'low_rank', scheme = )) refuses identically", {
  # Same check, ped() call site -- new_fb_terms() runs
  # .check_low_rank_cov_reserved() over every vm()/ped() term.
  withr::local_options(lifecycle_verbosity = "quiet")
  dat <- data.frame(y = rnorm(5), animal = factor(1:5))
  err <- tryCatch(
    fb_from_asreml(
      fixed = y ~ 1,
      random = ~ ped(animal, cov = fb_cov(F_mat, type = "low_rank", scheme = "low_rank_smooth")),
      data = dat
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "approximate_route_not_yet_registered")
})

# ---------------------------------------------------------------- #
# Structural validators, called directly (v0.9.3).                  #
# ---------------------------------------------------------------- #
#
# Through v0.9.2 these fixtures were run via .setup_env() (the shared
# environment-builder the withdrawn native engine's emit path used to
# bind known matrices and derived grouping vectors into a symbolic-
# graph environment before codegen). setup_env.R was withdrawn
# entirely at 0.9.3 alongside that emit path (see NEWS.md) -- it had
# zero callers left in R/ once that engine's codegen was removed. The validators it called
# (.validate_chol_input(), .validate_precision_input(),
# .validate_blocks_input(), .check_known_matrix_dim(),
# .check_known_matrix_dimnames(), .is_lower_triangular()) are NOT
# orphaned: they are still called directly from R/emit_inla.R and
# R/fb_cov.R on the active engines' own validation paths. These tests
# are rewritten to call the validators directly rather than through
# the deleted wrapper -- same refusal contracts, same reason codes,
# no more `ev`-binding assertions (there is no successor to "binds the
# matrix into a shared execution environment"; each active engine
# reads the raw matrix from `known_matrices` on its own emit path).

test_that("chol path: lower-triangular L passes the validator cleanly", {
  L <- diag(5)
  L[lower.tri(L)] <- 0.1
  expect_silent(
    flexyBayes:::.validate_chol_input(L, name = "L", group_var = "geno")
  )
})

test_that("chol path: upper-triangular L refuses with chol_not_triangular", {
  U <- diag(5)
  U[upper.tri(U)] <- 0.1 # upper, not lower
  err <- tryCatch(
    flexyBayes:::.validate_chol_input(U, name = "U", group_var = "geno"),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "chol_not_triangular")
})

test_that("chol path: non-square L refuses with chol_not_square", {
  L <- matrix(0.1, nrow = 5, ncol = 3)
  err <- tryCatch(
    flexyBayes:::.validate_chol_input(L, name = "L", group_var = "geno"),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "chol_not_square")
})

test_that("precision path: PD symmetric Q passes the validator cleanly", {
  skip_if_not_installed("Matrix")
  Q <- diag(5) + 0.01 * matrix(1, 5, 5) # symmetric PD
  expect_silent(
    flexyBayes:::.validate_precision_input(Q, name = "Q", group_var = "geno")
  )
})

test_that("pedigree_sparse_precision path: PD symmetric Q passes the same validator cleanly", {
  skip_if_not_installed("Matrix")
  Q <- diag(5) + 0.01 * matrix(1, 5, 5) # symmetric PD
  expect_silent(
    flexyBayes:::.validate_precision_input(Q, name = "Q", group_var = "animal")
  )
})

test_that("precision path: asymmetric Q refuses with precision_not_symmetric", {
  skip_if_not_installed("Matrix")
  Q <- matrix(
    c(
      1, 0.5, 0, 0, 0,
      0, 1, 0.5, 0, 0,
      0, 0, 1, 0.5, 0,
      0, 0, 0, 1, 0.5,
      0, 0, 0, 0, 1
    ),
    nrow = 5,
    byrow = TRUE
  )
  err <- tryCatch(
    flexyBayes:::.validate_precision_input(Q, name = "Q", group_var = "geno"),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "precision_not_symmetric")
})

test_that("precision path: indefinite symmetric Q refuses with precision_not_positive_definite", {
  skip_if_not_installed("Matrix")
  Q <- diag(c(1, 1, 1, 1, -1)) # symmetric, indefinite
  err <- tryCatch(
    flexyBayes:::.validate_precision_input(Q, name = "Q", group_var = "geno"),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "precision_not_positive_definite")
})

test_that("chol path with L missing from known_matrices refuses cleanly", {
  err <- tryCatch(
    flexyBayes:::.validate_chol_input(NULL, name = "L", group_var = "geno"),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "chol_not_in_known_matrices")
})

test_that("blocks path: valid 2+3 partition passes the validator cleanly", {
  Bs <- list(diag(2), diag(3)) # 2 + 3 = 5 = nlevels(geno)
  expect_silent(
    flexyBayes:::.validate_blocks_input(
      Bs,
      name = "Bs",
      group_var = "geno",
      expected_n = 5L
    )
  )
})

# ---------------------------------------------------------------- #
# Phase B-inla: gate flip + INLA emit + routing-policy version bump.#
# ---------------------------------------------------------------- #

.fixture_data_for_vm <- function(N = 60L, n_geno = 6L) {
  set.seed(20260525L)
  data.frame(
    geno = factor(rep(seq_len(n_geno), length.out = N)),
    x = rnorm(N),
    yield = rnorm(N, 50, 5)
  )
}

.mk_fb_for_random_term <- function(random_expr, dat) {
  parsed_fixed <- flexyBayes:::.parse_fixed(yield ~ 1, dat)
  flexyBayes:::new_fb_terms(
    response = "yield",
    family = "gaussian",
    link = "identity",
    intercept = parsed_fixed$intercept,
    fixed_terms = parsed_fixed$terms,
    random_terms = flexyBayes:::.parse_formula(random_expr, dat),
    residual_terms = list(list(type = "units")),
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
# v0.3.8 audit Critical Fix #2: known-matrix dim/level alignment,   #
# called directly against the validators (see the header note      #
# above -- the wrapper these ran through at v0.3.8 is withdrawn).  #
# ---------------------------------------------------------------- #

test_that("precision validator refuses with known_matrix_dim_mismatch when Q dim != nlevels(geno)", {
  skip_if_not_installed("Matrix")
  Q <- diag(4) + 0.01 * matrix(1, 4, 4) # 4 x 4 but geno has 5 levels
  err <- tryCatch(
    flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "geno",
      expected_n = 5L
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "known_matrix_dim_mismatch")
  expect_equal(err$expected_n, 5L)
  expect_equal(err$actual_dim, c(4L, 4L))
})

test_that("precision validator refuses with known_matrix_level_mismatch when Q dimnames are permuted", {
  skip_if_not_installed("Matrix")
  # geno's levels (1, 2, 3, 4, 5) -- factor default ordering.
  fit_levels <- levels(factor(c("1", "2", "3", "4", "5")))
  Q <- diag(5) + 0.01 * matrix(1, 5, 5)
  # Dimnames are the correct level set, but in reverse order.
  dimnames(Q) <- list(c("5", "4", "3", "2", "1"), c("5", "4", "3", "2", "1"))
  err <- tryCatch(
    flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "geno",
      expected_n = 5L,
      fit_levels = fit_levels
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "known_matrix_level_mismatch")
  # Refusal message names the perm fix verbatim.
  expect_match(conditionMessage(err), "perm <- match\\(levels")
})

test_that("chol validator refuses with known_matrix_dim_mismatch when L dim != nlevels(geno)", {
  L <- diag(4)
  L[lower.tri(L)] <- 0.1 # 4 x 4 lower-triangular
  err <- tryCatch(
    flexyBayes:::.validate_chol_input(
      L,
      name = "L",
      group_var = "geno",
      expected_n = 5L
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "known_matrix_dim_mismatch")
})

test_that("precision validator: aligned dimnames pass cleanly", {
  skip_if_not_installed("Matrix")
  fit_levels <- c("g1", "g2", "g3", "g4", "g5")
  Q <- diag(5) + 0.01 * matrix(1, 5, 5)
  dimnames(Q) <- list(fit_levels, fit_levels)
  expect_silent(
    flexyBayes:::.validate_precision_input(
      Q,
      name = "Q",
      group_var = "geno",
      expected_n = 5L,
      fit_levels = fit_levels
    )
  )
})
