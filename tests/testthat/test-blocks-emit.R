# test-blocks-emit.R --- ADR 0025 Decision 3 (v0.3.10) block-diagonal
# carrier on greta + INLA, plus the upgraded low_rank refusal-stub
# (Decision 4) and the three first-migration entries into the v0.3.8
# C7 refusal-registry scaffold.
#
# Activates ADR 0025 §8 subtests (f), (g), (h), (i), (j) end-to-end.
# Heavy MCMC paths are stress-gated via FLEXYBAYES_RUN_STRESS; the
# default-on subtests exercise validators, codegen, formula assembly,
# fb_plan() shape, and registry surface --- all deterministic and
# cheap.
#
# Fixture shape: two blocks of sizes 2 and 3 totalling 5 = nlevels(geno).
# The validator's block_sizes + total_n contract holds independently
# of K; the small fixture keeps the tests fast.

# ---------------------------------------------------------------- #
# (a) Validator happy path                                          #
# ---------------------------------------------------------------- #

test_that(".validate_blocks_input() returns block_sizes + total_n + canonical list on a valid 2+3 partition", {
  Bs <- list(diag(2), diag(3) + 0.1)
  meta <- flexyBayes:::.validate_blocks_input(
    Bs,
    name = "Bs",
    group_var = "geno",
    expected_n = 5L
  )
  expect_type(meta, "list")
  expect_equal(meta$block_sizes, c(2L, 3L))
  expect_equal(meta$total_n, 5L)
  expect_identical(meta$blocks, Bs)
})

# ---------------------------------------------------------------- #
# (b) block_partition_incomplete refusal                            #
# ---------------------------------------------------------------- #

test_that(".validate_blocks_input() refuses with block_partition_incomplete when sum(n_k) != expected_n", {
  Bs <- list(diag(2), diag(2)) # 2 + 2 = 4, not 5
  err <- tryCatch(
    flexyBayes:::.validate_blocks_input(
      Bs,
      name = "Bs",
      group_var = "geno",
      expected_n = 5L
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "block_partition_incomplete")
  expect_equal(err$expected_n, 5L)
  expect_equal(err$actual_n, 4L)
  expect_equal(err$block_sizes, c(2L, 2L))
  expect_match(err$message, "2 \\+ 2 = 4")
  expect_match(err$message, "5 levels|level count")
})

# ---------------------------------------------------------------- #
# (c) block_not_positive_definite refusal                           #
# ---------------------------------------------------------------- #

test_that(".validate_blocks_input() refuses with block_not_positive_definite when V_k is indefinite (block index named)", {
  Bs <- list(diag(2), diag(3) * -1) # block 2 is negative-definite
  err <- tryCatch(
    flexyBayes:::.validate_blocks_input(
      Bs,
      name = "Bs",
      group_var = "geno",
      expected_n = 5L
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_s3_class(err, "flexybayes_structured_cov_refusal")
  expect_equal(err$reason_code, "block_not_positive_definite")
  expect_equal(err$block_index, 2L)
  expect_match(err$message, "block 2")
})

# ---------------------------------------------------------------- #
# (c2/c3/c4) Structural refusals at the carrier level                #
# ---------------------------------------------------------------- #

test_that(".validate_blocks_input() refuses with blocks_not_in_known_matrices when NULL", {
  err <- tryCatch(
    flexyBayes:::.validate_blocks_input(
      NULL,
      name = "Bs",
      group_var = "geno",
      expected_n = 5L
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_equal(err$reason_code, "blocks_not_in_known_matrices")
  expect_match(err$message, "known_matrices = list\\(Bs = list\\(")
})

test_that(".validate_blocks_input() refuses with blocks_not_a_list when given a matrix", {
  err <- tryCatch(
    flexyBayes:::.validate_blocks_input(
      diag(5),
      name = "Bs",
      group_var = "geno",
      expected_n = 5L
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_equal(err$reason_code, "blocks_not_a_list")
})

test_that(".validate_blocks_input() refuses with blocks_empty_list", {
  err <- tryCatch(
    flexyBayes:::.validate_blocks_input(
      list(),
      name = "Bs",
      group_var = "geno",
      expected_n = 0L
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_equal(err$reason_code, "blocks_empty_list")
})

# ---------------------------------------------------------------- #
# (d) greta codegen: Matrix::bdiag() in the sqrt expression          #
# ---------------------------------------------------------------- #

test_that("codegen: blocks path emits t(chol(as.matrix(Matrix::bdiag(...)))) sqrt expression", {
  dat <- data.frame(
    geno = factor(rep(seq_len(5L), length.out = 60L)),
    yield = rnorm(60L, 50, 5)
  )
  Bs <- list(diag(2), diag(3) + 0.1)
  code <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
    data = dat,
    known_matrices = list(Bs = Bs),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl(
    "t(chol(as.matrix(Matrix::bdiag(Bs))))",
    code,
    fixed = TRUE
  ))
  expect_false(grepl("Bs %*%", code, fixed = TRUE))
})

# ---------------------------------------------------------------- #
# (e) INLA formula build: K f() calls with per-block Cmatrix         #
# ---------------------------------------------------------------- #

test_that(".build_inla_formula() blocks path emits K f() calls, one per block", {
  skip_if_not_installed("Matrix")
  dat <- data.frame(
    geno = factor(rep(seq_len(5L), length.out = 60L)),
    yield = rnorm(60L, 50, 5)
  )
  Bs <- list(diag(2), diag(3) + 0.1)
  fb <- flexyBayes:::fb_from_asreml(
    fixed = yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
    data = dat,
    family = "gaussian",
    known_matrices = list(Bs = Bs)
  )
  form_chr <- deparse(flexyBayes:::.build_inla_formula(
    fb,
    known_matrices = list(Bs = Bs)
  ))
  form_str <- paste(form_chr, collapse = " ")
  expect_true(grepl(
    "f(geno_id_block_1, model = \"generic0\", Cmatrix = Bs_Q_1)",
    form_str,
    fixed = TRUE
  ))
  expect_true(grepl(
    "f(geno_id_block_2, model = \"generic0\", Cmatrix = Bs_Q_2)",
    form_str,
    fixed = TRUE
  ))
  expect_false(grepl("Bs_Q_3", form_str, fixed = TRUE))
})

# ---------------------------------------------------------------- #
# (f) Representation print line: exact (block-diagonal, K blocks)    #
# ---------------------------------------------------------------- #

test_that("flexybayes(plan = TRUE) Representation label renders 'exact (block-diagonal, 2 blocks)' for blocks-format vm()", {
  dat <- data.frame(
    geno = factor(rep(seq_len(5L), length.out = 60L)),
    yield = rnorm(60L, 50, 5)
  )
  Bs <- list(diag(2), diag(3) + 0.1)
  plan <- flexybayes(
    fixed = yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
    data = dat,
    known_matrices = list(Bs = Bs),
    backend = "greta",
    plan = TRUE,
    verbose = FALSE
  )
  expect_match(
    plan$representation_label,
    "^exact \\(block-diagonal, 2 blocks\\)$"
  )
})

# ---------------------------------------------------------------- #
# (g) fb_plan() representation_plan = block_diagonal with K count   #
# ---------------------------------------------------------------- #

test_that("flexybayes(plan = TRUE) representation_plan carries block_diagonal class + block_count for blocks-format vm()", {
  dat <- data.frame(
    geno = factor(rep(seq_len(5L), length.out = 60L)),
    yield = rnorm(60L, 50, 5)
  )
  Bs <- list(diag(2), diag(3) + 0.1)
  plan <- flexybayes(
    fixed = yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
    data = dat,
    known_matrices = list(Bs = Bs),
    backend = "greta",
    plan = TRUE,
    verbose = FALSE
  )
  expect_true(length(plan$representation_plan) >= 1L)
  block_entry <- Filter(
    function(rp) {
      identical(rp$representation_class, "block_diagonal")
    },
    plan$representation_plan
  )
  expect_length(block_entry, 1L)
  expect_equal(block_entry[[1L]]$block_count, 2L)
})

# ---------------------------------------------------------------- #
# (h) lgm_gate accepts blocks format on vm/ped                      #
# ---------------------------------------------------------------- #

test_that("lgm_gate accepts vm() with cov_representation$format = 'blocks'", {
  dat <- data.frame(
    geno = factor(rep(seq_len(5L), length.out = 60L)),
    yield = rnorm(60L, 50, 5)
  )
  Bs <- list(diag(2), diag(3) + 0.1)
  fb <- flexyBayes:::fb_from_asreml(
    fixed = yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
    data = dat,
    family = "gaussian",
    known_matrices = list(Bs = Bs)
  )
  gated <- lgm_gate(fb)
  expect_false(inherits(gated, "lgm_refusal"))
  expect_true("lgm_compatible" %in% gated$capabilities)
})

# ---------------------------------------------------------------- #
# (i) low_rank refusal upgrade names the registry + workaround      #
# ---------------------------------------------------------------- #

test_that("low_rank refusal message names the reserved fb_cov() carrier, v0.4.0, and the dense materialisation workaround", {
  # Legacy keyword carrier (deprecated v0.4.0); quiet the lifecycle
  # warning -- the deprecation itself is asserted in
  # test-fb-cov-constructor.R.
  withr::local_options(lifecycle_verbosity = "quiet")
  dat <- data.frame(geno = factor(1:5))
  random_terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, low_rank_factor = U, low_rank_scheme = "pca"),
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
      list(U = matrix(0, 5, 2)),
      NULL
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_match(err$message, "reserved type")
  expect_match(err$message, "fb_cov\\(")
  expect_match(err$message, "U %\\*% t\\(U\\)")
})

# ---------------------------------------------------------------- #
# (j) Refusal registry: the three v0.3.10 first-migration entries    #
# ---------------------------------------------------------------- #

test_that(".refusal_registry post-.onLoad() carries the three v0.3.10 reason codes (block_partition_incomplete, block_not_positive_definite, approximate_route_not_yet_registered)", {
  reg <- flexyBayes:::.refusal_registry
  expect_true(exists(
    "block_partition_incomplete",
    envir = reg,
    inherits = FALSE
  ))
  expect_true(exists(
    "block_not_positive_definite",
    envir = reg,
    inherits = FALSE
  ))
  expect_true(exists(
    "approximate_route_not_yet_registered",
    envir = reg,
    inherits = FALSE
  ))
})

test_that(".refusal_registry remains locked post-registration (binding-locked, not just env-locked)", {
  reg <- flexyBayes:::.refusal_registry
  expect_true(environmentIsLocked(reg))
  expect_error(
    assign("block_partition_incomplete", list(), envir = reg),
    regexp = "locked|cannot"
  )
})
