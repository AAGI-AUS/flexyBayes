# test-blocks-emit.R --- ADR 0025 Decision 3 (v0.3.10) block-diagonal
# carrier on INLA, plus the upgraded low_rank refusal-stub
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
# (d) [removed] codegen: Matrix::bdiag() in the sqrt expression      #
# ---------------------------------------------------------------- #
#
# Deleted (not rewritten): this test checked the Cholesky-square-root
# code text the withdrawn native engine's codegen produced for a
# block-diagonal carrier (see NEWS.md, 0.9.3). No active engine has a
# successor -- brms refuses a block-diagonal vm() carrier outright
# (stan_cannot_represent_structured_cov via `.capability_brms()`), and
# `return_code = TRUE` on `backend = "auto"` resolves only to brms
# ("the only active code-producing engine"; verified live: the auto
# path raises `auto_no_active_route` for this exact call). INLA
# represents the carrier (see subtest (e) below) but has no
# code-generation stage to inspect via `return_code` at all -- it
# builds an R formula/call object directly, not a text blob.

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
    backend = "inla",
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
    backend = "inla",
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
# (i) [removed] low_rank refusal upgrade names the registry + workaround #
# ---------------------------------------------------------------- #
#
# Deleted (not rewritten): this test called `.setup_env()` directly to
# exercise a "reserved type" refusal on the legacy
# `vm(low_rank_factor =, low_rank_scheme =)` keyword spelling.
# `.setup_env()` was the withdrawn native engine's model-environment
# builder (see NEWS.md, 0.9.3) and had zero callers from any other
# emit path even before the withdrawal -- confirmed via a caller sweep
# of R/*.R during this session's orphan-detection pass -- so this
# refusal never fired for an INLA or brms fit and is not a capability
# either active engine ever had.

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
