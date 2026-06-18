# fb_cov() constructor noun + vm/ped migration + legacy deprecation
# (ADR 0030 C3; v0.4.0 Wave 2 Phase 2A).
#
# Coverage:
#   (a) standalone constructor: classed object + attribute schema per
#       carrier type.
#   (b) construction-time refusals (unknown type, missing matrix,
#       blocks-not-a-list, low_rank scheme rules).
#   (c) inline `cov = fb_cov(...)` formula parsing maps onto the
#       existing cov_representation IR slot for every carrier type.
#   (d) the inline path is codegen-equivalent to the legacy keyword
#       form (the migration is behaviour-preserving).
#   (e) the four legacy keyword forms fire a lifecycle deprecation
#       warning and continue to parse.

# Small SPD fixtures reused across blocks.
.fb_cov_spd <- function(n = 3L) {
  M <- matrix(stats::rnorm(n * n), n, n)
  crossprod(M) + diag(n)
}

.fb_cov_fixture_data <- function(n_geno = 6L, reps = 4L) {
  data.frame(
    geno = factor(rep(seq_len(n_geno), each = reps)),
    yield = stats::rnorm(n_geno * reps)
  )
}

# ---------------------------------------------------------------- #
# (a) Standalone constructor: classed object + attribute schema.    #
# ---------------------------------------------------------------- #

test_that("fb_cov() returns a classed object with the documented attribute schema", {
  G <- .fb_cov_spd(3L)
  cov <- fb_cov(G, type = "dense")

  expect_s3_class(cov, "fb_cov")
  expect_true(is_fb_cov(cov))
  expect_true(is.list(cov))
  expect_equal(cov$type, "dense")
  expect_identical(cov$M, G)
  expect_equal(attr(cov, "type"), "dense")
  expect_equal(attr(cov, "representation_class"), "dense_cov")
  expect_type(attr(cov, "validation_summary"), "character")
})

test_that("fb_cov() maps each carrier type to its locked representation class", {
  G <- .fb_cov_spd(3L)
  L <- t(chol(G))
  Q <- solve(G)
  expect_equal(
    attr(fb_cov(G, type = "dense"), "representation_class"),
    "dense_cov"
  )
  expect_equal(
    attr(fb_cov(L, type = "chol"), "representation_class"),
    "chol_cov"
  )
  expect_equal(
    attr(fb_cov(Q, type = "precision"), "representation_class"),
    "sparse_precision"
  )
  expect_equal(
    attr(fb_cov(list(G, G), type = "blocks"), "representation_class"),
    "block_diagonal"
  )
  expect_equal(
    attr(
      fb_cov(matrix(0, 5L, 2L), type = "low_rank", scheme = "low_rank_smooth"),
      "representation_class"
    ),
    "low_rank"
  )
})

test_that("fb_cov() validation_summary reflects the carrier's structure", {
  G <- .fb_cov_spd(3L)
  L <- t(chol(G))
  expect_match(
    attr(fb_cov(L, type = "chol"), "validation_summary"),
    "lower-triangular"
  )
  expect_match(
    attr(fb_cov(G, type = "precision"), "validation_summary"),
    "symmetric"
  )
  expect_match(
    attr(
      fb_cov(list(G, .fb_cov_spd(2L)), type = "blocks"),
      "validation_summary"
    ),
    "2 blocks"
  )
})

test_that("fb_cov() carries optional levels metadata", {
  G <- .fb_cov_spd(3L)
  cov <- fb_cov(G, type = "dense", levels = c("a", "b", "c"))
  expect_equal(attr(cov, "levels"), c("a", "b", "c"))
})

test_that("print.fb_cov() renders type, representation, and carrier summary", {
  G <- .fb_cov_spd(3L)
  out <- utils::capture.output(print(fb_cov(t(chol(G)), type = "chol")))
  expect_true(any(grepl("<fb_cov>", out, fixed = TRUE)))
  expect_true(any(grepl("type = \"chol\"", out, fixed = TRUE)))
  expect_true(any(grepl("chol_cov", out, fixed = TRUE)))
})

# ---------------------------------------------------------------- #
# (b) Construction-time refusals.                                   #
# ---------------------------------------------------------------- #

test_that("fb_cov() refuses an unknown carrier type", {
  err <- tryCatch(
    fb_cov(.fb_cov_spd(3L), type = "bogus"),
    flexybayes_refusal = identity
  )
  expect_s3_class(err, "flexybayes_refusal")
  expect_equal(err$reason_code, "fb_cov_type_unknown")
})

test_that("fb_cov() refuses a missing carrier matrix", {
  err <- tryCatch(fb_cov(type = "dense"), flexybayes_refusal = identity)
  expect_equal(err$reason_code, "fb_cov_missing_matrix")
})

test_that("fb_cov(type = 'blocks') refuses a non-list carrier", {
  err <- tryCatch(
    fb_cov(.fb_cov_spd(3L), type = "blocks"),
    flexybayes_refusal = identity
  )
  expect_equal(err$reason_code, "blocks_not_a_list")
})

test_that("fb_cov(type = 'low_rank') requires a scheme", {
  err <- tryCatch(
    fb_cov(matrix(0, 5L, 2L), type = "low_rank"),
    flexybayes_refusal = identity
  )
  expect_equal(err$reason_code, "low_rank_scheme_required")
})

test_that("fb_cov(type = 'low_rank') refuses an unregistered scheme", {
  err <- tryCatch(
    fb_cov(matrix(0, 5L, 2L), type = "low_rank", scheme = "not_a_scheme"),
    flexybayes_refusal = identity
  )
  expect_equal(err$reason_code, "approximation_scheme_unknown")
})

# ---------------------------------------------------------------- #
# (c) Inline `cov = fb_cov(...)` formula parsing -> IR slot.        #
# ---------------------------------------------------------------- #

test_that("vm(geno, cov = fb_cov(L, type = 'chol')) parses to the chol slot", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(L, type = "chol")),
    dat
  )
  expect_equal(terms[[1]]$type, "vm")
  expect_equal(terms[[1]]$cov_representation$format, "chol")
  expect_equal(terms[[1]]$cov_representation$data, "L")
  expect_true(is.na(terms[[1]]$mat))
})

test_that("vm(geno, cov = fb_cov(G, type = 'dense')) keeps the dense slot + mat", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Gmat, type = "dense")),
    dat
  )
  expect_equal(terms[[1]]$cov_representation$format, "dense")
  expect_equal(terms[[1]]$cov_representation$data, "Gmat")
  expect_equal(terms[[1]]$mat, "Gmat")
})

test_that("vm() cov = fb_cov() defaults type to dense", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(~ vm(geno, cov = fb_cov(Gmat)), dat)
  expect_equal(terms[[1]]$cov_representation$format, "dense")
  expect_equal(terms[[1]]$cov_representation$data, "Gmat")
})

test_that("vm(geno, cov = fb_cov(Q, type = 'precision')) parses to the precision slot", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Q, type = "precision")),
    dat
  )
  expect_equal(terms[[1]]$cov_representation$format, "precision")
  expect_equal(terms[[1]]$cov_representation$data, "Q")
})

test_that("vm(geno, cov = fb_cov(Bs, type = 'blocks')) parses to the blocks slot", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(Bs, type = "blocks")),
    dat
  )
  expect_equal(terms[[1]]$cov_representation$format, "blocks")
  expect_equal(terms[[1]]$cov_representation$data, "Bs")
})

test_that("vm() cov = fb_cov(type = 'low_rank') carries the scheme", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ vm(geno, cov = fb_cov(U, type = "low_rank", scheme = "low_rank_smooth")),
    dat
  )
  expect_equal(terms[[1]]$cov_representation$format, "low_rank")
  expect_equal(terms[[1]]$cov_representation$data, "U")
  expect_equal(terms[[1]]$cov_representation$scheme, "low_rank_smooth")
})

test_that("ped() cov = fb_cov(type = 'precision', sparse_precision = TRUE) folds to the pedigree route", {
  dat <- data.frame(animal = factor(1:5))
  terms <- flexyBayes:::.parse_formula(
    ~ ped(animal, cov = fb_cov(A, type = "precision", sparse_precision = TRUE)),
    dat
  )
  expect_equal(
    terms[[1]]$cov_representation$format,
    "pedigree_sparse_precision"
  )
  expect_equal(terms[[1]]$cov_representation$data, "A")
})

test_that("vm() cov = <bare symbol> refuses with cov_arg_not_fb_cov", {
  dat <- data.frame(geno = factor(1:5))
  err <- tryCatch(
    flexyBayes:::.parse_formula(~ vm(geno, cov = some_object), dat),
    flexybayes_structured_cov_refusal = identity
  )
  expect_equal(err$reason_code, "cov_arg_not_fb_cov")
})

test_that("vm() cov = fb_cov(type = '<bogus>') refuses with fb_cov_type_unknown", {
  dat <- data.frame(geno = factor(1:5))
  err <- tryCatch(
    flexyBayes:::.parse_formula(
      ~ vm(geno, cov = fb_cov(Gmat, type = "bogus")),
      dat
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_equal(err$reason_code, "fb_cov_type_unknown")
})

test_that("vm() cov = fb_cov(type = 'low_rank') without scheme refuses at parse", {
  dat <- data.frame(geno = factor(1:5))
  err <- tryCatch(
    flexyBayes:::.parse_formula(
      ~ vm(geno, cov = fb_cov(U, type = "low_rank")),
      dat
    ),
    flexybayes_structured_cov_refusal = identity
  )
  expect_equal(err$reason_code, "low_rank_scheme_required")
})

# ---------------------------------------------------------------- #
# (d) Codegen equivalence: fb_cov() form == legacy keyword form.    #
# ---------------------------------------------------------------- #

test_that("codegen: cov = fb_cov(L, type = 'chol') emits the same square root as chol = L", {
  dat <- .fb_cov_fixture_data()
  G <- diag(6) + 0.1
  L <- t(chol(G))

  code_new <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(L_chol, type = "chol")),
    data = dat,
    known_matrices = list(L_chol = L),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl("L_G_geno <- as.matrix(L_chol)", code_new, fixed = TRUE))
  expect_false(grepl("t(chol(", code_new, fixed = TRUE))
})

test_that("codegen: cov = fb_cov(Q, type = 'precision') emits solve(chol(Q))", {
  skip_if_not_installed("Matrix")
  dat <- .fb_cov_fixture_data()
  Q <- solve(diag(6) + 0.1)
  code_new <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Qprec, type = "precision")),
    data = dat,
    known_matrices = list(Qprec = Q),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl(
    "L_G_geno <- as.matrix(solve(chol(Qprec)))",
    code_new,
    fixed = TRUE
  ))
})

test_that("codegen: cov = fb_cov(G, type = 'dense') keeps the t(chol()) wrap", {
  dat <- .fb_cov_fixture_data()
  G <- diag(6) + 0.1
  code_new <- flexybayes(
    yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(Gmat, type = "dense")),
    data = dat,
    known_matrices = list(Gmat = G),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(grepl("L_G_geno <- t(chol(Gmat))", code_new, fixed = TRUE))
})

# ---------------------------------------------------------------- #
# (e) Legacy keyword deprecation (ADR 0030 §5).                     #
# ---------------------------------------------------------------- #

test_that("the four legacy vm() keyword carriers fire a lifecycle deprecation warning", {
  withr::local_options(lifecycle_verbosity = "warning")
  dat <- data.frame(geno = factor(1:5))

  lifecycle::expect_deprecated(
    flexyBayes:::.parse_formula(~ vm(geno, chol = L), dat)
  )
  lifecycle::expect_deprecated(
    flexyBayes:::.parse_formula(~ vm(geno, precision = Q), dat)
  )
  lifecycle::expect_deprecated(
    flexyBayes:::.parse_formula(~ vm(geno, blocks = Bs), dat)
  )
  lifecycle::expect_deprecated(
    flexyBayes:::.parse_formula(
      ~ vm(geno, low_rank_factor = U, low_rank_scheme = "pca"),
      dat
    )
  )
})

test_that("a deprecated legacy keyword still parses to the correct IR slot", {
  withr::local_options(lifecycle_verbosity = "warning")
  dat <- data.frame(geno = factor(1:5))
  terms <- suppressWarnings(
    flexyBayes:::.parse_formula(~ vm(geno, chol = L), dat)
  )
  expect_equal(terms[[1]]$cov_representation$format, "chol")
  expect_equal(terms[[1]]$cov_representation$data, "L")
})

test_that("the deprecation message names the fb_cov() migration without a 'report the issue' footer", {
  withr::local_options(lifecycle_verbosity = "warning")
  dat <- data.frame(geno = factor(1:5))
  w <- tryCatch(
    flexyBayes:::.parse_formula(~ vm(geno, chol = L), dat),
    warning = identity
  )
  expect_s3_class(w, "lifecycle_warning_deprecated")
  expect_match(conditionMessage(w), "fb_cov")
  expect_false(grepl("report the issue", conditionMessage(w)))
})

test_that("the new cov = fb_cov() form does NOT fire a deprecation warning", {
  withr::local_options(lifecycle_verbosity = "warning")
  dat <- data.frame(geno = factor(1:5))
  expect_no_warning(
    flexyBayes:::.parse_formula(
      ~ vm(geno, cov = fb_cov(L, type = "chol")),
      dat
    )
  )
})
