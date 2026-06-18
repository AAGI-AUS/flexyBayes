# Tests for the representation registry --- ADR 0030 C4 (v0.4.0 Wave 1
# Phase 1A). The registry formalises the representation-class
# vocabulary previously spread informally across fb_preflight.R,
# fb_plan.R, and methods_truth_display.R into a single locked
# .representation_registry container at .onLoad(). The container,
# the .register_representation() / .lookup_representation() /
# .representation_class() helpers, and the v0.4.0 populate +
# lock helpers are all internal.
#
# Twelve acceptance criteria:
#
#   (a) .representation_registry exists post-.onLoad() and is locked
#   (b) all sixteen v0.4.0 entries are registered
#   (c) every registered entry carries the documented schema
#   (d) .lookup_representation() returns the registered entry
#   (e) .lookup_representation() refuses unknown name
#   (f) .lookup_representation() refuses non-string input
#   (g) .representation_class() returns the validated name string
#   (h) .register_representation() validates inputs
#   (i) .register_representation() refuses duplicate name
#   (j) .register_representation() refuses when registry is locked
#   (k) user-side assign() into the locked registry raises
#   (l) snapshot test on the registered vocabulary (vocab lock gate)

# ---------------------------------------------------------------- #
# (a) registry exists + is locked                                   #
# ---------------------------------------------------------------- #

test_that(".representation_registry exists and is locked post-load", {
  reg <- flexyBayes:::.representation_registry
  expect_true(is.environment(reg))
  expect_true(environmentIsLocked(reg))
})

# ---------------------------------------------------------------- #
# (b) all sixteen v0.4.0 entries are registered                     #
# ---------------------------------------------------------------- #

test_that("all sixteen v0.4.0 representation classes are registered", {
  reg <- flexyBayes:::.representation_registry
  expected <- c(
    "dense_cov",
    "chol_cov",
    "sparse_precision",
    "pedigree_sparse_precision",
    "block_diagonal",
    "dense_smooth",
    "sparse_smooth",
    "banded_smooth",
    "indexed_structured_known",
    "indexed_structured_estimate",
    "indexed_random_intercept",
    "dense_baseline",
    "indexed_fixed_numeric",
    "indexed_fixed_factor",
    "indexed_fixed_factor_numeric",
    # v0.4.0 Wave 1 Phase 1B (ADR 0027): low-rank truncated smooth basis
    "low_rank"
  )
  actual <- ls(envir = reg, all.names = FALSE)
  expect_setequal(actual, expected)
  expect_length(actual, 16L)
})

# ---------------------------------------------------------------- #
# (c) registered-entry schema                                       #
# ---------------------------------------------------------------- #

test_that("every registered entry carries the documented schema", {
  reg <- flexyBayes:::.representation_registry
  required_fields <- c(
    "name",
    "description",
    "category",
    "registered_in_adr",
    "since_version"
  )
  for (nm in ls(envir = reg, all.names = FALSE)) {
    entry <- get(nm, envir = reg, inherits = FALSE)
    expect_true(
      is.list(entry),
      info = paste0("entry '", nm, "' must be a list")
    )
    expect_named(
      entry,
      required_fields,
      ignore.order = TRUE,
      info = paste0("entry '", nm, "' schema")
    )
    expect_identical(
      entry$name,
      nm,
      info = paste0("entry '", nm, "' name field")
    )
    expect_true(
      is.character(entry$description) &&
        length(entry$description) == 1L,
      info = paste0("entry '", nm, "' description")
    )
    expect_true(
      is.character(entry$category) &&
        length(entry$category) == 1L,
      info = paste0("entry '", nm, "' category")
    )
  }
})

# ---------------------------------------------------------------- #
# (d) lookup returns the entry                                      #
# ---------------------------------------------------------------- #

test_that(".lookup_representation() returns the registered entry", {
  entry <- flexyBayes:::.lookup_representation("block_diagonal")
  expect_true(is.list(entry))
  expect_identical(entry$name, "block_diagonal")
  expect_identical(entry$since_version, "0.3.10")
  expect_match(entry$registered_in_adr, "ADR 0025")
  expect_match(entry$registered_in_adr, "ADR 0030")
})

# ---------------------------------------------------------------- #
# (e) lookup refuses unknown name                                   #
# ---------------------------------------------------------------- #

test_that(".lookup_representation() refuses unknown name", {
  expect_error(
    flexyBayes:::.lookup_representation("definitely_not_a_class"),
    regexp = "not a registered representation class"
  )
})

# ---------------------------------------------------------------- #
# (f) lookup refuses non-string input                               #
# ---------------------------------------------------------------- #

test_that(".lookup_representation() refuses non-string input", {
  expect_error(
    flexyBayes:::.lookup_representation(NULL),
    regexp = "non-empty single string"
  )
  expect_error(
    flexyBayes:::.lookup_representation(""),
    regexp = "non-empty single string"
  )
  expect_error(
    flexyBayes:::.lookup_representation(c("a", "b")),
    regexp = "non-empty single string"
  )
  expect_error(
    flexyBayes:::.lookup_representation(1L),
    regexp = "non-empty single string"
  )
})

# ---------------------------------------------------------------- #
# (g) .representation_class() round-trips the name                  #
# ---------------------------------------------------------------- #

test_that(".representation_class() returns the validated name", {
  expect_identical(
    flexyBayes:::.representation_class("block_diagonal"),
    "block_diagonal"
  )
  expect_identical(
    flexyBayes:::.representation_class("indexed_structured_known"),
    "indexed_structured_known"
  )
  expect_identical(
    flexyBayes:::.representation_class("dense_baseline"),
    "dense_baseline"
  )
  expect_identical(
    flexyBayes:::.representation_class("indexed_fixed_numeric"),
    "indexed_fixed_numeric"
  )
  expect_identical(
    flexyBayes:::.representation_class("indexed_fixed_factor"),
    "indexed_fixed_factor"
  )
  expect_identical(
    flexyBayes:::.representation_class("indexed_fixed_factor_numeric"),
    "indexed_fixed_factor_numeric"
  )
  expect_error(
    flexyBayes:::.representation_class("not_a_class"),
    regexp = "not a registered representation class"
  )
})

# ---------------------------------------------------------------- #
# (h) .register_representation() validates inputs                   #
# ---------------------------------------------------------------- #
#
# The post-load registry is locked, so we exercise the validation
# logic by passing bad arguments. Validation fires before the lock
# check, so these errors surface even with a locked registry.

test_that(".register_representation() validates inputs", {
  expect_error(
    flexyBayes:::.register_representation(
      name = NULL,
      description = "x",
      category = "x",
      registered_in_adr = "ADR 0030",
      since_version = "0.4.0"
    ),
    regexp = "`name` must be a non-empty single string"
  )
  expect_error(
    flexyBayes:::.register_representation(
      name = "x",
      description = c("a", "b"),
      category = "x",
      registered_in_adr = "ADR 0030",
      since_version = "0.4.0"
    ),
    regexp = "`description` must be a single string"
  )
  expect_error(
    flexyBayes:::.register_representation(
      name = "x",
      description = "x",
      category = "",
      registered_in_adr = "ADR 0030",
      since_version = "0.4.0"
    ),
    regexp = "`category` must be a non-empty single string"
  )
  expect_error(
    flexyBayes:::.register_representation(
      name = "x",
      description = "x",
      category = "x",
      registered_in_adr = 42L,
      since_version = "0.4.0"
    ),
    regexp = "`registered_in_adr` must be a single string"
  )
  expect_error(
    flexyBayes:::.register_representation(
      name = "x",
      description = "x",
      category = "x",
      registered_in_adr = "ADR 0030",
      since_version = ""
    ),
    regexp = "`since_version` must be a non-empty single string"
  )
})

# ---------------------------------------------------------------- #
# (i) .register_representation() refuses duplicate name             #
# ---------------------------------------------------------------- #
#
# The duplicate-name check happens after the lock check, so on the
# post-load locked registry the duplicate path is reachable only
# via the lock-error path. We assert the registry refuses any attempt
# to re-register a known name --- exercising whichever of the two
# refusals fires first. Both surface the same load-bearing contract:
# vocabulary is append-only.

test_that(".register_representation() refuses duplicate name", {
  err <- tryCatch(
    flexyBayes:::.register_representation(
      name = "block_diagonal",
      description = "duplicate attempt",
      category = "dense_covariance",
      registered_in_adr = "ADR 0030",
      since_version = "0.4.0"
    ),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  # Either the lock check or the duplicate check refuses; both
  # protect vocabulary-lock semantics.
  expect_match(
    conditionMessage(err),
    "is locked|is already registered"
  )
})

# ---------------------------------------------------------------- #
# (j) .register_representation() refuses when registry is locked    #
# ---------------------------------------------------------------- #

test_that(".register_representation() refuses on locked registry", {
  expect_error(
    flexyBayes:::.register_representation(
      name = "brand_new_unregistered_class",
      description = "would never land",
      category = "dense_covariance",
      registered_in_adr = "ADR 9999",
      since_version = "0.4.0"
    ),
    regexp = "registry is locked"
  )
})

# ---------------------------------------------------------------- #
# (k) user-side assign() into the locked registry raises            #
# ---------------------------------------------------------------- #

test_that("user-side assign() into the registry raises", {
  reg <- flexyBayes:::.representation_registry
  expect_error(
    assign("user_injected_entry", list(name = "x"), envir = reg),
    regexp = "locked|cannot add bindings"
  )
})

# ---------------------------------------------------------------- #
# (l) snapshot test on the registered vocabulary (vocab-lock gate)  #
# ---------------------------------------------------------------- #
#
# Any vocabulary churn surfaces as a snapshot diff that the maintainer
# must consciously accept or reject; an unintentional addition or
# removal of a representation class will not pass silently.

test_that("registered representation vocabulary snapshot", {
  reg <- flexyBayes:::.representation_registry
  vocab <- sort(ls(envir = reg, all.names = FALSE))
  expect_snapshot(vocab)
})
