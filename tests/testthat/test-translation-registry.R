# =============================================================================
# The translation registry and its generated surfaces.
#
# Each ASReml spelling lowers in one of three ways: an exact structural
# translation, a covariance-family equivalent with different identification
# and priors, or an analogue built from a different basis. The registry,
# .fb_translation_registry(), says which. inst/TRANSLATION.md and
# inst/translation/translation_registry.csv are generated from it by
# tools/generate_translation_registry.R. These tests pin the vocabulary, pin
# the generated files to the table, and check that every `anchor` names a
# test file that exists.
# =============================================================================

.tr_root <- function() {
  root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE
  )
  if (file.exists(file.path(root, "DESCRIPTION"))) root else NA_character_
}


test_that("the registry's class and vintage vocabularies are closed", {
  reg <- .fb_translation_registry()
  expect_true(all(reg$class %in% .FB_TRANSLATION_CLASSES))
  expect_true(all(reg$vintage %in% .FB_TRANSLATION_VINTAGES))
  expect_true(all(nzchar(reg$structure)))
  expect_true(all(nzchar(reg$spelling)))
  expect_true(all(nzchar(reg$oracle)))
})

test_that("every registry row names a test file that exists", {
  root <- .tr_root()
  if (is.na(root)) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }
  reg <- .fb_translation_registry()
  anchors <- unique(reg$anchor)
  missing <- anchors[
    !file.exists(file.path(root, "tests", "testthat", anchors))
  ]
  expect_equal(missing, character(0))
})

test_that("a qualified row states its precondition", {
  # `conditional` means covariance-family equivalent only where a stated
  # precondition holds; `refused` means the request terminates and the row
  # has to say on what grounds. Both are promises carrying a qualification,
  # and a blank condition is that qualification removed.
  #
  # Asserted over both classes, since `conditional` alone may be empty.
  reg <- .fb_translation_registry()
  qualified <- reg[reg$class %in% c("conditional", "refused"), , drop = FALSE]
  expect_gt(nrow(qualified), 0L)
  expect_true(all(nzchar(qualified$condition)))
})

test_that("the ASReml 3 aliases are labelled as aliases, not as current syntax", {
  # ASReml 4 removed residual-side at() in favour of dsum(), and ped() in
  # favour of vm(). Both still work here, but neither may be presented as
  # current ASReml 4 syntax.
  reg <- .fb_translation_registry()
  aliases <- reg[reg$vintage == "ASReml 3 alias", , drop = FALSE]
  expect_true(any(grepl("at(env):units", aliases$spelling, fixed = TRUE)))
  expect_true(any(grepl("ped(g, A)", aliases$spelling, fixed = TRUE)))
  expect_true(all(grepl("ASReml 4", aliases$condition, fixed = TRUE)))
})

test_that("spl() is registered as an analogue, not as a translation", {
  reg <- .fb_translation_registry()
  row <- reg[grepl("spl(x)", reg$spelling, fixed = TRUE), , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_equal(row$class[[1]], "analogue")
})

test_that("inst/TRANSLATION.md matches a fresh render", {
  skip_if_covr()
  root <- .tr_root()
  if (is.na(root)) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }
  path <- file.path(root, "inst", "TRANSLATION.md")
  committed <- readLines(path, warn = FALSE)
  b <- which(committed == .FB_TRANSLATION_MARKER_BEGIN)
  e <- which(committed == .FB_TRANSLATION_MARKER_END)
  expect_length(b, 1L)
  expect_length(e, 1L)
  expect_equal(
    committed[seq.int(b, e)],
    .fb_translation_markdown(markers = TRUE)
  )
})

test_that("the shipped CSV matches the table", {
  skip_if_covr()
  root <- .tr_root()
  if (is.na(root)) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }
  path <- file.path(root, "inst", "translation", "translation_registry.csv")
  expect_true(file.exists(path))
  on_disk <- utils::read.csv(path, stringsAsFactors = FALSE)
  reg <- .fb_translation_registry()
  expect_equal(nrow(on_disk), nrow(reg))
  expect_equal(names(on_disk), names(reg))
  expect_equal(on_disk$class, reg$class)
  expect_equal(on_disk$structure, reg$structure)
})
