# =============================================================================
# The translation registry and its generated surfaces (audit 2026-09-11,
# F02 / F05 / F06 / F07, work package 0).
#
# The audit's central finding was not a bug but a category error: "ASReml
# syntax" was used as a blanket claim over lowerings of three different kinds
# -- exact structural translations, covariance-family equivalents with
# different identification and priors, and analogues built from a different
# basis. A user choosing between this package and ASReml cannot act on a claim
# that does not say which they are getting.
#
# .fb_translation_registry() is where that is now written down once.
# inst/TRANSLATION.md and inst/translation/translation_registry.csv are
# generated from it by tools/generate_translation_registry.R. These tests pin
# the vocabulary closed, pin the generated surfaces to the table, and pin the
# `anchor` column to test files that exist -- seven anchors in the first draft
# of this table named test files that did not, which is the same defect class
# the registry exists to close.
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
  missing <- anchors[!file.exists(file.path(root, "tests", "testthat", anchors))]
  expect_equal(missing, character(0))
})

test_that("a conditional row states its precondition", {
  # `conditional` means covariance-family equivalent only where a stated
  # precondition holds. A conditional row with an empty condition is a
  # promise with the qualification removed.
  reg <- .fb_translation_registry()
  conditional <- reg[reg$class == "conditional", , drop = FALSE]
  expect_gt(nrow(conditional), 0L)
  expect_true(all(nzchar(conditional$condition)))
})

test_that("the ASReml 3 aliases are labelled as aliases, not as current syntax", {
  # ASReml 4 removed residual-side at() in favour of dsum(), and ped() in
  # favour of vm(). Both still work here; neither may be presented as
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
