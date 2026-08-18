# =============================================================================
# fb_complete_grid() -- padding a field book by hand.
#
# The helper and `na_action = "augment"` run one implementation
# (.fb_grid_complete()), which is the property that keeps them from
# completing the same trial two different ways. Three things are asserted:
#
#   * the arithmetic -- a 4 x 3 array with two rows removed comes back with
#     twelve rows and two missing responses;
#   * the refusal -- a varying design factor is not invented, because a cell
#     absent from the file is a cell whose assignment nobody recorded;
#   * the hatch -- `unused_level` opens that door, warns with the name of
#     every column it wrote to, and still refuses a varying numeric column,
#     which has no level to name.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


.fcg_grid <- function(n_row = 4L, n_col = 3L, seed = 4242L) {
  set.seed(seed)
  g <- expand.grid(col = factor(seq_len(n_col)), row = factor(seq_len(n_row)))
  g$yield <- stats::rnorm(nrow(g), mean = 5)
  g
}


# ---------------------------------------------------------------- #
# 1. The arithmetic                                                 #
# ---------------------------------------------------------------- #

test_that("two removed cells come back with a missing response", {
  g <- .fcg_grid()
  holed <- g[-c(2L, 7L), , drop = FALSE]
  expect_identical(nrow(holed), 10L)

  out <- fb_complete_grid(holed, ~ row * col, response = "yield")
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 12L)
  expect_identical(sum(is.na(out$yield)), 2L)
  # Every node of the array is present exactly once.
  expect_identical(
    nrow(unique(out[, c("row", "col")])),
    12L
  )
  # The observed rows are unchanged and come first.
  expect_identical(out[seq_len(10L), , drop = FALSE], holed,
    ignore_attr = "row.names"
  )
  # The reinstated cells are the two that were removed.
  gone <- out[is.na(out$yield), c("row", "col")]
  expect_identical(
    sort(paste(gone$row, gone$col)),
    sort(paste(g$row[c(2L, 7L)], g$col[c(2L, 7L)]))
  )
})

test_that("a complete grid is returned unchanged", {
  g <- .fcg_grid()
  out <- fb_complete_grid(g, ~ row * col, response = "yield")
  expect_identical(out, g, ignore_attr = "row.names")
})

test_that("the index formula is read as a variable set", {
  g <- .fcg_grid()
  holed <- g[-2L, , drop = FALSE]
  # `*`, `+` and `:` name the same two variables; the operator is not a
  # model term here, so all three complete the same array.
  a <- fb_complete_grid(holed, ~ row * col, response = "yield")
  b <- fb_complete_grid(holed, ~ row + col, response = "yield")
  d <- fb_complete_grid(holed, ~ row:col, response = "yield")
  expect_identical(a, b)
  expect_identical(a, d)
  expect_identical(nrow(a), 12L)
})


# ---------------------------------------------------------------- #
# 2. The refusal                                                    #
# ---------------------------------------------------------------- #

test_that("a varying design factor refuses rather than being invented", {
  g <- .fcg_grid()
  g$geno <- factor(rep(c("A", "B", "C"), length.out = nrow(g)))
  holed <- g[-c(2L, 7L), , drop = FALSE]

  err <- tryCatch(
    fb_complete_grid(holed, ~ row * col, response = "yield"),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_augment_cell_not_determinable")
  msg <- conditionMessage(err)
  expect_match(msg, "geno", fixed = TRUE)
  expect_match(msg, "LANCER", fixed = TRUE)
  expect_match(msg, "nin89", fixed = TRUE)
  expect_match(msg, "unused_level", fixed = TRUE)
})

test_that("a constant column is determined and does not refuse", {
  g <- .fcg_grid()
  g$site <- factor("Roseworthy")
  holed <- g[-2L, , drop = FALSE]
  out <- fb_complete_grid(holed, ~ row * col, response = "yield")
  expect_identical(nrow(out), 12L)
  expect_true(all(out$site == "Roseworthy"))
})


# ---------------------------------------------------------------- #
# 3. The hatch                                                      #
# ---------------------------------------------------------------- #

test_that("unused_level fills the varying factor and names the column", {
  g <- .fcg_grid()
  g$geno <- factor(rep(c("A", "B", "C"), length.out = nrow(g)))
  holed <- g[-c(2L, 7L), , drop = FALSE]

  expect_warning(
    out <- fb_complete_grid(
      holed, ~ row * col,
      response = "yield", unused_level = "UNSOWN"
    ),
    "geno"
  )
  expect_identical(nrow(out), 12L)
  expect_identical(sum(is.na(out$yield)), 2L)
  expect_true(is.factor(out$geno))
  expect_true("UNSOWN" %in% levels(out$geno))
  expect_identical(sum(out$geno == "UNSOWN"), 2L)
  # The observed rows keep their own varieties.
  expect_identical(
    as.character(out$geno[!is.na(out$yield)]),
    as.character(holed$geno)
  )
})

test_that("a varying numeric column refuses whatever unused_level says", {
  # The hatch is for design factors. A measured covariate has no level to
  # name, so filling one would be a fabricated measurement rather than a
  # named gap in the design.
  g <- .fcg_grid()
  g$moisture <- stats::runif(nrow(g))
  holed <- g[-c(2L, 7L), , drop = FALSE]

  err <- tryCatch(
    fb_complete_grid(
      holed, ~ row * col,
      response = "yield", unused_level = "UNSOWN"
    ),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_augment_cell_not_determinable")
  expect_match(conditionMessage(err), "moisture", fixed = TRUE)
  expect_match(conditionMessage(err), "no level to name", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# 4. Argument contract                                              #
# ---------------------------------------------------------------- #

test_that("the arguments are checked by name", {
  g <- .fcg_grid()
  expect_error(fb_complete_grid(as.matrix(1:4), ~ row * col, "yield"),
    "must be a data frame"
  )
  expect_error(fb_complete_grid(g, "row * col", "yield"),
    "one-sided formula"
  )
  expect_error(fb_complete_grid(g, yield ~ row * col, "yield"),
    "one-sided formula"
  )
  expect_error(fb_complete_grid(g, ~ row * col), "single string")
  expect_error(fb_complete_grid(g, ~ row * plot, "yield"), "not a column")
  expect_error(fb_complete_grid(g, ~ row * col, "harvest"), "not a column")
  expect_error(
    fb_complete_grid(g, ~ row * col, "yield", unused_level = c("a", "b")),
    "single non-missing value"
  )
})


# ---------------------------------------------------------------- #
# 5. One implementation, two entry points (E-p4)                    #
# ---------------------------------------------------------------- #

test_that("the helper and na_action = 'augment' complete the same grid", {
  g <- .fcg_grid(n_row = 8L, n_col = 6L)
  holed <- g[-c(5L, 19L), , drop = FALSE]

  by_hand <- fb_complete_grid(holed, ~ row * col, response = "yield")
  ir <- flexyBayes:::new_fb_terms(
    response = "yield", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(list(
      type = "ar1_spatial", row_var = "row", col_var = "col"
    )),
    residual_terms = list(list(type = "units")),
    source = "asreml"
  )
  automatic <- flexyBayes:::.fb_complete_design_grid(
    ir, holed, c("row", "col")
  )

  expect_identical(automatic$added, 2L)
  expect_identical(nrow(automatic$data), nrow(by_hand))
  expect_identical(
    sum(is.na(automatic$data$yield)),
    sum(is.na(by_hand$yield))
  )
  expect_identical(
    automatic$data[, c("row", "col", "yield")],
    by_hand[, c("row", "col", "yield")]
  )
})
