# The seven penalised-selection exports stay callable until the removal is
# announced in a public release; each refuses with selection_removed.

removed_fns <- c(
  "fb_select", "fb_select_path", "fb_select_refit", "fb_select_fold_data",
  "fb_penalty", "fb_resamples", "fb_decision"
)

test_that("each removed function is still exported", {
  ns <- getNamespaceExports("flexyBayes")
  expect_true(all(removed_fns %in% ns))
})

test_that("each removed function signals selection_removed", {
  for (nm in removed_fns) {
    f <- getExportedValue("flexyBayes", nm)
    err <- tryCatch(f(), error = function(e) e)
    expect_s3_class(err, "flexybayes_refusal_selection_removed")
    expect_s3_class(err, "flexybayes_refusal")
    expect_match(conditionMessage(err), paste0(nm, "()"), fixed = TRUE)
  }
})

test_that("the message names the vignette section", {
  err <- tryCatch(fb_select(), error = function(e) e)
  expect_match(conditionMessage(err), "flexyBayes-06-two-stage", fixed = TRUE)
  expect_match(conditionMessage(err),
               "Penalised selection on the whitened design", fixed = TRUE)
})

test_that("the two-stage vignette has the section the message points at", {
  vig <- testthat::test_path("..", "..", "vignettes", "flexyBayes-06-two-stage.Rmd")
  skip_if_not(file.exists(vig))
  expect_true(any(grepl("Penalised selection on the whitened design",
                        readLines(vig), fixed = TRUE)))
})

test_that("the Rd and API_STABILITY name the same removal release", {
  rd <- testthat::test_path("..", "..", "man", "fb_selection_removed.Rd")
  api <- testthat::test_path("..", "..", "API_STABILITY.md")
  skip_if_not(file.exists(rd) && file.exists(api))
  phrase <- "release after the first public release that announces their removal"
  rd_text <- gsub("\\s+", " ", paste(readLines(rd), collapse = " "))
  expect_true(grepl(phrase, rd_text, fixed = TRUE))
  api_lines <- readLines(api)
  rows <- api_lines[grepl("deprecated in 0.16.0", api_lines, fixed = TRUE)]
  expect_length(rows, 7L)
  expect_true(all(grepl(phrase, rows, fixed = TRUE)))
})
