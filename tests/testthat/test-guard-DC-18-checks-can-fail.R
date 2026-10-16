# Class guard DC-18: a check that cannot fail.
#
# Reads every test_that() block in this directory. A block must hold at
# least one expectation that can fail: an expect_*() call, fail(), or a
# call to a helper of the same file that holds one. succeed() alone does
# not count, and a block that only calls skip() is reported as skipped.
# No expectation may be true by construction: expect_true(TRUE),
# expect_false(FALSE), the same expression on both sides of an equality,
# or a length, count or sum compared as at least zero. The release
# criteria script exited zero with every check skipped (DEF-090), and
# several checks passed whatever the code did (DEF-098, DEF-117).
#
# A check that runs but is too weak for its defect (DEF-110, DEF-151,
# DEF-160) is found only by planting the defect. That is what the
# register's falsification runs do.

.g18_calls <- function(e) {
  out <- list()
  walk <- function(x) {
    if (is.call(x)) {
      out[[length(out) + 1L]] <<- x
      for (a in as.list(x)) {
        if (!missing(a)) walk(a)
      }
    } else if (is.pairlist(x)) {
      for (a in as.list(x)) {
        if (!missing(a)) walk(a)
      }
    }
  }
  walk(e)
  out
}

.g18_name <- function(cl) {
  f <- cl[[1L]]
  if (is.symbol(f)) {
    return(as.character(f))
  }
  if (is.call(f) && as.character(f[[1L]])[1L] %in% c("::", ":::")) {
    return(as.character(f[[3L]]))
  }
  ""
}

.g18_can_fail <- function(fn, helpers) {
  (grepl("^expect_", fn) | fn == "fail" | fn %in% helpers) & fn != "succeed"
}

# Why an expectation is true by construction, or "".
.g18_tautology <- function(cl) {
  nm <- .g18_name(cl)
  a <- as.list(cl)[-1L]
  if (!length(a)) {
    return("")
  }
  if (nm == "expect_true" && isTRUE(a[[1L]])) {
    return("expect_true(TRUE)")
  }
  if (nm == "expect_false" && identical(a[[1L]], FALSE)) {
    return("expect_false(FALSE)")
  }
  if (
    nm %in%
      c("expect_equal", "expect_identical", "expect_equivalent") &&
      length(a) >= 2L &&
      identical(a[[1L]], a[[2L]])
  ) {
    return("the same expression on both sides")
  }
  counts <- c("length", "nrow", "ncol", "nchar", "sum", "NROW", "NCOL")
  if (
    nm == "expect_gte" &&
      length(a) >= 2L &&
      is.call(a[[1L]]) &&
      .g18_name(a[[1L]]) %in% counts &&
      is.numeric(a[[2L]]) &&
      a[[2L]] <= 0
  ) {
    return("a count compared as at least zero")
  }
  ""
}

test_that("DC-18 guard: every test block holds an expectation that can fail", {
  files <- list.files(
    testthat::test_path("."),
    "^test-.*\\.R$",
    full.names = TRUE
  )
  empty <- character(0)
  tautologies <- character(0)
  blocks <- 0L
  for (f in files) {
    exprs <- as.list(parse(f, keep.source = FALSE))
    helpers <- character(0)
    for (e in exprs) {
      if (
        is.call(e) &&
          .g18_name(e) %in% c("<-", "=") &&
          is.symbol(e[[2L]]) &&
          is.call(e[[3L]]) &&
          identical(e[[3L]][[1L]], quote(`function`))
      ) {
        fn <- vapply(.g18_calls(e[[3L]]), .g18_name, character(1L))
        if (any(.g18_can_fail(fn, character(0)))) {
          helpers <- c(helpers, as.character(e[[2L]]))
        }
      }
    }
    for (e in exprs) {
      if (!(is.call(e) && .g18_name(e) == "test_that")) {
        next
      }
      blocks <- blocks + 1L
      calls <- .g18_calls(e[[3L]])
      fn <- vapply(calls, .g18_name, character(1L))
      # A block that only calls skip() reports a skip, never a pass.
      skip_only <- identical(unique(fn[fn != "{"]), "skip")
      if (!any(.g18_can_fail(fn, helpers)) && !skip_only) {
        empty <- c(empty, paste0(basename(f), ": ", e[[2L]]))
      }
      for (cl in calls[grepl("^expect_", fn)]) {
        why <- .g18_tautology(cl)
        if (nzchar(why)) {
          tautologies <- c(
            tautologies,
            paste0(basename(f), ": ", e[[2L]], ": ", why)
          )
        }
      }
    }
  }
  expect_gt(blocks, 1000L)
  expect_identical(empty, character(0))
  expect_identical(tautologies, character(0))
})
