# Tests for the internal Stage 2 MVP `.fb_dataset()` wrapper
# (ADR 0021 / v0.3.0). Covers the three constructor paths
# (data.frame, data.table, metadata-only), the frozen-dictionary
# contract, the `origin` slot, the cached metadata slots, and the
# refusal surface on missing / malformed args.

test_that(".fb_dataset(data.frame): basic wrap + dictionary freeze", {
  df <- data.frame(
    y = c(1.5, 2.0, 3.1),
    g = factor(c("a", "b", "a"), levels = c("a", "b", "c")),
    x = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  ds <- .fb_dataset(df)

  expect_s3_class(ds, "fb_dataset")
  expect_identical(ds$n_rows, 3L)
  expect_identical(ds$origin, "data.frame")
  expect_true(inherits(ds$data, "data.table"))
  # col_types named character; one entry per column
  expect_identical(sort(names(ds$col_types)), sort(c("y", "g", "x")))
  expect_identical(unname(ds$col_types[["y"]]), "double")
  expect_identical(unname(ds$col_types[["g"]]), "factor")
  expect_identical(unname(ds$col_types[["x"]]), "double")
  # factor dictionary frozen at full level set, not the rows
  expect_identical(ds$dictionaries$g, c("a", "b", "c"))
})

test_that(".fb_dataset(data.table): identity origin + dictionary preserved", {
  dt <- data.table::data.table(
    y = c(1.5, 2.0),
    f = factor(c("x", "y"), levels = c("x", "y", "z"))
  )
  ds <- .fb_dataset(dt)

  expect_identical(ds$origin, "data.table")
  expect_true(inherits(ds$data, "data.table"))
  expect_identical(ds$dictionaries$f, c("x", "y", "z"))
  # No mutation of caller's data.table (copy taken)
  data.table::set(ds$data, j = "y", value = c(99, 99))
  expect_identical(dt$y, c(1.5, 2.0))
})

test_that(".fb_dataset(data = NULL, ...): metadata-only path", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e7,
    col_types = list(y = "double", x = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(100L)))
  )

  expect_s3_class(ds, "fb_dataset")
  expect_null(ds$data)
  expect_identical(ds$origin, "metadata-only")
  expect_identical(ds$n_rows, 10000000L)
  expect_identical(length(ds$dictionaries$g), 100L)
  expect_true(.fb_dataset_is_metadata(ds))
})

test_that(".fb_dataset(data = NULL): refuses without n_rows", {
  expect_error(
    .fb_dataset(data = NULL),
    regexp = "n_rows"
  )
})

test_that(".fb_dataset(data = NULL): refuses without col_types", {
  expect_error(
    .fb_dataset(data = NULL, n_rows = 100),
    regexp = "col_types"
  )
})

test_that(".fb_dataset(data = NULL): rejects malformed col_types / dictionaries", {
  # Unnamed col_types
  expect_error(
    .fb_dataset(data = NULL, n_rows = 10, col_types = c("double", "factor")),
    regexp = "named"
  )
  # Unnamed dictionaries (when non-empty)
  expect_error(
    .fb_dataset(
      data = NULL,
      n_rows = 10,
      col_types = list(g = "factor"),
      dictionaries = list(as.character(1:5))
    ),
    regexp = "named"
  )
  # Non-list dictionaries
  expect_error(
    .fb_dataset(
      data = NULL,
      n_rows = 10,
      col_types = list(g = "factor"),
      dictionaries = "not_a_list"
    ),
    regexp = "named list"
  )
})

test_that(".fb_dataset() rejects non-data.frame data input", {
  expect_error(
    .fb_dataset(data = list(x = 1:5)),
    regexp = "data\\.frame"
  )
  expect_error(
    .fb_dataset(data = 1:10),
    regexp = "data\\.frame"
  )
})

test_that(".fb_dataset(): character columns also get frozen dictionaries", {
  df <- data.frame(
    y = 1:3,
    site = c("PIRSA", "AGRF", "PIRSA"),
    stringsAsFactors = FALSE
  )
  ds <- .fb_dataset(df)

  expect_identical(unname(ds$col_types[["site"]]), "character")
  expect_identical(ds$dictionaries$site, c("AGRF", "PIRSA"))
})

test_that(".fb_dataset_levels() + .fb_dataset_type() reads via dictionary", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e6,
    col_types = list(y = "double", g = "factor", site = "character"),
    dictionaries = list(g = as.character(seq_len(50L)), site = c("a", "b", "c"))
  )

  expect_identical(.fb_dataset_levels(ds, "g"), 50L)
  expect_identical(.fb_dataset_levels(ds, "site"), 3L)
  expect_identical(.fb_dataset_levels(ds, "y"), NA_integer_)
  expect_identical(.fb_dataset_levels(ds, "absent"), NA_integer_)

  expect_identical(.fb_dataset_type(ds, "y"), "double")
  expect_identical(.fb_dataset_type(ds, "g"), "factor")
  expect_identical(.fb_dataset_type(ds, "absent"), NA_character_)
})

test_that("format.fb_dataset() / print.fb_dataset() produce diagnostic output", {
  df <- data.frame(
    y = 1:3,
    g = factor(letters[1:3])
  )
  ds <- .fb_dataset(df)

  fmt <- format(ds)
  expect_true(grepl("<fb_dataset>", fmt, fixed = TRUE))
  expect_true(grepl("n_rows = 3", fmt, fixed = TRUE))
  expect_true(grepl("origin = data.frame", fmt, fixed = TRUE))

  printed <- capture.output(print(ds))
  expect_true(any(grepl("col_types", printed, fixed = TRUE)))
  expect_true(any(grepl("g:factor", printed, fixed = TRUE)))
  expect_true(any(grepl("dicts", printed, fixed = TRUE)))
  expect_true(any(grepl("g(3)", printed, fixed = TRUE)))
})

test_that(".fb_dataset is internal -- no exported binding", {
  # The dotted-name internal contract: `.fb_dataset` is reachable
  # via `flexyBayes:::.fb_dataset` (triple-colon, internal access)
  # but NOT via `flexyBayes::.fb_dataset` (double-colon, exported
  # access). Asserts the export gate has not been opened.
  ns <- asNamespace("flexyBayes")
  exp <- getNamespaceExports(ns)
  expect_false(".fb_dataset" %in% exp)
  expect_false("fb_dataset" %in% exp)
  # Internal reachability holds
  expect_true(exists(".fb_dataset", envir = ns, inherits = FALSE))
})
