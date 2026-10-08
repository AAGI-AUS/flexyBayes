# Class guard DC-19: a test passes only on the development machine or in
# one order.
#
# Reads every test file in this directory and checks four habits that
# made tests pass here and fail in a clean build, with an engine missing,
# or in another file order:
#
# 1. every package a test reaches with `pkg::` is declared in DESCRIPTION
#    (DEF-106);
# 2. a test_that() block that reaches a suggested package with `pkg::`
#    skips when that package is absent, through skip_if_not_installed(),
#    requireNamespace(), a helper that calls one of them, or a package
#    that imports it (DEF-027, DEF-087, DEF-181);
# 3. no test sets an option with a bare options() call, which outlives
#    the test (DEF-016); setup.R and helper files set suite-wide options
#    on purpose and are not read for this rule;
# 4. no test writes into the installed library (DEF-074).
#
# DESCRIPTION must also not declare parallel testthat, which stalled the
# suite with orphaned workers (DEF-097).

.g19_test_dir <- function() normalizePath(testthat::test_path("."))

.g19_root <- function() {
  candidate <- testthat::test_path("..", "..")
  if (file.exists(file.path(candidate, "DESCRIPTION"))) {
    return(normalizePath(candidate))
  }
  NA_character_
}

.g19_declared <- function(dcf) {
  fields <- intersect(c("Depends", "Imports", "Suggests"), colnames(dcf))
  out <- lapply(fields, function(f) {
    x <- trimws(sub("\\(.*", "", unlist(strsplit(dcf[1L, f], ","))))
    x[nzchar(x)]
  })
  stats::setNames(out, fields)
}

# Every sub-call of an expression, depth first.
.g19_calls <- function(e) {
  out <- list()
  walk <- function(x) {
    if (is.call(x)) {
      out[[length(out) + 1L]] <<- x
      # The function position is walked too: `pkg::f(...)` holds `pkg::f`.
      for (a in as.list(x)) {
        if (!missing(a)) walk(a)
      }
    } else if (is.pairlist(x) || is.list(x)) {
      for (a in as.list(x)) {
        if (!missing(a)) walk(a)
      }
    }
  }
  walk(e)
  out
}

.g19_fun_name <- function(cl) {
  f <- cl[[1L]]
  if (is.symbol(f)) {
    return(as.character(f))
  }
  if (is.call(f) && as.character(f[[1L]]) %in% c("::", ":::")) {
    return(as.character(f[[3L]]))
  }
  ""
}

# Packages reached with `pkg::` or `pkg:::`.
.g19_namespaces_used <- function(calls) {
  pk <- vapply(calls, function(cl) {
    f <- cl[[1L]]
    if (is.symbol(f) && as.character(f) %in% c("::", ":::")) {
      return(as.character(cl[[2L]]))
    }
    ""
  }, character(1L))
  unique(pk[nzchar(pk)])
}

# Packages an expression checks for: skip_if_not_installed("p"),
# requireNamespace("p"), system.file(package = "p"), the string literals
# of a for loop whose body calls skip_if_not_installed(<loop variable>),
# and the packages of helpers it calls.
.g19_packages_checked <- function(calls, helpers) {
  out <- character(0)
  for (cl in calls) {
    fn <- .g19_fun_name(cl)
    if (fn %in% c("skip_if_not_installed", "requireNamespace")) {
      arg <- if (length(cl) >= 2L) cl[[2L]] else NULL
      if (is.character(arg)) out <- c(out, arg)
    } else if (fn == "system.file") {
      arg <- cl[["package"]]
      if (is.character(arg)) out <- c(out, arg)
    } else if (fn == "for") {
      body_calls <- .g19_calls(cl[[4L]])
      loops <- any(vapply(body_calls, function(b) {
        identical(.g19_fun_name(b), "skip_if_not_installed") &&
          length(b) >= 2L && identical(b[[2L]], cl[[2L]])
      }, logical(1L)))
      if (loops) {
        out <- c(out, unlist(Filter(is.character, as.list(cl[[3L]]))))
      }
    }
    if (fn %in% names(helpers)) out <- c(out, helpers[[fn]])
  }
  unique(out)
}

# Top-level function definitions whose body checks for a package.
.g19_helpers <- function(exprs, helpers = list()) {
  for (e in exprs) {
    if (
      is.call(e) && .g19_fun_name(e) %in% c("<-", "=") &&
        is.symbol(e[[2L]]) && is.call(e[[3L]]) &&
        identical(e[[3L]][[1L]], quote(`function`))
    ) {
      pk <- .g19_packages_checked(.g19_calls(e[[3L]][[3L]]), helpers)
      if (length(pk)) helpers[[as.character(e[[2L]])]] <- pk
    }
  }
  helpers
}

.g19_is_test_that <- function(e) {
  is.call(e) && identical(.g19_fun_name(e), "test_that")
}

.g19_label <- function(file, e) {
  desc <- if (.g19_is_test_that(e)) paste0(" \"", e[[2L]], "\"") else ""
  paste0(basename(file), desc)
}

.g19_parse <- function(file) {
  as.list(parse(file, keep.source = FALSE))
}

.g19_test_files <- function() {
  list.files(.g19_test_dir(), "^test-.*\\.R$", full.names = TRUE)
}

.g19_helper_files <- function() {
  list.files(.g19_test_dir(), "^helper-.*\\.R$", full.names = TRUE)
}

# Packages that a check for `checked` also guarantees: their hard
# dependencies, as installed here.
.g19_implied <- function(checked) {
  if (!length(checked)) {
    return(character(0))
  }
  db <- utils::installed.packages(fields = c("Depends", "Imports"))
  have <- intersect(checked, rownames(db))
  deps <- tools::package_dependencies(
    have,
    db = db,
    which = c("Depends", "Imports"),
    recursive = TRUE
  )
  unique(c(checked, unlist(deps, use.names = FALSE)))
}

test_that("DC-19 guard: every package a test reaches is declared", {
  root <- .g19_root()
  skip_if(is.na(root), "package source tree not available")
  dcf <- read.dcf(file.path(root, "DESCRIPTION"))
  declared <- c(unlist(.g19_declared(dcf)), "flexyBayes")
  base_pkgs <- rownames(utils::installed.packages(priority = "base"))
  hits <- character(0)
  for (f in c(.g19_test_files(), .g19_helper_files())) {
    used <- .g19_namespaces_used(.g19_calls(.g19_parse(f)))
    for (p in setdiff(used, c(declared, base_pkgs))) {
      hits <- c(hits, paste0(basename(f), ": ", p))
    }
  }
  expect_identical(hits, character(0))
})

test_that("DC-19 guard: a block that reaches a suggested package skips without it", {
  root <- .g19_root()
  skip_if(is.na(root), "package source tree not available")
  dcf <- read.dcf(file.path(root, "DESCRIPTION"))
  suggested <- .g19_declared(dcf)$Suggests
  always_there <- rownames(utils::installed.packages(
    priority = c("base", "recommended")
  ))
  # testthat runs the suite and withr is used by every file; neither can
  # be absent where the tests run.
  optional <- setdiff(suggested, c(always_there, "testthat", "withr"))

  shared <- list()
  for (f in .g19_helper_files()) {
    shared <- .g19_helpers(.g19_parse(f), shared)
  }

  hits <- character(0)
  blocks <- 0L
  for (f in .g19_test_files()) {
    exprs <- .g19_parse(f)
    helpers <- .g19_helpers(exprs, shared)
    top <- exprs[!vapply(exprs, .g19_is_test_that, logical(1L))]
    file_checked <- .g19_packages_checked(
      unlist(lapply(top, .g19_calls), recursive = FALSE),
      list()
    )
    for (e in exprs[vapply(exprs, .g19_is_test_that, logical(1L))]) {
      blocks <- blocks + 1L
      calls <- .g19_calls(e[[3L]])
      used <- intersect(.g19_namespaces_used(calls), optional)
      if (!length(used)) next
      checked <- .g19_implied(c(
        .g19_packages_checked(calls, helpers),
        file_checked
      ))
      missing_skip <- setdiff(used, checked)
      if (length(missing_skip)) {
        hits <- c(hits, paste0(
          .g19_label(f, e), ": ", paste(missing_skip, collapse = ", ")
        ))
      }
    }
  }
  expect_gt(blocks, 1000L)
  expect_identical(hits, character(0))
})

test_that("DC-19 guard: no test sets an option it does not restore", {
  hits <- character(0)
  for (f in .g19_test_files()) {
    for (e in .g19_parse(f)) {
      calls <- .g19_calls(e)
      # `old <- options(...)` keeps what it replaced, to restore it.
      kept <- unlist(lapply(calls, function(cl) {
        if (
          .g19_fun_name(cl) %in% c("<-", "=") &&
            is.call(cl[[3L]]) && identical(.g19_fun_name(cl[[3L]]), "options")
        ) {
          list(cl[[3L]])
        }
      }), recursive = FALSE)
      restored <- unlist(lapply(calls, function(cl) {
        if (identical(.g19_fun_name(cl), "on.exit")) .g19_calls(cl)
      }), recursive = FALSE)
      for (cl in calls) {
        if (!identical(.g19_fun_name(cl), "options")) next
        if (!length(names(cl)) || !any(nzchar(names(cl)[-1L]))) next
        if (any(vapply(c(kept, restored), identical, logical(1L), cl))) next
        hits <- c(hits, paste0(
          .g19_label(f, e), ": ",
          paste(deparse(cl, width.cutoff = 200L), collapse = "")
        ))
      }
    }
  }
  expect_identical(hits, character(0))
})

test_that("DC-19 guard: no test writes into the installed library", {
  writers <- c(
    "saveRDS", "save", "writeLines", "write", "write.csv", "write.table",
    "fwrite", "cat", "file.create", "dir.create", "file.copy", "sink",
    "png", "pdf"
  )
  hits <- character(0)
  for (f in .g19_test_files()) {
    for (e in .g19_parse(f)) {
      calls <- .g19_calls(e)
      # Names bound to a system.file() path in the same block.
      lib_vars <- unlist(lapply(calls, function(cl) {
        if (
          .g19_fun_name(cl) %in% c("<-", "=") && is.symbol(cl[[2L]]) &&
            any(vapply(.g19_calls(cl[[3L]]), function(x) {
              identical(.g19_fun_name(x), "system.file")
            }, logical(1L)))
        ) {
          as.character(cl[[2L]])
        }
      }))
      for (cl in calls) {
        if (!.g19_fun_name(cl) %in% writers) next
        args <- as.list(cl)[-1L]
        # file.copy() may read from the library; only its destination counts.
        if (identical(.g19_fun_name(cl), "file.copy")) {
          args <- if (!is.null(cl[["to"]])) list(cl[["to"]]) else args[2L]
        }
        reads_lib <- any(vapply(args, function(a) {
          if (is.symbol(a)) {
            return(as.character(a) %in% lib_vars)
          }
          is.call(a) && (
            identical(.g19_fun_name(a), "system.file") ||
              any(vapply(.g19_calls(a), function(x) {
                identical(.g19_fun_name(x), "system.file")
              }, logical(1L))) ||
              any(all.names(a) %in% lib_vars)
          )
        }, logical(1L)))
        if (reads_lib) {
          hits <- c(hits, paste0(
            .g19_label(f, e), ": ",
            paste(deparse(cl, width.cutoff = 200L), collapse = "")
          ))
        }
      }
    }
  }
  expect_identical(hits, character(0))
})

test_that("DC-19 guard: DESCRIPTION does not declare parallel testthat", {
  root <- .g19_root()
  skip_if(is.na(root), "package source tree not available")
  dcf <- read.dcf(file.path(root, "DESCRIPTION"))
  value <- if ("Config/testthat/parallel" %in% colnames(dcf)) {
    tolower(trimws(unname(dcf[1L, "Config/testthat/parallel"])))
  } else {
    "false"
  }
  expect_identical(value, "false")
})
