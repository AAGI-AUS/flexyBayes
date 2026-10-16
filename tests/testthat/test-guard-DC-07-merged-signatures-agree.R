# Class guard DC-07: a merge drops an argument or keeps a stale check.
#
# Reads the package source:
#
# 1. every call to a package function names only arguments that function
#    takes, matched exactly, so a signature that lost an argument its
#    callers still pass is found before the call runs (DEF-178: a merge
#    dropped decision_resolution from the optimiser's signature);
# 2. every refusal code is registered once (DEF-122: two merged units
#    registered the same codes).
#
# A function with `...` in its signature is not checked for 1.

.g07_root <- function() {
  candidate <- testthat::test_path("..", "..")
  if (dir.exists(file.path(candidate, "R"))) {
    return(normalizePath(candidate))
  }
  NA_character_
}

.g07_walk <- function(x, visit) {
  if (is.call(x)) {
    visit(x)
    for (a in as.list(x)) {
      if (!missing(a)) .g07_walk(a, visit)
    }
  } else if (is.pairlist(x)) {
    for (a in as.list(x)) {
      if (!missing(a)) .g07_walk(a, visit)
    }
  }
}

test_that("DC-07 guard: every call names only arguments its callee takes", {
  root <- .g07_root()
  skip_if(is.na(root), "package source tree not available")
  ns <- asNamespace("flexyBayes")
  funs <- Filter(
    function(n) {
      is.function(get(n, envir = ns)) && !is.primitive(get(n, envir = ns))
    },
    ls(ns, all.names = TRUE)
  )
  formals_of <- lapply(stats::setNames(funs, funs), function(n) {
    names(formals(get(n, envir = ns)))
  })
  hits <- character(0)
  checked <- 0L
  for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) {
    for (e in parse(f, keep.source = FALSE)) {
      where <- if (is.call(e) && length(e) >= 2L && is.symbol(e[[2L]])) {
        as.character(e[[2L]])
      } else {
        ""
      }
      .g07_walk(e, function(x) {
        fn <- x[[1L]]
        if (!is.symbol(fn) || !as.character(fn) %in% funs) {
          return(invisible())
        }
        fm <- formals_of[[as.character(fn)]]
        if (is.null(fm) || "..." %in% fm) {
          return(invisible())
        }
        given <- names(as.list(x)[-1L])
        given <- given[!is.na(given) & nzchar(given)]
        checked <<- checked + 1L
        bad <- setdiff(given, fm)
        if (length(bad)) {
          hits <<- c(
            hits,
            paste0(
              basename(f),
              " (",
              where,
              "): ",
              as.character(fn),
              "() has no ",
              paste0("`", bad, "`", collapse = ", ")
            )
          )
        }
      })
    }
  }
  expect_gt(checked, 1000L)
  expect_identical(hits, character(0))
})

test_that("DC-07 guard: every refusal code is registered once", {
  root <- .g07_root()
  skip_if(is.na(root), "package source tree not available")
  codes <- character(0)
  for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) {
    for (e in parse(f, keep.source = FALSE)) {
      .g07_walk(e, function(x) {
        if (identical(x[[1L]], quote(.register_refusal))) {
          code <- x[["reason_code"]]
          if (is.character(code)) codes <<- c(codes, code)
        }
      })
    }
  }
  expect_gt(length(codes), 100L)
  expect_identical(unique(codes[duplicated(codes)]), character(0))
})
