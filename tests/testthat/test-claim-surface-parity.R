# =============================================================================
# Claim-surface parity (audit 2026-09-11, F09 + F13).
#
# The audit found four numbers and names that had drifted apart: the package
# help said eight vignettes ship, `cran-comments.md` said eleven, four
# actually ship; the package help and the startup banner still called the
# cross-engine comparison the signature feature after it was unexported at
# 0.10.0.9000; and `vpredict()` was exported but absent from both the pkgdown
# reference index and API_STABILITY.md.
#
# Correcting the four was the easy half. This file is the other half: the
# project's own standing rule is that every claim surface gets a generator, a
# parity check or a grep-gate, and these had none. A count written into prose
# is the first thing to drift.
#
# These gates read the source tree, which is not present when the tests run
# from an installed package (`_pkgdown.yml` and `cran-comments.md` are
# .Rbuildignore'd). They skip with a named reason there rather than passing
# vacuously -- a gate's silence is not a pass.
# =============================================================================

.cs_root <- function() {
  root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE
  )
  if (file.exists(file.path(root, "DESCRIPTION")) &&
        file.exists(file.path(root, "NAMESPACE"))) {
    return(root)
  }
  NA_character_
}

.cs_skip_if_no_source <- function() {
  root <- .cs_root()
  if (is.na(root)) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }
  root
}

.cs_read <- function(root, path) {
  full <- file.path(root, path)
  if (!file.exists(full)) {
    testthat::skip(paste0(path, " not present in this tree"))
  }
  readLines(full, warn = FALSE)
}

# "Four" -> 4. The prose spells the count, so the gate has to read the word.
.cs_word_to_int <- function(w) {
  words <- c(
    one = 1L, two = 2L, three = 3L, four = 4L, five = 5L, six = 6L,
    seven = 7L, eight = 8L, nine = 9L, ten = 10L, eleven = 11L, twelve = 12L
  )
  key <- tolower(w)
  if (!key %in% names(words)) {
    return(NA_integer_)
  }
  words[[key]]
}

.cs_shipped_vignettes <- function(root) {
  # The .orig sources are the precompile inputs, not shipped vignettes.
  vs <- list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$")
  vs[!grepl("\\.orig$", vs)]
}


test_that("the vignette count in the package help matches the vignettes that ship", {
  root <- .cs_skip_if_no_source()
  n_actual <- length(.cs_shipped_vignettes(root))
  txt <- paste(.cs_read(root, "R/flexyBayes-package.R"), collapse = " ")
  m <- regmatches(
    txt,
    regexpr("([A-Za-z]+) vignettes ship with the package", txt)
  )
  expect_length(m, 1L)
  stated <- .cs_word_to_int(sub(" vignettes ship with the package", "", m))
  expect_equal(stated, n_actual)
})

test_that("the vignette count in cran-comments matches the vignettes that ship", {
  root <- .cs_skip_if_no_source()
  n_actual <- length(.cs_shipped_vignettes(root))
  txt <- paste(.cs_read(root, "cran-comments.md"), collapse = " ")
  m <- regmatches(txt, regexpr("The ([A-Za-z]+) vignettes ship", txt))
  expect_length(m, 1L)
  stated <- .cs_word_to_int(sub("^The ", "", sub(" vignettes ship", "", m)))
  expect_equal(stated, n_actual)
})

test_that("the landing surfaces do not offer a verb the user cannot call", {
  # The cross-engine comparison was unexported at 0.10.0.9000 on the R
  # Journal review's objection. It kept its help page, so a cross-reference
  # to it still resolved -- but the package-level help went on calling it
  # "the package's signature feature" and listing it among the entry points,
  # and a reader
  # following that finds a function they cannot call. The defect is a
  # landing surface offering an uncallable verb, not a dangling link.
  #
  # Scope is deliberate: the package landing page and the README are what a
  # new user reads first. Internal documentation may name internal functions.
  root <- .cs_skip_if_no_source()
  ns <- .cs_read(root, "NAMESPACE")
  exported <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))

  files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  defined <- unique(unlist(lapply(files, function(f) {
    lines <- readLines(f, warn = FALSE)
    sub("^([A-Za-z._][A-Za-z0-9._]*)\\s*<-\\s*function.*$", "\\1",
        grep("^[A-Za-z._][A-Za-z0-9._]*\\s*<-\\s*function", lines, value = TRUE))
  })))
  uncallable <- setdiff(defined, exported)
  uncallable <- uncallable[!startsWith(uncallable, ".")]

  offered <- character(0)
  for (path in c("R/flexyBayes-package.R", "README.md")) {
    if (!file.exists(file.path(root, path))) next
    lines <- .cs_read(root, path)
    prose <- if (endsWith(path, ".R")) grep("^#'", lines, value = TRUE) else lines
    text <- paste(prose, collapse = "\n")
    for (fn in uncallable) {
      if (grepl(paste0("[", fn, "()]"), text, fixed = TRUE) ||
            grepl(paste0("`", fn, "()`"), text, fixed = TRUE)) {
        offered <- c(offered, paste0(path, ": ", fn))
      }
    }
  }
  expect_equal(offered, character(0))
})

test_that("the startup banner does not advertise a withdrawn feature", {
  root <- .cs_skip_if_no_source()
  zzz <- paste(.cs_read(root, "R/zzz.R"), collapse = "\n")
  banner <- regmatches(
    zzz,
    regexpr("packageStartupMessage\\((?:[^()]|\\([^()]*\\))*\\)", zzz, perl = TRUE)
  )
  expect_length(banner, 1L)
  # The banner once closed by advertising the cross-engine comparison,
  # which 0.14.0 removes.
  expect_false(grepl("cross-engine", banner, fixed = TRUE))
})


test_that("every export appears in the pkgdown reference index", {
  root <- .cs_skip_if_no_source()
  ns <- .cs_read(root, "NAMESPACE")
  exported <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))

  # The reference index entries are YAML list items: "      - fb_prior".
  pkgdown_lines <- trimws(.cs_read(root, "_pkgdown.yml"))
  entries <- sub("^-\\s*", "", grep("^-\\s*[A-Za-z._]+$", pkgdown_lines, value = TRUE))

  missing <- Filter(
    function(x) {
      # An export is indexed either by its own name or, for a generic, by
      # the method entries that carry its flexyBayes implementations.
      !(x %in% entries) && !any(startsWith(entries, paste0(x, ".")))
    },
    exported
  )
  expect_equal(missing, character(0))
})


test_that("every export appears in API_STABILITY.md", {
  root <- .cs_skip_if_no_source()
  ns <- .cs_read(root, "NAMESPACE")
  exported <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
  api <- paste(.cs_read(root, "API_STABILITY.md"), collapse = "\n")
  missing <- Filter(function(x) !grepl(paste0("`", x, "("), api, fixed = TRUE), exported)
  expect_equal(missing, character(0))
})


test_that("the landing surfaces carry no withdrawn capability claim", {
  # The 0.12.0 package audit found claims the code does not support: a
  # "full MET" label on the small-sample baseline model, unstructured and
  # factor-analytic covariance said to fit on INLA "in seconds", a
  # disagreement under heavy missingness called a "property of the
  # estimand", pinning advice naming 0.10.x, and a known-issues sentence
  # saying INLA refuses the MET halves it fits.
  root <- .cs_skip_if_no_source()
  withdrawn <- c(
    "full MET", "full-MET", "in seconds", "property of the estimand",
    "0.10.x patch", "refuses both halves"
  )
  surfaces <- c(
    "README.md", "inst/KNOWN_ISSUES.md", "API_STABILITY.md",
    "R/na_action.R", "R/capability_matrix.R"
  )
  found <- character(0)
  for (path in surfaces) {
    text <- paste(.cs_read(root, path), collapse = "\n")
    hits <- withdrawn[vapply(
      withdrawn,
      function(p) grepl(p, text, fixed = TRUE),
      logical(1L)
    )]
    if (length(hits) > 0L) {
      found <- c(found, paste0(path, ": ", hits))
    }
  }
  expect_equal(found, character(0))
})

test_that("the 0.12.0 bound wording in NEWS carries its 0.12.1 correction", {
  root <- .cs_skip_if_no_source()
  news <- gsub("\\s+", " ", paste(.cs_read(root, "NEWS.md"), collapse = " "))

  expect_false(grepl(
    "standard deviation binds(?! \\(corrected in 0\\.12\\.1)",
    news,
    perl = TRUE
  ))
  # The 0.12.0 heading ends its bold span with a full stop, so the pattern
  # allows one; without it the check could not see the heading it guards.
  expect_false(grepl(
    "\\*\\*A tight uniform bound binds on INLA\\.?\\*\\*(?! \\(corrected in 0\\.12\\.1)",
    news,
    perl = TRUE
  ))
})
