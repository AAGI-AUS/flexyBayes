# =============================================================================
# Version coherence across the metadata surfaces.
#
# The package already had a coherence check, and it reported PASS while
# cran-comments.md, SECURITY.md, API_STABILITY.md, .zenodo.json and
# SUPPORT.md all still described the 0.8.x line. The check was looking at
# CITATION.cff, codemeta.json and inst/CITATION only, so five surfaces a
# reviewer reads before any source file went stale unobserved.
#
# This gate closes that. Two assertions per surface:
#
#   (a) the DESCRIPTION version string appears in the file;
#   (b) no superseded version is described as the *current* one -- a line
#       naming an older version may stand as history, but not alongside
#       "current", "development line", or "supported".
#
# The files are source-tree artefacts (SUPPORT.md is .Rbuildignore'd, and
# neither README.md nor cran-comments.md is installed), so each assertion
# skips when its file is absent rather than failing on an installed tree.
# =============================================================================

.vc_root <- function() {
  candidate <- testthat::test_path("..", "..")
  if (file.exists(file.path(candidate, "DESCRIPTION"))) {
    return(normalizePath(candidate))
  }
  NA_character_
}

.vc_version <- function() {
  root <- .vc_root()
  if (is.na(root)) {
    return(NA_character_)
  }
  as.character(read.dcf(file.path(root, "DESCRIPTION"), fields = "Version"))
}

# Surfaces that must name the current version. Kept as one vector so a new
# metadata file is added in exactly one place.
.VC_SURFACES <- c(
  "cran-comments.md",
  "SECURITY.md",
  "API_STABILITY.md",
  "SUPPORT.md",
  ".zenodo.json"
)

# A superseded version described as live. The marker set is deliberately
# narrow: a changelog line such as "0.8.3 -- added glance()" is history and
# must keep working.
.VC_CURRENCY_MARKERS <- paste(
  "current", "development line", "development release", "supported version",
  sep = "|"
)

test_that("the DESCRIPTION version is a plain three-part string", {
  version <- .vc_version()
  skip_if(is.na(version), "not running from the source tree")
  expect_match(version, "^[0-9]+\\.[0-9]+\\.[0-9]+$")
})

test_that("every metadata surface names the DESCRIPTION version", {
  root <- .vc_root()
  skip_if(is.na(root), "not running from the source tree")
  version <- .vc_version()

  for (surface in .VC_SURFACES) {
    path <- file.path(root, surface)
    skip_if_not(file.exists(path), paste(surface, "is not in this tree"))
    lines <- readLines(path, warn = FALSE)
    expect_true(
      any(grepl(version, lines, fixed = TRUE)),
      label = paste0(surface, " names version ", version)
    )
  }
})

# Version strings on a line, compared numerically against the current one.
# A hand-written regex for "any line older than this one" was wrong here:
# it read "0\\.[0-8]\\." and so could not see the 0.9.x line at all once
# DESCRIPTION reached 0.10.0 -- the blind spot covered exactly the
# transition being made. Parsing each version and comparing has no such
# gap, and needs no edit at the next bump.
.vc_stale_on <- function(line, current) {
  found <- regmatches(line, gregexpr("[0-9]+\\.[0-9]+\\.[0-9x]+", line))[[1]]
  if (!length(found)) {
    return(FALSE)
  }
  any(vapply(
    found,
    function(v) {
      isTRUE(tryCatch(
        package_version(sub("x$", "0", v)) < current,
        error = function(e) FALSE
      ))
    },
    logical(1)
  ))
}

test_that("no superseded version is described as the current one", {
  root <- .vc_root()
  skip_if(is.na(root), "not running from the source tree")
  version <- .vc_version()
  current <- package_version(version)

  for (surface in .VC_SURFACES) {
    path <- file.path(root, surface)
    skip_if_not(file.exists(path), paste(surface, "is not in this tree"))
    lines <- readLines(path, warn = FALSE)
    marked <- lines[grepl(.VC_CURRENCY_MARKERS, lines, ignore.case = TRUE)]
    offenders <- marked[vapply(marked, .vc_stale_on, logical(1), current = current)]
    expect_identical(
      offenders,
      character(0),
      label = paste0(
        surface, " describes a superseded version as current (",
        paste(trimws(offenders), collapse = " / "), ")"
      )
    )
  }
})

test_that(".zenodo.json's version field equals the DESCRIPTION version", {
  root <- .vc_root()
  skip_if(is.na(root), "not running from the source tree")
  path <- file.path(root, ".zenodo.json")
  skip_if_not(file.exists(path), ".zenodo.json is not in this tree")

  version <- .vc_version()
  lines <- readLines(path, warn = FALSE)
  field <- grep("^[[:space:]]*\"version\"[[:space:]]*:", lines, value = TRUE)
  expect_length(field, 1L)
  expect_match(field, paste0("\"", version, "\""), fixed = TRUE)
})

test_that("the citation surfaces still agree with DESCRIPTION", {
  root <- .vc_root()
  skip_if(is.na(root), "not running from the source tree")
  version <- .vc_version()

  for (surface in c("CITATION.cff", "codemeta.json", "inst/CITATION")) {
    path <- file.path(root, surface)
    skip_if_not(file.exists(path), paste(surface, "is not in this tree"))
    lines <- readLines(path, warn = FALSE)
    expect_true(
      any(grepl(version, lines, fixed = TRUE)),
      label = paste0(surface, " names version ", version)
    )
  }
})
