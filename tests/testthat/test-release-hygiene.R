# These tests read source files that are not installed, and skip without them.

.rh_root <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "DESCRIPTION"))) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }
  root
}

.rh_read <- function(root, path) {
  full <- file.path(root, path)
  if (!file.exists(full)) {
    testthat::skip(paste0(path, " not present in this tree"))
  }
  readLines(full, warn = FALSE, encoding = "UTF-8")
}

# "owner/repo" from the GitHub URL in DESCRIPTION's BugReports field.
.rh_github_repo <- function(root) {
  bug <- unname(read.dcf(file.path(root, "DESCRIPTION"), fields = "BugReports")[
    1L,
    1L
  ])
  sub("^https://github\\.com/([^/]+/[^/]+)/issues/?$", "\\1", bug)
}

# The lines of the fenced code block that holds `install_github(`.
.rh_install_block <- function(lines) {
  hit <- grep("install_github(", lines, fixed = TRUE)
  if (length(hit) != 1L) {
    return(NULL)
  }
  fences <- grep("^```", lines)
  open <- max(fences[fences < hit])
  close <- min(fences[fences > hit])
  lines[seq.int(open + 1L, close - 1L)]
}

test_that("the README install line names the package's own repository", {
  root <- .rh_root()
  repo <- .rh_github_repo(root)
  expect_identical(repo, "AAGI-AUS/flexyBayes")

  for (readme in c("README.Rmd", "README.md")) {
    block <- .rh_install_block(.rh_read(root, readme))
    expect_false(is.null(block), label = paste(readme, "has one install line"))
    call <- grep("install_github(", block, fixed = TRUE, value = TRUE)
    expect_match(
      call,
      paste0("install_github(\"", repo, "\""),
      fixed = TRUE,
      label = readme
    )
  }
})

test_that("the README install block installs what building the vignettes needs", {
  root <- .rh_root()
  rmd <- .rh_install_block(.rh_read(root, "README.Rmd"))
  md <- .rh_install_block(.rh_read(root, "README.md"))
  expect_identical(md, rmd)

  # remotes installs Imports only, and R CMD build stops without the
  # vignette builder named in DESCRIPTION.
  if (any(grepl("build_vignettes = TRUE", rmd, fixed = TRUE))) {
    first <- paste(
      rmd[!grepl("install_github(", rmd, fixed = TRUE)],
      collapse = " "
    )
    for (pkg in c("remotes", "knitr", "rmarkdown")) {
      expect_match(first, paste0("\"", pkg, "\""), fixed = TRUE, label = pkg)
    }
  }
})

# Family name and ORCID of each person, as "Family|0000-0000-0000-0000".
.rh_key <- function(family, orcid) {
  paste(family, sub("^https://orcid\\.org/", "", orcid), sep = "|")
}

.rh_description_people <- function(root, role) {
  field <- read.dcf(file.path(root, "DESCRIPTION"), fields = "Authors@R")[
    1L,
    1L
  ]
  people <- eval(parse(text = field))
  people <- people[vapply(
    people,
    function(p) role %in% p$role && length(p$family) == 1L,
    logical(1L)
  )]
  vapply(
    people,
    function(p) .rh_key(p$family, p$comment[["ORCID"]]),
    character(1L)
  )
}

# The values of `field` inside the JSON array `key` of codemeta.json.
.rh_codemeta_values <- function(text, key, field) {
  start <- regexpr(paste0("\"", key, "\"\\s*:\\s*\\["), text)
  if (start < 0L) {
    return(character(0))
  }
  chars <- strsplit(substring(text, start), "")[[1L]]
  depth <- cumsum((chars == "[") - (chars == "]"))
  block <- paste(
    chars[seq_len(which(depth == 0L & chars == "]")[1L])],
    collapse = ""
  )
  hits <- regmatches(
    block,
    gregexpr(paste0("\"", field, "\"\\s*:\\s*\"[^\"]*\""), block)
  )[[1L]]
  sub("^.*:\\s*\"([^\"]*)\"$", "\\1", hits)
}

test_that("every citation surface names the DESCRIPTION authors", {
  root <- .rh_root()
  authors <- .rh_description_people(root, "aut")
  contributors <- .rh_description_people(root, "ctb")
  expect_gt(length(authors), 0L)

  # inst/CITATION
  path <- file.path(root, "inst", "CITATION")
  skip_if_not(file.exists(path), "inst/CITATION not present")
  bib <- utils::readCitationFile(path, meta = list(Encoding = "UTF-8"))
  cited <- vapply(
    unclass(bib)[[1L]]$author,
    function(p) .rh_key(p$family, p$comment[["ORCID"]]),
    character(1L)
  )
  expect_identical(cited, authors, label = "inst/CITATION authors")

  # CITATION.cff: the authors block runs to the next top-level key.
  cff <- .rh_read(root, "CITATION.cff")
  from <- grep("^authors:", cff)
  to <- grep("^[a-z-]+:", cff)
  block <- cff[seq.int(from + 1L, min(to[to > from]) - 1L)]
  family <- sub(
    "^.*family-names:\\s*\"(.*)\"\\s*$",
    "\\1",
    grep("family-names:", block, value = TRUE)
  )
  orcid <- sub(
    "^.*orcid:\\s*\"(.*)\"\\s*$",
    "\\1",
    grep("orcid:", block, value = TRUE)
  )
  expect_identical(
    .rh_key(family, orcid),
    authors,
    label = "CITATION.cff authors"
  )

  # codemeta.json: authors and contributors, in DESCRIPTION's order.
  meta <- paste(.rh_read(root, "codemeta.json"), collapse = "\n")
  for (key in c("author", "contributor")) {
    expected <- if (key == "author") authors else contributors
    got <- .rh_key(
      .rh_codemeta_values(meta, key, "familyName"),
      .rh_codemeta_values(meta, key, "@id")
    )
    expect_identical(got, unname(expected), label = paste("codemeta.json", key))
  }
})

# Does any .Rbuildignore pattern exclude this top-level path from the build?
.rh_build_ignored <- function(root, path) {
  patterns <- readLines(file.path(root, ".Rbuildignore"), warn = FALSE)
  patterns <- patterns[nzchar(patterns) & !startsWith(patterns, "#")]
  any(vapply(
    patterns,
    function(p) grepl(p, path, perl = TRUE, ignore.case = TRUE),
    logical(1L)
  ))
}

test_that("cran-comments.md sits at the root, outside the build", {
  root <- .rh_root()
  # devtools::submit_cran() reads it from the package root.
  lines <- .rh_read(root, "cran-comments.md")
  expect_true(.rh_build_ignored(root, "cran-comments.md"))
  expect_false(file.exists(file.path(root, "inst", "cran-comments.md")))

  # A CRAN reviewer cannot open the maintainer's working folders.
  private <- grep("WS/|review/|/Users/", lines, value = TRUE)
  expect_identical(private, character(0))
})

test_that("every working folder at the root is kept out of the build", {
  root <- .rh_root()
  package <- c(
    "DESCRIPTION",
    "NAMESPACE",
    "LICENSE",
    "NEWS.md",
    "README.md",
    "R",
    "man",
    "tests",
    "vignettes",
    "inst",
    "data",
    "src",
    "tools",
    ".Rbuildignore",
    ".gitignore"
  )
  # Named even when absent, since they exist only in some checkouts.
  working <- c(
    "review",
    "oracle_work",
    "validation",
    "vignettes-reserve",
    "data-raw",
    "ledger",
    "paper",
    "dev-narrative",
    "docs",
    "flexyBayes.Rcheck"
  )
  present <- list.files(root, all.files = TRUE, no.. = TRUE)
  candidates <- union(setdiff(present, package), working)
  shipped <- candidates[
    !vapply(
      candidates,
      function(path) .rh_build_ignored(root, path),
      logical(1L)
    )
  ]
  expect_identical(shipped, character(0))
})
