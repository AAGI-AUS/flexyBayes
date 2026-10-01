# Release hygiene: the README install line, the citation authors, the place
# of cran-comments.md and the working folders kept out of the tarball. These
# read source-tree files that are not installed, so each test skips when the
# source tree is not reachable.

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
  bug <- unname(read.dcf(file.path(root, "DESCRIPTION"), fields = "BugReports")[1L, 1L])
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
    first <- paste(rmd[!grepl("install_github(", rmd, fixed = TRUE)],
                   collapse = " ")
    for (pkg in c("remotes", "knitr", "rmarkdown")) {
      expect_match(first, paste0("\"", pkg, "\""), fixed = TRUE, label = pkg)
    }
  }
})
