# tools/check_build_parity.R -- the tarball must be reproducible from
# tracked source.
#
# Builds two source tarballs of the same commit -- one from the working
# tree, one from a clean `git archive` of HEAD -- and compares their file
# lists. Any difference fails.
#
# The 2026-08-16 adversarial review found the failure this gate exists to
# catch: `inst/extdata/inla-verification/` was gitignored as host-local
# verification material but was not build-excluded, so the locally built
# tarball carried a file a clean clone of the same commit did not, and
# the runtime gate read that file to decide what INLA would fit. A public
# reviewer building from a clone would not have built the same package.
#
# Run from the package root:
#
#   Rscript tools/check_build_parity.R
#
# Exit status 0 on parity, 1 on any difference or on a build failure.
# Uncommitted tracked changes are reported as a note, not a failure: the
# gate is about untracked and ignored files entering the build, and a
# working tree mid-edit is the normal state while it runs.

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.", call. = FALSE)
}

pkg_root <- normalizePath(".")
work_dir <- file.path(tempdir(), "fb-build-parity")
unlink(work_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ------------------------------------------------------- #

.fail <- function(...) {
  cat("[build-parity] FAIL: ", ..., "\n", sep = "")
  quit(status = 1L, save = "no")
}

# Build a source tarball of `src` into `dest` and return the file list
# the tarball carries, with the leading package directory stripped so
# the two builds are comparable.
.build_and_list <- function(src, dest, label) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "build", "--no-build-vignettes", "--no-manual", shQuote(src)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  tarballs <- list.files(dest, pattern = "\\.tar\\.gz$", full.names = TRUE)
  if (!is.null(status) && status != 0L || length(tarballs) != 1L) {
    cat(paste(out, collapse = "\n"), "\n")
    .fail("R CMD build failed for the ", label, " source.")
  }
  files <- utils::untar(tarballs[[1L]], list = TRUE)
  files <- sub("^[^/]+/", "", files)
  files <- files[nzchar(files)]
  sort(unique(files))
}

# ---- (0) every .Rbuildignore line must compile as a regex ---------- #
#
# `.Rbuildignore` has no comment syntax: `R CMD build` reads every
# non-empty line as a PCRE pattern. A line beginning with `#` is
# therefore a pattern too, and one containing an unbalanced parenthesis
# aborts the build with `invalid regular expression` rather than being
# skipped. Check the file before building, so the failure names the
# offending line instead of arriving as a PCRE compilation error.

.check_rbuildignore <- function(path = ".Rbuildignore") {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }
  lines <- readLines(path, warn = FALSE)
  keep <- which(nzchar(trimws(lines)))
  bad <- character(0L)
  for (i in keep) {
    ok <- tryCatch(
      {
        grepl(lines[[i]], "probe", perl = TRUE)
        TRUE
      },
      error = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (!ok) {
      bad <- c(bad, sprintf("  line %d: %s", i, lines[[i]]))
    }
  }
  if (length(bad)) {
    cat("[build-parity] .Rbuildignore lines that are not valid regexes:\n")
    cat(paste(bad, collapse = "\n"), "\n", sep = "")
    .fail(
      ".Rbuildignore has no comment syntax -- every line is a PCRE ",
      "pattern, and the lines above do not compile."
    )
  }
  cat(
    "[build-parity] .Rbuildignore: ", length(keep),
    " pattern(s), all compile.\n",
    sep = ""
  )
  invisible(NULL)
}

.check_rbuildignore()


# ---- (1) the working-tree build ------------------------------------ #

wt_dir <- file.path(work_dir, "worktree-build")
dir.create(wt_dir, recursive = TRUE, showWarnings = FALSE)
owd <- setwd(wt_dir)
wt_files <- .build_and_list(pkg_root, wt_dir, "working-tree")
setwd(owd)

# ---- (2) the clean git-archive build ------------------------------- #

export_dir <- file.path(work_dir, "head-export")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
archive <- file.path(work_dir, "head.tar")
ga <- suppressWarnings(system2(
  "git",
  c("archive", "--format=tar", "-o", shQuote(archive), "HEAD"),
  stdout = TRUE,
  stderr = TRUE
))
if (!is.null(attr(ga, "status")) && attr(ga, "status") != 0L) {
  cat(paste(ga, collapse = "\n"), "\n")
  .fail("git archive HEAD failed.")
}
utils::untar(archive, exdir = export_dir)

ga_dir <- file.path(work_dir, "archive-build")
dir.create(ga_dir, recursive = TRUE, showWarnings = FALSE)
owd <- setwd(ga_dir)
ga_files <- .build_and_list(export_dir, ga_dir, "git-archive")
setwd(owd)

# ---- (3) compare ---------------------------------------------------- #

head_sha <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
cat("[build-parity] HEAD ", head_sha, "\n", sep = "")
cat(
  "[build-parity] working tree: ", length(wt_files), " files; ",
  "git archive: ", length(ga_files), " files\n",
  sep = ""
)

dirty <- system2("git", c("status", "--porcelain", "--untracked-files=no"),
                 stdout = TRUE)
if (length(dirty)) {
  cat(
    "[build-parity] note: ", length(dirty), " tracked path(s) modified in ",
    "the working tree; the gate compares file LISTS, so this is not a ",
    "failure.\n",
    sep = ""
  )
}

only_wt <- setdiff(wt_files, ga_files)
only_ga <- setdiff(ga_files, wt_files)

if (length(only_wt)) {
  cat("[build-parity] only in the working-tree tarball:\n")
  cat(paste0("  ", only_wt, collapse = "\n"), "\n", sep = "")
}
if (length(only_ga)) {
  cat("[build-parity] only in the git-archive tarball:\n")
  cat(paste0("  ", only_ga, collapse = "\n"), "\n", sep = "")
}

if (length(only_wt) || length(only_ga)) {
  .fail(
    "the two tarballs carry different files. An untracked or ignored ",
    "file is entering the build; add it to .Rbuildignore, or track it."
  )
}

cat("[build-parity] PASS: both tarballs carry the same file list.\n")
quit(status = 0L, save = "no")
