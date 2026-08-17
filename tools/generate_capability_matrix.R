# tools/generate_capability_matrix.R -- regenerate the capability block.
#
# One R-level table (.fb_capability_matrix(), R/capability_matrix.R) is the
# source of truth for what each active engine fits, emits, or refuses. This
# script renders it to Markdown and splices it between the
# <!-- capability-matrix:begin --> / <!-- capability-matrix:end --> markers
# in every surface that shows the table.
#
# Run from the package root after editing the R table:
#
#   Rscript tools/generate_capability_matrix.R
#
# tests/testthat/test-capability-matrix.R fails when a committed block
# differs from a fresh render, so the generator must be re-run whenever the
# table changes. It also re-derives every verdict from the gate and emit
# code, so a table edit that the code does not support fails there.

targets <- c("README.md", "inst/KNOWN_ISSUES.md")

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.", call. = FALSE)
}

# load_all() rather than library(): the generator reads internal
# functions, and the source tree is the version being documented.
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop(
    "tools/generate_capability_matrix.R needs pkgload (a devtools ",
    "dependency) to load the source tree.",
    call. = FALSE
  )
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

splice <- get(".fb_capability_splice", envir = asNamespace("flexyBayes"))

for (target in targets) {
  before <- readLines(target, warn = FALSE)
  splice(target, write = TRUE)
  after <- readLines(target, warn = FALSE)
  status <- if (identical(before, after)) "unchanged" else "updated"
  cat(sprintf("[capability-matrix] %s: %s\n", target, status))
}

cat("[capability-matrix] done\n")
