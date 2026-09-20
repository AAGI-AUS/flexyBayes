# tools/generate_capability_matrix.R -- regenerate the capability block.
#
# One R-level table (.fb_capability_matrix(), R/capability_matrix.R) is the
# source of truth for what each active engine fits, emits, or refuses. This
# script renders it to Markdown and splices it between the
# <!-- capability-matrix:begin --> / <!-- capability-matrix:end --> markers
# in every surface that shows the table (README.md with the rows at least one
# engine fits, inst/KNOWN_ISSUES.md with every row and its note; the layouts
# come from .fb_capability_targets()), and writes the same table as
# inst/validation/execution_grid/capability_matrix.csv -- the row builder in
# tools/capability_grid_matrix.R is shared with tools/execution_grid.R, so
# this CSV and the grid's own copy cannot drift apart for the same registry.
#
# Run from the package root after editing the R table:
#
#   Rscript tools/generate_capability_matrix.R
#
# tests/testthat/test-capability-matrix.R fails when a committed Markdown
# block differs from a fresh render, and
# tests/testthat/test-installed-evidence-parity.R fails when the installed
# CSV differs from the live registry, so the generator must be re-run
# whenever the table changes. It also re-derives every verdict from the
# gate and emit code, so a table edit that the code does not support fails
# there.

csv_target <- "inst/validation/execution_grid/capability_matrix.csv"

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
md_targets <- get(".fb_capability_targets", envir = asNamespace("flexyBayes"))()
capability_matrix <- get(
  ".fb_capability_matrix",
  envir = asNamespace("flexyBayes")
)

# The CSV row builder is shared with tools/execution_grid.R (see that
# file's M section) so the two writers cannot drift apart for the same
# registry. The drift check fails loudly, before either surface is
# touched, if a registry row has no M-section call and no recorded
# elsewhere-execution -- the same guard the grid itself runs before
# driving a single cell.
source("tools/capability_grid_matrix.R")

cap <- capability_matrix()
m_calls <- .fb_capability_grid_m_calls()
m_elsewhere <- .fb_capability_grid_m_executed_elsewhere()
.fb_capability_grid_check_drift(cap, m_calls, m_elsewhere)
cells <- .fb_capability_grid_cells(cap, m_calls)

for (target in names(md_targets)) {
  before <- readLines(target, warn = FALSE)
  splice(target, write = TRUE, layout = md_targets[[target]])
  after <- readLines(target, warn = FALSE)
  status <- if (identical(before, after)) "unchanged" else "updated"
  cat(sprintf("[capability-matrix] %s: %s\n", target, status))
}

before_csv <- if (file.exists(csv_target)) {
  readLines(csv_target, warn = FALSE)
} else {
  NA_character_
}
utils::write.csv(cells, csv_target, row.names = FALSE)
after_csv <- readLines(csv_target, warn = FALSE)
csv_status <- if (identical(before_csv, after_csv)) "unchanged" else "updated"
cat(sprintf("[capability-matrix] %s: %s\n", csv_target, csv_status))

cat("[capability-matrix] done\n")
