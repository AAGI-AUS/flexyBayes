# tools/generate_translation_registry.R: regenerate the translation registry.
#
# One R-level table (.fb_translation_registry(), R/translation_registry.R) is
# the source of truth for what each ASReml structure is lowered to, which
# grammar vintage the spelling belongs to, what the lowering promises, and
# what independently grounds it. This script renders it to the Markdown block
# in inst/TRANSLATION.md and to the machine-readable CSV that ships at
# inst/translation/translation_registry.csv.
#
# Run from the package root after editing the R table:
#
#   Rscript tools/generate_translation_registry.R
#
# tests/testthat/test-translation-registry.R fails when either committed
# artefact differs from a fresh render, so the generator must be re-run
# whenever the table changes.

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.", call. = FALSE)
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop(
    "tools/generate_translation_registry.R needs pkgload (a devtools ",
    "dependency) to load the source tree.",
    call. = FALSE
  )
}
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)

ns <- asNamespace("flexyBayes")
registry <- get(".fb_translation_registry", envir = ns)
splice <- get(".fb_translation_splice", envir = ns)

target <- "inst/TRANSLATION.md"
before <- readLines(target, warn = FALSE)
splice(target, write = TRUE)
after <- readLines(target, warn = FALSE)
message(
  if (identical(before, after)) "unchanged: " else "rewritten: ",
  target
)

csv_path <- "inst/translation/translation_registry.csv"
dir.create(dirname(csv_path), showWarnings = FALSE, recursive = TRUE)
reg <- registry()
before_csv <- if (file.exists(csv_path)) {
  readLines(csv_path, warn = FALSE)
} else {
  character(0)
}
utils::write.csv(reg, csv_path, row.names = FALSE)
after_csv <- readLines(csv_path, warn = FALSE)
message(
  if (identical(before_csv, after_csv)) "unchanged: " else "rewritten: ",
  csv_path
)
