# release.R -- canonical release-ceremony gate for flexyBayes.
#
# Three sequential steps:
#   1. Full tally via tools/tally.R (workspace root).
#   2. Skip-ledger refresh (inst/skip-ledger.md).
#   3. R CMD build + R CMD check --as-cran --no-manual on a fresh tarball.
#
# Target wall-time: 5 -- 10 minutes on a laptop.  Run from the workspace
# root:
#   Rscript flexyBayes/tools/profiles/release.R

cat("=== release profile ===\n")

# ---- Step 1/3: full tally --------------------------------------------

cat("\n-- Step 1/3: full tally (Rscript tools/tally.R) --\n")
Sys.setenv(NOT_CRAN = "true")
source("tools/tally.R", chdir = FALSE)

# ---- Step 2/3: skip-ledger refresh -----------------------------------

cat("\n-- Step 2/3: refresh inst/skip-ledger.md --\n")
source("flexyBayes/tools/skip_ledger.R")
out <- build_skip_ledger(
  test_dir = "flexyBayes/tests/testthat",
  output_path = "flexyBayes/inst/skip-ledger.md"
)
cat(sprintf(
  "Ledger: %d sites across %d files.\n",
  nrow(out),
  length(unique(out$file))
))

# ---- Step 3/3: R CMD build + check --as-cran --no-manual -------------

cat("\n-- Step 3/3: R CMD build + R CMD check --as-cran --no-manual --\n")
if (
  !requireNamespace("pkgbuild", quietly = TRUE) ||
    !requireNamespace("rcmdcheck", quietly = TRUE)
) {
  stop(
    call. = FALSE,
    "release profile requires `pkgbuild` and `rcmdcheck`; install both"
  )
}
tarball <- pkgbuild::build(
  path = "flexyBayes",
  dest_path = tempdir(),
  args = "--no-build-vignettes",
  quiet = TRUE
)
chk <- rcmdcheck::rcmdcheck(
  path = tarball,
  args = c("--as-cran", "--no-manual"),
  error_on = "never",
  quiet = FALSE
)
cat(sprintf(
  "\nR CMD check: %d errors, %d warnings, %d notes.\n",
  length(chk$errors),
  length(chk$warnings),
  length(chk$notes)
))

if (length(chk$errors) > 0L || length(chk$warnings) > 0L) {
  cat("\nRelease profile FAILED on R CMD check.\n")
  quit(save = "no", status = 1L)
}
cat("\nRelease profile PASSED.\n")
