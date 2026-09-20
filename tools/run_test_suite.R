# tools/run_test_suite.R - run the whole test suite as the release does, and
# record every block with the reason it was skipped, if it was.
#
# Usage, from the package root:
#   FLEXYBAYES_RUN_LONG=true Rscript tools/run_test_suite.R suite_record.tsv
#
# NOT_CRAN is set, so skip_on_cran() blocks run, and INLA runs on one thread.
# FLEXYBAYES_RUN_LONG is read from the environment and never set here: a
# release run that leaves it unset records every long-tier block as skipped,
# and criterion D11 of tools/definition_of_done.R fails on that record
#   FLEXYBAYES_SUITE_RECORD=suite_record.tsv \
#     Rscript tools/definition_of_done.R --full --section D11
# Without the variable the same command times the default tier.
#
# Not shipped (`^tools$` in .Rbuildignore).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript tools/run_test_suite.R <record.tsv>", call. = FALSE)
}
if (!file.exists("DESCRIPTION")) {
  stop("run from the package root (no DESCRIPTION here)", call. = FALSE)
}
record_path <- args[[1L]]

Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
if (requireNamespace("INLA", quietly = TRUE)) {
  INLA::inla.setOption(num.threads = "1:1")
}
testthat::set_max_fails(Inf)

t0 <- Sys.time()
res <- testthat::test_dir(
  "tests/testthat",
  reporter = testthat::SummaryReporter$new(
    show_praise = FALSE,
    max_reports = Inf
  ),
  stop_on_failure = FALSE,
  package = "flexyBayes",
  load_package = "none"
)
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# One row per block. The skip reason is the message of the block's first skip
# condition; tabs and newlines are flattened so the record stays one line per
# block.
tab <- as.data.frame(res)
tab$skip_reason <- vapply(res, function(r) {
  sk <- Filter(function(e) inherits(e, "expectation_skip"), r$results)
  if (length(sk)) conditionMessage(sk[[1L]]) else ""
}, character(1))
tab$test <- gsub("[\t\r\n]", " ", tab$test)
tab$skip_reason <- gsub("[\t\r\n]", " ", tab$skip_reason)
cols <- c("file", "test", "nb", "failed", "skipped", "error", "warning",
          "passed", "real", "skip_reason")
utils::write.table(tab[, cols], record_path, sep = "\t", row.names = FALSE,
                   quote = FALSE)

long_held <- grepl("FLEXYBAYES_RUN_LONG", tab$skip_reason, fixed = TRUE)
cat(sprintf(
  paste0("SUITE files=%d blocks=%d blocks_run=%d expectations=%d ",
         "failed=%d errored_blocks=%d warned=%d skipped_blocks=%d ",
         "long_tier_skipped=%d wall_s=%.0f FLEXYBAYES_RUN_LONG=%s\n"),
  length(unique(tab$file)), nrow(tab), sum(!tab$skipped), sum(tab$nb),
  sum(tab$failed), sum(tab$error), sum(tab$warning), sum(tab$skipped),
  sum(long_held), wall, Sys.getenv("FLEXYBAYES_RUN_LONG", "(unset)")
))
cat("record: ", record_path, "\n", sep = "")
