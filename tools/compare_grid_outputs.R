# tools/compare_grid_outputs.R -- normalised grid-output comparator.
#
# The release gate confirms a candidate-tarball grid run reproduces a
# reference run (recipe 62 / U-EVID-1's acceptance: "the final tarball's
# grid must reproduce them under a declared normalisation"). This script
# compares ledger.csv and misses.csv between two
# inst/validation/execution_grid/-shaped directories on that declared
# normalisation, which drops run-specific noise (commit sha, date,
# wall_sec, engine/package versions, free-text detail/message columns) so
# a re-run on a different day, against a different INLA point release,
# compares clean when the outcomes themselves agree:
#
#   ledger.csv  -- cell, verdict, class. `class` is the typed outcome
#                  (fit, refuse_typed, error_untyped, crash, timeout,
#                  silent_drop, translated, construct_accept, value);
#                  `verdict` is as_claimed / DIVERGENT.
#   misses.csv  -- cell, observed, class. `class` is the typed miss code
#                  (M1-M6); misses.csv carries no `verdict` column, so
#                  `observed` stands in for it -- every row already
#                  records a divergence from `expected`.
#
# Usage:
#
#   Rscript tools/compare_grid_outputs.R <reference-dir> <candidate-dir>
#
# Prints every row present, at its normalised key, in one directory's file
# and not the other's (both directions), and exits 1 if either file
# carries such a difference. A file absent from both directories is
# skipped, not a difference; a file present in one and absent from the
# other is. Exits 0 only when every declared file compares clean.

args_vec <- commandArgs(trailingOnly = TRUE)
if (length(args_vec) != 2L) {
  stop(
    "usage: Rscript tools/compare_grid_outputs.R <reference-dir> ",
    "<candidate-dir>",
    call. = FALSE
  )
}
ref_dir <- args_vec[[1L]]
cand_dir <- args_vec[[2L]]
for (d1 in c(ref_dir, cand_dir)) {
  if (!dir.exists(d1)) {
    stop("no such directory: ", d1, call. = FALSE)
  }
}

# --- the declared normalisation ----------------------------------------- #

.NORM_COLS <- list(
  "ledger.csv" = c("cell", "verdict", "class"),
  "misses.csv" = c("cell", "observed", "class")
)

.read_normalised <- function(dir, file, cols) {
  path <- file.path(dir, file)
  if (!file.exists(path)) {
    return(NULL)
  }
  tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  absent <- setdiff(cols, names(tab))
  if (length(absent) > 0L) {
    stop(
      file, " at ", dir, " is missing declared column(s): ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  tab[, cols, drop = FALSE]
}

.row_key <- function(tab) {
  do.call(paste, c(as.list(tab), sep = ""))
}

# --- one file's comparison ------------------------------------------------ #

.compare_file <- function(file) {
  cols <- .NORM_COLS[[file]]
  ref <- .read_normalised(ref_dir, file, cols)
  cand <- .read_normalised(cand_dir, file, cols)

  if (is.null(ref) && is.null(cand)) {
    cat(sprintf(
      "[compare-grid] %s: absent from both directories, skipped\n", file
    ))
    return(TRUE)
  }
  if (is.null(ref) || is.null(cand)) {
    cat(sprintf(
      "[compare-grid] %s: present in %s, absent from %s\n",
      file,
      if (is.null(ref)) cand_dir else ref_dir,
      if (is.null(ref)) ref_dir else cand_dir
    ))
    return(FALSE)
  }

  ref_key <- .row_key(ref)
  cand_key <- .row_key(cand)
  only_ref <- ref[!ref_key %in% cand_key, , drop = FALSE]
  only_cand <- cand[!cand_key %in% ref_key, , drop = FALSE]

  if (nrow(only_ref) == 0L && nrow(only_cand) == 0L) {
    cat(sprintf("[compare-grid] %s: %d row(s), clean\n", file, nrow(ref)))
    return(TRUE)
  }

  cat(sprintf(
    "[compare-grid] %s: %d row(s) only in reference (%s), %d row(s) only in candidate (%s)\n",
    file, nrow(only_ref), ref_dir, nrow(only_cand), cand_dir
  ))
  if (nrow(only_ref) > 0L) {
    cat("  -- reference only --\n")
    for (i1 in seq_len(nrow(only_ref))) {
      cat("  ", paste(only_ref[i1, ], collapse = " | "), "\n")
    }
  }
  if (nrow(only_cand) > 0L) {
    cat("  -- candidate only --\n")
    for (i1 in seq_len(nrow(only_cand))) {
      cat("  ", paste(only_cand[i1, ], collapse = " | "), "\n")
    }
  }
  FALSE
}

.clean_vec <- vapply(names(.NORM_COLS), .compare_file, logical(1L))

if (all(.clean_vec)) {
  cat("[compare-grid] clean: reference and candidate agree.\n")
  quit(save = "no", status = 0L)
}
cat("[compare-grid] NOT CLEAN: see the differences above.\n")
quit(save = "no", status = 1L)
