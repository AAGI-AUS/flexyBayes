# vignettes/_precompile.R
#
# Precompile every `flexyBayes-*.Rmd.orig` into its sibling
# `flexyBayes-*.Rmd` so that `R CMD build` and `R CMD check` render
# the static (already-evaluated) `.Rmd` in seconds rather than
# re-running the live MCMC chunks every time.
#
# Run interactively whenever a `.Rmd.orig` source is touched:
#
#   Rscript vignettes/_precompile.R                # all .Rmd.orig
#   Rscript vignettes/_precompile.R --only 05,08   # selective
#
# The `--only` flag accepts a comma-separated list of vignette
# basenames or basename suffixes. `--only 08` matches
# `flexyBayes-08-downstream-analysis.Rmd.orig`; `--only 05,08`
# precompiles both. Added 2026-05-26 (v0.3.8 Phase C) for the
# selective-refresh workflow ratified by v038-plan-2026-05-25
# §12.3 -- avoids paying the full 12-vignette MCMC cost when only
# a subset needs refreshing.
#
# Pre-requisites (the precompile is a one-time payment that must
# succeed in your interactive session, not in CI):
#
# - flexyBayes installed at the version matching DESCRIPTION
#   (`devtools::install()` or `R CMD INSTALL` of a recent tarball).
# - greta installed and `greta::install_greta_deps()` completed.
# - INLA installed (Additional_repositories binary path).
#
# Outputs:
#
# - `vignettes/<name>.Rmd` — static, ships in the package tarball.
# - `vignettes/<name>-figs/*.png` — figures, ship alongside the
#   static `.Rmd`. Per-vignette directory keeps figure namespaces
#   separate.
#
# `.Rbuildignore` excludes:
#
# - `^vignettes/_precompile\.R$` — this driver.
# - `^vignettes/.*\.Rmd\.orig$` — the live-MCMC sources.
#
# Audit recipe 09 (`/rpkg`).

vignettes_dir <- if (basename(getwd()) == "vignettes") {
  getwd()
} else if (dir.exists("vignettes")) {
  normalizePath("vignettes")
} else {
  stop(
    "Run from the package root or from the vignettes/ directory.",
    call. = FALSE
  )
}

orig_files <- list.files(
  vignettes_dir,
  pattern = "\\.Rmd\\.orig$",
  full.names = TRUE
)

if (length(orig_files) == 0L) {
  stop(
    "No .Rmd.orig files found in ",
    vignettes_dir,
    ". Did you cp the .Rmd files first?",
    call. = FALSE
  )
}

# Selective filter: --only <comma-list> restricts the precompile
# set to vignettes whose basename (without .Rmd.orig) contains any
# of the comma-separated tokens. Substring match -- "08" matches
# "flexyBayes-08-downstream-analysis"; "downstream" matches the
# same vignette. Tokens that match no vignette raise; the empty
# filter (no --only flag) keeps the full set.
.cli_args <- commandArgs(trailingOnly = TRUE)
.only_idx <- which(.cli_args == "--only")
if (length(.only_idx) == 1L && .only_idx < length(.cli_args)) {
  tokens <- strsplit(.cli_args[.only_idx + 1L], ",", fixed = TRUE)[[1L]]
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0L) {
    stop("--only requires at least one comma-separated token.", call. = FALSE)
  }
  bases <- sub("\\.Rmd\\.orig$", "", basename(orig_files))
  keep <- vapply(
    bases,
    function(b) {
      any(vapply(
        tokens,
        function(tok) grepl(tok, b, fixed = TRUE),
        logical(1L)
      ))
    },
    logical(1L)
  )
  unmatched <- tokens[
    !vapply(
      tokens,
      function(tok) {
        any(grepl(tok, bases, fixed = TRUE))
      },
      logical(1L)
    )
  ]
  if (length(unmatched) > 0L) {
    stop(
      "--only token(s) matched no vignette: ",
      paste(unmatched, collapse = ", "),
      ". Available basenames: ",
      paste(bases, collapse = ", "),
      call. = FALSE
    )
  }
  orig_files <- orig_files[keep]
  cat(sprintf(
    "[precompile] --only filter: keeping %d of %d vignette(s)\n",
    length(orig_files),
    length(bases)
  ))
}

cat(sprintf(
  "[precompile] %d vignette source(s) found in %s\n",
  length(orig_files),
  vignettes_dir
))

# Knit each .Rmd.orig in its own working directory so figure paths
# stay relative to vignettes/.
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(vignettes_dir)

n_ok <- 0L
n_fail <- 0L
failures <- character()

for (f in basename(orig_files)) {
  out <- sub("\\.orig$", "", f)
  base <- sub("\\.Rmd\\.orig$", "", f)

  cat(sprintf("[precompile] %s -> %s\n", f, out))

  # Per-vignette fig.path keeps figure files namespaced.
  knitr::opts_chunk$set(fig.path = paste0(base, "-figs/"))

  t0 <- Sys.time()
  result <- tryCatch(
    knitr::knit(input = f, output = out, envir = new.env(), quiet = FALSE),
    error = function(e) {
      cat(sprintf(
        "[precompile] FAILED on %s: %s\n",
        f,
        conditionMessage(e)
      ))
      structure("error", condition_message = conditionMessage(e))
    }
  )
  dt <- Sys.time() - t0

  if (identical(unclass(result), "error")) {
    n_fail <- n_fail + 1L
    failures <- c(failures, f)
  } else {
    n_ok <- n_ok + 1L
    cat(sprintf("[precompile] OK %s (%.1f s)\n", out, as.numeric(dt)))
  }
}

cat(sprintf(
  "\n[precompile] done: %d OK, %d FAILED\n",
  n_ok,
  n_fail
))

if (n_fail > 0L) {
  cat("[precompile] failures:\n")
  for (fl in failures) {
    cat(sprintf("  - %s\n", fl))
  }
  quit(status = 1L)
}
