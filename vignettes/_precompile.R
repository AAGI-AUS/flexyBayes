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

# The knit runs against the INSTALLED flexyBayes, not this tree. A bake
# after source edits without a reinstall ships stale output that reads
# as fresh -- it happened twice in the v0.9.1 documentation pass (the
# tree had moved, the installed build had not, and only a diff caught
# it). The build stamp of the installed package is printed here so the
# operator can refuse a stale bake; reinstall from a tarball built OFF
# the working tree (never from a cloud-synced source directory) before
# any bake that follows a source change.
.installed_stamp <- utils::packageDescription("flexyBayes")[["Packaged"]]
cat(sprintf(
  paste0(
    "[precompile] knitting against the INSTALLED flexyBayes ",
    "(Packaged: %s).\n[precompile] If R/ has changed since that ",
    "install, STOP and reinstall first -- the bake would be stale.\n"
  ),
  if (is.null(.installed_stamp)) "unknown" else .installed_stamp
))
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

# Chunk errors are failures, not output.
#
# `knitr::knit()` defaults to `error = TRUE`, which renders an error into
# the document and carries on. The v0.9.1 documentation pass reported
# "10 OK, 0 FAILED" on a batch in which one chunk had printed
# `Error in round(): non-numeric argument to mathematical function` into
# the page, and the only thing that caught it was reading the diff. An
# undeclared error now aborts that vignette's knit and reaches the
# FAILED path below. A chunk that declares `error = TRUE` in its own
# header still renders its error, which is what the refusal pages want.
knitr::opts_chunk$set(error = FALSE)

# .chunk_error_lines() --- error output the knit put in the page.
#
# The backstop to the option above, for the shapes it cannot reach: a
# chunk that declares `error = TRUE` and errors somewhere its author did
# not intend. Two things are not errors and are excluded:
#
#   * a line that appears verbatim in the .Rmd.orig -- prose passes
#     through the knit unchanged, and two shipped pages write a `#>
#     Error:` line by hand inside a plain fenced block to show a refusal
#     without paying for the fit that raises it;
#   * every hit in a source that declares `error = TRUE` on a chunk of
#     its own, because that page's subject is the refusals themselves.
#
# The exclusions are deliberately coarse. This is the second gate, and
# the first one is the option above.
.chunk_error_lines <- function(baked, source_file) {
  if (!file.exists(baked)) {
    return(character())
  }
  src <- readLines(source_file, warn = FALSE)
  if (any(grepl("error\\s*=\\s*TRUE", src))) {
    return(character())
  }
  out <- readLines(baked, warn = FALSE)
  hits <- out[grepl("^#>\\s*Error|Error in ", out)]
  setdiff(hits, src)
}

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
    next
  }

  bad <- .chunk_error_lines(out, f)
  if (length(bad) > 0L) {
    n_fail <- n_fail + 1L
    failures <- c(failures, f)
    cat(sprintf(
      "[precompile] FAILED on %s: %d chunk-error line(s) in %s\n",
      f,
      length(bad),
      out
    ))
    for (b in utils::head(bad, 5L)) {
      cat(sprintf("             %s\n", b))
    }
    next
  }

  n_ok <- n_ok + 1L
  # Fix the units: difftime picks its own (secs / mins / hours), so a
  # bare as.numeric() reported a 1.3-minute knit as "1.3 s".
  cat(sprintf(
    "[precompile] OK %s (%.1f s)\n",
    out,
    as.numeric(dt, units = "secs")
  ))
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
