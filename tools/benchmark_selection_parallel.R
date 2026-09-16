#!/usr/bin/env Rscript

# benchmark_selection_parallel.R -- what concurrency buys an fb_select()
# call, and what it changes.
#
# The fixture, the shapes and the memory instrument are P1-FIX's
# (tools/fixture_selection_cost.R, review/units/P1-FIX_2026-09-16): the
# same package-owned synthetic data at the pilot's recorded shape, and the
# same process-tree sampler, because `/usr/bin/time -l` reports a maximum
# over the tree and never a sum, and a concurrent run is exactly the case
# where those two differ. What this script adds is the worker axis and the
# comparison between a serial run and a concurrent one.
#
# Two things are measured and they answer different questions.
#
#   Cost. Wall time and aggregate memory per arm, repeated, with the spread
#   across repeats reported. A speed-up inside the repeat spread is not a
#   speed-up.
#
#   Identity. Each run writes a digest of everything the envelope claims is
#   exactly identical -- the fold allocation, the candidate space, the score
#   contract hashes, the certificate, the selected support -- beside the
#   numbers the engine produced. The orchestrator compares digests between
#   arms rather than eyeballing scores.
#
# Repeats are INTERLEAVED across arms (S, P, S, P, ...) rather than run in
# blocks. A machine that drifts -- another job starting, thermal
# throttling -- biases a blocked schedule entirely into whichever arm ran
# during the drift, and the arms here are compared against each other.
#
# Usage:
#   Rscript tools/benchmark_selection_parallel.R probe <outdir> <rows> \
#       <candidates> <clusters> <trials> <folds> <workers> <threads>
#   Rscript tools/benchmark_selection_parallel.R worker <outdir> <tag> \
#       <rows> <candidates> <clusters> <trials> <folds> <workers> <threads>
#   Rscript tools/benchmark_selection_parallel.R run <outdir> <reps> <cell>
#
# A <cell> is rows:candidates:clusters:trials:folds, and an arm is
# workers:threads, where `threads` is the INLA thread count the SERIAL arm
# is pinned to (0 leaves INLA at its own default). On a concurrent arm the
# package divides the threads itself and the value is ignored.

suppressWarnings(suppressMessages({
  this_file <- normalizePath(sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  ))
}))
pkg_dir <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(dirname(this_file), "fixture_selection_cost.R"))

args_vec <- commandArgs(trailingOnly = TRUE)
if (length(args_vec) < 2L) {
  stop("usage: benchmark_selection_parallel.R probe|worker|run <outdir> ...")
}
mode_arg <- args_vec[[1L]]
out_dir <- args_vec[[2L]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------------ #
# The sampler: aggregate memory over the whole process tree           #
# ------------------------------------------------------------------ #

# As P1-FIX's, and for its reason: a concurrent run holds several R
# processes and several INLA binaries at once, so the figure that says what
# the machine needs is the SUM across the processes alive at one instant,
# which no single-process instrument reports.
fb_par_sampler_start <- function(pid, path_sample, interval = "0.25") {
  awk_prog <- paste0(
    "{ p[NR]=$1; q[NR]=$2; r[NR]=$3; n=NR } ",
    "END { m[root]=1; ch=1; ",
    "while (ch) { ch=0; for (i=1;i<=n;i++) ",
    "if (!m[p[i]] && m[q[i]]) { m[p[i]]=1; ch=1 } } ",
    "t=0; x=0; c=0; ",
    "for (i=1;i<=n;i++) if (m[p[i]]) { t+=r[i]; if (r[i]>x) x=r[i]; c++ } ",
    "print t, x, c }"
  )
  cmd <- paste0(
    "while kill -0 ", pid, " 2>/dev/null; do ",
    "ps -A -o pid=,ppid=,rss= | awk -v root=", pid, " '", awk_prog, "' ",
    ">> ", shQuote(path_sample), "; ",
    "sleep ", interval, "; done"
  )
  system(cmd, wait = FALSE)
  invisible(NULL)
}


fb_par_sampler_read <- function(path_sample) {
  if (!file.exists(path_sample)) {
    return(list(aggregate_mb = NA_real_, process_mb = NA_real_, n_sample = 0L))
  }
  tab_s <- utils::read.table(path_sample, header = FALSE, col.names = c(
    "total_kb", "max_kb", "n_proc"
  ))
  if (nrow(tab_s) == 0L) {
    return(list(aggregate_mb = NA_real_, process_mb = NA_real_, n_sample = 0L))
  }
  list(
    aggregate_mb = max(tab_s$total_kb) / 1024,
    process_mb = max(tab_s$max_kb) / 1024,
    n_sample = nrow(tab_s),
    max_n_proc = max(tab_s$n_proc)
  )
}


# ------------------------------------------------------------------ #
# One run, measured, with its identity digest                         #
# ------------------------------------------------------------------ #

fb_par_load <- function() {
  suppressPackageStartupMessages(pkgload::load_all(pkg_dir, quiet = TRUE))
  options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE
  )
}


# fb_par_digest() -- everything the envelope claims is exactly identical,
# and the engine numbers it does not.
#
# The certificate is taken whole: the worker count is deliberately not on
# it (R/selection_parallel.R), so a serial and a concurrent certificate are
# comparable field for field, and taking it whole means a field added later
# is compared without this script being edited.
fb_par_digest <- function(sel) {
  cert <- sel$certificate
  list(
    folds = sel$folds,
    candidate_space = cert$design$candidate_space,
    certificate = cert,
    reference_contract = cert$reference_predictive$contract_hash,
    decision_contract = cert$decision_predictive$contract_hash,
    chosen_terms = sel$selection$chosen_terms,
    scores_model = as.character(sel$scores$model),
    scores = as.numeric(sel$scores$score),
    se = as.numeric(sel$scores$se),
    n_units = as.integer(sel$scores$n_units),
    increment = sel$increment$improvement,
    clusters = sel$clusters
  )
}


fb_par_one <- function(tag, n_rows, n_cand, n_clu, n_tri, n_fold, workers,
                       threads) {
  fb_par_load()
  if (threads > 0L) {
    INLA::inla.setOption(num.threads = paste0(threads, ":1"))
  }
  options(flexyBayes.selection_workers = as.integer(workers))
  print(paste0(
    "arm ", tag, ": workers=", workers,
    " serial_threads=", threads,
    " INLA num.threads=", INLA::inla.getOption("num.threads")
  ))

  tab_d <- fb_cost_fixture(n_rows, n_cand, n_clu, n_tri, seed = FB_COST_SEED)
  t0 <- proc.time()
  base_fit <- flexybayes(
    fb_cost_base_fixed(),
    random = fb_cost_base_random(),
    data = tab_d,
    backend = "inla",
    verbose = FALSE
  )
  s_base <- unname((proc.time() - t0)["elapsed"])

  res <- fb_resamples(unit = ~cluster, v = as.integer(n_fold))
  t0 <- proc.time()
  sel <- fb_select(
    base_fit,
    candidates = fb_cost_candidate_formula(n_cand),
    penalty = NULL,
    resamples = res,
    decision = fb_decision(method = "full", metric = "elpd"),
    seed = 1L
  )
  s_select <- unname((proc.time() - t0)["elapsed"])

  saveRDS(
    fb_par_digest(sel),
    file.path(out_dir, paste0("digest_", tag, ".rds"))
  )
  print(paste0(
    "RESULT|tag=", tag, "|rows=", n_rows, "|candidates=", n_cand,
    "|folds=", n_fold, "|workers=", workers, "|threads=", threads,
    "|base_s=", format(s_base, digits = 6),
    "|select_s=", format(s_select, digits = 6),
    "|folds_ran=", length(sel$folds)
  ))
  invisible(NULL)
}


fb_par_worker <- function(tag, n_rows, n_cand, n_clu, n_tri, n_fold, workers,
                          threads) {
  path_sample <- file.path(out_dir, paste0("sample_", tag, ".txt"))
  unlink(path_sample)
  fb_par_sampler_start(Sys.getpid(), path_sample)
  fb_par_one(tag, n_rows, n_cand, n_clu, n_tri, n_fold, workers, threads)
  samp <- fb_par_sampler_read(path_sample)
  print(paste0(
    "SAMPLER|aggregate_mb=", format(samp$aggregate_mb, digits = 6),
    "|process_mb=", format(samp$process_mb, digits = 6),
    "|n_sample=", samp$n_sample, "|max_n_proc=", samp$max_n_proc
  ))
  invisible(NULL)
}


# ------------------------------------------------------------------ #
# The orchestrator                                                    #
# ------------------------------------------------------------------ #

fb_par_field <- function(line, key) {
  m <- regmatches(line, regexpr(paste0(key, "=[^|\"]*"), line))
  if (length(m) == 0L) {
    return(NA_character_)
  }
  sub(paste0(key, "="), "", m)
}


fb_par_spread <- function(x) {
  x <- x[!is.na(x)]
  c(
    n = length(x),
    median = if (length(x)) stats::median(x) else NA_real_,
    min = if (length(x)) min(x) else NA_real_,
    max = if (length(x)) max(x) else NA_real_,
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_
  )
}


fb_par_cell <- function(tag, num, workers, threads) {
  rscript <- file.path(R.home("bin"), "Rscript")
  path_log <- file.path(out_dir, paste0("cell_", tag, ".log"))
  print(paste0(
    "running ", tag, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ))
  cmd <- paste(
    "/usr/bin/time -l", shQuote(rscript), "--no-init-file",
    shQuote(this_file), "worker", shQuote(out_dir), shQuote(tag),
    num[[1L]], num[[2L]], num[[3L]], num[[4L]], num[[5L]],
    workers, threads,
    ">", shQuote(path_log), "2>&1"
  )
  t0 <- proc.time()
  status <- system(cmd)
  s_outer <- unname((proc.time() - t0)["elapsed"])
  txt <- readLines(path_log, warn = FALSE)
  res_line <- grep("RESULT|", txt, value = TRUE, fixed = TRUE)
  samp_line <- grep("SAMPLER|", txt, value = TRUE, fixed = TRUE)
  rss_line <- grep("maximum resident set size", txt, value = TRUE)
  if (status != 0L || length(res_line) == 0L) {
    print(paste0(
      "  cell did not complete: status ", status, ", see ", basename(path_log)
    ))
    return(NULL)
  }
  data.frame(
    tag = tag,
    workers = as.integer(workers),
    threads = as.integer(threads),
    folds = as.integer(num[[5L]]),
    base_s = as.numeric(fb_par_field(res_line[[1L]], "base_s")),
    select_s = as.numeric(fb_par_field(res_line[[1L]], "select_s")),
    outer_s = s_outer,
    peak_process_mb = if (length(rss_line)) {
      as.numeric(sub("^ *([0-9]+).*", "\\1", rss_line[[1L]])) / 1048576
    } else {
      NA_real_
    },
    aggregate_mb = if (length(samp_line)) {
      as.numeric(fb_par_field(samp_line[[1L]], "aggregate_mb"))
    } else {
      NA_real_
    },
    max_n_proc = if (length(samp_line)) {
      as.integer(fb_par_field(samp_line[[1L]], "max_n_proc"))
    } else {
      NA_integer_
    },
    stringsAsFactors = FALSE
  )
}


# fb_par_compare() -- the envelope, read off two runs' digests.
fb_par_compare <- function(tag_a, tag_b) {
  path_a <- file.path(out_dir, paste0("digest_", tag_a, ".rds"))
  path_b <- file.path(out_dir, paste0("digest_", tag_b, ".rds"))
  if (!file.exists(path_a) || !file.exists(path_b)) {
    print(paste0("COMPARE|", tag_a, "|", tag_b, "|missing digest"))
    return(NULL)
  }
  a <- readRDS(path_a)
  b <- readRDS(path_b)
  exact <- c(
    "folds", "candidate_space", "certificate", "reference_contract",
    "decision_contract", "chosen_terms", "scores_model", "n_units"
  )
  for (nm in exact) {
    print(paste0(
      "COMPARE|", tag_a, "|", tag_b, "|", nm, "|identical=",
      identical(a[[nm]], b[[nm]])
    ))
  }
  d_score <- max(abs(a$scores - b$scores))
  d_se <- max(abs(a$se - b$se))
  d_inc <- abs((a$increment %||% NA_real_) - (b$increment %||% NA_real_))
  print(paste0(
    "COMPARE|", tag_a, "|", tag_b, "|numeric|max_abs_score=",
    format(d_score, digits = 6), "|max_abs_se=", format(d_se, digits = 6),
    "|abs_increment=", format(d_inc, digits = 6),
    "|same_order=", identical(order(a$scores), order(b$scores))
  ))
  invisible(list(score = d_score, se = d_se, increment = d_inc))
}


`%||%` <- function(x, y) if (is.null(x)) y else x


fb_par_run <- function(reps, cell, arms) {
  print(R.version.string)
  print(paste0("platform ", utils::sessionInfo()$platform))
  print(paste0("cores ", parallel::detectCores()))
  print(paste0("load average at start: ", paste(
    readLines(textConnection(system("uptime", intern = TRUE))), collapse = " "
  )))
  for (p1 in c("INLA", "brms", "glmnet", "posterior")) {
    v1 <- tryCatch(
      as.character(utils::packageVersion(p1)),
      error = function(e) "absent"
    )
    print(paste0("package ", p1, " ", v1))
  }

  num <- as.integer(strsplit(cell, ":", fixed = TRUE)[[1L]])
  if (length(num) != 5L) {
    stop("a cell is rows:candidates:clusters:trials:folds")
  }

  tab_all <- NULL
  # Interleaved: every arm runs its first repeat before any arm runs its
  # second, so a machine that drifts drifts through all of them.
  for (r1 in seq_len(reps)) {
    for (a1 in seq_along(arms)) {
      arm <- as.integer(strsplit(arms[[a1]], ":", fixed = TRUE)[[1L]])
      tag <- paste0("w", arm[[1L]], "t", arm[[2L]], "_rep", r1)
      row_local <- fb_par_cell(tag, num, arm[[1L]], arm[[2L]])
      if (!is.null(row_local)) {
        row_local$arm <- arms[[a1]]
        row_local$repeat_id <- r1
        tab_all <- rbind(tab_all, row_local)
        print(paste0(
          "  select_s ", format(row_local$select_s, digits = 6),
          "  aggregate_mb ", format(row_local$aggregate_mb, digits = 6),
          "  max_n_proc ", row_local$max_n_proc
        ))
      }
    }
  }

  if (is.null(tab_all)) {
    print("no cell completed")
    return(invisible(NULL))
  }

  for (a1 in unique(tab_all$arm)) {
    sub_tab <- tab_all[tab_all$arm == a1, , drop = FALSE]
    for (meas in c("select_s", "aggregate_mb")) {
      sp <- fb_par_spread(sub_tab[[meas]])
      print(paste0(
        "SPREAD|arm=", a1, "|measure=", meas, "|n=", sp[["n"]],
        "|median=", format(sp[["median"]], digits = 6),
        "|min=", format(sp[["min"]], digits = 6),
        "|max=", format(sp[["max"]], digits = 6),
        "|sd=", format(sp[["sd"]], digits = 4),
        "|rsd_pct=", format(100 * sp[["sd"]] / sp[["median"]], digits = 3)
      ))
    }
  }

  path_csv <- file.path(out_dir, "parallel_cells.csv")
  utils::write.csv(tab_all, path_csv, row.names = FALSE)
  print(paste0("wrote ", path_csv))
  print(tab_all, digits = 5)
  invisible(tab_all)
}


if (identical(mode_arg, "run")) {
  a <- args_vec[-(1:2)]
  fb_par_run(as.integer(a[[1L]]), a[[2L]], a[-(1:2)])
} else if (identical(mode_arg, "worker")) {
  a <- args_vec[-(1:2)]
  fb_par_worker(
    a[[1L]], as.integer(a[[2L]]), as.integer(a[[3L]]), as.integer(a[[4L]]),
    as.integer(a[[5L]]), as.integer(a[[6L]]), as.integer(a[[7L]]),
    as.integer(a[[8L]])
  )
} else if (identical(mode_arg, "probe")) {
  a <- args_vec[-(1:2)]
  fb_par_worker(
    "probe", as.integer(a[[1L]]), as.integer(a[[2L]]), as.integer(a[[3L]]),
    as.integer(a[[4L]]), as.integer(a[[5L]]), as.integer(a[[6L]]),
    as.integer(a[[7L]])
  )
} else if (identical(mode_arg, "compare")) {
  a <- args_vec[-(1:2)]
  fb_par_compare(a[[1L]], a[[2L]])
} else {
  stop("unknown mode: ", mode_arg)
}
