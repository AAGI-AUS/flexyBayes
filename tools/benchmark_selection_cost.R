#!/usr/bin/env Rscript

# benchmark_selection_cost.R -- what one fb_select() call costs on the
# package-owned fixture at the pilot's recorded shapes.
#
# Every figure this script reports is measured. Nothing is scaled from a
# smaller cell and reported as if it had been run.
#
# How a cell is measured. Each (cell, repeat) runs in a fresh R
# subprocess under `/usr/bin/time -l`, so the wall time and the peak
# resident set are read from the operating system rather than from inside
# a long-lived session that has already allocated. Two memory figures are
# reported and they are not the same quantity:
#
#   * peak process memory -- the largest resident set any single process
#     in the tree reached. This is what `/usr/bin/time -l` reports:
#     measured on this machine, its "maximum resident set size" follows a
#     waited-for child rather than the parent, so it is a maximum over the
#     tree and never a sum.
#   * aggregate memory -- the largest total resident set across every
#     process alive at one instant. INLA fits by handing the problem to a
#     separate binary, so R and that binary hold memory at the same time
#     and the sum is the figure that says what the machine needs. It is
#     sampled by polling the worker's process tree while the worker runs.
#
# The fold count is read back from the returned object rather than from
# the request, so a cell reports the folds that ran.
#
# Usage:
#   Rscript tools/benchmark_selection_cost.R run <outdir> <cell> [<cell> ...]
#   Rscript tools/benchmark_selection_cost.R worker <outdir> <tag> <rows> \
#       <candidates> <clusters> <trials> <folds> <seed>
#   Rscript tools/benchmark_selection_cost.R decompose <outdir> <rows> \
#       <candidates> <clusters> <trials> <folds> <seed>
#
# A <cell> is name:rows:candidates:clusters:trials:folds:repeats, for
# example  barley:10890:82:22:48:3:3
#
# `decompose` attributes one cell's cost between the engine refits and the
# posterior draws that score them, by running the fold loop through the
# same internals fb_select() calls. It exists because the per-fit timing
# the selection route records is not carried onto the returned object, so
# the split cannot be read back from a completed call.

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
  stop("usage: benchmark_selection_cost.R run|worker|decompose <outdir> ...")
}
mode_arg <- args_vec[[1L]]
out_dir <- args_vec[[2L]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------------ #
# The sampler: aggregate memory over the worker's process tree        #
# ------------------------------------------------------------------ #

# fb_cost_sampler_start() -- poll `ps` for every process descended from
# `pid` and append one line per sample: total resident kilobytes, the
# largest single process's kilobytes, and the process count. The loop
# exits on its own when the worker goes away, so a worker that dies still
# leaves the samples it produced.
#
# The descendant set is recomputed at every poll because the INLA binary
# appears and disappears within a call.
fb_cost_sampler_start <- function(pid, path_sample, interval = "0.25") {
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


# fb_cost_sampler_read() -- the peaks the sampler saw, in megabytes.
fb_cost_sampler_read <- function(path_sample) {
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
# The worker: one fb_select() call, measured                          #
# ------------------------------------------------------------------ #

fb_cost_load <- function() {
  suppressPackageStartupMessages(pkgload::load_all(pkg_dir, quiet = TRUE))
  options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE
  )
}


fb_cost_worker <- function(tag, n_rows, n_cand, n_clu, n_tri, n_fold, seed) {
  path_sample <- file.path(out_dir, paste0("sample_", tag, ".txt"))
  unlink(path_sample)
  fb_cost_sampler_start(Sys.getpid(), path_sample)

  fb_cost_load()
  print(paste0("INLA num.threads: ", INLA::inla.getOption("num.threads")))
  print(paste0(
    "shape rows=", n_rows, " candidates=", n_cand, " clusters=", n_clu,
    " trials=", n_tri, " folds=", n_fold, " seed=", seed
  ))

  tab_d <- fb_cost_fixture(n_rows, n_cand, n_clu, n_tri, seed = seed)
  print(paste0(
    "fixture rows ", nrow(tab_d), " clusters ", nlevels(tab_d$cluster),
    " trials ", nlevels(tab_d$trial), " years ", nlevels(tab_d$year)
  ))

  # No `seed` on the base fit: INLA's nested Laplace approximation draws no
  # random numbers, so the argument has no effect and the package says so.
  t0 <- proc.time()
  base_fit <- flexybayes(
    fb_cost_base_fixed(),
    random = fb_cost_base_random(),
    data = tab_d,
    backend = "inla",
    verbose = FALSE
  )
  s_base <- unname((proc.time() - t0)["elapsed"])
  print(paste0("base fit seconds ", format(s_base, digits = 6)))

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

  # Three readbacks of the fold count, none of them the request: the folds
  # the object carries, the certificate's resampling record, and the
  # candidate-space record, which compiles one entry per fold.
  n_fold_ran <- length(sel$folds)
  n_fold_cert <- sel$certificate$resampling$n_folds
  space <- sel$certificate$design$candidate_space
  n_fold_space <- if (is.null(space)) NA_integer_ else space$n_folds
  n_dropped <- if (is.null(space) || is.null(space$dropped)) {
    NA_integer_
  } else {
    nrow(as.data.frame(space$dropped))
  }
  gc_mat <- gc()
  heap_mb <- sum(gc_mat[, which(colnames(gc_mat) == "max used") + 1L])

  print(paste0("fb_select seconds ", format(s_select, digits = 6)))
  print(paste0(
    "folds: object ", n_fold_ran, ", certificate ", n_fold_cert,
    ", candidate space ", n_fold_space, "; dropped rows ", n_dropped
  ))
  print(paste0("R heap peak Mb ", format(heap_mb, digits = 6)))

  # The separator is "|" rather than a tab: print() escapes a tab into a
  # backslash and a t, and the orchestrator would not find the field.
  print(paste0(
    "RESULT|tag=", tag, "|rows=", n_rows, "|candidates=", n_cand,
    "|clusters=", n_clu, "|trials=", n_tri, "|folds_requested=", n_fold,
    "|folds_ran=", n_fold_ran, "|folds_cert=", n_fold_cert,
    "|folds_space=", n_fold_space, "|dropped_rows=", n_dropped,
    "|base_s=", format(s_base, digits = 6),
    "|select_s=", format(s_select, digits = 6),
    "|heap_mb=", format(heap_mb, digits = 6)
  ))

  samp <- fb_cost_sampler_read(path_sample)
  print(paste0(
    "SAMPLER|aggregate_mb=", format(samp$aggregate_mb, digits = 6),
    "|process_mb=", format(samp$process_mb, digits = 6),
    "|n_sample=", samp$n_sample, "|max_n_proc=", samp$max_n_proc
  ))
  invisible(NULL)
}


# ------------------------------------------------------------------ #
# The decomposition: engine refits against posterior draws            #
# ------------------------------------------------------------------ #

fb_cost_decompose <- function(n_rows, n_cand, n_clu, n_tri, n_fold, seed) {
  fb_cost_load()
  ns <- asNamespace("flexyBayes")
  print(paste0("INLA num.threads: ", INLA::inla.getOption("num.threads")))
  print(paste0(
    "decompose rows=", n_rows, " candidates=", n_cand, " clusters=", n_clu,
    " trials=", n_tri, " folds=", n_fold
  ))

  tab_d <- fb_cost_fixture(n_rows, n_cand, n_clu, n_tri, seed = seed)
  cand_f <- fb_cost_candidate_formula(n_cand)
  base_fit <- ns$.fb_sel_bayes_normalise_base(flexybayes(
    fb_cost_base_fixed(),
    random = fb_cost_base_random(),
    data = tab_d,
    backend = "inla",
    verbose = FALSE,
    seed = 1L
  ))
  spec <- ns$.fb_sel_bayes_base_spec(base_fit)
  provider <- ns$.fb_sel_bayes_provider(base_fit)
  folds <- ns$.fb_sel_folds(
    fb_resamples(unit = ~cluster, v = as.integer(n_fold)), tab_d, seed = 1L
  )
  space <- ns$.fb_sel_candidate_space(
    spec$protected, cand_f, tab_d, folds,
    groupings = c("cluster", "year", "trial"), route = "bayes"
  )
  cand_lab <- attr(stats::terms(cand_f, data = tab_d), "term.labels")
  sets <- list(protected_only = integer(0), full = seq_along(cand_lab))

  for (f1 in seq_along(folds)) {
    split <- folds[[f1]]
    inputs <- ns$.fb_sel_bayes_fold_inputs(
      tab_d, split,
      response = spec$response, protected = spec$protected,
      candidates = cand_f, random = spec$random_vars,
      prediction = "new_cluster"
    )
    dropped <- space$folds[[split$label]]$dropped_terms
    for (s1 in seq_along(sets)) {
      name <- names(sets)[[s1]]
      fixed <- ns$.fb_sel_bayes_fixed(
        spec$response, spec$protected_labels,
        setdiff(cand_lab[sets[[name]]], dropped), spec$offsets
      )
      fit_seed <- 1000L * f1 + s1
      t0 <- proc.time()
      fit <- ns$.fb_sel_bayes_refit(
        base_fit, inputs$engine_data, fixed, provider,
        seed = fit_seed, na_action = "augment"
      )
      s_fit <- unname((proc.time() - t0)["elapsed"])
      t0 <- proc.time()
      ns$.fb_sel_posterior_predict_fold(
        fit, split$assess_rows, inputs$engine_data,
        prediction = "new_cluster", n_draws = 4000L, seed = fit_seed
      )
      s_draw <- unname((proc.time() - t0)["elapsed"])
      print(paste0(
        "DECOMP|fold=", f1, "|model=", name,
        "|terms=", length(setdiff(cand_lab[sets[[name]]], dropped)),
        "|dropped=", length(dropped),
        "|refit_s=", format(s_fit, digits = 6),
        "|draws_s=", format(s_draw, digits = 6)
      ))
    }
  }
  invisible(NULL)
}


# ------------------------------------------------------------------ #
# The orchestrator                                                    #
# ------------------------------------------------------------------ #

fb_cost_spread <- function(x) {
  c(
    n = length(x),
    median = stats::median(x),
    min = min(x),
    max = max(x),
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_
  )
}


fb_cost_run <- function(cell_vec) {
  rscript <- file.path(R.home("bin"), "Rscript")
  print(R.version.string)
  print(paste0("platform ", utils::sessionInfo()$platform))
  print(paste0("cores ", parallel::detectCores()))
  for (p1 in c("INLA", "brms", "glmnet", "cmdstanr", "posterior")) {
    v1 <- tryCatch(
      as.character(utils::packageVersion(p1)),
      error = function(e) "absent"
    )
    print(paste0("package ", p1, " ", v1))
  }

  tab_all <- NULL
  for (cell in cell_vec) {
    part <- strsplit(cell, ":", fixed = TRUE)[[1L]]
    if (length(part) != 7L) {
      stop("a cell is name:rows:candidates:clusters:trials:folds:repeats")
    }
    nm <- part[[1L]]
    num <- as.integer(part[-1L])
    reps <- num[[6L]]
    row_local <- NULL

    for (r1 in seq_len(reps)) {
      tag <- paste0(nm, "_rep", r1)
      path_log <- file.path(out_dir, paste0("cell_", tag, ".log"))
      print(paste0(
        "running ", tag, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ))
      cmd <- paste(
        "/usr/bin/time -l", shQuote(rscript), "--no-init-file",
        shQuote(this_file), "worker", shQuote(out_dir), shQuote(tag),
        num[[1L]], num[[2L]], num[[3L]], num[[4L]], num[[5L]], 20260916L,
        ">", shQuote(path_log), "2>&1"
      )
      t0 <- proc.time()
      status <- system(cmd)
      s_wall_outer <- unname((proc.time() - t0)["elapsed"])
      txt <- readLines(path_log, warn = FALSE)

      # `/usr/bin/time -l` writes its report on stderr, which is in the
      # same log; a failed cell is reported rather than skipped.
      real_line <- grep("real .* user .* sys", txt, value = TRUE)
      rss_line <- grep("maximum resident set size", txt, value = TRUE)
      res_line <- grep("RESULT|", txt, value = TRUE, fixed = TRUE)
      samp_line <- grep("SAMPLER|", txt, value = TRUE, fixed = TRUE)
      if (status != 0L || length(res_line) == 0L) {
        print(paste0("  cell did not complete: status ", status,
                     ", see ", basename(path_log)))
        next
      }
      field <- function(line, key) {
        m <- regmatches(line, regexpr(paste0(key, "=[^|\"]*"), line))
        if (length(m) == 0L) return(NA_character_)
        sub(paste0(key, "="), "", m)
      }
      wall_s <- if (length(real_line)) {
        as.numeric(sub("^ *([0-9.]+) real.*", "\\1", real_line[[1L]]))
      } else {
        s_wall_outer
      }
      proc_mb <- if (length(rss_line)) {
        as.numeric(sub("^ *([0-9]+).*", "\\1", rss_line[[1L]])) / 1048576
      } else {
        NA_real_
      }
      row_local <- rbind(row_local, data.frame(
        cell = nm, repeat_id = r1,
        rows = num[[1L]], candidates = num[[2L]],
        clusters = num[[3L]], trials = num[[4L]],
        folds_ran = as.integer(field(res_line[[1L]], "folds_ran")),
        folds_cert = as.integer(field(res_line[[1L]], "folds_cert")),
        folds_space = as.integer(field(res_line[[1L]], "folds_space")),
        dropped_rows = as.integer(field(res_line[[1L]], "dropped_rows")),
        base_s = as.numeric(field(res_line[[1L]], "base_s")),
        select_s = as.numeric(field(res_line[[1L]], "select_s")),
        wall_s = wall_s,
        peak_process_mb = proc_mb,
        aggregate_mb = if (length(samp_line)) {
          as.numeric(field(samp_line[[1L]], "aggregate_mb"))
        } else {
          NA_real_
        },
        sampler_process_mb = if (length(samp_line)) {
          as.numeric(field(samp_line[[1L]], "process_mb"))
        } else {
          NA_real_
        },
        heap_mb = as.numeric(field(res_line[[1L]], "heap_mb")),
        stringsAsFactors = FALSE
      ))
      print(paste0(
        "  select_s ", field(res_line[[1L]], "select_s"),
        "  wall_s ", format(wall_s, digits = 6),
        "  peak_process_mb ", format(proc_mb, digits = 5),
        "  aggregate_mb ",
        if (length(samp_line)) field(samp_line[[1L]], "aggregate_mb") else "not sampled"
      ))
    } # ends r1, over repeats

    if (!is.null(row_local)) {
      tab_all <- rbind(tab_all, row_local)
      for (meas in c("select_s", "wall_s", "peak_process_mb", "aggregate_mb")) {
        sp <- fb_cost_spread(row_local[[meas]])
        print(paste0(
          "SPREAD|cell=", nm, "|measure=", meas,
          "|n=", sp[["n"]], "|median=", format(sp[["median"]], digits = 6),
          "|min=", format(sp[["min"]], digits = 6),
          "|max=", format(sp[["max"]], digits = 6),
          "|sd=", format(sp[["sd"]], digits = 4)
        ))
      }
    }
  } # ends cell, over cells

  if (!is.null(tab_all)) {
    path_csv <- file.path(out_dir, "cost_cells.csv")
    utils::write.csv(tab_all, path_csv, row.names = FALSE)
    print(paste0("wrote ", path_csv))
    print(tab_all, digits = 5)
  }
  invisible(tab_all)
}


if (identical(mode_arg, "run")) {
  fb_cost_run(args_vec[-(1:2)])
} else if (identical(mode_arg, "worker")) {
  a <- args_vec[-(1:2)]
  fb_cost_worker(
    a[[1L]], as.integer(a[[2L]]), as.integer(a[[3L]]), as.integer(a[[4L]]),
    as.integer(a[[5L]]), as.integer(a[[6L]]), as.integer(a[[7L]])
  )
} else if (identical(mode_arg, "decompose")) {
  a <- args_vec[-(1:2)]
  fb_cost_decompose(
    as.integer(a[[1L]]), as.integer(a[[2L]]), as.integer(a[[3L]]),
    as.integer(a[[4L]]), as.integer(a[[5L]]), as.integer(a[[6L]])
  )
} else {
  stop("unknown mode: ", mode_arg)
}
