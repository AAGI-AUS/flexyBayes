#!/usr/bin/env Rscript

# stress_bigdata_extreme.R -- Push the streaming exact-aggregation path
# to extreme row counts under the realistic scenario where a collaborator
# shares a large dataset already partitioned into `.fst` shards on disk.
#
# For each target N the script (1) writes the shards one at a time so
# generation memory stays bounded, then (2) in a fresh subprocess under
# `/usr/bin/time -l`, streams every shard, accumulates per-cell
# sufficient statistics, and fits the K-cell model with INLA -- measuring
# the wall-clock and the OS peak resident memory of the whole fit. The
# headline evidence is that peak memory stays flat as N grows from tens
# of millions to a billion rows, because only one chunk plus the K-cell
# accumulator is ever resident.
#
# Usage:
#   Rscript tools/stress_bigdata_extreme.R                 # default ladder
#   Rscript tools/stress_bigdata_extreme.R run <N1,N2,...> <shard_n>
#   Rscript tools/stress_bigdata_extreme.R shards <N> <shard_n> <dir>
#   Rscript tools/stress_bigdata_extreme.R fit <dir> <chunk_rows>

suppressWarnings(suppressMessages({
  this_file <- normalizePath(sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  ))
}))
pkg_dir <- normalizePath(file.path(dirname(this_file), ".."))

N_ENV <- 6L
N_GENO <- 200L # 6 x 200 = 1200 cells, independent of N
ENV_EFF <- NULL
GENO_EFF <- NULL
.init_truth <- function() {
  set.seed(2024L)
  ENV_EFF <<- stats::rnorm(N_ENV, 0, 1)
  GENO_EFF <<- stats::rnorm(N_GENO, 0, 0.7)
}

.write_shards <- function(total_n, shard_n, dir) {
  if (!requireNamespace("fst", quietly = TRUE)) {
    stop("fst is required.", call. = FALSE)
  }
  .init_truth()
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  n_shards <- as.integer(ceiling(total_n / shard_n))
  written <- 0
  for (s in seq_len(n_shards)) {
    ns <- as.integer(min(shard_n, total_n - written))
    set.seed(1000L + s)
    env <- sample.int(N_ENV, ns, replace = TRUE)
    geno <- sample.int(N_GENO, ns, replace = TRUE)
    y <- ENV_EFF[env] + GENO_EFF[geno] + stats::rnorm(ns, 0, 1.5)
    shard <- data.frame(
      env = factor(env, levels = seq_len(N_ENV)),
      geno = factor(geno, levels = seq_len(N_GENO)),
      y = y
    )
    fst::write_fst(
      shard,
      file.path(dir, sprintf("part-%04d.fst", s)),
      compress = 60
    )
    written <- written + ns
    rm(shard, env, geno, y)
    gc(FALSE)
  }
  cat(sprintf("SHARDS n=%d total=%.0f dir=%s\n", n_shards, written, dir))
}

.fit_shards <- function(dir, chunk_rows) {
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE))
  options(flexyBayes.silence_uniform_inla_approx = TRUE)
  paths <- sort(list.files(dir, pattern = "\\.fst$", full.names = TRUE))
  pr <- flexyBayes::fb_prior(
    sigma ~ uniform(lower = 0, upper = 10),
    sd(group = "geno") ~ uniform(lower = 0, upper = 10)
  )
  t_all <- system.time(
    fit <- flexyBayes::flexybayes_stream(
      y ~ env,
      random = ~geno,
      source = paths,
      family = "gaussian",
      chunk_rows = chunk_rows,
      prior = pr,
      fit = TRUE,
      verbose = FALSE
    )
  )[["elapsed"]]
  meta <- fit$extras$aggregation_meta
  b <- stats::coef(fit)
  cat(sprintf(
    "RESULT N=%.0f K=%d compression=%.6f t_all=%.2f intercept=%.5f\n",
    meta$N,
    meta$K,
    meta$compression,
    t_all,
    unname(b[1L])
  ))
}

parse_rss_mb <- function(lines) {
  hit <- grep("maximum resident set size", lines, value = TRUE)
  if (!length(hit)) {
    return(NA_real_)
  }
  as.numeric(sub("^\\s*([0-9]+).*$", "\\1", hit[[1L]])) / 1024^2
}
dir_size_gb <- function(dir) {
  fs <- list.files(dir, pattern = "\\.fst$", full.names = TRUE)
  sum(file.info(fs)$size) / 1024^3
}

# ---- dispatch ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "run"

if (identical(mode, "shards")) {
  .write_shards(as.numeric(args[[2L]]), as.numeric(args[[3L]]), args[[4L]])
  quit(save = "no")
}
if (identical(mode, "fit")) {
  .fit_shards(args[[2L]], as.numeric(args[[3L]]))
  quit(save = "no")
}

# ---- orchestrate the ladder -------------------------------------------
targets <- if (length(args) >= 2L) {
  as.numeric(strsplit(args[[2L]], ",")[[1L]])
} else {
  c(1e7, 1e8, 1e9)
}
shard_n <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 5e7
chunk_rows <- 5e6
work_dir <- file.path(tempdir(), "fb_shards")

cat("flexyBayes extreme streaming demo (partitioned .fst, INLA)\n")
cat(
  "design: ",
  N_ENV,
  " env x ",
  N_GENO,
  " geno = ",
  N_ENV * N_GENO,
  " cells\n\n",
  sep = ""
)

rows <- list()
for (tn in targets) {
  cat(sprintf(
    "== target N = %s ==\n",
    format(tn, big.mark = ",", scientific = FALSE)
  ))
  unlink(work_dir, recursive = TRUE)
  dir.create(work_dir, recursive = TRUE)

  t_gen <- system.time(
    system2(
      "Rscript",
      c(
        "--vanilla",
        this_file,
        "shards",
        format(tn, scientific = FALSE),
        format(shard_n, scientific = FALSE),
        work_dir
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  )[["elapsed"]]
  disk_gb <- dir_size_gb(work_dir)
  cat(sprintf("  shards written in %.1fs, on-disk %.2f GB\n", t_gen, disk_gb))

  out <- system2(
    "/usr/bin/time",
    c(
      "-l",
      "Rscript",
      "--vanilla",
      this_file,
      "fit",
      work_dir,
      format(chunk_rows, scientific = FALSE)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  res_line <- grep("^RESULT", out, value = TRUE)
  if (length(res_line) != 1L) {
    cat("  FIT FAILED; tail of output:\n")
    cat(paste0("    ", utils::tail(out, 8)), sep = "\n")
    next
  }
  grab <- function(k) {
    as.numeric(sub(paste0(".*", k, "=([-0-9.]+).*"), "\\1", res_line))
  }
  rows[[length(rows) + 1L]] <- data.frame(
    N = grab("N"),
    K = grab("K"),
    compression = grab("compression"),
    gen_s = t_gen,
    fit_s = grab("t_all"),
    peak_mb = parse_rss_mb(out),
    disk_gb = disk_gb,
    intercept = grab("intercept")
  )
  print(rows[[length(rows)]], row.names = FALSE)
  cat("\n")
}
unlink(work_dir, recursive = TRUE)

res <- do.call(rbind, rows)
cat("\n===== extreme streaming summary =====\n")
print(res, row.names = FALSE)

out_dir <- file.path(pkg_dir, "..", "benchmark_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date())
write.csv(
  res,
  file.path(out_dir, paste0("bigdata_extreme_", stamp, ".csv")),
  row.names = FALSE
)
saveRDS(res, file.path(out_dir, paste0("bigdata_extreme_", stamp, ".rds")))
cat("artefacts written to ", normalizePath(out_dir), "\n", sep = "")
