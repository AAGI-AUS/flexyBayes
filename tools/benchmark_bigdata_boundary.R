#!/usr/bin/env Rscript

# benchmark_bigdata_boundary.R -- Honest performance comparison of the
# per-row and the streaming-aggregated fit paths at sample sizes where
# both are feasible, plus the crossover where the per-row path becomes
# infeasible while the aggregated path stays cheap.
#
# Each (method, N) cell is run in a fresh R subprocess under
# `/usr/bin/time -l` so wall-clock and peak resident memory are measured
# by the OS, not estimated inside a long-lived session. The orchestrator
# collects the rows, prints a table, and writes a CSV + RDS artefact to
# benchmark_results/.
#
# Usage:
#   Rscript tools/benchmark_bigdata_boundary.R              # orchestrate
#   Rscript tools/benchmark_bigdata_boundary.R worker <method> <N>

suppressWarnings(suppressMessages({
  this_file <- normalizePath(sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  ))
}))
pkg_dir <- normalizePath(file.path(dirname(this_file), ".."))

# Fixed factorial design: 6 environments x 60 genotypes = 360 cells,
# regardless of N. This is the regime aggregation is built for.
N_ENV <- 6L
N_GENO <- 60L

.make_data <- function(n, seed = 42L) {
  set.seed(seed)
  env_eff <- stats::rnorm(N_ENV, 0, 1)
  geno_eff <- stats::rnorm(N_GENO, 0, 0.7)
  env <- sample.int(N_ENV, n, replace = TRUE)
  geno <- sample.int(N_GENO, n, replace = TRUE)
  y <- env_eff[env] + geno_eff[geno] + stats::rnorm(n, 0, 1.5)
  data.frame(env = factor(env), geno = factor(geno), y = y)
}


# ---- worker mode -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && identical(args[[1L]], "worker")) {
  method <- args[[2L]]
  n <- as.numeric(args[[3L]])
  suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE))
  options(flexyBayes.silence_uniform_inla_approx = TRUE)
  df <- .make_data(n)
  pr <- flexyBayes::fb_prior(
    sigma ~ uniform(lower = 0, upper = 10),
    sd(group = "geno") ~ uniform(lower = 0, upper = 10)
  )

  t <- system.time({
    if (identical(method, "per_row")) {
      fit <- flexyBayes::flexybayes(
        y ~ env,
        random = ~geno,
        data = df,
        backend = "inla",
        aggregate = FALSE,
        prior = pr,
        verbose = FALSE
      )
      kk <- NA_integer_
    } else {
      fit <- flexyBayes::flexybayes_stream(
        y ~ env,
        random = ~geno,
        source = df,
        backend = "inla",
        chunk_rows = 2e6,
        prior = pr,
        verbose = FALSE
      )
      kk <- fit$extras$aggregation_meta$K
    }
  })[["elapsed"]]
  b <- stats::coef(fit)
  cat(sprintf(
    "RESULT elapsed=%.3f K=%s intercept=%.5f\n",
    t,
    as.character(kk),
    unname(b[1L])
  ))
  quit(save = "no", status = 0L)
}


# ---- orchestrator ------------------------------------------------------
parse_rss_bytes <- function(lines) {
  # macOS `/usr/bin/time -l` reports "<bytes>  maximum resident set size".
  hit <- grep("maximum resident set size", lines, value = TRUE)
  if (!length(hit)) {
    return(NA_real_)
  }
  as.numeric(sub("^\\s*([0-9]+).*$", "\\1", hit[[1L]]))
}

run_cell <- function(method, n, timeout_s = 600) {
  out <- tryCatch(
    system2(
      "/usr/bin/time",
      args = c(
        "-l",
        "Rscript",
        "--vanilla",
        this_file,
        "worker",
        method,
        format(n, scientific = FALSE)
      ),
      stdout = TRUE,
      stderr = TRUE,
      timeout = timeout_s
    ),
    error = function(e) attr(e, "result") %||% character(0)
  )
  res_line <- grep("^RESULT", out, value = TRUE)
  status <- attr(out, "status")
  ok <- length(res_line) == 1L && (is.null(status) || status == 0L)
  if (!ok) {
    return(data.frame(
      method = method,
      n = n,
      ok = FALSE,
      elapsed_s = NA_real_,
      peak_mb = NA_real_,
      K = NA_integer_,
      intercept = NA_real_
    ))
  }
  el <- as.numeric(sub(".*elapsed=([0-9.]+).*", "\\1", res_line))
  kk <- suppressWarnings(as.integer(sub(".*K=([0-9NA]+).*", "\\1", res_line)))
  ic <- as.numeric(sub(".*intercept=([-0-9.]+).*", "\\1", res_line))
  data.frame(
    method = method,
    n = n,
    ok = TRUE,
    elapsed_s = el,
    peak_mb = parse_rss_bytes(out) / 1024^2,
    K = kk,
    intercept = ic
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

grid <- list(
  list(method = "per_row", n = 1e5),
  list(method = "per_row", n = 5e5),
  list(method = "per_row", n = 1e6),
  list(method = "per_row", n = 5e6),
  list(method = "streamed", n = 1e5),
  list(method = "streamed", n = 1e6),
  list(method = "streamed", n = 5e6),
  list(method = "streamed", n = 1e7),
  list(method = "streamed", n = 5e7)
)

cat("flexyBayes big-data boundary benchmark\n")
cat(
  "design: ",
  N_ENV,
  "env x ",
  N_GENO,
  "geno = ",
  N_ENV * N_GENO,
  " cells; backend INLA\n\n",
  sep = ""
)

rows <- lapply(grid, function(g) {
  cat(sprintf(
    "  running %-9s N = %s ...\n",
    g$method,
    format(g$n, big.mark = ",", scientific = FALSE)
  ))
  run_cell(g$method, g$n)
})
res <- do.call(rbind, rows)
res$compression <- ifelse(
  res$method == "streamed",
  (N_ENV * N_GENO) / res$n,
  NA_real_
)

print(res, row.names = FALSE)

out_dir <- file.path(pkg_dir, "..", "benchmark_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date())
write.csv(
  res,
  file.path(out_dir, paste0("bigdata_boundary_", stamp, ".csv")),
  row.names = FALSE
)
saveRDS(res, file.path(out_dir, paste0("bigdata_boundary_", stamp, ".rds")))
cat("\nartefacts written to ", normalizePath(out_dir), "\n", sep = "")
