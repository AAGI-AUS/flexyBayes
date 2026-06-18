#!/usr/bin/env Rscript

# study_bigdata_equivalence.R -- Side-by-side numerical comparison of the
# per-row and the exact-aggregated fitting paths, on the same data and
# matched priors, across moderate-to-border sample sizes. For each
# parameter the script records the per-row estimate, the aggregated
# estimate, and their difference, so the two paths can be assessed
# directly. INLA backend (deterministic Laplace), so any difference is
# numerical, not Monte-Carlo.
#
# Output: benchmark_results/bigdata_study_estimates_<date>.csv (long
# format: family, N, K, param, scale, per_row, aggregated, abs_diff).

suppressWarnings(suppressMessages({
  this_file <- normalizePath(sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  ))
}))
pkg_dir <- normalizePath(file.path(dirname(this_file), ".."))
suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE))
options(flexyBayes.silence_uniform_inla_approx = TRUE)

N_ENV <- 4L
N_GENO <- 30L

make_data <- function(family, n, seed = 2024L) {
  set.seed(seed)
  env_eff <- c(0, 1.2, -0.7, 0.5)
  geno_eff <- stats::rnorm(N_GENO, 0, 0.6)
  env <- sample.int(N_ENV, n, replace = TRUE)
  geno <- sample.int(N_GENO, n, replace = TRUE)
  eta <- env_eff[env] + geno_eff[geno]
  y <- switch(
    family,
    gaussian = eta + stats::rnorm(n, 0, 1.5),
    poisson = stats::rpois(n, exp(eta - 0.5)),
    binomial = stats::rbinom(n, 1, stats::plogis(eta - 1))
  )
  data.frame(env = factor(env), geno = factor(geno), y = y)
}

prior_for <- function(family) {
  if (identical(family, "gaussian")) {
    fb_prior(
      sigma ~ uniform(lower = 0, upper = 10),
      sd(group = "geno") ~ uniform(lower = 0, upper = 10)
    )
  } else {
    fb_prior(sd(group = "geno") ~ uniform(lower = 0, upper = 10))
  }
}

# Pull a uniform parameter table off a fitted INLA-backed flexybayes fit.
extract <- function(fit) {
  sf <- fit$inla$summary.fixed
  h <- fit$inla$summary.hyperpar
  out <- data.frame(
    param = rownames(sf),
    scale = "coef",
    est = sf$mean,
    sd = sf$sd,
    stringsAsFactors = FALSE
  )
  prec_resid <- h$mean[grepl("Gaussian observations", rownames(h))]
  if (length(prec_resid)) {
    out <- rbind(
      out,
      data.frame(
        param = "sigma (residual SD)",
        scale = "sd",
        est = 1 / sqrt(prec_resid),
        sd = NA
      )
    )
  }
  prec_geno <- h$mean[grepl("^Precision for geno$", rownames(h))]
  if (length(prec_geno)) {
    out <- rbind(
      out,
      data.frame(
        param = "tau (geno SD)",
        scale = "sd",
        est = 1 / sqrt(prec_geno),
        sd = NA
      )
    )
  }
  out
}

families <- c("gaussian", "binomial", "poisson")
grid <- c(1e4, 1e5, 1e6)
rows <- list()

for (fam in families) {
  for (n in grid) {
    df <- make_data(fam, n)
    pr <- prior_for(fam)
    fr <- flexybayes(
      y ~ env,
      random = ~geno,
      data = df,
      family = fam,
      backend = "inla",
      aggregate = FALSE,
      prior = pr,
      verbose = FALSE
    )
    fa <- flexybayes(
      y ~ env,
      random = ~geno,
      data = df,
      family = fam,
      backend = "inla",
      aggregate = TRUE,
      prior = pr,
      verbose = FALSE
    )
    er <- extract(fr)
    ea <- extract(fa)
    m <- merge(er, ea, by = c("param", "scale"), suffixes = c("_row", "_agg"))
    m$abs_diff <- abs(m$est_row - m$est_agg)
    k <- fa$extras$aggregation_meta$K
    rows[[length(rows) + 1L]] <- data.frame(
      family = fam,
      N = n,
      K = k,
      param = m$param,
      scale = m$scale,
      per_row = round(m$est_row, 6),
      aggregated = round(m$est_agg, 6),
      abs_diff = signif(m$abs_diff, 3)
    )
    cat(sprintf(
      "%-9s N=%.0e K=%d  max|diff|=%.2e\n",
      fam,
      n,
      k,
      max(m$abs_diff)
    ))
  }
}

res <- do.call(rbind, rows)
out_dir <- file.path(pkg_dir, "..", "benchmark_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
f <- file.path(out_dir, "bigdata_study_estimates_2026-06-01.csv")
write.csv(res, f, row.names = FALSE)
cat("\nwrote ", normalizePath(f), "\n", sep = "")
