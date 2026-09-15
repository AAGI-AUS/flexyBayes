# tools/calibrate_boundary_collapse.R -- measure the floor the
# boundary-collapse detector is calibrated against.
#
# The detector in R/boundary_collapse.R warns when a variance component's
# upper credible bound is a negligible fraction of the residual SD. The
# question that fixes the threshold is: how small does that ratio get
# when the component is GENUINELY zero? Anything at or above that floor
# must not warn, or honest nulls are reported as failed modes.
#
# Design: group SD set to exactly zero, crossed over n and the number of
# groups, 8 seeds per cell, on the package's default prior. A non-null
# arm (group SD = 0.8) runs alongside as a positive control that the
# ratio separates at all.
#
# Run from the package root; needs INLA. Writes
# inst/validation/boundary_calibration.csv, which
# tests/testthat/test-boundary-collapse.R reads so that the threshold is
# asserted against measured data rather than a number typed into a mock.
#
#   Rscript tools/calibrate_boundary_collapse.R

suppressMessages(library(INLA))
INLA::inla.setOption(num.threads = "2:1")

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.", call. = FALSE)
}

fit_ratio <- function(n, k, sd_g, seed) {
  set.seed(seed)
  g <- factor(rep(seq_len(k), length.out = n))
  u <- rnorm(k, 0, sd_g)
  y <- as.numeric(u[as.integer(g)]) + rnorm(n)
  d <- data.frame(y = y, g = as.integer(g))
  r <- try(
    INLA::inla(y ~ 1 + f(g, model = "iid"), family = "gaussian",
               data = d, control.compute = list(return.marginals = FALSE)),
    silent = TRUE
  )
  if (inherits(r, "try-error")) {
    return(NA_real_)
  }
  h <- r$summary.hyperpar
  ig <- grep("^Precision for g$", rownames(h))
  ie <- grep("Gaussian observations", rownames(h))
  if (!length(ig) || !length(ie)) {
    return(NA_real_)
  }
  # on the SD scale: q97.5(sd) = 1 / sqrt(q0.025(precision))
  (1 / sqrt(h[ig, "0.025quant"])) / (1 / sqrt(h[ie, "0.5quant"]))
}

grid <- expand.grid(
  n = c(30L, 60L, 120L, 240L, 480L),
  k = c(5L, 10L, 20L),
  arm = c("null", "real"),
  seed = 1:8,
  stringsAsFactors = FALSE
)
grid <- grid[grid$k * 3L <= grid$n, ]
grid$ratio <- NA_real_
for (i in seq_len(nrow(grid))) {
  grid$ratio[i] <- fit_ratio(
    grid$n[i], grid$k[i],
    if (grid$arm[i] == "null") 0 else 0.8,
    grid$seed[i] * 1000L + i
  )
}

grid <- grid[order(grid$arm, grid$n, grid$k, grid$seed),
             c("arm", "n", "k", "seed", "ratio")]
grid$ratio <- signif(grid$ratio, 6)
out <- "inst/validation/boundary_calibration.csv"
write.csv(grid, out, row.names = FALSE)

nu <- grid[grid$arm == "null" & is.finite(grid$ratio), ]
re <- grid[grid$arm == "real" & is.finite(grid$ratio), ]
cat(sprintf("wrote %s: %d rows (%d null, %d real)\n",
            out, nrow(grid), nrow(nu), nrow(re)))
cat(sprintf("null floor: min %.4g, 5th pct %.4g, median %.4g\n",
            min(nu$ratio), quantile(nu$ratio, 0.05), median(nu$ratio)))
cat(sprintf("real arm  : min %.4g, median %.4g\n",
            min(re$ratio), median(re$ratio)))
