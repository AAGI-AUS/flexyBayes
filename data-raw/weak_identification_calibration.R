# data-raw/weak_identification_calibration.R
#
# The interval-to-median ratio, (q97.5 - q2.5) / median on the standard
# deviation, for the ten-seed designs that calibrate the weak-identification
# flag (R/weak_identification.R). Each design is `y ~ 1 + (1 | g)` with group
# and residual SD 1 unless stated, fitted on INLA under the package defaults
# (best of 3, one thread, default half-normal prior), seeds 1 to 10. The
# random copy of a fixed factor is the design of the DEF-004 test.
#
# Run from the package root with INLA installed:
#
#   Rscript --no-init-file data-raw/weak_identification_calibration.R
#
# Writes inst/extdata/weak_identification_calibration.csv, which
# tests/testthat/test-weak-variance-flag.R reads to check the threshold.

devtools::load_all(".", quiet = TRUE)
options(
  flexyBayes.inla_threads = "1:1",
  flexyBayes.silence_default_prior_note = TRUE,
  flexyBayes.silence_weak_identification_warning = TRUE,
  flexyBayes.silence_boundary_collapse_warning = TRUE
)

designs <- data.frame(
  design = c("3x10", "4x10", "5x10", "8x10", "30x8", "5x6_sd0", "20x8_sd0",
             "copy_of_fixed_6"),
  levels = c(3L, 4L, 5L, 8L, 30L, 5L, 20L, 6L),
  rows = c(10L, 10L, 10L, 10L, 8L, 6L, 8L, 10L),
  sd_g = c(1, 1, 1, 1, 1, 0, 0, 1),
  stringsAsFactors = FALSE
)

one_fit <- function(design, levels, rows, sd_g, seed) {
  set.seed(seed)
  if (identical(design, "copy_of_fixed_6")) {
    d <- data.frame(blk = factor(rep(seq_len(levels), each = rows)))
    d$plot <- factor(paste0("p", d$blk))
    d$y <- stats::rnorm(levels)[d$blk] + stats::rnorm(nrow(d))
    fit <- suppressWarnings(flexybayes(
      y ~ blk, random = ~ plot, data = d, backend = "inla", verbose = FALSE
    ))
    component <- "sd_plot"
  } else {
    d <- data.frame(g = factor(rep(seq_len(levels), each = rows)))
    d$y <- stats::rnorm(levels, sd = sd_g)[d$g] + stats::rnorm(nrow(d))
    fit <- suppressWarnings(fb(
      y ~ 1 + (1 | g), data = d, backend = "inla", verbose = FALSE
    ))
    component <- "sd_g"
  }
  vc <- fit$extras$variance_comps
  i <- match(component, vc$component)
  med <- attr(vc, "posterior_median")[[component]]
  share <- attr(vc, "prior_share")[component]
  data.frame(
    design = design,
    levels = levels,
    rows = rows,
    sd_g = sd_g,
    seed = seed,
    component = component,
    q2.5 = vc$q2.5[[i]],
    median = med,
    q97.5 = vc$q97.5[[i]],
    ratio = (vc$q97.5[[i]] - vc$q2.5[[i]]) / med,
    prior_share = unname(share),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (i1 in seq_len(nrow(designs))) {
  for (seed in seq_len(10L)) {
    rows[[length(rows) + 1L]] <- one_fit(
      designs$design[[i1]], designs$levels[[i1]], designs$rows[[i1]],
      designs$sd_g[[i1]], seed
    )
    message(designs$design[[i1]], " seed ", seed, " ratio ",
            signif(rows[[length(rows)]]$ratio, 4))
  }
}
out <- do.call(rbind, rows)
utils::write.csv(out, "inst/extdata/weak_identification_calibration.csv",
                 row.names = FALSE)

# The threshold: the 99th percentile over the identified designs with five or
# more levels. The random copy of a fixed factor is unidentified by
# construction and is a positive case, not a calibration row.
calib <- out$levels >= 5L & out$design != "copy_of_fixed_6"
message("threshold (q99, type 7) = ",
        format(stats::quantile(out$ratio[calib], 0.99, names = FALSE),
               digits = 6))
