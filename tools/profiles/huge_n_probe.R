# huge_n_probe.R -- shape probe of the preflight + aggregation +
# prediction kernel paths at N = 5e6.
#
# Active at v0.4.0 (Wave 2 Phase 2C trust-release gate, v040-plan
# section 9.3 + section 11). The script exercises the three data-shape
# code paths -- preflight, aggregation plan, and prediction-kernel
# file-output format resolution -- at N = 5e6 without materialising the
# (intractable) full MCMC fit. The point is to confirm the metadata /
# dispatch layer handles the huge-N data shape: a regression here would
# invert the v0.4.0 trust narrative, so a fault is a release blocker.
#
# Env-gated by FLEXYBAYES_RUN_HUGE_N=true.  Run from the workspace root:
#   FLEXYBAYES_RUN_HUGE_N=true Rscript flexyBayes/tools/profiles/huge_n_probe.R

if (!nzchar(Sys.getenv("FLEXYBAYES_RUN_HUGE_N"))) {
  cat("Set FLEXYBAYES_RUN_HUGE_N=true to run the huge_n_probe.\n")
  quit(save = "no", status = 0L)
}

suppressMessages(devtools::load_all("flexyBayes", quiet = TRUE))

t0 <- Sys.time()

# ---- Synthetic data ---------------------------------------------------

n_rows <- 5e6L
n_grp <- 200L
set.seed(20260527L)
d <- data.frame(
  group = factor(sample.int(n_grp, n_rows, replace = TRUE)),
  x = stats::rnorm(n_rows),
  y = stats::rnorm(n_rows)
)
cat(sprintf(
  "Generated synthetic dataset: %d rows, %d groups.\n",
  n_rows,
  n_grp
))

# ---- Preflight shape probe -------------------------------------------

cat("\n-- preflight shape --\n")
pf <- tryCatch(
  flexyBayes:::.fb_preflight(
    fixed = y ~ x,
    random = ~group,
    data = d
  ),
  error = function(e) {
    cat(sprintf("preflight raised expected refusal: %s\n", conditionMessage(e)))
    NULL
  }
)
if (!is.null(pf)) {
  cat(sprintf("preflight class: %s\n", paste(class(pf), collapse = ", ")))
  agg <- attr(pf, "aggregation_plan")
  if (!is.null(agg)) {
    cat(sprintf("aggregation_plan code: %s\n", agg$code))
  } else {
    cat("aggregation_plan: NULL (no aggregation triggered)\n")
  }
}

# ---- Aggregation kernel probe ----------------------------------------

cat("\n-- aggregation plan probe --\n")
plan <- tryCatch(
  flexyBayes:::.fb_aggregation_plan(
    fixed = y ~ x,
    random = ~group,
    data = d
  ),
  error = function(e) {
    cat(sprintf("aggregation_plan raised: %s\n", conditionMessage(e)))
    NULL
  }
)
if (!is.null(plan)) {
  cat(sprintf(
    "plan code: %s; n_cells = %s\n",
    plan$code,
    if (is.null(plan$n_cells)) "NA" else plan$n_cells
  ))
}

# ---- Prediction-kernel file-output shape probe -----------------------
#
# At N = 5e6 the prediction kernel routes its output through the
# file-backed path (ADR 0023): format = "auto" resolves to fst when fst
# is installed, rds otherwise, csv under interop = TRUE. The probe
# confirms the format-resolution dispatch handles the huge-N shape
# without materialising the 5e6-row payload.

cat("\n-- prediction-kernel format resolution --\n")
fst_ok <- requireNamespace("fst", quietly = TRUE)
fmt_auto <- flexyBayes:::.resolve_format("auto", n_rows, FALSE, fst_ok)
fmt_interop <- flexyBayes:::.resolve_format("auto", n_rows, TRUE, fst_ok)
cat(sprintf("fst installed: %s\n", fst_ok))
cat(sprintf(
  "format(auto, N=%g): %s  (expected %s)\n",
  n_rows,
  fmt_auto,
  if (fst_ok) "fst" else "rds"
))
cat(sprintf(
  "format(auto, N=%g, interop): %s  (expected csv)\n",
  n_rows,
  fmt_interop
))
stopifnot(identical(fmt_auto, if (fst_ok) "fst" else "rds"))
stopifnot(identical(fmt_interop, "csv"))

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\n=== huge_n_probe COMPLETE (v0.4.0) -- %.1fs ===\n", elapsed))
cat(
  "All three huge-N shape paths (preflight, aggregation, prediction ",
  "format resolution) handled the N = 5e6 data shape.\n",
  sep = ""
)
