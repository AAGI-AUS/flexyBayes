# data-raw/at_ar1_reference.R
#
# Simulated field trials fitted once in ASReml and stored, so the tests can
# check flexyBayes against ASReml without an ASReml licence. Four fits:
#
#   one_trial     residual = ~ ar1(ROW):ar1(RANGE), one trial, no nugget
#   at_field      random = ~ at(TRIAL):ar1(ROW):ar1(RANGE), residual units,
#                 three trials with their own field each
#   at_field_lev  random = ~ at(TRIAL, "T2"):ar1(ROW):ar1(RANGE), residual
#                 units, three trials of which only T2 carries a field
#   dsum_field    residual = ~ dsum(~ ar1(ROW):ar1(RANGE) | TRIAL), three
#                 trials, no nugget
#   at_level      random = ~ TRIAL:REP + at(DESIGN, "split"):MAINPLOT, six
#                 trials, main plots in the three split-plot trials only
#
# Run from the package root with a licensed ASReml on the library path:
#
#   FLEXYBAYES_ASREML_LIB=/path/to/library \
#     Rscript --vanilla data-raw/at_ar1_reference.R
#
# Unset, the default library is used. Every number the tests attribute to
# ASReml comes from this file.
asreml_lib <- Sys.getenv("FLEXYBAYES_ASREML_LIB", unset = "")
if (nzchar(asreml_lib)) {
  if (!dir.exists(asreml_lib)) {
    stop(
      "FLEXYBAYES_ASREML_LIB is set to a directory that does not exist: ",
      asreml_lib,
      call. = FALSE
    )
  }
  .libPaths(c(asreml_lib, .libPaths()))
}
suppressMessages(library(asreml))

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261008L)

n_row <- 16L
n_range <- 10L
n_variety <- 40L

# --- 1. simulators -------------------------------------------------------- #

ar1_chol <- function(n, rho) {
  t(chol(rho^abs(outer(seq_len(n), seq_len(n), "-"))))
}

# A separable AR1 field on the n_row x n_range grid, laid out as the rows of
# expand.grid(ROW, RANGE): ROW varies fastest.
sim_field <- function(sd, rho_row, rho_range) {
  z <- matrix(stats::rnorm(n_row * n_range), n_row, n_range)
  f <- ar1_chol(n_row, rho_row) %*% z %*% t(ar1_chol(n_range, rho_range))
  sd * as.vector(f)
}

# One trial: four complete replicates of 40 varieties, each replicate a
# strip of four rows across the ten ranges.
sim_trial <- function(trial, sd, rho_row, rho_range, sd_e, var_eff) {
  d <- expand.grid(
    ROW = seq_len(n_row),
    RANGE = seq_len(n_range),
    KEEP.OUT.ATTRS = FALSE
  )
  d$REP <- paste0("R", (d$ROW - 1L) %/% 4L + 1L)
  d$VARIETY <- NA_character_
  for (r in unique(d$REP)) {
    i <- which(d$REP == r)
    d$VARIETY[i] <- sample(sprintf("V%02d", seq_len(n_variety)))
  }
  d$TRIAL <- trial
  d$y <- 5 +
    var_eff[d$VARIETY] +
    sim_field(sd, rho_row, rho_range) +
    stats::rnorm(nrow(d), 0, sd_e)
  d
}

as_design <- function(d) {
  d$ROW <- factor(d$ROW)
  d$RANGE <- factor(d$RANGE)
  d$TRIAL <- factor(d$TRIAL)
  d$REP <- factor(d$REP)
  d$VARIETY <- factor(d$VARIETY)
  d[order(d$TRIAL, d$ROW, d$RANGE), ]
}

var_eff <- stats::setNames(
  stats::rnorm(n_variety, 0, 0.5),
  sprintf("V%02d", seq_len(n_variety))
)

truth <- list(
  one_trial = c(sd = 0.8, rho_row = 0.6, rho_range = 0.4),
  at_field = rbind(
    T1 = c(sd = 0.8, rho_row = 0.7, rho_range = 0.3, sd_e = 0.3),
    T2 = c(sd = 0.5, rho_row = 0.4, rho_range = 0.6, sd_e = 0.3),
    T3 = c(sd = 1.0, rho_row = 0.8, rho_range = 0.5, sd_e = 0.3)
  ),
  at_field_lev = rbind(
    T1 = c(sd = 0, rho_row = 0, rho_range = 0, sd_e = 0.5),
    T2 = c(sd = 0.7, rho_row = 0.6, rho_range = 0.5, sd_e = 0.5),
    T3 = c(sd = 0, rho_row = 0, rho_range = 0, sd_e = 0.5)
  ),
  dsum_field = rbind(
    T1 = c(sd = 0.8, rho_row = 0.7, rho_range = 0.3),
    T2 = c(sd = 0.5, rho_row = 0.4, rho_range = 0.6),
    T3 = c(sd = 1.0, rho_row = 0.8, rho_range = 0.5)
  )
)

one_trial <- as_design(with(
  as.list(truth$one_trial),
  sim_trial("T1", sd, rho_row, rho_range, 0, var_eff)
))

at_data <- as_design(do.call(
  rbind,
  lapply(rownames(truth$at_field), function(k) {
    p <- truth$at_field[k, ]
    sim_trial(
      k,
      p[["sd"]],
      p[["rho_row"]],
      p[["rho_range"]],
      p[["sd_e"]],
      var_eff
    )
  })
))

at_lev_data <- as_design(do.call(
  rbind,
  lapply(rownames(truth$at_field_lev), function(k) {
    p <- truth$at_field_lev[k, ]
    sim_trial(
      k,
      p[["sd"]],
      p[["rho_row"]],
      p[["rho_range"]],
      p[["sd_e"]],
      var_eff
    )
  })
))

dsum_data <- as_design(do.call(
  rbind,
  lapply(rownames(truth$dsum_field), function(k) {
    p <- truth$dsum_field[k, ]
    sim_trial(k, p[["sd"]], p[["rho_row"]], p[["rho_range"]], 0, var_eff)
  })
))

# Six trials, three of them split plots: within each replicate, ten main
# plots of four subplots. Main plots carry their own effect only in the
# split-plot trials.
truth$at_level <- c(sd_rep = 0.4, sd_mainplot = 0.5, sd_e = 0.4)
level_data <- do.call(
  rbind,
  lapply(seq_len(6L), function(k) {
    d <- expand.grid(
      PLOT = seq_len(40L),
      REP = paste0("R", seq_len(3L)),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    d$TRIAL <- sprintf("T%d", k)
    d$DESIGN <- if (k <= 3L) "split" else "rcb"
    d$MAINPLOT <- paste0(d$TRIAL, "_", d$REP, "_M", (d$PLOT - 1L) %/% 4L + 1L)
    d$VARIETY <- NA_character_
    for (r in unique(d$REP)) {
      i <- which(d$REP == r)
      d$VARIETY[i] <- sample(sprintf("V%02d", seq_len(n_variety)))
    }
    rep_eff <- stats::rnorm(3L, 0, truth$at_level[["sd_rep"]])
    names(rep_eff) <- paste0("R", seq_len(3L))
    mp <- unique(d$MAINPLOT)
    mp_eff <- if (k <= 3L) {
      stats::setNames(
        stats::rnorm(length(mp), 0, truth$at_level[["sd_mainplot"]]),
        mp
      )
    } else {
      stats::setNames(numeric(length(mp)), mp)
    }
    d$y <- 5 +
      k / 4 +
      var_eff[d$VARIETY] +
      rep_eff[d$REP] +
      mp_eff[d$MAINPLOT] +
      stats::rnorm(nrow(d), 0, truth$at_level[["sd_e"]])
    d
  })
)
level_data$PLOT <- NULL
for (v in c("TRIAL", "DESIGN", "REP", "MAINPLOT", "VARIETY")) {
  level_data[[v]] <- factor(level_data[[v]])
}

# --- 2. the ASReml fits --------------------------------------------------- #

varcomp_table <- function(m) {
  vc <- summary(m)$varcomp
  data.frame(
    term = rownames(vc),
    component = vc$component,
    std.error = vc$std.error,
    bound = as.character(vc$bound),
    stringsAsFactors = FALSE
  )
}

fit <- function(...) {
  m <- asreml(..., trace = FALSE, maxit = 60)
  stopifnot(isTRUE(m$converge))
  list(varcomp = varcomp_table(m), loglik = m$loglik)
}

fits <- list(
  one_trial = fit(
    fixed = y ~ VARIETY,
    residual = ~ ar1(ROW):ar1(RANGE),
    data = one_trial
  ),
  at_field = fit(
    fixed = y ~ TRIAL + VARIETY,
    random = ~ at(TRIAL):ar1(ROW):ar1(RANGE),
    residual = ~units,
    data = at_data
  ),
  at_field_lev = fit(
    fixed = y ~ TRIAL + VARIETY,
    random = ~ at(TRIAL, "T2"):ar1(ROW):ar1(RANGE),
    residual = ~units,
    data = at_lev_data
  ),
  dsum_field = fit(
    fixed = y ~ TRIAL + VARIETY,
    residual = ~ dsum(~ ar1(ROW):ar1(RANGE) | TRIAL),
    data = dsum_data
  ),
  at_level = fit(
    fixed = y ~ TRIAL + VARIETY,
    random = ~ TRIAL:REP + at(DESIGN, "split"):MAINPLOT,
    residual = ~units,
    data = level_data
  )
)

at_ar1_reference <- list(
  stamp = list(
    date = as.character(Sys.Date()),
    asreml = as.character(utils::packageVersion("asreml")),
    r = R.version.string,
    script = "data-raw/at_ar1_reference.R"
  ),
  truth = truth,
  data = list(
    one_trial = one_trial,
    at_field = at_data,
    at_field_lev = at_lev_data,
    dsum_field = dsum_data,
    at_level = level_data
  ),
  fits = fits
)

saveRDS(at_ar1_reference, "inst/extdata/at_ar1_reference.rds")
cat("written inst/extdata/at_ar1_reference.rds\n")
for (nm in names(fits)) {
  cat("\n==", nm, "\n")
  print(fits[[nm]]$varcomp, row.names = FALSE)
}
