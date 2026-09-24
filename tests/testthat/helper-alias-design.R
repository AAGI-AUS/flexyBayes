# helper-alias-design.R -- a synthetic multi-trial barley design with empty
# cells of VARIETY x NAPPLIED x SEEDRATE, for the aliased-protected-column
# tests. Every value is simulated.
#
# `shape = "study"`: 48 trials at 22 localities; N at 0/20/80 kg/ha in 8
# trials, 0/30/90 in 33 and 0/30/90/150 in 7; seed rate 75/150/300 and 3 reps
# everywhere; 8 of 22 varieties per trial (6 at one trial, 10 at another).
# `shape = "small"`: 12 trials at 6 localities, 8 varieties, 4 per trial, the
# same three N designs and seed rates, 2 reps.
#
# Varieties and N rates never sown together leave most empty cells. Three
# cells are lost at every trial that sowed them, so some columns are aliased
# without being zero. The recorded asreml 4.2 fits of both shapes are in
# fixtures/asreml_alias_<shape>.rds, written by the evidence script
# review/units/U-ALIAS_2026-09-24/asreml_alias_oracle.R.
alias_design <- function(shape = c("study", "small"), seed = 20260924L) {
  shape <- match.arg(shape)
  set.seed(seed)
  if (identical(shape, "study")) {
    n_trial <- 48L
    per_locality <- c(rep(4L, 3L), rep(3L, 6L), rep(2L, 5L), rep(1L, 8L))
    design_a <- c(3L, 8L, 14L, 19L, 25L, 31L, 38L, 45L)
    design_c <- c(5L, 11L, 22L, 28L, 34L, 41L, 47L)
    common <- sprintf("V%02d", 1:5)
    pool <- list(
      A = sprintf("V%02d", 6:16),
      B = sprintf("V%02d", 6:22),
      C = sprintf("V%02d", 12:22)
    )
    n_var <- rep(8L, n_trial)
    n_var[[17L]] <- 6L
    n_var[[36L]] <- 10L
    k_common <- 4:5
    reps <- 3L
    lost <- list(c("V13", 0, 300), c("V07", 30, 75), c("V19", 90, 150))
  } else {
    n_trial <- 12L
    per_locality <- rep(2L, 6L)
    design_a <- c(1L, 6L, 11L)
    design_c <- c(4L, 7L, 10L)
    common <- sprintf("V%02d", 1:3)
    pool <- list(
      A = sprintf("V%02d", 4:6),
      B = sprintf("V%02d", 4:8),
      C = sprintf("V%02d", 6:8)
    )
    n_var <- rep(4L, n_trial)
    k_common <- 2L
    reps <- 2L
    lost <- list(c("V05", 0, 300), c("V04", 30, 75), c("V07", 90, 150))
  }
  locality <- rep(sprintf("L%02d", seq_along(per_locality)), per_locality)
  design <- rep("B", n_trial)
  design[design_a] <- "A"
  design[design_c] <- "C"
  n_levels <- list(A = c(0, 20, 80), B = c(0, 30, 90), C = c(0, 30, 90, 150))
  rows <- list()
  for (t in seq_len(n_trial)) {
    k <- if (length(k_common) > 1L) sample(k_common, 1L) else k_common
    varieties <- c(
      sample(common, k),
      sample(pool[[design[[t]]]], n_var[[t]] - k)
    )
    g <- expand.grid(
      REP = seq_len(reps),
      SEEDRATE = c(75, 150, 300),
      NAPPLIED = n_levels[[design[[t]]]],
      VARIETY = varieties,
      stringsAsFactors = FALSE
    )
    g$TRIAL <- sprintf("T%02d", t)
    g$LOCALITY <- locality[[t]]
    rows[[t]] <- g
  }
  d <- do.call(rbind, rows)
  for (cell in lost) {
    d <- d[
      !(d$VARIETY == cell[[1L]] &
        d$NAPPLIED == as.numeric(cell[[2L]]) &
        d$SEEDRATE == as.numeric(cell[[3L]])),
      ,
      drop = FALSE
    ]
  }
  d$TRIAL <- factor(d$TRIAL)
  d$LOCALITY <- factor(d$LOCALITY)
  d$REP <- factor(d$REP)
  d$VARIETY <- factor(d$VARIETY)
  d$NAPPLIED <- factor(d$NAPPLIED, levels = c(0, 20, 30, 80, 90, 150))
  d$SEEDRATE <- factor(d$SEEDRATE, levels = c(75, 150, 300))
  n <- nrow(d)
  trial <- as.integer(d$TRIAL)
  trial_rep <- as.integer(interaction(d$TRIAL, d$REP, drop = TRUE))
  variety_effect <- stats::rnorm(nlevels(d$VARIETY), sd = 0.3)
  n_num <- as.numeric(as.character(d$NAPPLIED))
  s_num <- as.numeric(as.character(d$SEEDRATE))
  d$x1 <- stats::rnorm(n)
  d$x2 <- stats::rnorm(n)
  d$x3 <- stats::rnorm(n)
  d$x4 <- stats::rnorm(n)
  d$YIELD <- 3 +
    variety_effect[as.integer(d$VARIETY)] +
    0.012 * n_num -
    0.00004 * n_num^2 +
    0.002 * s_num +
    0.00002 * n_num * s_num * (as.integer(d$VARIETY) %% 3L) +
    0.3 * d$x1 +
    stats::rnorm(n_trial, sd = 0.6)[trial] +
    stats::rnorm(max(trial_rep), sd = 0.2)[trial_rep] +
    stats::rnorm(n, sd = 0.35)
  rownames(d) <- NULL
  d
}
