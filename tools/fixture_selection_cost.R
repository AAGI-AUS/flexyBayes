#!/usr/bin/env Rscript

# fixture_selection_cost.R -- the package-owned synthetic fixture the
# selection-cost measurements run on.
#
# The shapes are the ones recorded for the selection module's pricing
# pilot: 10,890 rows at the first shape and 73,341 rows at the second, 82
# candidate columns, 22 clusters, 3 folds. Nothing here reads a study
# table; every column is drawn from a seeded generator, so a cost figure
# measured on it can be reproduced from this file alone.
#
# Structure, and why it is the one the pilot priced. The base model
# carries a small biometric fixed part and three grouped random terms
# (cluster, year, trial), which is the site, year and site-by-year
# decomposition the pilot's base used. The candidate columns are constant
# within a trial, as an environmental covariate block is: that is what
# puts more candidate columns than trial degrees of freedom into the fixed
# part, and it is the condition under which a training fold leaves part of
# the candidate block unidentified.
#
# The random structure is held fixed across both row counts so that the
# row count is the only axis that moves between the two shapes. A cell
# that also grew its number of trials would change the latent field at the
# same time and the two effects could not be told apart.
#
# Sourced by tools/benchmark_selection_cost.R. Defines functions and
# constants only; it runs nothing on being sourced.

FB_COST_N_CANDIDATES <- 82L
FB_COST_N_CLUSTERS <- 22L
FB_COST_N_TRIALS <- 48L
FB_COST_N_YEARS <- 5L
FB_COST_N_VARIETIES <- 4L
FB_COST_N_FOLDS <- 3L
FB_COST_SEED <- 20260916L

# The pilot's two recorded row counts.
FB_COST_ROWS_BARLEY <- 10890L
FB_COST_ROWS_DEPOSIT <- 73341L


# fb_cost_candidate_names() -- the declared candidate block, in order.
# The order is the declared priority the candidate-space contract reads
# when a fold cannot identify every column, so it is fixed here rather
# than left to the caller.
fb_cost_candidate_names <- function(n_candidates = FB_COST_N_CANDIDATES) {
  sprintf("z%02d", seq_len(n_candidates))
}


# fb_cost_candidate_formula() -- the same block as the one-sided formula
# fb_select() takes.
fb_cost_candidate_formula <- function(n_candidates = FB_COST_N_CANDIDATES) {
  stats::as.formula(
    paste("~", paste(fb_cost_candidate_names(n_candidates), collapse = " + "))
  )
}


# fb_cost_base_fixed() / fb_cost_base_random() -- the base model. The
# fixed part is biometric and stays protected through the selection; the
# candidates are the block above.
fb_cost_base_fixed <- function() {
  y ~ variety + napplied + seedrate
}

fb_cost_base_random <- function() {
  ~ cluster + year + trial
}


# fb_cost_fixture() -- the data frame.
#
# `n_rows` rows are dealt round-robin over `n_trials` trials, each trial
# sits in one cluster, and each candidate column takes one value per
# trial. The response carries two protected signals, one candidate signal
# and the three random effects, so the model is identified and the
# selection has something to find; none of that affects cost, which is
# what this fixture exists to measure.
fb_cost_fixture <- function(n_rows,
                            n_candidates = FB_COST_N_CANDIDATES,
                            n_clusters = FB_COST_N_CLUSTERS,
                            n_trials = FB_COST_N_TRIALS,
                            n_years = FB_COST_N_YEARS,
                            seed = FB_COST_SEED) {
  stopifnot(n_rows >= n_trials, n_trials >= n_clusters)
  set.seed(seed)

  trial_of_row <- rep_len(seq_len(n_trials), n_rows)
  cluster_of_trial <- rep_len(seq_len(n_clusters), n_trials)

  tab_d <- data.frame(
    cluster = factor(sprintf("C%02d", cluster_of_trial[trial_of_row])),
    trial = factor(sprintf("T%03d", trial_of_row)),
    year = factor(sprintf("Y%d", rep_len(seq_len(n_years), n_rows))),
    variety = factor(sprintf("V%d", rep_len(seq_len(FB_COST_N_VARIETIES), n_rows))),
    napplied = stats::rnorm(n_rows),
    seedrate = stats::rnorm(n_rows)
  )

  # Trial-constant candidates: one draw per trial, spread over the rows of
  # that trial.
  cand_vec <- fb_cost_candidate_names(n_candidates)
  for (j1 in seq_along(cand_vec)) {
    tab_d[[cand_vec[[j1]]]] <- stats::rnorm(n_trials)[trial_of_row]
  }

  eff_trial <- stats::rnorm(n_trials, sd = 0.5)
  eff_cluster <- stats::rnorm(n_clusters, sd = 0.4)
  eff_year <- stats::rnorm(n_years, sd = 0.2)

  tab_d$y <- 1 +
    0.4 * tab_d$napplied +
    0.3 * tab_d$seedrate +
    0.5 * tab_d[[cand_vec[[1L]]]] +
    eff_trial[trial_of_row] +
    eff_cluster[as.integer(tab_d$cluster)] +
    eff_year[as.integer(tab_d$year)] +
    stats::rnorm(n_rows, sd = 0.6)

  tab_d
}
