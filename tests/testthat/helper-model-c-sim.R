# helper-model-c-sim.R -- the model C simulation of data-raw/model_c_reference.R.
#
# Section 1 of that script as a function of the seed. At seed 20261008 it
# reproduces the stored fixture inst/extdata/model_c_reference.rds; the
# recovery study in test-model-c-pooled.R draws its other fixtures from it.

.mc_simulate <- function(seed) {
  withr::with_seed(
    seed,
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection",
    {
      n_trial <- 12L
      n_variety <- 30L
      n_rep <- 3L
      trial_lev <- sprintf("T%02d", seq_len(n_trial))
      variety_lev <- sprintf("V%02d", seq_len(n_variety))
      sd_gxe <- exp(stats::rnorm(n_trial, log(0.35), 0.4))
      sd_res <- exp(stats::rnorm(n_trial, log(0.50), 0.35))
      names(sd_gxe) <- trial_lev
      names(sd_res) <- trial_lev
      d <- expand.grid(
        VARIETY = variety_lev,
        REP = sprintf("R%d", seq_len(n_rep)),
        TRIAL = trial_lev,
        stringsAsFactors = FALSE
      )
      d <- d[, c("TRIAL", "REP", "VARIETY")]
      d$TRIAL <- factor(d$TRIAL, levels = trial_lev)
      d$REP <- factor(d$REP)
      d$VARIETY <- factor(d$VARIETY, levels = variety_lev)
      t_i <- as.integer(d$TRIAL)
      v_i <- as.integer(d$VARIETY)
      eff_variety <- stats::rnorm(n_variety, 0, 0.5)
      eff_trial <- stats::rnorm(n_trial, 0, 1.0)
      eff_rep <- matrix(stats::rnorm(n_trial * n_rep, 0, 0.3), n_trial, n_rep)
      eff_gxe <- matrix(
        stats::rnorm(n_trial * n_variety, 0, rep(sd_gxe, n_variety)),
        n_trial,
        n_variety
      )
      d$y <- 5 +
        eff_variety[v_i] +
        eff_trial[t_i] +
        eff_rep[cbind(t_i, as.integer(d$REP))] +
        eff_gxe[cbind(t_i, v_i)] +
        stats::rnorm(nrow(d), 0, sd_res[t_i])
      d <- d[order(d$TRIAL, d$REP, d$VARIETY), ]
      rownames(d) <- NULL
      # The SDs of the log per-trial SDs are the generator's 0.4 and 0.35.
      list(
        data = d,
        truth = list(
          sd_gxe = sd_gxe,
          sd_res = sd_res,
          pool_gxe = 0.4,
          pool_res = 0.35
        )
      )
    }
  )
}
