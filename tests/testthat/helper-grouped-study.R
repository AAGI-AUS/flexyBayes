# A synthetic multi-trial design shaped like a grouped covariate selection:
# trials of about 70 plots, a random trial effect, one protected numeric
# column, and declared blocks of strongly correlated candidates that vary
# mostly between trials. `sizes` gives the columns per block and `signal`
# the non-zero candidate effects.
fbgs_study <- function(n_trial = 48L,
                       sizes = c(2L, 3L, 4L, 5L, 6L, 8L, 4L, 3L),
                       signal = c(g1_1 = 0.6, g3_2 = -0.4, g6_5 = 0.3),
                       seed = 240924L) {
  .fb_sel_with_seed(seed, {
    per <- sample(60:82, n_trial, replace = TRUE)
    trial <- factor(rep(sprintf("t%02d", seq_len(n_trial)), times = per))
    n <- length(trial)
    index <- as.integer(trial)
    d <- data.frame(trial = trial, x = stats::rnorm(n))
    for (g in seq_along(sizes)) {
      # Between-trial variation dominates, as for site-level covariates, and
      # members of a block correlate at about 0.97 with each other.
      latent <- stats::rnorm(n_trial)[index] + stats::rnorm(n, sd = 0.3)
      for (j in seq_len(sizes[[g]])) {
        d[[sprintf("g%d_%d", g, j)]] <- latent + stats::rnorm(n, sd = 0.18)
      }
    }
    d$y <- 2 +
      0.5 * d$x +
      drop(as.matrix(d[names(signal)]) %*% signal) +
      stats::rnorm(n_trial, sd = 1)[index] +
      stats::rnorm(n, sd = 1)
    d
  })
}


# One `group = "block"` specification per declared block, at one scale.
fbgs_penalty <- function(sizes, family_call) {
  specs <- lapply(seq_along(sizes), function(g) {
    lhs <- call("term", call("matches", sprintf("^g%d_", g)))
    stats::as.formula(call("~", lhs, family_call))
  })
  do.call(fb_penalty, specs)
}


fbgs_candidates <- function(d) {
  columns <- grep("^g[0-9]+_", names(d), value = TRUE)
  stats::as.formula(paste("~", paste(columns, collapse = " + ")))
}


fbgs_base <- function(d) {
  fb <- fb_from_asreml(y ~ x, random = ~trial, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = y ~ x, random = ~trial, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}
