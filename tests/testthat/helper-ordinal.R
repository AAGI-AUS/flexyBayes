# helper-ordinal.R -- simulated ordinal field trials for the ordinal tests.
#
# The layout is that of Cao, Megirian, Nguyen and Sparks (2026, arXiv
# 2610.37384): a 5 row by 18 column array, three replicates of six columns,
# sowing depth applied to strips of two columns within each replicate (the
# main plots), ten sub-plot treatments within each main plot, and a number of
# plants scored on an ordinal scale in every plot.

# Category probabilities, written out independently of the package: the
# cumulative family gives P(Y <= k) = F(tau_k - eta), the continuation ratio
# P(Y > k | Y >= k) = F(eta - tau_k).
.ord_probs_reference <- function(eta, tau, family, inv_link = stats::plogis) {
  K <- length(tau) + 1L
  n <- length(eta)
  if (identical(family, "cumulative")) {
    cum <- cbind(inv_link(outer(-eta, tau, "+")), 1)
    return(cum - cbind(0, cum[, -K, drop = FALSE]))
  }
  continue <- matrix(inv_link(outer(eta, tau, "-")), n, K - 1L)
  reach <- matrix(1, n, K)
  for (k in seq_len(K - 1L)) {
    reach[, k + 1L] <- reach[, k] * continue[, k]
  }
  reach * cbind(1 - continue, 1)
}

# One simulated trial. `truth` holds the depth contrasts (d80, d120), the
# standard deviations of the replicate, main-plot, field and nugget terms,
# the two field correlations and the thresholds `tau`.
.ord_sim_trial <- function(
  seed,
  family = "cumulative",
  truth = .ord_truth(family),
  n_plant = 10L,
  field = TRUE
) {
  set.seed(seed)
  n_row <- 5L
  n_col <- 18L
  plots <- expand.grid(row = seq_len(n_row), col = seq_len(n_col))
  plots$rep <- (plots$col - 1L) %/% 6L + 1L
  strip <- (plots$col - 1L) %/% 2L + 1L
  depth_of_strip <- unlist(lapply(seq_len(3L), function(r) {
    sample(c(30L, 80L, 120L))
  }))
  plots$depth <- factor(depth_of_strip[strip], levels = c(30L, 80L, 120L))
  plots$mp <- interaction(plots$rep, plots$depth, drop = TRUE)
  plots$plot <- factor(seq_len(nrow(plots)))
  plots$vh <- factor(stats::ave(
    seq_len(nrow(plots)),
    plots$mp,
    FUN = function(i) sample(seq_along(i))
  ))
  ar1_chol <- function(n, rho) {
    t(chol(rho^abs(outer(seq_len(n), seq_len(n), "-"))))
  }
  z <- matrix(stats::rnorm(n_row * n_col), n_row, n_col)
  s <- truth$s_field * ar1_chol(n_row, truth$rho_row) %*% z %*%
    t(ar1_chol(n_col, truth$rho_col))
  eta <- truth$d80 * (plots$depth == "80") +
    truth$d120 * (plots$depth == "120") +
    stats::rnorm(3L, 0, truth$s_rep)[plots$rep] +
    stats::rnorm(nlevels(plots$mp), 0, truth$s_mp)[as.integer(plots$mp)] +
    stats::rnorm(nrow(plots), 0, truth$s_nugget)
  if (field) {
    eta <- eta + s[cbind(plots$row, plots$col)]
  }
  d <- plots[rep(seq_len(nrow(plots)), each = n_plant), ]
  rownames(d) <- NULL
  p <- .ord_probs_reference(eta[as.integer(d$plot)], truth$tau, family)
  u <- stats::runif(nrow(d))
  y <- 1L + rowSums(u > t(apply(p, 1L, cumsum))[, -ncol(p), drop = FALSE])
  d$stage <- factor(
    y,
    levels = seq_len(length(truth$tau) + 1L),
    labels = paste0("Z", seq_len(length(truth$tau) + 1L)),
    ordered = TRUE
  )
  d$row <- factor(d$row)
  d$col <- factor(d$col)
  d$rep <- factor(d$rep)
  d
}

# The sealed recovery truth of PLAN_v0.1 (thresholds in brms's convention).
.ord_truth <- function(family) {
  list(
    d80 = 0.8,
    d120 = 1.6,
    s_rep = 0.3,
    s_mp = 0.35,
    s_field = 0.7,
    s_nugget = 0.5,
    rho_row = 0.3,
    rho_col = 0.6,
    tau = if (identical(family, "cumulative")) {
      c(-2.0, -0.8, 0.2, 1.2, 2.3)
    } else {
      c(-1.5, -0.5, 0.3, 1.0, 1.8)
    }
  )
}
