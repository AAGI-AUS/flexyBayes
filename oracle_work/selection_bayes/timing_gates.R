suppressPackageStartupMessages(library(brms)); options(brms.backend = "cmdstanr")
run <- function(label, n, p, n_groups, p0 = 10, chains = 2, iter = 2000) {
  set.seed(11); g <- factor(sample(seq_len(n_groups), n, replace = TRUE))
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("c", seq_len(p))
  beta <- numeric(p); beta[seq_len(p0)] <- rnorm(p0, 0, 1)
  y <- 1 + X %*% beta + rnorm(n_groups, 0, 0.7)[g] + rnorm(n, 0, 1)
  d <- data.frame(y = as.numeric(y), X, g)
  f <- as.formula(paste("y ~", paste(colnames(X), collapse = " + "), "+ (1|g)"))
  t0 <- Sys.time()
  m <- brm(f, d, prior = prior(horseshoe(par_ratio = p0 / (p - p0)), class = "b"),
           chains = chains, cores = chains, iter = iter, refresh = 0, seed = 1, silent = 2,
           control = list(adapt_delta = 0.95))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  np <- nuts_params(m); div <- sum(np$Value[np$Parameter == "divergent__"])
  td <- max(np$Value[np$Parameter == "treedepth__"])
  s <- posterior::summarise_draws(as_draws_array(m, variable = "^b_", regex = TRUE))
  cat(sprintf("\n=== %s: n=%d p=%d groups=%d | wall %.0f s | div %d | max treedepth %d | max rhat %.3f | min bulk ESS %.0f | min tail ESS %.0f\n",
      label, n, p, n_groups, el, div, td, max(s$rhat, na.rm=TRUE), min(s$ess_bulk, na.rm=TRUE), min(s$ess_tail, na.rm=TRUE)))
  cat("recovered top-", p0, " |b| are true nonzero:", mean(order(-abs(s$mean))[1:p0] %in% seq_len(p0)), "\n")
}
run("F1 barley shape", n = 48, p = 176, n_groups = 48)
run("F2 deposit shape (reduced trial level)", n = 3427, p = 758, n_groups = 400)
