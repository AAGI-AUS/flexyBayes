# F0b -- a hand-written regularised horseshoe over b[cand_idx], against the
# native class-wide horseshoe on the same data.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 4.4. The
# class-wide route cannot be aimed at a block, so a model with protected
# fixed effects has to give them up to the global-local prior or use a
# coefficient-row family. A stanvar-carried horseshoe density over an index
# vector would give both exact targeting and a linear formula that projpred
# accepts. This script decides whether it does.
#
# Exit rule, written before the run: the route ships only if divergences,
# R-hat, ESS and posterior means match the native class-wide horseshoe on
# the same simulated data within Monte-Carlo error. Otherwise it stays here
# as an artefact and is not wired into the compiler.
#
# Design: n = 100, p = 20, three true non-zero coefficients, one random
# intercept over 10 groups. Both arms see identical data and identical
# sampler settings.

suppressPackageStartupMessages(library(brms))
options(brms.backend = "cmdstanr", mc.cores = 4)

set.seed(20260912L)

n <- 100L
p <- 20L
p0 <- 3L
n_groups <- 10L

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", seq_len(p))
beta_true <- numeric(p)
beta_true[seq_len(p0)] <- c(2.0, -1.5, 1.0)
g <- factor(rep(seq_len(n_groups), length.out = n))
u <- rnorm(n_groups, 0, 0.8)
y <- as.numeric(1.0 + X %*% beta_true + u[g] + rnorm(n, 0, 1))
d <- data.frame(y = y, X, g = g)

fml <- as.formula(
  paste("y ~", paste(colnames(X), collapse = " + "), "+ (1 | g)")
)

par_ratio <- p0 / (p - p0)
scale_global <- par_ratio / sqrt(n)
df <- 1
df_global <- 1
scale_slab <- 2
df_slab <- 4

chains <- 4L
iter <- 2000L
adapt_delta <- 0.95
seed <- 4242L

report <- function(label, fit, elapsed) {
  np <- brms::nuts_params(fit)
  div <- sum(np$Value[np$Parameter == "divergent__"])
  td <- max(np$Value[np$Parameter == "treedepth__"])
  s <- posterior::summarise_draws(
    posterior::as_draws_array(fit, variable = "^b_x", regex = TRUE)
  )
  cat(sprintf(
    "\n=== %s | wall %.0f s | divergences %d | max treedepth %d | max rhat %.4f | min bulk ESS %.0f | min tail ESS %.0f\n",
    label, elapsed, div, td, max(s$rhat), min(s$ess_bulk), min(s$ess_tail)
  ))
  means <- stats::setNames(s$mean, sub("^b_", "", s$variable))
  means <- means[colnames(X)]
  cat("posterior means of the three true non-zero coefficients: ",
    paste(sprintf("%.3f", means[seq_len(p0)]), collapse = ", "),
    "  (truth ", paste(sprintf("%.3f", beta_true[seq_len(p0)]),
      collapse = ", "
    ), ")\n",
    sep = ""
  )
  cat(sprintf(
    "max |posterior mean| among the %d true zeros: %.3f\n",
    p - p0, max(abs(means[(p0 + 1L):p]))
  ))
  cat(sprintf(
    "top-%d by |posterior mean| are the true non-zeros: %s\n",
    p0, all(order(-abs(means))[seq_len(p0)] %in% seq_len(p0))
  ))
  list(
    divergences = div, max_rhat = max(s$rhat), min_bulk = min(s$ess_bulk),
    min_tail = min(s$ess_tail), means = means, sd = s$sd[match(
      colnames(X), sub("^b_", "", s$variable)
    )], elapsed = elapsed
  )
}


# --- arm A: the native class-wide horseshoe -------------------------------- #

prior_native <- brms::set_prior(
  sprintf(
    "horseshoe(df = %g, df_global = %g, scale_slab = %g, df_slab = %g, par_ratio = %.10g)",
    df, df_global, scale_slab, df_slab, par_ratio
  ),
  class = "b"
)

t0 <- Sys.time()
fit_native <- brms::brm(
  fml, d,
  prior = prior_native, chains = chains, cores = chains, iter = iter,
  control = list(adapt_delta = adapt_delta), seed = seed, refresh = 0,
  silent = 2
)
A <- report(
  "A native class-wide horseshoe", fit_native,
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
)


# --- arm B: the hand-written density over b[cand_idx] ----------------------- #
#
# Centred: b[cand_idx] is given a normal whose scale is the global-local
# product, rather than brms's non-centred b = zb .* sdb. The point of the
# route is that b keeps its own rows, so a protected coefficient can carry an
# ordinary prior beside the shrunk block -- and a centred b is what makes
# that possible.

cand_idx <- seq_len(p)

sv <- brms::stanvar(
  x = as.integer(cand_idx), name = "cand_idx", block = "data",
  scode = sprintf("array[%d] int cand_idx;", length(cand_idx))
) +
  brms::stanvar(
    x = length(cand_idx), name = "M_cand", block = "data",
    scode = "int<lower=1> M_cand;"
  ) +
  brms::stanvar(
    x = scale_global, name = "fbhs_scale_global", block = "data",
    scode = "real<lower=0> fbhs_scale_global;"
  ) +
  brms::stanvar(
    scode = paste(
      "real<lower=0> fbhs_global;",
      "real<lower=0> fbhs_slab;",
      "vector<lower=0>[M_cand] fbhs_local;",
      sep = "\n  "
    ),
    block = "parameters"
  ) +
  brms::stanvar(
    scode = paste(
      sprintf(
        "target += student_t_lpdf(fbhs_global | %g, 0, fbhs_scale_global * sigma) - log(0.5);",
        df_global
      ),
      sprintf(
        "target += inv_gamma_lpdf(fbhs_slab | %g, %g);",
        0.5 * df_slab, 0.5 * df_slab
      ),
      sprintf(
        "target += student_t_lpdf(fbhs_local | %g, 0, 1) - M_cand * log(0.5);",
        df
      ),
      "{",
      sprintf("  real fbhs_c2 = %g * fbhs_slab;", scale_slab^2),
      "  vector[M_cand] fbhs_l2 = square(fbhs_local);",
      "  vector[M_cand] fbhs_lt = sqrt(fbhs_c2 * fbhs_l2 ./ (fbhs_c2 + square(fbhs_global) * fbhs_l2));",
      "  target += normal_lpdf(b[cand_idx] | 0, fbhs_global * fbhs_lt);",
      "}",
      sep = "\n  "
    ),
    block = "model"
  )

cat("\n--- arm B generated program, shrinkage lines ---\n")
codeB <- brms::make_stancode(fml, d, stanvars = sv)
cat(
  grep(
    "fbhs|cand_idx|M_cand", strsplit(codeB, "\n", fixed = TRUE)[[1L]],
    value = TRUE
  ),
  sep = "\n"
)

t0 <- Sys.time()
fit_stanvar <- brms::brm(
  fml, d,
  stanvars = sv, chains = chains, cores = chains, iter = iter,
  control = list(adapt_delta = adapt_delta), seed = seed, refresh = 0,
  silent = 2
)
B <- report(
  "B stanvar centred horseshoe over b[cand_idx]", fit_stanvar,
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
)


# --- the comparison the exit rule is stated on ------------------------------ #
#
# Monte-Carlo error on a posterior mean is roughly sd / sqrt(ESS) per arm, so
# the two arms agree within Monte-Carlo error when the difference in means is
# small against the pooled figure. Two standard errors is the bar.

mcse <- sqrt(A$sd^2 / A$min_bulk + B$sd^2 / B$min_bulk)
delta <- abs(A$means - B$means)
z <- delta / mcse

cat("\n=== comparison ===\n")
cat(sprintf(
  "divergences: native %d, stanvar %d\n", A$divergences, B$divergences
))
cat(sprintf(
  "max rhat: native %.4f, stanvar %.4f\n", A$max_rhat, B$max_rhat
))
cat(sprintf(
  "min bulk ESS: native %.0f, stanvar %.0f\n", A$min_bulk, B$min_bulk
))
cat(sprintf(
  "min tail ESS: native %.0f, stanvar %.0f\n", A$min_tail, B$min_tail
))
cat(sprintf("wall: native %.0f s, stanvar %.0f s\n", A$elapsed, B$elapsed))
cat(sprintf(
  "posterior means: max |difference| %.4f, max |difference| / MCSE %.2f\n",
  max(delta), max(z)
))
print(data.frame(
  coefficient = colnames(X),
  native = round(A$means, 3),
  stanvar = round(B$means, 3),
  difference = round(A$means - B$means, 3),
  z = round(z, 2),
  row.names = NULL
))

verdict <- B$divergences <= A$divergences &&
  B$max_rhat <= 1.01 &&
  B$min_bulk >= 0.5 * A$min_bulk &&
  max(z) < 2
cat(sprintf(
  "\nEXIT RULE (divergences no worse, rhat <= 1.01, bulk ESS within a factor of two, every mean within 2 MCSE): %s\n",
  if (verdict) "MET -- the route is a candidate for wiring" else
    "NOT MET -- the route is not wired"
))

cat("\nsessionInfo:\n")
print(utils::sessionInfo()$otherPkgs$brms$Version)
print(cmdstanr::cmdstan_version())
