# data-raw/asreml_predict_new_level.R
#
# The ASReml 4.2 side of the new-level prediction test
# (tests/testthat/test-predict-new-levels-asreml.R), computed once and stored
# in tests/testthat/fixtures/asreml_predict_new_level.rds, so the test runs
# without an ASReml licence.
#
# A simulated multi-environment trial: 16 localities with 2 trials each, each
# trial a complete 4 x 6 lattice of plots holding 6 varieties x 2 N rates in
# 2 replicates. The model is the comparison study's base model in small:
#
#   y ~ VARIETY * NRATE,
#   random = ~ loc + TRIAL + TRIAL:REP + TRIAL:VARIETY + diag(design):mainplot
#              + id(TRIAL):ar1(RANGE):ar1(ROW)
#
# ASReml is fitted without locality L16, and the trials of L16 are predicted:
# each trial mean, three cell means, and the locality mean over both trials.
# For a set P of plots at the new locality the plug-in prediction is
#
#   mean      the mean over P of ASReml's fixed-effect cell predictions
#   variance  w' V w, V the covariance of the plots of P under the REML
#             variances with every level of P new, w = 1 / |P|; plus w' C w,
#             C the covariance of the cell predictions (predict(vcov = TRUE))
#
# w' V w is computed twice, as a dense matrix and by the sum-of-shares
# formula of the study's r_code_04, and the two must agree.
#
# Run with a licensed ASReml from the package root:
#   FLEXYBAYES_ASREML_LIB=/path/to/library \
#     Rscript --vanilla data-raw/asreml_predict_new_level.R
# Unset, the default library is used.

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
asreml.options(trace = FALSE, threads = 1L)

out_file <- file.path(
  "tests",
  "testthat",
  "fixtures",
  "asreml_predict_new_level.rds"
)
stopifnot(dir.exists(dirname(out_file)))


# --- 1. the data --------------------------------------------------------- #

sim_met <- function(seed, n_loc, n_trial) {
  set.seed(seed)
  n_range <- 4L
  n_row <- 6L
  var_lev <- paste0("V", 1:6)
  n_lev <- c("lo", "hi")
  sd <- c(
    loc = 0.6,
    trial = 0.5,
    rep = 0.15,
    tv = 0.25,
    mp = 0.12,
    field = 0.3,
    e = 0.2
  )
  beta_v <- stats::setNames(seq(-0.5, 0.5, length.out = 6L), var_lev)
  rr <- 0.5^abs(outer(seq_len(n_range), seq_len(n_range), "-"))
  rc <- 0.7^abs(outer(seq_len(n_row), seq_len(n_row), "-"))
  chol_field <- t(chol(kronecker(rr, rc)))
  out <- list()
  k <- 0L
  for (l in seq_len(n_loc)) {
    u_loc <- stats::rnorm(1L, 0, sd[["loc"]])
    for (t in seq_len(n_trial)) {
      k <- k + 1L
      g <- expand.grid(ROW = seq_len(n_row), RANGE = seq_len(n_range))
      g$REP <- ifelse(g$RANGE <= 2L, 1L, 2L)
      g$cell <- NA_integer_
      for (r in 1:2) {
        g$cell[g$REP == r] <- sample(12L)
      }
      g$VARIETY <- var_lev[(g$cell - 1L) %/% 2L + 1L]
      g$NRATE <- n_lev[(g$cell - 1L) %% 2L + 1L]
      mp <- paste(g$REP, g$VARIETY)
      u_mp <- stats::rnorm(length(unique(mp)), 0, sd[["mp"]])
      u_rep <- stats::rnorm(2L, 0, sd[["rep"]])
      u_tv <- stats::setNames(stats::rnorm(6L, 0, sd[["tv"]]), var_lev)
      # kronecker(rr, rc) orders the nodes range by range
      field <- drop(chol_field %*% stats::rnorm(n_range * n_row)) *
        sd[["field"]]
      node <- (g$RANGE - 1L) * n_row + g$ROW
      y <- 3 +
        beta_v[g$VARIETY] +
        ifelse(g$NRATE == "hi", 0.4, 0) +
        u_loc +
        stats::rnorm(1L, 0, sd[["trial"]]) +
        u_rep[g$REP] +
        u_tv[g$VARIETY] +
        u_mp[match(mp, unique(mp))] +
        field[node] +
        stats::rnorm(nrow(g), 0, sd[["e"]])
      out[[k]] <- data.frame(
        loc = sprintf("L%02d", l),
        TRIAL = sprintf("T%02d", k),
        design = if (k %% 2L == 0L) "split" else "rcb",
        REP = g$REP,
        RANGE = g$RANGE,
        ROW = g$ROW,
        VARIETY = g$VARIETY,
        NRATE = g$NRATE,
        y = unname(y),
        stringsAsFactors = FALSE
      )
    }
  }
  d <- do.call(rbind, out)
  d$mainplot <- paste(d$TRIAL, d$REP, d$VARIETY, sep = "_")
  for (v in c(
    "loc",
    "TRIAL",
    "design",
    "REP",
    "RANGE",
    "ROW",
    "VARIETY",
    "NRATE",
    "mainplot"
  )) {
    d[[v]] <- factor(d[[v]], levels = sort(unique(d[[v]])))
  }
  d
}

seed <- 20261008L
data <- sim_met(seed, n_loc = 16L, n_trial = 2L)
held_out <- "L16"
train <- droplevels(data[data$loc != held_out, ])
test <- droplevels(data[data$loc == held_out, ])
rownames(test) <- NULL


# --- 2. the ASReml fit --------------------------------------------------- #

fit <- asreml(
  fixed = y ~ VARIETY * NRATE,
  random = ~ loc +
    TRIAL +
    TRIAL:REP +
    TRIAL:VARIETY +
    diag(design):mainplot +
    id(TRIAL):ar1(RANGE):ar1(ROW),
  residual = ~units,
  data = train
)
n_update <- 0L
while (!isTRUE(fit$converge) && n_update < 10L) {
  fit <- update(fit)
  n_update <- n_update + 1L
}
stopifnot(isTRUE(fit$converge))
vc <- summary(fit)$varcomp

# ASReml's names to the package's variance-component names; variances, and
# the two correlations of the field as estimated
pick <- function(name) {
  stopifnot(name %in% rownames(vc))
  vc[name, "component"]
}
variance <- c(
  sd_loc = pick("loc"),
  sd_TRIAL = pick("TRIAL"),
  `sd_TRIAL:REP` = pick("TRIAL:REP"),
  `sd_TRIAL:VARIETY` = pick("TRIAL:VARIETY"),
  sd_mainplot__designrcb = pick("design:mainplot!design_rcb"),
  sd_mainplot__designsplit = pick("design:mainplot!design_split"),
  sd_spatial = pick("TRIAL:RANGE:ROW"),
  sigma = pick("units!R")
)
correlation <- c(
  rho_row = pick("TRIAL:RANGE:ROW!RANGE!cor"),
  rho_col = pick("TRIAL:RANGE:ROW!ROW!cor")
)
stopifnot(length(variance) + length(correlation) == nrow(vc))


# --- 3. the plug-in prediction of each set ------------------------------ #

pv <- predict(fit, classify = "VARIETY:NRATE", vcov = TRUE)
cell <- pv$pvals
cell_key <- paste(cell$VARIETY, cell$NRATE)
cell_vcov <- as.matrix(pv$vcov)

sets <- list(
  trial_T31 = which(test$TRIAL == "T31"),
  trial_T32 = which(test$TRIAL == "T32"),
  cell_T31_V1_lo = which(
    test$TRIAL == "T31" & test$VARIETY == "V1" & test$NRATE == "lo"
  ),
  cell_T31_V6_hi = which(
    test$TRIAL == "T31" & test$VARIETY == "V6" & test$NRATE == "hi"
  ),
  cell_T32_V3_hi = which(
    test$TRIAL == "T32" & test$VARIETY == "V3" & test$NRATE == "hi"
  ),
  locality_L16 = seq_len(nrow(test))
)
stopifnot(all(lengths(sets) > 0L))

# the covariance of the plots of `p`, every level new, as a dense matrix
dense_v <- function(p) {
  d <- test[p, ]
  same <- function(...) outer(paste(...), paste(...), "==") * 1
  r <- as.integer(as.character(d$RANGE))
  w <- as.integer(as.character(d$ROW))
  field <- same(d$TRIAL) *
    correlation[["rho_row"]]^abs(outer(r, r, "-")) *
    correlation[["rho_col"]]^abs(outer(w, w, "-"))
  mp <- same(d$mainplot) * (outer(d$design, d$design, "==") * 1)
  mp_var <- ifelse(
    d$design == "rcb",
    variance[["sd_mainplot__designrcb"]],
    variance[["sd_mainplot__designsplit"]]
  )
  variance[["sd_loc"]] *
    same(d$loc) +
    variance[["sd_TRIAL"]] * same(d$TRIAL) +
    variance[["sd_TRIAL:REP"]] * same(d$TRIAL, d$REP) +
    variance[["sd_TRIAL:VARIETY"]] * same(d$TRIAL, d$VARIETY) +
    mp * sqrt(outer(mp_var, mp_var)) +
    variance[["sd_spatial"]] * field +
    diag(variance[["sigma"]], nrow(d))
}

# the variance of the mean over the plots of `p`, one trial, by shares
share_v <- function(p) {
  d <- test[p, ]
  n <- nrow(d)
  share2 <- function(g) sum((tabulate(match(g, unique(g))) / n)^2)
  r <- as.integer(as.character(d$RANGE))
  w <- as.integer(as.character(d$ROW))
  mp_name <- paste0("sd_mainplot__design", as.character(d$design[1L]))
  variance[["sd_loc"]] +
    variance[["sd_TRIAL"]] +
    variance[["sd_TRIAL:REP"]] * share2(d$REP) +
    variance[["sd_TRIAL:VARIETY"]] * share2(d$VARIETY) +
    variance[[mp_name]] * share2(as.character(d$mainplot)) +
    variance[["sd_spatial"]] *
      mean(
        outer(r, r, function(a, b) correlation[["rho_row"]]^abs(a - b)) *
          outer(w, w, function(a, b) correlation[["rho_col"]]^abs(a - b))
      ) +
    variance[["sigma"]] / n
}

pred <- do.call(
  rbind,
  lapply(names(sets), function(nm) {
    p <- sets[[nm]]
    k <- match(paste(test$VARIETY[p], test$NRATE[p]), cell_key)
    stopifnot(!anyNA(k))
    wc <- tabulate(k, nrow(cell)) / length(p)
    v_dense <- sum(dense_v(p)) / length(p)^2
    v_share <- if (length(unique(test$TRIAL[p])) == 1L) share_v(p) else NA_real_
    data.frame(
      set = nm,
      n_rows = length(p),
      mean = sum(wc * cell$predicted.value),
      var_fixed = drop(t(wc) %*% cell_vcov %*% wc),
      var_random = v_dense,
      var_random_share = v_share,
      stringsAsFactors = FALSE
    )
  })
)
pred$sd <- sqrt(pred$var_fixed + pred$var_random)
ok <- is.na(pred$var_random_share) |
  abs(pred$var_random_share / pred$var_random - 1) < 1e-12
stopifnot(all(ok))


# --- 4. the stored fixture ------------------------------------------------ #

fixture <- list(
  data = data,
  held_out = held_out,
  sets = sets,
  prediction = pred,
  variance = variance,
  correlation = correlation,
  varcomp = vc,
  call = list(
    fixed = "y ~ VARIETY * NRATE",
    random = paste(
      "~ loc + TRIAL + TRIAL:REP + TRIAL:VARIETY + diag(design):mainplot +",
      "id(TRIAL):ar1(RANGE):ar1(ROW)"
    ),
    residual = "~ units"
  ),
  stamp = list(
    date = as.character(Sys.Date()),
    asreml = as.character(utils::packageVersion("asreml")),
    r = R.version.string,
    seed = seed,
    n_update = n_update,
    loglik = fit$loglik
  )
)
saveRDS(fixture, out_file, version = 2)
print(vc)
print(pred)
