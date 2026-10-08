# data-raw/model_c_reference.R
#
# A simulated series of trials with its own genotype-by-trial variance and
# its own residual variance in every trial, fitted once in ASReml and stored,
# so the tests can check flexyBayes against ASReml without an ASReml licence.
# The model is model C of the barley comparison study:
#
#   fixed    = y ~ VARIETY
#   random   = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY
#   residual = ~ dsum(~ units | TRIAL)
#
# Run from the package root with a licensed ASReml on the library path:
#
#   FLEXYBAYES_ASREML_LIB=/path/to/library \
#     Rscript --vanilla data-raw/model_c_reference.R
#
# Unset, the default library is used. Every number the tests attribute to
# ASReml comes from this file; edit the stored numbers by hand and the
# provenance is gone.
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

# --- 1. the data --------------------------------------------------------- #

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261008L)

n_trial <- 12L
n_variety <- 30L
n_rep <- 3L
trial_lev <- sprintf("T%02d", seq_len(n_trial))
variety_lev <- sprintf("V%02d", seq_len(n_variety))

# Per-trial standard deviations drawn on the log scale, so the trials differ
# by a factor of about three in each.
sd_gxe <- exp(stats::rnorm(n_trial, log(0.35), 0.4))
sd_res <- exp(stats::rnorm(n_trial, log(0.50), 0.35))
names(sd_gxe) <- trial_lev
names(sd_res) <- trial_lev
sd_trial <- 1.0
sd_rep <- 0.3

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
eff_trial <- stats::rnorm(n_trial, 0, sd_trial)
eff_rep <- matrix(stats::rnorm(n_trial * n_rep, 0, sd_rep), n_trial, n_rep)
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

# --- 2. the ASReml fit ---------------------------------------------------- #

t0 <- proc.time()
m <- asreml(
  fixed = y ~ VARIETY,
  random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
  residual = ~ dsum(~ units | TRIAL),
  data = d,
  trace = FALSE,
  maxit = 50
)
elapsed <- unname((proc.time() - t0)[["elapsed"]])
stopifnot(isTRUE(m$converge))

vc <- summary(m)$varcomp
varcomp <- data.frame(
  term = rownames(vc),
  component = vc$component,
  se = vc$std.error,
  bound = as.character(vc$bound),
  stringsAsFactors = FALSE
)

# ASReml names the genotype-by-trial variances `TRIAL:VARIETY!TRIAL_T01` and
# the residuals `TRIAL_T01!R`; read them by trial so the order is explicit.
pick <- function(pattern) {
  out <- vapply(
    trial_lev,
    function(l) {
      i <- grep(sprintf(pattern, l), varcomp$term, fixed = TRUE)
      stopifnot(length(i) == 1L)
      varcomp$component[i]
    },
    numeric(1)
  )
  out
}
pick_se <- function(pattern) {
  vapply(
    trial_lev,
    function(l) {
      i <- grep(sprintf(pattern, l), varcomp$term, fixed = TRUE)
      varcomp$se[i]
    },
    numeric(1)
  )
}

data_sha256 <- digest::digest(d, algo = "sha256")
script_path <- "data-raw/model_c_reference.R"

model_c_reference <- list(
  data = d,
  truth = list(
    sd_gxe = sd_gxe,
    sd_res = sd_res,
    sd_trial = sd_trial,
    sd_rep = sd_rep
  ),
  asreml = list(
    varcomp = varcomp,
    gxe_var = pick("TRIAL:VARIETY!TRIAL_%s"),
    gxe_se = pick_se("TRIAL:VARIETY!TRIAL_%s"),
    res_var = pick("TRIAL_%s!R"),
    res_se = pick_se("TRIAL_%s!R"),
    trial_var = varcomp$component[varcomp$term == "TRIAL"],
    rep_var = varcomp$component[varcomp$term == "TRIAL:REP"],
    loglik = m$loglik,
    converged = isTRUE(m$converge),
    elapsed_sec = elapsed
  ),
  stamp = list(
    date = as.character(Sys.Date()),
    asreml = as.character(utils::packageVersion("asreml")),
    r = R.version.string,
    platform = R.version$platform,
    rng_kind = paste(RNGkind(), collapse = " / "),
    seed = 20261008L,
    data_sha256 = data_sha256,
    script_sha256 = digest::digest(file = script_path, algo = "sha256"),
    run_command = paste0("Rscript --vanilla ", script_path)
  )
)

saveRDS(model_c_reference, "inst/extdata/model_c_reference.rds")
cat("written inst/extdata/model_c_reference.rds\n")
print(model_c_reference$stamp)
print(round(cbind(
  gxe_true = sd_gxe^2,
  gxe_asreml = model_c_reference$asreml$gxe_var,
  res_true = sd_res^2,
  res_asreml = model_c_reference$asreml$res_var
), 4))
