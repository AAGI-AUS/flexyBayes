# Model C with its per-trial SDs partially pooled, on brms only.
#
# `sd(group = "TRIAL:VARIETY", by = "TRIAL") ~ pooled(scale = ...)` gives the
# log genotype-by-trial SDs a normal distribution with an estimated mean and
# SD; `sigma(by = "TRIAL") ~ pooled(scale = ...)` does the same for the log
# residual SDs. The reference is the hand-written brms call of route B1 in
# review/units/c2_2026-10-08/spike_brms.R, with the package's parameter names.
# Recovery fixtures come from .mc_simulate() (helper-model-c-sim.R); the
# coverage check was declared in review/units/p1-pooling_2026-10-09/
# RECOVERY_CRITERION.md before any recovery fit ran.

.mcp_reference <- function() {
  file <- "model_c_reference.rds"
  path <- system.file("extdata", file, package = "flexyBayes")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "extdata", file)
  }
  readRDS(path)
}

.mcp_prior <- function(scale = 0.5) {
  fb_prior(
    sd(group = "TRIAL:VARIETY", by = "TRIAL") ~
      pooled(scale = half_normal(scale = scale)),
    sigma(by = "TRIAL") ~ pooled(scale = half_normal(scale = scale))
  )
}

# The priors of the B1 spike: brms's student_t(3, 0, 2.5) on the TRIAL and
# TRIAL:REP SDs, written out, and a unit half-normal on each pooling SD.
.mcp_b1_prior <- function() {
  fb_prior(
    sd(group = "TRIAL:VARIETY", by = "TRIAL") ~
      pooled(scale = half_normal(scale = 1)),
    sigma(by = "TRIAL") ~ pooled(scale = half_normal(scale = 1)),
    sd(group = "TRIAL") ~ student_t(df = 3, scale = 2.5),
    sd(group = "TRIAL:REP") ~ student_t(df = 3, scale = 2.5)
  )
}

.mcp_fit <- function(d, prior, backend = "brms", seed = 20261008L, ...) {
  suppressMessages(flexybayes(
    fixed = y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = d,
    prior = prior,
    backend = backend,
    chains = 2,
    warmup = 1000,
    n_samples = 1000,
    seed = seed,
    control = list(adapt_delta = 0.95),
    verbose = FALSE,
    mcmc_verbose = FALSE,
    ...
  ))
}

.mcp_terms <- function(d, prior) {
  .fb_from_asreml_impl(
    fixed = y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = d,
    family = "gaussian",
    prior = prior
  )
}

# Route B1 written by hand, with the names the package gives the pooling
# parameters and the location.
.mcp_b1_hand <- function(d) {
  f <- brms::bf(
    y ~ VARIETY + (1 | TRIAL) + (1 | TRIAL:REP) + (0 + TRIAL || VARIETY),
    sigma ~ 1 + (1 | TRIAL)
  )
  pr <- c(
    brms::set_prior("", class = "sd"),
    brms::set_prior("student_t(3, 0, 2.5)", class = "sd", group = "TRIAL"),
    brms::set_prior("student_t(3, 0, 2.5)", class = "sd", group = "TRIAL:REP"),
    brms::set_prior("normal(0, 1)", class = "sd", group = "TRIAL",
                    dpar = "sigma"),
    brms::set_prior("normal(pool_location, 2)", class = "Intercept",
                    dpar = "sigma")
  )
  sv <- brms::stanvar(x = log(stats::sd(d$y)), name = "pool_location") +
    brms::stanvar(
      scode = "real pool_mean_VARIETY;\n  real<lower=0> pool_sd_VARIETY;",
      block = "parameters"
    ) +
    brms::stanvar(scode = paste(
      "target += lognormal_lpdf(sd_3 | pool_mean_VARIETY, pool_sd_VARIETY);",
      "target += normal_lpdf(pool_mean_VARIETY | pool_location, 2);",
      "target += normal_lpdf(pool_sd_VARIETY | 0, 1);",
      sep = "\n  "
    ), block = "model")
  list(formula = f, prior = pr, stanvars = sv)
}

# The 26 per-trial SDs and the two pooling SDs, from a brmsfit's draws.
.mcp_sd_draws <- function(brmsfit, trials) {
  dr <- as.data.frame(posterior::as_draws_df(brmsfit))
  b0 <- dr[["b_sigma_Intercept"]]
  gxe <- vapply(trials, function(l) {
    dr[[paste0("sd_VARIETY__TRIAL", l)]]
  }, numeric(nrow(dr)))
  res <- vapply(trials, function(l) {
    exp(b0 + dr[[paste0("r_TRIAL__sigma[", l, ",Intercept]")]])
  }, numeric(nrow(dr)))
  colnames(gxe) <- paste0("gxe_", trials)
  colnames(res) <- paste0("res_", trials)
  cbind(
    gxe, res,
    pool_gxe = dr[["pool_sd_VARIETY"]],
    pool_res = dr[["sd_TRIAL__sigma_Intercept"]],
    mean_gxe = dr[["pool_mean_VARIETY"]],
    mean_res = b0
  )
}

.mcp_cache <- new.env()


test_that("fb_prior() reads pooled() with by = on sd() and on sigma", {
  p <- .mcp_prior()
  expect_length(p$specs, 2L)
  sd_spec <- p$specs[[1L]]
  expect_identical(sd_spec$target$type, "sd")
  expect_identical(sd_spec$target$group, "TRIAL:VARIETY")
  expect_identical(sd_spec$target$by, "TRIAL")
  expect_identical(sd_spec$spec$family, "pooled")
  expect_identical(sd_spec$spec$args$scale$family, "half_normal")
  expect_identical(sd_spec$spec$args$scale$args$scale, 0.5)
  expect_identical(p$specs[[2L]]$target, list(type = "sigma", by = "TRIAL"))
  printed <- utils::capture.output(print(p))
  expect_true(any(grepl(
    "sd(group = \"TRIAL:VARIETY\", by = \"TRIAL\") ~ pooled(scale = half_normal(scale = 0.5))",
    printed,
    fixed = TRUE
  )))
  expect_true(any(grepl("sigma(by = \"TRIAL\") ~ pooled(", printed,
                        fixed = TRUE)))
})

test_that("pooled() without its scale prior is refused, naming the missing half", {
  err <- expect_error(
    fb_prior(sigma(by = "TRIAL") ~ pooled()),
    class = "flexybayes_refusal_prior_argument_missing"
  )
  expect_match(conditionMessage(err), "pooling SD", fixed = TRUE)
  expect_match(conditionMessage(err), "`scale`", fixed = TRUE)
})

test_that("pooled() needs by =, and by = needs pooled()", {
  expect_error(
    fb_prior(sigma ~ pooled(scale = half_normal(scale = 0.5))),
    class = "flexybayes_refusal_prior_target_argument_missing"
  )
  expect_error(
    fb_prior(sd(group = "TRIAL:VARIETY", by = "TRIAL") ~
      half_normal(scale = 0.5)),
    class = "flexybayes_refusal_prior_target_unsupported"
  )
  expect_error(
    fb_prior(sigma(by = "TRIAL") ~
      pooled(scale = pooled(scale = half_normal(scale = 1)))),
    class = "flexybayes_refusal_prior_distribution_unknown"
  )
})

test_that("the INLA gate refuses pooled() by name and names the brms route", {
  ref <- .mcp_reference()
  gated <- lgm_gate(.mcp_terms(ref$data, .mcp_prior()))
  expect_true(is_lgm_refusal(gated))
  expect_identical(gated$primary_rule, "pooled_prior_requires_brms")
  expect_match(gated$failures[[1L]]$reason, "backend = \"brms\"",
               fixed = TRUE)
})

test_that("backend = \"inla\" refuses a pooled prior with a typed refusal", {
  skip_if_not_installed("INLA")
  ref <- .mcp_reference()
  err <- expect_error(
    .mcp_fit(ref$data, .mcp_prior(), backend = "inla"),
    class = "flexybayes_refusal_pooled_prior_requires_brms"
  )
  expect_identical(err$reason_code, "pooled_prior_requires_brms")
  expect_match(conditionMessage(err), "backend = \"brms\"", fixed = TRUE)
})

test_that("fb_plan() routes a pooled prior to brms under auto and says so", {
  skip_if_not_installed("brms")
  ref <- .mcp_reference()
  p <- fb_plan(
    y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = ref$data,
    priors = .mcp_prior()
  )
  expect_identical(p$backend_chosen, "brms")
  expect_identical(p$gate_primary_rule, "pooled_prior_requires_brms")
  expect_identical(p$pooled, "per-trial SDs partially pooled")
  printed <- utils::capture.output(print(p))
  expect_true(any(grepl("per-trial SDs partially pooled", printed,
                        fixed = TRUE)))
  expect_false(any(grepl("as ASReml fits it", printed, fixed = TRUE)))
})

test_that("a pooled prior on a term the model does not have is refused", {
  skip_if_not_installed("brms")
  ref <- .mcp_reference()
  expect_error(
    .mcp_fit(ref$data, fb_prior(
      sd(group = "TRIAL:VARIETY", by = "REP") ~
        pooled(scale = half_normal(scale = 0.5))
    ), return_code = TRUE),
    class = "flexybayes_refusal_prior_target_not_in_model"
  )
  expect_error(
    .mcp_fit(ref$data, fb_prior(
      sigma(by = "REP") ~ pooled(scale = half_normal(scale = 0.5))
    ), return_code = TRUE),
    class = "flexybayes_refusal_prior_target_not_in_model"
  )
})

test_that("a pooled residual is written sigma ~ 1 + (1 | TRIAL) on brms", {
  skip_if_not_installed("brms")
  ref <- .mcp_reference()
  f <- .fb_to_brms_formula(.mcp_terms(ref$data, .mcp_prior()))
  expect_identical(
    paste(deparse(f$pforms$sigma), collapse = ""),
    "sigma ~ 1 + (1 | TRIAL)"
  )
  gxe_only <- fb_prior(
    sd(group = "TRIAL:VARIETY", by = "TRIAL") ~
      pooled(scale = half_normal(scale = 0.5))
  )
  f <- .fb_to_brms_formula(.mcp_terms(ref$data, gxe_only))
  expect_identical(paste(deparse(f$pforms$sigma), collapse = ""),
                   "sigma ~ 0 + TRIAL")
})

test_that("the package writes the Stan program of the hand-written B1 call", {
  skip_if_not_installed("brms")
  skip_on_cran()
  ref <- .mcp_reference()
  d <- ref$data
  code <- .mcp_fit(d, .mcp_b1_prior(), return_code = TRUE)
  hand <- .mcp_b1_hand(d)
  hand_code <- brms::make_stancode(
    hand$formula, data = d, prior = hand$prior, stanvars = hand$stanvars
  )
  expect_identical(
    strsplit(code, "\n")[[1L]],
    strsplit(as.character(hand_code), "\n")[[1L]]
  )
})

# Code inspection under auto resolves to brms (see ?flexybayes).
test_that("auto returns the brms program of a pooled model", {
  skip_if_not_installed("brms")
  skip_on_cran()
  ref <- .mcp_reference()
  code <- flexybayes(
    fixed = y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = ref$data,
    prior = .mcp_prior(),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_match(code, "pool_sd_VARIETY", fixed = TRUE)
})

# The risk the plan names: the hierarchical prior is written on brms's own
# name for the per-level SDs, and a brms release can rename it.
test_that("the pooled prior lands on a per-level SD vector brms declares", {
  skip_if_not_installed("brms")
  skip_on_cran()
  ref <- .mcp_reference()
  d <- ref$data
  code <- .mcp_fit(d, .mcp_prior(), return_code = TRUE)
  line <- grep("lognormal_lpdf(sd_", strsplit(code, "\n")[[1L]],
               fixed = TRUE, value = TRUE)
  expect_length(line, 1L)
  sd_name <- sub("^.*lognormal_lpdf\\((sd_[0-9]+) .*$", "\\1", line)
  expect_true(.fb_brms_stancode_declares(code, sd_name))
  hand <- .mcp_b1_hand(d)
  sdata <- brms::make_standata(hand$formula, data = d, prior = hand$prior,
                               stanvars = hand$stanvars)
  k <- sub("^sd_", "", sd_name)
  expect_identical(as.integer(sdata[[paste0("M_", k)]]), 12L)
  expect_identical(as.integer(sdata[[paste0("N_", k)]]), 30L)
  # The check fails on a program in which the name is gone.
  planted <- gsub(paste0("\\b", sd_name, "\\b"), "sd_99", code, perl = TRUE)
  planted <- sub("sd_99 \\|", paste0(sd_name, " |"), planted)
  expect_false(.fb_brms_stancode_declares(planted, sd_name))
})

test_that("the recovery helper reproduces the stored model C fixture", {
  ref <- .mcp_reference()
  sim <- .mc_simulate(20261008L)
  expect_identical(sim$data, ref$data)
  expect_identical(sim$truth$sd_gxe, ref$truth$sd_gxe)
  expect_identical(sim$truth$sd_res, ref$truth$sd_res)
})

test_that("the package reproduces the B1 posterior to Monte Carlo error", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(mc.cores = 2L)
  ref <- .mcp_reference()
  d <- ref$data
  trials <- levels(d$TRIAL)
  fit <- .mcp_fit(d, .mcp_b1_prior())
  .mcp_cache$fit <- fit
  hand <- .mcp_b1_hand(d)
  b1 <- brms::brm(
    hand$formula, data = d, prior = hand$prior, stanvars = hand$stanvars,
    chains = 2, iter = 2000, warmup = 1000, seed = 20261009L,
    control = list(adapt_delta = 0.95), refresh = 0, silent = 2
  )
  a <- .mcp_sd_draws(fit$brms, trials)
  b <- .mcp_sd_draws(b1, trials)
  mcse <- function(m) {
    apply(m, 2L, function(x) posterior::mcse_mean(x))
  }
  z <- (colMeans(a) - colMeans(b)) / sqrt(mcse(a)^2 + mcse(b)^2)
  expect_true(all(abs(z) < 4), info = paste(names(z)[abs(z) >= 4],
                                            collapse = ", "))
})

test_that("every reader of a pooled fit says the SDs are partially pooled", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(mc.cores = 2L)
  fit <- .mcp_cache$fit %||% .mcp_fit(.mcp_reference()$data, .mcp_b1_prior())
  label <- "per-trial SDs partially pooled"
  readers <- list(
    print = utils::capture.output(print(fit)),
    summary = utils::capture.output(print(summary(fit))),
    backend_decision = backend_decision(fit)$reason,
    prior_summary = utils::capture.output(print(prior_summary(fit)))
  )
  for (nm in names(readers)) {
    expect_true(any(grepl(label, readers[[nm]], fixed = TRUE)), info = nm)
    expect_false(any(grepl("as ASReml fits it", readers[[nm]], fixed = TRUE)),
                 info = nm)
  }
  # The prior column names the pooled specification and the pooling SD's
  # own prior.
  vc <- summary(fit)$varcomp
  prior_of <- stats::setNames(vc$prior, vc$component)
  expect_identical(prior_of[["sd_VARIETY__TRIALT01"]],
                   "pooled(scale=half_normal(scale=1))")
  expect_identical(prior_of[["sigma_T06"]],
                   "pooled(scale=half_normal(scale=1))")
  expect_identical(prior_of[["sd_TRIAL__sigma_Intercept"]],
                   "half_normal(scale=1)")
  expect_identical(prior_of[["sd_TRIAL"]], "student_t(df=3, scale=2.5)")
  # The per-trial residual SDs are read from the pooled lowering.
  tab <- .brms_residual_by_level_table(fit)
  expect_identical(tab$level, levels(fit$data$TRIAL))
  # Every pooled block reads back from the program brms compiled.
  m <- fit$extras$manifest
  pooled_slots <- m$slots[grepl("^stan:", m$slots$slot_id), , drop = FALSE]
  expect_identical(nrow(pooled_slots), 3L)
  expect_true(all(pooled_slots$readback_agrees))
  expect_true(m$readback_complete)
})

test_that("auto fits a pooled model on brms and its decision says so", {
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(mc.cores = 2L)
  ref <- .mcp_reference()
  d <- droplevels(ref$data[ref$data$TRIAL %in% c("T01", "T02", "T03"), ])
  # One chain of 200 draws checks the route, not the posterior; Stan's
  # convergence warnings on it are expected.
  fit <- suppressWarnings(suppressMessages(flexybayes(
    fixed = y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = d, prior = .mcp_prior(), chains = 1, warmup = 200,
    n_samples = 200, seed = 1L, verbose = FALSE, mcmc_verbose = FALSE
  )))
  expect_s3_class(fit, "flexybayes_brms")
  bd <- backend_decision(fit)
  expect_identical(bd$backend, "brms")
  expect_match(bd$reason, "pooled_prior_requires_brms", fixed = TRUE)
  expect_match(bd$reason, "per-trial SDs partially pooled", fixed = TRUE)
})

# Declared in RECOVERY_CRITERION.md before the runs: central 90% intervals,
# 20 fixtures, each count inside the central 99% range of Binomial(n, 0.9).
.mcp_recovery_one <- function(seed) {
  sim <- .mc_simulate(seed)
  trials <- levels(sim$data$TRIAL)
  t0 <- proc.time()
  fit <- .mcp_fit(sim$data, .mcp_prior(), seed = seed)
  elapsed <- unname((proc.time() - t0)[["elapsed"]])
  dr <- .mcp_sd_draws(fit$brms, trials)
  q <- apply(dr, 2L, stats::quantile, probs = c(0.05, 0.95), names = FALSE)
  truth <- c(
    stats::setNames(sim$truth$sd_gxe, paste0("gxe_", trials)),
    stats::setNames(sim$truth$sd_res, paste0("res_", trials)),
    pool_gxe = sim$truth$pool_gxe,
    pool_res = sim$truth$pool_res
  )
  nm <- names(truth)
  covered <- truth >= q[1L, nm] & truth <= q[2L, nm]
  data.frame(
    seed = seed,
    elapsed = elapsed,
    max_rhat = max(brms::rhat(fit$brms), na.rm = TRUE),
    divergent = sum(brms::nuts_params(fit$brms, pars = "divergent__")$Value),
    pool_gxe_covered = covered[["pool_gxe"]],
    pool_res_covered = covered[["pool_res"]],
    pool_gxe_median = stats::median(dr[, "pool_gxe"]),
    pool_res_median = stats::median(dr[, "pool_res"]),
    gxe_covered = sum(covered[startsWith(nm, "gxe_")]),
    res_covered = sum(covered[startsWith(nm, "res_")])
  )
}

test_that("pooled intervals cover the simulated truth at the 90% rate", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(mc.cores = 2L)
  out_dir <- Sys.getenv("FLEXYBAYES_RECOVERY_DIR")
  rows <- list()
  for (seed in 20261008L + 0:19) {
    row <- .mcp_recovery_one(seed)
    if (nzchar(out_dir)) {
      utils::write.csv(row, file.path(out_dir, paste0("recovery_", seed, ".csv")),
                       row.names = FALSE)
    }
    rows[[length(rows) + 1L]] <- row
  }
  tab <- do.call(rbind, rows)
  expect_gte(sum(tab$pool_gxe_covered), 14L)
  expect_gte(sum(tab$pool_res_covered), 14L)
  expect_true(sum(tab$gxe_covered) >= 203L && sum(tab$gxe_covered) <= 227L)
  expect_true(sum(tab$res_covered) >= 203L && sum(tab$res_covered) <= 227L)
})
