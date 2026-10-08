# Model C: a genotype-by-trial variance and a residual variance for every
# trial, `diag(TRIAL):VARIETY` with `dsum(~ units | TRIAL)`. With 12 trials it
# carries 26 variance parameters, past the 15 INLA integrates over
# numerically. Above that limit INLA fits it by empirical Bayes
# (`int.strategy = "eb"`), and the output says so.
#
# The oracle is ASReml 4.2 on the same data, stored by
# data-raw/model_c_reference.R. The tolerance is half of ASReml's standard
# error on every component ASReml identifies (std.error / component < 0.5),
# as in data-raw/asreml_oracle/SPEC.md. The largest deviation measured on
# this fixture was 0.17 of a standard error on INLA and 0.28 on brms.

.mc_reference <- function() {
  file <- "model_c_reference.rds"
  path <- system.file("extdata", file, package = "flexyBayes")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "extdata", file)
  }
  readRDS(path)
}

.mc_terms <- function(d) {
  fb_from_asreml(
    fixed = y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = d,
    family = "gaussian"
  )
}

.mc_cache <- new.env()

# One INLA fit of the full fixture, shared by the blocks that read it.
.mc_inla_fit <- function() {
  if (is.null(.mc_cache$fit)) {
    withr::local_options(
      flexyBayes.inla_threads = "1:1",
      flexyBayes.silence_default_prior_note = TRUE
    )
    .mc_cache$fit <- .mc_fit(.mc_reference()$data, backend = "inla")
  }
  .mc_cache$fit
}

.mc_fit <- function(d, backend = "inla") {
  suppressMessages(flexybayes(
    fixed = y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = d,
    backend = backend,
    verbose = FALSE
  ))
}

# Posterior medians of the per-trial variances from a fit's variance draws.
.mc_medians <- function(fit, trials) {
  dr <- as.data.frame(fb_variance_draws(fit, n = 4000L, seed = 1L))
  list(
    gxe = vapply(trials, function(l) {
      stats::median(dr[[paste0("sd_VARIETY__TRIAL", l)]]^2)
    }, numeric(1)),
    res = vapply(trials, function(l) {
      stats::median(dr[[paste0("sigma_", l)]]^2)
    }, numeric(1))
  )
}

# The oracle comparison. Returns one message per failure.
.mc_against_asreml <- function(med, a) {
  ident_g <- a$gxe_se / a$gxe_var < 0.5
  ident_r <- a$res_se / a$res_var < 0.5
  dev_g <- abs(med$gxe - a$gxe_var) / a$gxe_se
  dev_r <- abs(med$res - a$res_var) / a$res_se
  ratio_dev <- abs(med$gxe / (med$gxe + med$res) -
    a$gxe_var / (a$gxe_var + a$res_var))
  both <- ident_g & ident_r
  off <- function(what, bad) {
    if (length(bad)) paste0(what, " off: ", paste(bad, collapse = ", "))
  }
  c(
    off("genotype-by-trial variance", names(dev_g)[ident_g & dev_g > 0.5]),
    off("residual variance", names(dev_r)[ident_r & dev_r > 0.5]),
    off("variance ratio", names(ratio_dev)[both & ratio_dev > 0.05])
  ) %||% character(0)
}


test_that("DEF-207: the count has one precision per dsum() section", {
  ref <- .mc_reference()
  fb <- .mc_terms(ref$data)
  # INLA reports 26 hyperparameters for this model: 12 observation
  # precisions, 12 genotype-by-trial precisions, TRIAL and TRIAL:REP.
  expect_identical(.lgm_count_hyperparams(fb), 26L)
})

test_that("model C passes the INLA gate with empirical-Bayes integration", {
  ref <- .mc_reference()
  gated <- lgm_gate(.mc_terms(ref$data))
  expect_false(is_lgm_refusal(gated))
  expect_true("inla_integration:eb" %in% gated$capabilities)
})

test_that("the gate still refuses when the structured part is over the limit", {
  set.seed(1L)
  d <- data.frame(y = stats::rnorm(64L))
  for (k in seq_len(15L)) {
    d[[paste0("f", k)]] <- factor(sample(letters[1:4], 64L, replace = TRUE))
  }
  fb <- fb_from_asreml(
    fixed = y ~ 1,
    random = stats::as.formula(paste("~", paste0("f", 1:15, collapse = " + "))),
    data = d,
    family = "gaussian"
  )
  gated <- lgm_gate(fb)
  expect_true(is_lgm_refusal(gated))
  ids <- vapply(gated$failures, function(r) r$rule_id, character(1))
  expect_true("hyperparam_budget" %in% ids)
})

test_that("a model within the budget keeps INLA's own integration", {
  ref <- .mc_reference()
  d <- droplevels(ref$data[ref$data$TRIAL %in% c("T01", "T02", "T03"), ])
  gated <- lgm_gate(.mc_terms(d))
  expect_false(is_lgm_refusal(gated))
  expect_false("inla_integration:eb" %in% gated$capabilities)
})

test_that("model C fits on INLA by empirical Bayes and agrees with ASReml", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE
  )
  ref <- .mc_reference()
  fit <- .mc_inla_fit()

  expect_identical(fit$inla$.args$control.inla$int.strategy, "eb")
  expect_identical(fit$extras$inla_integration$strategy, "eb")
  expect_identical(
    nrow(fit$inla$summary.hyperpar),
    .lgm_count_hyperparams(fit$fb)
  )
  expect_match(backend_decision(fit)$reason, "empirical Bayes", fixed = TRUE)
  printed <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "empirical Bayes", fixed = TRUE)
  printed <- paste(utils::capture.output(print(summary(fit))), collapse = "\n")
  expect_match(printed, "empirical Bayes", fixed = TRUE)

  med <- .mc_medians(fit, levels(ref$data$TRIAL))
  expect_identical(.mc_against_asreml(med, ref$asreml), character(0))

  # The comparison is shown to reject a fit whose trials are mislabelled.
  swapped <- med
  swapped$gxe[c("T05", "T10")] <- swapped$gxe[c("T10", "T05")]
  expect_gt(length(.mc_against_asreml(swapped, ref$asreml)), 0L)
})

test_that("DEF-208: the print of a dsum() fit on INLA names its family once", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  d <- .mc_reference()$data
  d <- droplevels(d[d$TRIAL %in% c("T01", "T02", "T03"), ])
  printed <- utils::capture.output(print(.mc_fit(d)))
  expect_true(any(grepl("^  Family   : gaussian / identity $", printed)))
})

test_that("DEF-209: the soft-limit warning does not claim a ccd switch", {
  re <- replicate(10, list(type = "simple", var = "g"), simplify = FALSE)
  fb <- new_fb_terms(
    response = "y", family = "gaussian", fixed_terms = list(),
    random_terms = re, residual_terms = list(list(type = "units")),
    addition_terms = list(), priors = list(legacy = TRUE), source = "asreml"
  )
  r <- .lgm_check_hyperparam_budget(fb)
  expect_match(r$warning, "exceeds soft limit", fixed = TRUE)
  expect_false(grepl("ccd", r$warning, fixed = TRUE))
})

test_that("fb_plan() says model C fits on INLA by empirical Bayes", {
  skip_if_not_installed("INLA")
  ref <- .mc_reference()
  p <- fb_plan(
    y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = ref$data
  )
  expect_true(p$will_fit)
  expect_identical(p$backend_chosen, "inla")
  expect_identical(p$inla_integration, "eb")
  printed <- utils::capture.output(print(p))
  expect_true(any(grepl("empirical Bayes", printed, fixed = TRUE)))
})

test_that("backend = auto fits model C on INLA by empirical Bayes", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE
  )
  ref <- .mc_reference()
  fit <- .mc_fit(ref$data, backend = "auto")
  expect_s3_class(fit, "flexybayes_inla")
  expect_identical(fit$extras$inla_integration$strategy, "eb")
})

test_that("model C on brms agrees with ASReml", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_if_not_long_tier()
  ref <- .mc_reference()
  lev <- levels(ref$data$TRIAL)
  fit <- suppressMessages(flexybayes(
    fixed = y ~ VARIETY,
    random = ~ TRIAL + TRIAL:REP + diag(TRIAL):VARIETY,
    residual = ~ dsum(~ units | TRIAL),
    data = ref$data, backend = "brms", chains = 2, warmup = 1000,
    n_samples = 1000, seed = 20261008L, control = list(adapt_delta = 0.95),
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  dr <- posterior::as_draws_df(fit$brms)
  med <- list(
    gxe = vapply(lev, function(l) {
      stats::median(dr[[paste0("sd_VARIETY__TRIAL", l)]]^2)
    }, numeric(1)),
    res = vapply(lev, function(l) {
      stats::median(exp(2 * dr[[paste0("b_sigma_TRIAL", l)]]))
    }, numeric(1))
  )
  expect_identical(.mc_against_asreml(med, ref$asreml), character(0))
})
