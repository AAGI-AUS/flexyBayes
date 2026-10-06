# fb_plan() must agree with the fit. A multi-environment trial model with
# a per-design main-plot variance, diag(design):mainplot, and a per-trial
# AR1 x AR1 field was planned as "Will fit: no  (preflight refused)" and
# then fitted on INLA by flexybayes(). The fit runs the design-memory
# preflight only from 100,000 rows; the plan applied it to every model.
# The model with one variety-by-trial variance per trial is still refused
# by the hyperparameter budget, by both the plan and the fit.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

.met_field_data <- function(n_trial = 6L, seed = 20261006L) {
  set.seed(seed)
  d <- expand.grid(
    ROW = factor(seq_len(6L)),
    RANGE = factor(seq_len(4L)),
    TRIAL = factor(sprintf("T%d", seq_len(n_trial)))
  )
  d$REP <- factor(ifelse(as.integer(d$RANGE) <= 2L, "1", "2"))
  d$VARIETY <- factor(sample(c("V1", "V2", "V3"), nrow(d), replace = TRUE))
  d$NAPPLIED <- factor(sample(c("N0", "N1"), nrow(d), replace = TRUE))
  d$SEEDRATE <- factor(sample(c("S1", "S2"), nrow(d), replace = TRUE))
  d$loc_f <- factor(sprintf("L%d", (as.integer(d$TRIAL) + 1L) %/% 2L))
  d$design_f <- factor(ifelse(as.integer(d$TRIAL) %% 2L == 0L, "split", "rcb"))
  d$state_f <- factor(ifelse(as.integer(d$TRIAL) <= n_trial / 2L, "A", "B"))
  d$mainplot_f <- factor(paste(d$TRIAL, d$REP, d$VARIETY, d$NAPPLIED, sep = "_"))
  d$tv_f <- factor(paste(d$TRIAL, d$VARIETY, sep = "_"))
  d$y <- stats::rnorm(nrow(d), 3 + as.integer(d$TRIAL) / 2)
  d$w <- rep(seq(0.5, 1.5, length.out = n_trial), each = 24L)
  d
}

.met_field_random <- function(gxe) {
  stats::as.formula(paste(
    "~ loc_f + TRIAL + TRIAL:REP + TRIAL:RANGE + TRIAL:ROW +",
    gxe,
    "+ diag(design_f):mainplot_f + id(TRIAL):ar1(RANGE):ar1(ROW)"
  ))
}

.met_field_plan <- function(d, gxe) {
  fb_plan(
    y ~ VARIETY * NAPPLIED * SEEDRATE,
    random = .met_field_random(gxe),
    residual = ~units,
    data = d,
    weights = d$w
  )
}

.met_field_fit <- function(d, gxe) {
  withr::local_options(list(
    flexyBayes.silence_spatial_collapse_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE
  ))
  flexybayes(
    y ~ VARIETY * NAPPLIED * SEEDRATE,
    random = .met_field_random(gxe),
    residual = ~units,
    data = d,
    weights = d$w,
    fixed_alias = "drop",
    backend = "inla",
    verbose = FALSE
  )
}

test_that("the plan says a main-plot model with one GxE variance fits, and it does", {
  skip_if_no_inla()
  d <- .met_field_data()
  p <- .met_field_plan(d, "TRIAL:VARIETY")

  expect_identical(p$gate_outcome, "accept")
  expect_false(p$preflight_refused)
  expect_true(p$will_fit)
  expect_identical(p$backend_chosen, "inla")

  fit <- .met_field_fit(d, "TRIAL:VARIETY")
  expect_s3_class(fit, "flexybayes")
})

test_that("the plan says a main-plot model with GxE variance per state fits, and it does", {
  skip_if_no_inla()
  d <- .met_field_data()
  p <- .met_field_plan(d, "diag(state_f):tv_f")

  expect_identical(p$gate_outcome, "accept")
  expect_true(p$will_fit)
  expect_identical(p$backend_chosen, "inla")

  fit <- .met_field_fit(d, "diag(state_f):tv_f")
  expect_s3_class(fit, "flexybayes")
})

test_that("the plan and the fit both refuse GxE variance per trial on the hyperparameter budget", {
  skip_if_no_inla()
  d <- .met_field_data()
  p <- .met_field_plan(d, "diag(TRIAL):VARIETY")
  expect_identical(p$gate_outcome, "refuse_structural")
  expect_identical(p$gate_primary_rule, "hyperparam_budget")
  expect_false(p$will_fit)
  # Six per-trial variances take the count from 12 to 17, above the limit of 15.
  expect_match(p$gate_checks[[1L]]$reason, "count = 17 ", fixed = TRUE)

  expect_error(
    .met_field_fit(d, "diag(TRIAL):VARIETY"),
    "hyperparameter count = 17 ",
    fixed = TRUE
  )
})
