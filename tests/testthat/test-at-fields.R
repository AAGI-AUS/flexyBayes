# test-at-fields.R -- per-level at() terms and ASReml's nugget-free residual.
#
# at(f, level):g gives the random effect g only where f takes the named
# levels, each level with its own variance. at(trial):ar1(row):ar1(col)
# fits a separate AR1 x AR1 field, with its own variance and correlations,
# on each trial; at(trial, level) keeps the named trials. residual =
# ~ ar1(row):ar1(col) and ~ dsum(~ ar1(row):ar1(col) | trial) are the same
# fields with no nugget.
#
# The oracle is ASReml-R 4.2 on simulated trials, fitted once by
# data-raw/at_ar1_reference.R and stored in inst/extdata/. A posterior median
# is compared with the REML estimate, on the SD scale for a variance, and
# passes when the two are within 0.3 posterior SDs. Across the 29 compared
# parameters the largest gap is 0.25 posterior SDs (the 18-level TRIAL:REP
# SD); the planted defects below move a parameter by more than 0.3 posterior
# SDs or change the set of parameters.

.at_ref <- function() {
  path <- system.file("extdata", "at_ar1_reference.rds", package = "flexyBayes")
  if (!nzchar(path)) {
    path <- testthat::test_path(
      "..",
      "..",
      "inst",
      "extdata",
      "at_ar1_reference.rds"
    )
  }
  readRDS(path)
}

# ASReml's component names, as flexyBayes names the same parameters.
.at_asreml_names <- function(term) {
  sfx <- character(length(term))
  at <- grepl("^at\\(", term)
  sfx[at] <- sub("^at\\(([A-Za-z]+), '([^']+)'\\).*$", "__\\1\\2", term[at])
  ds <- grepl("^TRIAL_", term)
  sfx[ds] <- sub("^TRIAL_([^!]+)!.*$", "__TRIAL\\1", term[ds])
  out <- ifelse(
    grepl("!ROW!cor$", term),
    paste0("rho_row", sfx),
    ifelse(
      grepl("!RANGE!cor$", term),
      paste0("rho_col", sfx),
      paste0("sd_spatial", sfx)
    )
  )
  out[term == "units!R"] <- "sigma"
  out[term == "TRIAL:REP"] <- "sd_TRIAL:REP"
  out[term == "at(DESIGN, 'split'):MAINPLOT"] <- "sd_MAINPLOT__DESIGNsplit"
  out
}

# Posterior median against ASReml, in posterior SDs, per parameter. NA for
# a parameter one side has and the other lacks.
.at_oracle_gap <- function(fit, asreml_varcomp) {
  vc <- fit$extras$variance_comps
  med <- attr(vc, "posterior_median")
  cmp <- .at_asreml_names(asreml_varcomp$term)
  est <- ifelse(
    startsWith(cmp, "rho"),
    asreml_varcomp$component,
    sqrt(asreml_varcomp$component)
  )
  all_cmp <- union(cmp, as.character(vc$component))
  gap <- stats::setNames(rep(NA_real_, length(all_cmp)), all_cmp)
  both <- intersect(cmp, as.character(vc$component))
  i <- match(both, as.character(vc$component))
  gap[both] <- (med[both] - est[match(both, cmp)]) / vc$sd[i]
  gap
}

.expect_matches_asreml <- function(fit, asreml_varcomp, tol = 0.3) {
  gap <- .at_oracle_gap(fit, asreml_varcomp)
  expect_false(
    anyNA(gap),
    label = paste(names(gap)[is.na(gap)], collapse = ", ")
  )
  expect_true(
    all(abs(gap) <= tol, na.rm = TRUE),
    label = paste(sprintf("%s %.2f", names(gap), gap), collapse = "; ")
  )
}

.at_fails_oracle <- function(fit, asreml_varcomp, tol = 0.3) {
  gap <- .at_oracle_gap(fit, asreml_varcomp)
  anyNA(gap) || any(abs(gap) > tol)
}

.at_fit <- function(...) {
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE
  )
  suppressMessages(flexybayes(..., backend = "inla"))
}


# --- parsing --------------------------------------------------------------- #

test_that("at(f, level):g keeps the named levels, by value or by position", {
  d <- .at_ref()$data$at_level
  by_name <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ at(DESIGN, "split"):MAINPLOT,
    data = d
  )$random_terms[[1L]]
  by_pos <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ at(DESIGN, 2):MAINPLOT,
    data = d
  )$random_terms[[1L]]
  expect_identical(by_name$outer_levels, "split")
  expect_identical(by_pos$outer_levels, "split")
  expect_identical(by_name$n_outer, 1L)
  every <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ at(DESIGN):MAINPLOT,
    data = d
  )$random_terms[[1L]]
  expect_identical(every$outer_levels, c("rcb", "split"))
})

test_that("a level the data lack, or one not written as a literal, refuses", {
  d <- .at_ref()$data$at_field
  for (random in list(
    ~ at(TRIAL, "T9"):REP,
    ~ at(TRIAL, 4):REP,
    ~ at(TRIAL, paste0("T", 1)):REP,
    ~ at(TRIAL, "T9"):ar1(ROW):ar1(RANGE)
  )) {
    err <- tryCatch(
      fb_from_asreml(fixed = y ~ 1, random = random, data = d),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_at_level_not_in_data")
    expect_match(conditionMessage(err), "T1, T2, T3", fixed = TRUE)
  }
})

test_that("at(trial):ar1(row):ar1(col) is one field per trial", {
  d <- .at_ref()$data$at_field
  term <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ at(TRIAL):ar1(ROW):ar1(RANGE),
    data = d
  )$random_terms[[1L]]
  expect_identical(term$type, "ar1_spatial_at")
  expect_identical(term$by_levels, c("T1", "T2", "T3"))
  expect_true(term$col_ar1)
  one <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ at(TRIAL, "T2"):ar1(ROW):ar1(RANGE),
    data = d
  )$random_terms[[1L]]
  expect_identical(one$by_levels, "T2")
})

test_that("the residual spellings become a field with no nugget", {
  ref <- .at_ref()
  fb <- fb_from_asreml(
    fixed = y ~ 1,
    residual = ~ ar1(ROW):ar1(RANGE),
    data = ref$data$one_trial
  )
  expect_identical(fb$random_terms[[1L]]$type, "ar1_spatial")
  expect_true(fb$random_terms[[1L]]$residual)
  expect_true(fb$residual_terms[[1L]]$nugget_free)
  fb <- fb_from_asreml(
    fixed = y ~ 1,
    residual = ~ dsum(~ ar1(ROW):ar1(RANGE) | TRIAL),
    data = ref$data$dsum_field
  )
  expect_identical(fb$random_terms[[1L]]$type, "ar1_spatial_at")
  expect_identical(fb$random_terms[[1L]]$by_levels, c("T1", "T2", "T3"))
  expect_true(fb$random_terms[[1L]]$residual)
  # The random-side spelling keeps its nugget.
  fb <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ ar1(ROW):ar1(RANGE),
    data = ref$data$one_trial
  )
  expect_null(fb$random_terms[[1L]]$residual)
  expect_null(fb$residual_terms[[1L]]$nugget_free)
})


# --- what reaches INLA ----------------------------------------------------- #

test_that("each trial's field is its own f() term on that trial's plots", {
  skip_if_not_installed("INLA")
  d <- .at_ref()$data$at_field
  code <- suppressMessages(flexybayes(
    y ~ TRIAL,
    random = ~ at(TRIAL, c("T1", "T3")):ar1(ROW):ar1(RANGE),
    data = d,
    backend = "inla",
    return_code = TRUE
  ))
  txt <- gsub("[[:space:]]+", " ", paste(deparse(code$formula), collapse = ""))
  expect_match(
    txt,
    "f(ROW_at_TRIAL_T1_id, model = \"ar1\", group = RANGE_at_TRIAL_T1_id",
    fixed = TRUE
  )
  expect_match(
    txt,
    "f(ROW_at_TRIAL_T3_id, model = \"ar1\", group = RANGE_at_TRIAL_T3_id",
    fixed = TRUE
  )
  expect_false(grepl("TRIAL_T2", txt, fixed = TRUE))
  expect_false(grepl("replicate", txt, fixed = TRUE))
})

test_that("the nugget-free residual fixes the observation precision", {
  skip_if_not_installed("INLA")
  ref <- .at_ref()
  code <- suppressMessages(flexybayes(
    y ~ 1,
    residual = ~ ar1(ROW):ar1(RANGE),
    data = ref$data$one_trial,
    backend = "inla",
    return_code = TRUE
  ))
  expect_true(code$control_family$hyper$prec$fixed)
  expect_gte(code$control_family$hyper$prec$initial, 10)
  nugget <- suppressMessages(flexybayes(
    y ~ 1,
    random = ~ ar1(ROW):ar1(RANGE),
    data = ref$data$one_trial,
    backend = "inla",
    return_code = TRUE
  ))
  expect_null(nugget$control_family$hyper$prec$fixed)
})


# --- ASReml oracle --------------------------------------------------------- #

test_that("the nugget-free residual on one trial matches ASReml", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  ref <- .at_ref()
  fit <- .at_fit(
    y ~ VARIETY,
    residual = ~ ar1(ROW):ar1(RANGE),
    data = ref$data$one_trial
  )
  .expect_matches_asreml(fit, ref$fits$one_trial$varcomp)
  expect_false("sigma" %in% summary(fit)$varcomp$component)
  printed <- paste(utils::capture.output(print(summary(fit))), collapse = "\n")
  expect_match(
    printed,
    "R: ar1(ROW):ar1(RANGE) residual, no nugget",
    fixed = TRUE
  )
  expect_match(
    printed,
    "The residual is the autoregressive field",
    fixed = TRUE
  )
  expect_false(grepl("field + nugget", printed, fixed = TRUE))
  # Planted defect: the random-side field with its nugget is a different
  # model, and the comparison refuses it.
  nugget <- .at_fit(
    y ~ VARIETY,
    random = ~ ar1(ROW):ar1(RANGE),
    data = ref$data$one_trial
  )
  expect_true(.at_fails_oracle(nugget, ref$fits$one_trial$varcomp))
})

test_that("at(TRIAL):ar1(ROW):ar1(RANGE) matches ASReml per trial", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  ref <- .at_ref()
  fit <- .at_fit(
    y ~ TRIAL + VARIETY,
    random = ~ at(TRIAL):ar1(ROW):ar1(RANGE),
    data = ref$data$at_field
  )
  .expect_matches_asreml(fit, ref$fits$at_field$varcomp)
  # Planted defect: one field shared by the three trials, its parameters
  # copied to each trial's name.
  shared <- .at_fit(
    y ~ TRIAL + VARIETY,
    random = ~ id(TRIAL):ar1(ROW):ar1(RANGE),
    data = ref$data$at_field
  )
  vc <- shared$extras$variance_comps
  med <- attr(vc, "posterior_median")
  rows <- lapply(c("T1", "T2", "T3"), function(k) {
    r <- vc[vc$component != "sigma", ]
    r$component <- paste0(r$component, "__TRIAL", k)
    r
  })
  planted <- rbind(vc[vc$component == "sigma", ], do.call(rbind, rows))
  attr(planted, "posterior_median") <- stats::setNames(
    c(med[["sigma"]], rep(med[c("sd_spatial", "rho_row", "rho_col")], 3L)),
    planted$component
  )
  shared$extras$variance_comps <- planted
  expect_true(.at_fails_oracle(shared, ref$fits$at_field$varcomp))
})

test_that("at(TRIAL, \"T2\"):ar1(ROW):ar1(RANGE) matches ASReml", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  ref <- .at_ref()
  fit <- .at_fit(
    y ~ TRIAL + VARIETY,
    random = ~ at(TRIAL, "T2"):ar1(ROW):ar1(RANGE),
    data = ref$data$at_field_lev
  )
  .expect_matches_asreml(fit, ref$fits$at_field_lev$varcomp)
  # Planted defect: the field on every trial.
  every <- .at_fit(
    y ~ TRIAL + VARIETY,
    random = ~ at(TRIAL):ar1(ROW):ar1(RANGE),
    data = ref$data$at_field_lev
  )
  expect_true(.at_fails_oracle(every, ref$fits$at_field_lev$varcomp))
})

test_that("dsum(~ ar1(ROW):ar1(RANGE) | TRIAL) matches ASReml per trial", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  ref <- .at_ref()
  fit <- .at_fit(
    y ~ TRIAL + VARIETY,
    residual = ~ dsum(~ ar1(ROW):ar1(RANGE) | TRIAL),
    data = ref$data$dsum_field
  )
  .expect_matches_asreml(fit, ref$fits$dsum_field$varcomp)
  expect_false("sigma" %in% summary(fit)$varcomp$component)
})

test_that("at(DESIGN, \"split\"):MAINPLOT matches ASReml on INLA", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  ref <- .at_ref()
  fit <- .at_fit(
    y ~ TRIAL + VARIETY,
    random = ~ TRIAL:REP + at(DESIGN, "split"):MAINPLOT,
    data = ref$data$at_level
  )
  .expect_matches_asreml(fit, ref$fits$at_level$varcomp)
  # Planted defect: a main-plot variance in every design.
  every <- .at_fit(
    y ~ TRIAL + VARIETY,
    random = ~ TRIAL:REP + at(DESIGN):MAINPLOT,
    data = ref$data$at_level
  )
  expect_true(.at_fails_oracle(every, ref$fits$at_level$varcomp))
})

test_that("at(DESIGN, \"split\"):MAINPLOT matches ASReml on brms", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  ref <- .at_ref()
  # posterior caps an effective sample size above the number of draws; that
  # note is the only warning allowed through.
  fit <- withCallingHandlers(
    suppressMessages(flexybayes(
      y ~ TRIAL + VARIETY,
      random = ~ TRIAL:REP + at(DESIGN, "split"):MAINPLOT,
      data = ref$data$at_level,
      backend = "brms",
      seed = 1L,
      chains = 2L,
      n_samples = 1000L,
      warmup = 1000L
    )),
    warning = function(w) {
      if (grepl("ESS has been capped", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  .expect_matches_asreml(fit, ref$fits$at_level$varcomp)
  pred <- predict(fit, newdata = ref$data$at_level[1:4, ])
  expect_length(pred, 4L)
})


# --- refusals -------------------------------------------------------------- #

test_that("fields refuse on brms; a residual field refuses beside another", {
  skip_if_not_installed("brms")
  ref <- .at_ref()
  err <- tryCatch(
    flexybayes(
      y ~ 1,
      random = ~ at(TRIAL):ar1(ROW):ar1(RANGE),
      data = ref$data$at_field,
      backend = "brms"
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_stan_cannot_represent_ar1_field")
  err <- tryCatch(
    flexybayes(
      y ~ 1,
      residual = ~ dsum(~ ar1(ROW):ar1(RANGE) | TRIAL),
      data = ref$data$dsum_field,
      backend = "brms"
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_stan_cannot_represent_ar1_residual")
  err <- tryCatch(
    fb_from_asreml(
      fixed = y ~ 1,
      random = ~ ar1(ROW):ar1(RANGE),
      residual = ~ ar1(ROW):ar1(RANGE),
      data = ref$data$one_trial
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_ar1_residual_with_random_field")
})

test_that("the nugget-free residual refuses off the Gaussian identity link", {
  skip_if_not_installed("INLA")
  d <- .at_ref()$data$one_trial
  d$count <- stats::rpois(nrow(d), 3)
  err <- tryCatch(
    suppressMessages(flexybayes(
      count ~ 1,
      residual = ~ ar1(ROW):ar1(RANGE),
      family = "poisson",
      data = d,
      backend = "inla"
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_ar1_residual_not_representable")
  expect_match(conditionMessage(err), "Gaussian response", fixed = TRUE)
})

test_that("a residual SD prior refuses on the nugget-free residual", {
  skip_if_not_installed("INLA")
  d <- .at_ref()$data$one_trial
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ 1,
      residual = ~ ar1(ROW):ar1(RANGE),
      prior = fb_prior(sigma() ~ half_normal(1)),
      data = d,
      backend = "inla"
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_prior_target_not_in_model")
})

test_that("a trial whose plots do not fill its grid refuses", {
  skip_if_not_installed("INLA")
  d <- .at_ref()$data$at_field
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ TRIAL,
      random = ~ at(TRIAL):ar1(ROW):ar1(RANGE),
      data = d[-5L, ],
      backend = "inla",
      na_action = "omit"
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_ar1_spatial_requires_complete_grid")
  expect_match(conditionMessage(err), "level(s) T1", fixed = TRUE)
})
