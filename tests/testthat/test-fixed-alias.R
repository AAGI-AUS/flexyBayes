# =============================================================================
# flexybayes(fixed_alias = ): a fixed-effect column the rows with a response
# do not estimate refuses by default, and with "drop" is fixed at zero as
# asreml 4.2 drops it, with predictions of the cells it leaves unestimated
# returned as NA and status "Aliased".
#
# Fixture: the unbalanced multi-trial design of helper-alias-design.R,
# VARIETY x NAPPLIED x SEEDRATE with empty and lost cells, random TRIAL and
# TRIAL:REP. Oracle: asreml 4.2 on the same design, recorded in
# fixtures/asreml_alias_<shape>.rds (review/units/U-ALIAS_2026-09-24).
# =============================================================================

fa_fixed <- YIELD ~ VARIETY * NAPPLIED * SEEDRATE
fa_random <- ~ TRIAL + TRIAL:REP

fa_oracle <- function(shape) {
  readRDS(test_path("fixtures", paste0("asreml_alias_", shape, ".rds")))
}

fa_cells <- function(d) paste(d$VARIETY, d$NAPPLIED, d$SEEDRATE, sep = "/")

fa_cache <- new.env(parent = emptyenv())

fa_fit <- function(...) {
  key <- paste0("fit", paste(names(list(...)), unlist(list(...)), collapse = "/"))
  if (!exists(key, envir = fa_cache, inherits = FALSE)) {
    withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
    assign(key, flexybayes(
      fa_fixed,
      random = fa_random,
      data = alias_design("small"),
      backend = "inla",
      verbose = FALSE,
      fixed_alias = "drop",
      ...
    ), envir = fa_cache)
  }
  get(key, envir = fa_cache, inherits = FALSE)
}


# --- the default --------------------------------------------------------- #

test_that("flexybayes() refuses an aliased fixed part by default and names asreml's aliased columns", {
  for (shape in c("small", "study")) {
    err <- tryCatch(
      flexybayes(
        fa_fixed,
        random = fa_random,
        data = alias_design(shape),
        backend = "inla",
        verbose = FALSE
      ),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_fixed_effects_aliased")
    expect_identical(sort(err$columns), fa_oracle(shape)$subsets$all$aliased)
    expect_match(conditionMessage(err), "fixed_alias = \"drop\"", fixed = TRUE)
    expect_match(conditionMessage(err), err$columns[[1L]], fixed = TRUE)
  }
})

test_that("a fixed part the rows estimate fits as before, with no alias record", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- alias_design("small")
  fit <- flexybayes(
    YIELD ~ VARIETY + NAPPLIED + SEEDRATE,
    random = fa_random,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  expect_null(fit$extras$fixed_alias)
  expect_null(fit$extras$call_info$fixed_alias)
  expect_null(attributes(predict(fit)))
  expect_false("status" %in% names(predict(fit, classify = "VARIETY")))
})

test_that("fixed_alias = \"drop\" refuses brms and the aggregated route before any fit", {
  d <- alias_design("small")
  expect_error(
    flexybayes(fa_fixed, random = fa_random, data = d, backend = "brms",
               fixed_alias = "drop"),
    "fixed_alias = \"drop\" is fitted on INLA",
    fixed = TRUE
  )
  expect_error(
    flexybayes(fa_fixed, random = fa_random, data = d, backend = "inla",
               aggregate = TRUE, fixed_alias = "drop"),
    "aggregate = TRUE asks for the aggregated one",
    fixed = TRUE
  )
  expect_error(
    flexybayes(fa_fixed, random = fa_random, data = d, fixed_alias = "keep"),
    "should be one of"
  )
  # A model routed to brms after the check refuses in the emit.
  fb <- .build_ir_polymorphic(
    fixed = fa_fixed, random = fa_random, residual = NULL, data = d,
    family = "gaussian", link = NULL, weights = NULL, known_matrices = list(),
    prior = NULL, prior_fixed_sd = 100, prior_vc_sd = 1, syntax = "auto"
  )
  fb$fixed_alias <- .fb_fixed_alias(fb, d)
  expect_error(
    emit_brms(fb = fb, data = d, return_code = TRUE),
    "this model was routed to brms",
    fixed = TRUE
  )
})


# --- "drop", against asreml ---------------------------------------------- #

test_that("fixed_alias = \"drop\" drops asreml's aliased columns and agrees with its estimates and cell predictions", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  a <- fa_oracle("small")$subsets$all
  fit <- fa_fit()
  expect_identical(sort(fit$extras$fixed_alias$dropped$column), a$aliased)
  expect_identical(fit$extras$call_info$fixed_alias, "drop")

  fx <- summary(fit)$fixed
  aliased <- fx[fx$term %in% a$aliased, ]
  expect_true(all(aliased$estimate == 0))
  expect_true(all(is.na(aliased$std.error)))
  expect_true(all(is.na(aliased$conf.low)))
  # Measured: the largest gap is 0.0027 posterior standard deviations.
  est <- fx[match(names(a$fixed), fx$term), ]
  expect_false(anyNA(est$std.error))
  expect_true(all(abs(est$estimate - a$fixed) <= 0.01 * est$std.error))

  p <- predict(fit, classify = "VARIETY:NAPPLIED:SEEDRATE")
  m <- a$predictions[match(fa_cells(p), a$predictions$cell), ]
  expect_identical(p$status, ifelse(m$status == "Estimable", "Estimable", "Aliased"))
  expect_identical(sum(p$status == "Aliased"), 21L)
  expect_true(all(is.na(p$estimate[p$status == "Aliased"])))
  ok <- p$status == "Estimable"
  # Measured: the largest gap is 0.0025 standard errors.
  expect_true(all(abs(p$estimate[ok] - m$predicted[ok]) <= 0.01 * p$std.error[ok]))
})

test_that("predict() on new rows returns NA with status \"Aliased\" for the cells asreml cannot predict", {
  skip_if_not_installed("INLA")
  d <- alias_design("small")
  a <- fa_oracle("small")$subsets$all
  fit <- fa_fit()
  grid <- expand.grid(
    VARIETY = levels(d$VARIETY),
    NAPPLIED = levels(d$NAPPLIED),
    SEEDRATE = levels(d$SEEDRATE)
  )
  pr <- predict(fit, newdata = grid)
  status <- a$predictions$status[match(fa_cells(grid), a$predictions$cell)]
  expect_identical(attr(pr, "status"), ifelse(status == "Estimable", "Estimable", "Aliased"))
  expect_identical(is.na(as.numeric(pr)), status != "Estimable")
  se <- predict(fit, newdata = grid, se.fit = TRUE)
  expect_identical(se$status, attr(pr, "status"))
  expect_identical(unname(is.na(se$se.fit)), status != "Estimable")

  # Every fitted row is estimable.
  fitted_rows <- predict(fit)
  expect_true(all(attr(fitted_rows, "status") == "Estimable"))
  expect_false(anyNA(fitted_rows))
})

test_that("a dropped column is reported as fixed at zero by confint(), vcov() and summary()", {
  skip_if_not_installed("INLA")
  fit <- fa_fit()
  dropped <- fit$extras$fixed_alias$dropped$column
  ci <- confint(fit)
  expect_true(all(is.na(ci[dropped, ])))
  expect_false(anyNA(ci[setdiff(rownames(ci), dropped), ]))
  v <- vcov(fit, n_samples = 200L)
  expect_true(all(v[dropped, ] == 0))
  expect_match(
    paste(utils::capture.output(print(summary(fit))), collapse = "\n"),
    "21 aliased columns fixed at zero",
    fixed = TRUE
  )
})

test_that("fixing a column at zero fits the model without it", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- alias_design("small")
  fit <- fa_fit(prior_fixed_sd = 100)
  dropped <- fit$extras$fixed_alias$dropped$column
  x <- stats::model.matrix(fa_fixed, d)
  kept <- setdiff(colnames(x), c("(Intercept)", dropped))
  for (j in seq_along(kept)) {
    d[[paste0("k", j)]] <- x[, kept[[j]]]
  }
  without <- flexybayes(
    stats::as.formula(paste(
      "YIELD ~", paste0("k", seq_along(kept), collapse = " + ")
    )),
    random = fa_random,
    data = d,
    backend = "inla",
    verbose = FALSE,
    prior_fixed_sd = 100
  )
  # Measured: 2.0e-7.
  expect_equal(
    unname(coef(fit)[c("(Intercept)", kept)]),
    unname(coef(without)[c("(Intercept)", paste0("k", seq_along(kept)))]),
    tolerance = 1e-5
  )
})
