# A shared field written id(g):ar1(r):ar1(c) under prediction = "new_cluster".

field_term <- function(at_var = "TRIAL") {
  list(random_terms = list(list(
    type = "ar1_spatial", row_var = "RANGE", col_var = "ROW",
    col_ar1 = TRUE, at_var = at_var
  )))
}

test_that("id(g) makes the field relabelable by g; a field without it still refuses", {
  cols <- flexyBayes:::.fb_sel_grouping_columns(field_term())
  expect_identical(cols$columns, "TRIAL")
  expect_identical(cols$unrelabelled, character(0))

  plain <- field_term(at_var = NULL)
  expect_identical(
    flexyBayes:::.fb_sel_grouping_columns(plain)$unrelabelled,
    "ar1_spatial(RANGE, ROW)"
  )
  expect_error(
    flexyBayes:::.fb_sel_require_relabelable(plain, "new_cluster"),
    class = "flexybayes_refusal_selection_new_cluster_not_representable"
  )
  expect_null(flexyBayes:::.fb_sel_require_relabelable(field_term(), "new_cluster"))
})

test_that("a relabelled trial's field is independent of the trial it came from", {
  d <- expand.grid(ROW = factor(1:4), RANGE = factor(1:3), TRIAL = c("a", "b"))
  d$TRIAL <- as.character(d$TRIAL)
  # Half of trial b held out and relabelled as a new level.
  held <- which(d$TRIAL == "b" & as.integer(d$RANGE) <= 2L)
  d$TRIAL[held] <- "b::fb_heldout"
  d$TRIAL <- factor(d$TRIAL)
  rho <- c(row = 0.6, col = -0.3)

  slot <- flexyBayes:::.fb_sel_slots_ar1("RANGE", "ROW", d, "TRIAL")[[1L]]
  z <- flexyBayes:::.fb_sel_slot_dense(
    slot, nrow(d), list(rho = list(spatial = rho))
  )
  got <- tcrossprod(z)

  # Hand-built: rho_row^|dr| * rho_col^|dc| within one TRIAL level, 0 across.
  r <- as.integer(d$RANGE)
  cc <- as.integer(d$ROW)
  want <- outer(seq_len(nrow(d)), seq_len(nrow(d)), function(i, j) {
    ifelse(
      d$TRIAL[i] == d$TRIAL[j],
      rho[["row"]]^abs(r[i] - r[j]) * rho[["col"]]^abs(cc[i] - cc[j]),
      0
    )
  })
  expect_equal(unname(got), want, tolerance = 1e-12)
  expect_true(all(got[held, d$TRIAL == "b"] == 0))
})

test_that("the study route runs under new_cluster on a base with the shared field", {
  skip_if_not_installed("INLA")
  set.seed(5L)
  d <- expand.grid(ROW = factor(1:6), RANGE = factor(1:4), TRIAL = factor(1:5))
  d$x1 <- stats::rnorm(nrow(d))
  d$x2 <- stats::rnorm(nrow(d))
  d$y <- 0.5 * d$x1 + stats::rnorm(5L)[as.integer(d$TRIAL)] +
    stats::rnorm(nrow(d))
  base <- suppressWarnings(suppressMessages(flexybayes(
    y ~ 1,
    random = ~ TRIAL + id(TRIAL):ar1(RANGE):ar1(ROW),
    data = d,
    backend = "inla",
    verbose = FALSE
  )))
  sel <- suppressWarnings(fb_select(
    base = base,
    candidates = ~ x1 + x2,
    penalty = fb_penalty(term(everything()) ~ elastic_net(lambda = 1, alpha = 1)),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~ TRIAL, v = 5L),
    prediction = "new_cluster",
    decision = fb_decision(
      method = "full", metric = "rmse", rule = "one_se",
      stability_subsamples = 0L
    ),
    covariance = "base_model",
    refit = list(method = "bayes"),
    seed = 1L
  ))
  expect_s3_class(sel, "fb_select")
  expect_length(sel$folds, 5L)
  expect_identical(sel$certificate$covariance$builder, "block_cholesky")
  expect_identical(nrow(sel$scores), 2L)
})
