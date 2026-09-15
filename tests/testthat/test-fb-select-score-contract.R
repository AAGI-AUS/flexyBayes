# =============================================================================
# fb_select(): an out-of-fold comparison exists only under one score contract,
# and at one model stage.
#
# The defects this closes. The side-by-side stanza compared the SET of fold
# labels with setequal(), so two different assignments of units to folds
# labelled 1..K passed; it received no unit-to-fold map, no row identity, no
# prediction target and no design signature, yet certified one partition and
# one target; and it placed the Bayesian reference beside the optimiser's
# SELECTED model, two different model stages, in one table.
#
# A score contract now carries the canonical unit-to-fold map, the scored
# rows, the prediction target, the response, the metric, the scoring unit,
# the candidate design, the preprocessing boundary and the model stage, and
# an arm enters a stage table only when its contract hash equals the primary
# arm's. No sampler is involved.
# =============================================================================

fbsc_data <- function() {
  set.seed(31L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(12L)), each = 4L)))
  d$p1 <- stats::rnorm(nrow(d))
  d$x1 <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d))
  d
}

fbsc_contract <- function(folds,
                          data = fbsc_data(),
                          stage = "reference",
                          prediction = "new_cluster",
                          metric = "nlpd",
                          design_columns = c("(Intercept)", "p1", "x1"),
                          learned_on = "analysis_rows",
                          scoring_unit = NULL,
                          candidate_space = NULL) {
  .fb_select_score_contract(
    folds = folds,
    data = data,
    response = "y",
    prediction = prediction,
    metric = metric,
    design_columns = design_columns,
    learned_on = learned_on,
    stage = stage,
    scoring_unit = scoring_unit,
    candidate_space = candidate_space
  )
}

# A per-unit score table of the shape both scorers write, one row per unit.
fbsc_clusters <- function(folds, model = "full", shift = 0) {
  key <- attr(folds, "unit_key")
  do.call(rbind, lapply(folds, function(s) {
    units <- unique(key[s$assess_rows])
    data.frame(
      model = model, fold = s$label, repeat_id = s$repeat_id,
      cluster = units, n = 4L,
      mse = 1 + shift, rmse = 1 + shift,
      elpd = -(seq_along(units) + shift),
      nlpd = seq_along(units) / 10 + shift,
      stringsAsFactors = FALSE
    )
  }))
}

fbsc_folds <- function(seed = 3L, data = fbsc_data()) {
  .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), data, seed = seed)
}


test_that("one computation gives one hash, and the hash is label-invariant", {
  folds <- fbsc_folds()
  a <- fbsc_contract(folds)
  b <- fbsc_contract(folds)
  expect_s3_class(a, "fb_select_score_contract")
  expect_identical(a$hash, b$hash)
  expect_match(a$hash, "^[0-9a-f]{8}-[0-9a-f]{8}$")

  # The same partition under different fold labels and in a different order
  # is the same partition: the map is canonical, not a copy of the labels.
  relabelled <- folds
  for (i in seq_along(relabelled)) {
    relabelled[[i]]$label <- paste0("F", 4L - i)
    relabelled[[i]]$fold_id <- paste0("F", 4L - i)
  }
  relabelled <- structure(
    rev(unclass(relabelled)),
    class = class(folds),
    unit = attr(folds, "unit"), unit_key = attr(folds, "unit_key"),
    n_rows = attr(folds, "n_rows")
  )
  expect_identical(fbsc_contract(relabelled)$hash, a$hash)
})


test_that("two assignments with the same fold labels have different contracts", {
  # The defect, restated as data. Both partitions deal twelve sites into
  # folds labelled 1..3, so a comparison of label sets passes them as one.
  a <- fbsc_folds(seed = 3L)
  b <- fbsc_folds(seed = 4L)
  labels <- function(folds) vapply(folds, function(s) s$label, character(1L))
  expect_true(setequal(labels(a), labels(b)))
  expect_false(identical(
    lapply(a, function(s) s$assess_rows),
    lapply(b, function(s) s$assess_rows)
  ))

  ca <- fbsc_contract(a)
  cb <- fbsc_contract(b)
  expect_false(identical(ca$hash, cb$hash))
  expect_true("unit_to_fold" %in% .fb_select_contract_differences(ca, cb))

  # And the stage table refuses to place them side by side, naming the field.
  table <- .fb_select_stage_table(
    "reference",
    list(
      list(arm = "bayes", model = "full", clusters = fbsc_clusters(a),
           contract = ca),
      list(arm = "gls", model = "full", clusters = fbsc_clusters(b),
           contract = cb)
    ),
    metric = "nlpd"
  )
  expect_false(table$available)
  expect_identical(table$excluded$arm, "gls")
  expect_match(table$reason, "unit_to_fold", fixed = TRUE)
})


test_that("every contract field moves the hash", {
  d <- fbsc_data()
  folds <- fbsc_folds(data = d)
  base <- fbsc_contract(folds, data = d)
  shifted <- d
  shifted$y <- shifted$y + 1

  variants <- list(
    prediction = fbsc_contract(folds, data = d, prediction = "known_cluster"),
    response = fbsc_contract(folds, data = shifted),
    metric = fbsc_contract(folds, data = d, metric = "rmse"),
    design_columns = fbsc_contract(
      folds, data = d, design_columns = c("(Intercept)", "p1")
    ),
    learned_on = fbsc_contract(folds, data = d, learned_on = "all_rows"),
    stage = fbsc_contract(folds, data = d, stage = "decision"),
    scoring_unit = fbsc_contract(folds, data = d, scoring_unit = "fold"),
    candidate_space = fbsc_contract(
      folds, data = d,
      candidate_space = list(dropped = data.frame(
        fold = folds[[1L]]$label, column = "x1", stringsAsFactors = FALSE
      ))
    )
  )
  for (field in names(variants)) {
    expect_false(identical(variants[[field]]$hash, base$hash), info = field)
    expect_identical(
      .fb_select_contract_differences(base, variants[[field]]), field,
      info = field
    )
  }
})


test_that("a stage table is written only when every arm shares the primary's contract", {
  folds <- fbsc_folds()
  contract <- fbsc_contract(folds)
  entries <- list(
    list(arm = "bayes", model = "full", clusters = fbsc_clusters(folds),
         contract = contract),
    list(arm = "gls", model = "full",
         clusters = fbsc_clusters(folds, shift = 0.2), contract = contract)
  )
  table <- .fb_select_stage_table("reference", entries, metric = "nlpd")
  expect_true(table$available)
  expect_setequal(table$arms$arm, c("bayes", "gls"))
  expect_true(all(table$arms$stage == "reference"))
  expect_identical(unique(table$arms$contract_hash), contract$hash)
  expect_true(all(table$arms$n_units == 12L))
  # The scores are re-aggregated from the per-unit tables by one function.
  expect_equal(
    table$arms$score[table$arms$arm == "gls"],
    .fb_sel_aggregate_bayes_scores(
      entries[[2L]]$clusters, metric = "nlpd", unit = "cluster"
    )$score
  )

  # A stage mismatch is a contract mismatch: the reference posterior is not
  # placed beside a selected model.
  selected <- entries
  selected[[2L]]$model <- "benchmark_7"
  selected[[2L]]$clusters <- fbsc_clusters(folds, model = "benchmark_7")
  selected[[2L]]$contract <- fbsc_contract(folds, stage = "decision")
  mixed <- .fb_select_stage_table("reference", selected, metric = "nlpd")
  expect_false(mixed$available)
  expect_match(mixed$reason, "stage", fixed = TRUE)

  # Equal unit COUNTS are not the same units. One arm scores a unit its
  # contract does not name in place of one it does.
  swapped <- entries
  swapped[[2L]]$clusters$cluster[[1L]] <- "S99"
  expect_identical(
    length(unique(swapped[[2L]]$clusters$cluster)), 12L
  )
  swapped_table <- .fb_select_stage_table("reference", swapped, metric = "nlpd")
  expect_false(swapped_table$available)
  expect_match(swapped_table$excluded$reason, "units", fixed = TRUE)

  # A table with no primary score, or with one arm, says why and holds no
  # rows.
  alone <- .fb_select_stage_table("decision", entries[1L], metric = "nlpd")
  expect_false(alone$available)
  expect_null(alone$arms)
  named <- .fb_select_stage_table(
    "decision", entries, metric = "nlpd",
    unavailable = "a posterior threshold defines no predictive distribution"
  )
  expect_false(named$available)
  expect_match(named$reason, "posterior threshold", fixed = TRUE)
})


test_that("the hash is a content checksum that moves with one character", {
  expect_identical(.fb_select_hash("abc"), .fb_select_hash("abc"))
  expect_false(identical(.fb_select_hash("abc"), .fb_select_hash("abd")))
  expect_false(identical(.fb_select_hash("ab"), .fb_select_hash("ba")))
  expect_match(.fb_select_hash(""), "^[0-9a-f]{8}-[0-9a-f]{8}$")
})
