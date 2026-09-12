# =============================================================================
# The projection predictive decision.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 6.1.
#
# The arithmetic and the provider are tested apart, because they fail
# differently and only one of them costs a sampler.
#
#   the arithmetic  fold conversion, unit-level aggregation and the three
#                   decision rules, on a constructed table whose answer is
#                   known in closed form. The rules are shown to DISAGREE
#                   on the same table -- a one-standard-error rule that
#                   silently returned the best size would pass a test that
#                   only ever saw one rule
#   the provider    one validated grouped k-fold run, which is where the
#                   forced protected block, the rank stability, the
#                   labelled projected draws and the unit-level standard
#                   error are read off a real search
# =============================================================================


# --- a table whose answer is known -------------------------------------- #

# Ten units of three rows each. The unit totals are chosen so that the best
# size is 4, size 2 is inside one standard error of the paired difference
# against it, and size 3 is not (its differences are constant, so their
# spread is zero and any loss is outside it).
fb_sel_fake_summaries <- function() {
  n_units <- 10L
  rows_per_unit <- 3L
  unit_key <- rep(sprintf("u%02d", seq_len(n_units)), each = rows_per_unit)
  base_vals <- seq(-6, -4, length.out = n_units)
  diffs <- c(-1.1, -0.6, -0.1, 0.2, 0.4, 0.6, -0.4, 0.9, 0.7, 0.4)

  as_lppd <- function(unit_values) rep(unit_values / rows_per_unit, each = rows_per_unit)
  list(
    summaries = list(
      sub = list(
        list(lppd = as_lppd(base_vals - 5.0)),
        list(lppd = as_lppd(base_vals - 3.0)),
        list(lppd = as_lppd(base_vals - diffs)),
        list(lppd = as_lppd(base_vals - 0.05)),
        list(lppd = as_lppd(base_vals))
      ),
      ref = list(lppd = as_lppd(base_vals + 0.2))
    ),
    unit_key = unit_key,
    diffs = diffs,
    n_units = n_units
  )
}


test_that("elpd and its standard error are aggregated at the unit, not the row", {
  fake <- fb_sel_fake_summaries()
  agg <- .fb_sel_unit_elpd(fake$summaries, fake$unit_key)

  expect_identical(agg$n_units, fake$n_units)
  expect_identical(nrow(agg$table), 5L)
  expect_identical(agg$table$size, 0:4)
  expect_identical(agg$n_nonfinite, 0L)

  # The unit totals, and the standard error of their sum: sd across units
  # times the square root of their number.
  best <- agg$unit_matrix[, 5L]
  expect_equal(agg$table$elpd[[5L]], sum(best), tolerance = 1e-10)
  expect_equal(
    agg$table$se[[5L]],
    stats::sd(best) * sqrt(fake$n_units),
    tolerance = 1e-10
  )
  # The difference to the reference is PAIRED before its spread is taken.
  paired <- best - agg$reference_units
  expect_equal(agg$table$delta_to_ref[[5L]], sum(paired), tolerance = 1e-10)
  expect_equal(
    agg$table$delta_se[[5L]],
    stats::sd(paired) * sqrt(fake$n_units),
    tolerance = 1e-10
  )
})


test_that("a non-finite pointwise density is counted, not dropped", {
  fake <- fb_sel_fake_summaries()
  fake$summaries$sub[[3L]]$lppd[[2L]] <- -Inf
  agg <- .fb_sel_unit_elpd(fake$summaries, fake$unit_key)
  expect_identical(agg$n_nonfinite, 1L)
  expect_false(is.finite(agg$table$elpd[[3L]]))
})


test_that("the three decision rules disagree, each in its own direction", {
  fake <- fb_sel_fake_summaries()
  agg <- .fb_sel_unit_elpd(fake$summaries, fake$unit_key)

  min_rule <- .fb_sel_decide(agg, .fb_sel_decision_spec(list(rule = "min")))
  expect_identical(min_rule$size, 4L)

  one_se <- .fb_sel_decide(agg, .fb_sel_decision_spec(list(rule = "one_se")))
  expect_identical(one_se$size, 2L)
  # And the rule is the paired one: the loss against the best size is
  # inside the spread of that paired difference, and outside the raw
  # standard error the size's own elpd carries.
  loss <- sum(agg$unit_matrix[, 5L] - agg$unit_matrix[, 3L])
  se_paired <- stats::sd(agg$unit_matrix[, 5L] - agg$unit_matrix[, 3L]) *
    sqrt(fake$n_units)
  expect_lte(loss, se_paired)

  tol <- .fb_sel_decide(
    agg,
    .fb_sel_decision_spec(list(rule = "tolerance", tolerance = 0.28))
  )
  expect_identical(tol$size, 3L)
  expect_false(identical(tol$size, one_se$size))
})


test_that("a size below the forced block is not an admissible answer", {
  # The injected defect: the same table, decided with no floor. The rule
  # then picks a submodel that has dropped a protected term.
  fake <- fb_sel_fake_summaries()
  agg <- .fb_sel_unit_elpd(fake$summaries, fake$unit_key)
  agg$table$elpd[[1L]] <- max(agg$table$elpd) + 100

  unfloored <- .fb_sel_decide(agg, .fb_sel_decision_spec(list(rule = "min")))
  expect_identical(unfloored$size, 0L)

  floored <- .fb_sel_decide(
    agg,
    .fb_sel_decision_spec(list(rule = "min")),
    floor_size = 1L
  )
  expect_gte(floored$size, 1L)
})


test_that("the tolerance rule refuses to be defaulted", {
  expect_error(
    .fb_sel_decision_spec(list(rule = "tolerance")),
    "tolerance"
  )
  expect_error(
    .fb_sel_decision_spec(list(rule = "tolerance", tolerance = -1)),
    "tolerance"
  )
})


# --- folds --------------------------------------------------------------- #

test_that("an <fb_sel_folds> becomes one fold number per row", {
  d <- fb_sel_bayes_data()
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), d, seed = 5L)
  fold_vector <- .fb_sel_projpred_folds(folds)

  expect_length(fold_vector, nrow(d))
  expect_identical(sort(unique(fold_vector)), 1:3)
  # Every row of a cluster is in the same fold: the unit is the cluster,
  # and a fold that cut through one would leak.
  by_site <- tapply(fold_vector, d$site, function(f) length(unique(f)))
  expect_true(all(by_site == 1L))
  expect_identical(attr(fold_vector, "unit"), "site")
  expect_length(attr(fold_vector, "unit_key"), nrow(d))
})


test_that("a repeated resampling scheme refuses rather than being flattened", {
  d <- fb_sel_bayes_data()
  folds <- .fb_sel_folds(
    fb_resamples(unit = ~site, v = 3L, repeats = 2L),
    d,
    seed = 5L
  )
  expect_error(.fb_sel_projpred_folds(folds), "one partition")
})


test_that("the projection refuses without a seed", {
  skip_if_not_installed("projpred")
  expect_error(
    .fb_sel_project(structure(list(), class = "list"), folds = 1L, seed = NULL),
    class = "flexybayes_refusal_selection_seed_required"
  )
})


# --- one real search ------------------------------------------------------ #

fb_sel_proj_bits <- function() {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("projpred")
  d <- fb_sel_bayes_data()
  cand <- c("x1", "x2", "z1", "z2")
  comp <- .fb_sel_test_blocks(cand, family = "regularized_horseshoe", p0 = 1)
  base <- fb_sel_bayes_memo("base", fb_sel_bayes_base(d))
  # The sampler's own warnings about a handful of divergent transitions at
  # this deliberately small budget are what the convergence gate reads and
  # reports; they are not the subject of these tests.
  # The memo cache is shared across every file that sources
  # helper-selection-bayes.R, so the key names the FIT rather than its role:
  # test-selection-reference.R builds a different reference (2500 draws,
  # 1000 warm-up) and under a plain "reference" key whichever file ran first
  # handed its fit to the other. The adequacy assertions there then read a
  # posterior they had not asked for.
  ref <- fb_sel_bayes_memo("reference_projection", suppressWarnings(.fb_sel_reference_fit(
    base, ~ x1 + x2 + z1 + z2,
    prior = comp$prior,
    column_map = .fb_sel_test_column_map(cand),
    data = d,
    engine = fb_engine("brms", chains = 2L, n_samples = 1500L, warmup = 750L),
    seed = 20260912L
  )))
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), d, seed = 5L)
  list(d = d, cand = cand, ref = ref, folds = folds)
}

fb_sel_proj_run <- function(bits) {
  fb_sel_bayes_memo("projection_search", suppressWarnings(.fb_sel_project(
    bits$ref, bits$folds,
    decision = list(rule = "one_se", validate_search = TRUE),
    seed = 20260912L,
    nclusters = 8L, nclusters_pred = 20L, ndraws_project = 100L
  )))
}


test_that("a validated grouped k-fold search runs on the caller's own folds", {
  bits <- fb_sel_proj_bits()
  projection <- fb_sel_proj_run(bits)

  expect_true(projection$settings$validate_search)
  expect_false(projection$settings$optimistic)
  expect_identical(projection$settings$cv_method, "kfold")
  expect_identical(projection$settings$K, 3L)
  expect_identical(projection$settings$search, "forward")
  expect_identical(projection$settings$rule, "one_se")
  expect_match(projection$settings$projpred_version, "^2\\.")

  # The scoring unit is the resampling unit, not the row and not the fold.
  expect_identical(projection$scoring$unit, "site")
  expect_identical(projection$scoring$n_units, 12L)
  expect_identical(projection$scoring$n_nonfinite_lppd, 0L)
})


test_that("the protected block is forced into every submodel", {
  bits <- fb_sel_proj_bits()
  projection <- fb_sel_proj_run(bits)

  expect_identical(projection$forced_terms, c("p1", "(1 | site)"))
  block <- paste(projection$forced_terms, collapse = " + ")
  # Every search path starts from the block: the full-data ranking, and
  # every fold's own ranking.
  expect_identical(projection$ranking[[1L]], block)
  expect_true(all(projection$foldwise_rankings[, 1L] == block))
  # And the chosen submodel carries the protected terms whatever its size.
  expect_true(all(c("p1", "(1 | site)") %in% projection$chosen_terms))
  expect_gte(projection$chosen_size, 1L)
})


test_that("the elpd table carries both rulers and the unit-level one is wider", {
  bits <- fb_sel_proj_bits()
  projection <- fb_sel_proj_run(bits)
  tab <- projection$elpd_by_size

  expect_true(all(
    c("size", "elpd", "se", "delta_to_ref", "delta_se", "se_row_level") %in%
      names(tab)
  ))
  expect_identical(tab$size, seq_len(nrow(tab)) - 1L)
  expect_true(all(is.finite(tab$elpd)))
  # The point of aggregating at the cluster: treating six correlated rows
  # as six independent replicates understates the spread. The chosen row
  # is the one a decision is actually read off.
  row_at_choice <- tab[tab$size == projection$chosen_size, ]
  expect_gt(row_at_choice$se, row_at_choice$se_row_level)
})


test_that("the projected draws are labelled and are not the reference posterior", {
  bits <- fb_sel_proj_bits()
  projection <- fb_sel_proj_run(bits)
  draws <- projection$projected_draws

  expect_identical(draws$label, "projected_parameters")
  expect_identical(attr(draws$draws, "fb_draw_label"), "projected_parameters")
  expect_true(is.matrix(draws$draws))
  expect_identical(nrow(draws$draws), 100L)
  # A projected submodel carries the terms it kept and no others.
  expect_true("b_p1" %in% colnames(draws$draws))
  dropped <- setdiff(bits$cand, projection$chosen_terms)
  if (length(dropped) > 0L) {
    expect_false(any(paste0("b_", dropped) %in% colnames(draws$draws)))
  }
})


test_that("the rank stability table reports where the folds disagreed", {
  bits <- fb_sel_proj_bits()
  projection <- fb_sel_proj_run(bits)
  stability <- projection$rank_stability

  expect_false(is.null(stability))
  expect_true(is.matrix(unclass(stability)))
  # Proportions over the folds: each row sums to one.
  expect_equal(unname(rowSums(unclass(stability))), rep(1, nrow(stability)),
    tolerance = 1e-8
  )
  expect_identical(nrow(projection$foldwise_rankings), 3L)
})


test_that("the search is deterministic under its seed", {
  # Run with the search unvalidated, which is the cheap mode and is
  # labelled optimistic on the object; determinism is a property of the
  # seed, not of the validation.
  #
  # This is the test that found the defect it now guards. The reference
  # model was constructed OUTSIDE the seeded region, and `get_refmodel()`
  # draws from R's stream itself: three successive runs of the same seeded
  # call returned the same ranking and elpd values that differed by up to a
  # tenth of a nat. A decision rule reading a one-standard-error margin off
  # that column is reading noise.
  bits <- fb_sel_proj_bits()
  once <- suppressWarnings(.fb_sel_project(
    bits$ref, bits$folds,
    decision = list(rule = "min", validate_search = FALSE),
    seed = 4242L, nclusters = 5L, nclusters_pred = 10L, ndraws_project = 50L
  ))
  twice <- suppressWarnings(.fb_sel_project(
    bits$ref, bits$folds,
    decision = list(rule = "min", validate_search = FALSE),
    seed = 4242L, nclusters = 5L, nclusters_pred = 10L, ndraws_project = 50L
  ))

  expect_identical(once$ranking, twice$ranking)
  expect_equal(once$elpd_by_size$elpd, twice$elpd_by_size$elpd)
  expect_identical(once$chosen_size, twice$chosen_size)
  expect_equal(once$projected_draws$draws, twice$projected_draws$draws)
  # The fast mode says what it is.
  expect_true(once$settings$optimistic)
  expect_false(once$settings$validate_search)
})


test_that("an unseen cluster is predicted with more spread than a known one", {
  # Section 6.3's last paragraph, executed: the new-cluster target draws
  # the held-out cluster's effect from the estimated population, and the
  # known-cluster target conditions on the fitted level.
  bits <- fb_sel_proj_bits()
  rows <- bits$d[bits$d$site == levels(bits$d$site)[[1L]], , drop = FALSE]

  known <- .fb_sel_predict_cluster(
    bits$ref, rows, prediction = "known_cluster", ndraws = 400L, seed = 11L
  )
  unseen <- rows
  unseen$site <- factor("S99")
  new <- .fb_sel_predict_cluster(
    bits$ref, unseen, prediction = "new_cluster", ndraws = 400L, seed = 11L
  )

  expect_identical(dim(known), dim(new))
  expect_gt(mean(apply(new, 2L, stats::sd)), mean(apply(known, 2L, stats::sd)))
  # A new level without the flag is refused by brms, which is why the
  # target has to be declared rather than inferred.
  expect_error(
    .fb_sel_predict_cluster(bits$ref, unseen, prediction = "known_cluster"),
    "S99"
  )
})
