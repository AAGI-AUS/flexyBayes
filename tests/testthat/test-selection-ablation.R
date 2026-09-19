# =============================================================================
# fb_select(ablation = TRUE): grouped ablation of the penalised optimiser.
# Each declared block of penalised columns is dropped in turn, the path is
# walked again on the same folds and tuned by the same rule, and each held-out
# unit's loss is paired with its loss under the full path.
#
# The oracles are written in this file and read stored values only:
#
#   The paired loss. Each held-out unit's loss is recomputed from the stored
#   per-row predictions, the full path's in `assessment_predictions` and each
#   ablated path's in `ablation$rows`: the root mean squared residual over
#   the unit's rows in each repeat, averaged over repeats. The paired
#   difference, its mean and its standard error over units are recomputed
#   from those, and compared with `identical()`, because both sides do the
#   same arithmetic on the same stored values in the same order.
#
#   Planted blocks. A block carrying the signal raises the held-out loss by
#   more than two standard errors when it is dropped; a block of pure noise
#   never raises it by more than one.
#
#   The same folds and contract. Every ablated path holds out each row in
#   the same fold as the full path, read row by row from the two stored
#   tables, and carries the full path's ablation contract hash.
#
#   The protected model. Dropping the only declared block leaves the
#   protected columns; the ablated loss of every unit equals the
#   generalised-least-squares protected-only model the benchmark scored on
#   the same folds (`clusters`, model `protected_only`), a separate route.
#
#   The log predictive density. With `refit = list(method = "bayes")`, the
#   paired difference of held-out log predictive density is recomputed from
#   the stored per-row densities of the full and the ablated refits.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

fba_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

# Sites of rows: one protected covariate, six candidates of which the first
# three carry signal and the last three are noise, a site effect and a
# residual.
fba_data <- function(seed = 1L, n_site = 24L, per = 5L) {
  .fb_sel_with_seed(seed, {
    n <- n_site * per
    site <- factor(rep(sprintf("s%02d", seq_len(n_site)), each = per))
    d <- data.frame(site = site, x = stats::rnorm(n))
    for (j in 1:6) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    u <- stats::rnorm(n_site)
    d$y <- 2 + 0.5 * d$x + 1.0 * d$z1 - 0.8 * d$z2 + 0.6 * d$z3 +
      u[as.integer(site)] + stats::rnorm(n)
    d
  })
}

fba_candidates <- ~ z1 + z2 + z3 + z4 + z5 + z6

# Two declared blocks at one set of hyperparameters: the signal and the
# noise.
fba_two_blocks <- function() {
  fb_penalty(
    term(z1, z2, z3) ~ laplace(scale = 0.15),
    term(z4, z5, z6) ~ laplace(scale = 0.15)
  )
}

# A base carrying the parsed model, the data and the call record the fold
# routes read, with no engine behind it.
fba_stub_base <- function(d) {
  fb <- suppressWarnings(fb_from_asreml(y ~ x, random = ~site, data = d))
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = y ~ x, random = ~site, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

fba_select <- function(d,
                       penalty = fba_two_blocks(),
                       candidates = fba_candidates,
                       ablation = TRUE,
                       v = 6L,
                       plan = FALSE,
                       repeats = 1L) {
  fb_select(
    fba_stub_base(d), candidates,
    penalty = penalty,
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = v, repeats = repeats),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 1L,
    plan = plan,
    ablation = ablation
  )
}

# --- the hand oracle ----------------------------------------------------- #

# Each held-out unit's loss from stored per-row predictions: the root mean
# squared residual over the unit's rows in each repeat, rows in row order,
# then the mean over repeats. Units in radix order.
fba_unit_rmse <- function(repeat_id, unit, row, observed, predicted) {
  o <- order(repeat_id, row, method = "radix")
  repeat_id <- repeat_id[o]
  unit <- unit[o]
  residual <- observed[o] - predicted[o]
  cell <- paste(repeat_id, unit, sep = "\r")
  per_cell <- vapply(
    split(residual, factor(cell, levels = unique(cell))),
    function(r) sqrt(mean(r^2)),
    numeric(1L)
  )
  cell_unit <- unit[match(names(per_cell), cell)]
  units <- sort(unique(unit), method = "radix")
  stats::setNames(
    vapply(units, function(u) mean(per_cell[cell_unit == u]), numeric(1L)),
    units
  )
}

fba_full_units <- function(sel) {
  ap <- sel$assessment_predictions
  fba_unit_rmse(
    ap$repeat_id, ap$unit, ap$row,
    ap$observed_model_scale, ap$predicted_model_scale
  )
}

fba_ablated_units <- function(sel, block) {
  rows <- sel$ablation$rows[sel$ablation$rows$block == block, , drop = FALSE]
  fba_unit_rmse(rows$repeat_id, rows$unit, rows$row, rows$observed,
                rows$predicted)
}

# --- (i) planted blocks ---------------------------------------------------- #

# Measured on data seeds 1 to 8 under both rules (probe 02 of the unit's
# evidence): the signal block's difference is 3.6 to 8.0 standard errors. The
# noise block's is at most +0.34 standard errors, and down to -3.2: without
# the noise the path is tuned to less shrinkage and predicts better, so an
# irrelevant block is never read as important but can be read as harmful. The
# test asserts that, on the first two seeds, under the default rule and under
# `rule = "min"`, where the noise block is in the full support and its
# difference is not zero by construction.
test_that("dropping a planted signal block raises the held-out loss beyond two standard errors, and a noise block never raises it beyond one", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  for (rule in c("one_se", "min")) {
    for (seed in 1:2) {
      d <- fba_data(seed = seed)
      sel <- fb_select(
        fba_stub_base(d), fba_candidates,
        penalty = fba_two_blocks(),
        workflow = "penalized_optimizer",
        resamples = fb_resamples(unit = ~site, v = 6L),
        decision = fb_decision(method = "full", metric = "rmse", rule = rule),
        seed = 1L,
        ablation = TRUE
      )
      label <- paste0("rule ", rule, ", seed ", seed)
      tg <- tidy(sel, level = "group")
      expect_setequal(tg$group, c("block_1", "block_2"))
      signal <- tg[tg$group == "block_1", ]
      noise <- tg[tg$group == "block_2", ]
      expect_gt(
        signal$ablation_loss_difference,
        2 * signal$ablation_loss_se,
        label = paste0(label, ": the signal block's loss difference")
      )
      expect_lte(
        noise$ablation_loss_difference,
        noise$ablation_loss_se,
        label = paste0(label, ": the noise block's loss difference")
      )
      if (identical(rule, "min")) {
        full <- sel$assessment_predictions
        held <- as.data.frame(sel$benchmark$glmnet$support$by_lambda)
        held <- held[held$arm == "benchmark" &
          held$lambda_index == full$selected_lambda_index[[1L]], ]
        expect_true(any(grepl("z[456]", held$terms)), label = label)
        expect_gt(noise$ablation_loss_se, 0)
      }
      # The dropped block's columns are in no ablated support.
      sup <- sel$ablation$supports
      expect_false(any(
        sup$column[sup$block == "block_1"] %in% c("z1", "z2", "z3")
      ))
      expect_false(any(
        sup$column[sup$block == "block_2"] %in% c("z4", "z5", "z6")
      ))
    }
  }
})

# --- (ii) the paired differences, by hand ---------------------------------- #

test_that("the paired loss differences equal a hand computation from the stored per-row predictions", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  for (repeats in c(1L, 2L)) {
    sel <- fba_select(fba_data(), repeats = repeats)
    full <- fba_full_units(sel)
    tg <- tidy(sel, level = "group")
    for (block in c("block_1", "block_2")) {
      ablated <- fba_ablated_units(sel, block)
      expect_identical(names(ablated), names(full))
      difference <- ablated - full
      stored <- sel$ablation$units[sel$ablation$units$block == block, ]
      expect_identical(stored$unit, names(full))
      expect_identical(stored$loss_full, unname(full))
      expect_identical(stored$loss_ablated, unname(ablated))
      expect_identical(stored$loss_difference, unname(difference))
      row <- tg[tg$group == block, ]
      expect_identical(row$ablation_loss_difference, mean(difference))
      expect_identical(
        row$ablation_loss_se,
        stats::sd(difference) / sqrt(length(difference))
      )
      expect_identical(row$ablation_n_units, length(difference))
      # The ablated path is tuned by the one-standard-error rule on its own
      # curve, and its loss is read at the point the rule chose.
      curve <- sel$ablation$paths[sel$ablation$paths$block == block, ]
      best <- which.min(curve$score)
      inside <- which(curve$score <= curve$score[[best]] + curve$se[[best]])
      chosen <- inside[which.min(curve$complexity[inside])]
      blocks <- sel$ablation$blocks
      expect_identical(blocks$selected_s[blocks$block == block],
                       curve$s[[chosen]])
      expect_identical(mean(ablated), curve$score[[chosen]])
    }
  }
})

# --- (iii) the same folds and contract ------------------------------------- #

test_that("every ablated path holds out each row in the full path's fold and carries its contract hash", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  sel <- fba_select(fba_data(), repeats = 2L)
  ap <- sel$assessment_predictions
  full_key <- paste(ap$repeat_id, ap$row)
  for (block in c("block_1", "block_2")) {
    rows <- sel$ablation$rows[sel$ablation$rows$block == block, ]
    key <- paste(rows$repeat_id, rows$row)
    expect_setequal(key, full_key)
    expect_identical(rows$fold, ap$fold[match(key, full_key)])
    expect_identical(rows$unit, ap$unit[match(key, full_key)])
  }
  blocks <- sel$ablation$blocks
  expect_true(all(blocks$ran))
  expect_identical(
    blocks$contract_hash,
    rep(sel$ablation$contract$hash, nrow(blocks))
  )
  expect_true(grepl("^[0-9a-f]{32}$", sel$ablation$contract$hash))
  expect_identical(sel$ablation$contract$stage, "ablation")
  expect_identical(sel$ablation$contract$metric, "rmse")
  expect_identical(glance(sel)$ablation_contract_hash,
                   sel$ablation$contract$hash)
})

# --- (iv) the only block, against the protected-only model ----------------- #

test_that("dropping the only declared block leaves the protected model the benchmark scored", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  sel <- fba_select(
    fba_data(),
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15))
  )
  blocks <- sel$ablation$blocks
  expect_identical(blocks$block, "block_1")
  expect_true(blocks$ran)
  expect_true(is.na(blocks$selected_s))
  protected <- as.data.frame(sel$clusters)
  protected <- protected[protected$model == "protected_only", ]
  hand <- vapply(
    split(protected$rmse, protected$cluster),
    mean,
    numeric(1L)
  )
  units <- sel$ablation$units
  expect_identical(units$loss_ablated, unname(hand[units$unit]))
  expect_identical(nrow(sel$ablation$supports), 0L)
})

test_that("a column aliased only with the dropped block is estimable on the path without it", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  d <- fba_data()
  # z7 is exactly the signal of the first block, so every fold's training
  # rows alias it with z1, z2 and z3, written before it.
  d$z7 <- d$z1 - 0.8 * d$z2 + 0.6 * d$z3
  sel <- fba_select(
    d,
    penalty = fb_penalty(
      term(z1, z2, z3) ~ laplace(scale = 0.15),
      term(z4, z5, z6, z7) ~ laplace(scale = 0.15)
    ),
    candidates = ~ z1 + z2 + z3 + z4 + z5 + z6 + z7
  )
  labels <- vapply(sel$folds, function(f) f$label, character(1L))
  dropped <- sel$certificate$design$candidate_space$dropped
  expect_setequal(dropped$fold[dropped$column == "z7"], labels)
  # Without the first block, z7 is identified on every fold and carries the
  # signal, so every fold's chosen support holds it.
  sup <- sel$ablation$supports
  expect_setequal(sup$fold[sup$block == "block_1" & sup$column == "z7"],
                  labels)
})

# --- (v) accessors, print and plan ----------------------------------------- #

test_that("glance, print and tidy state the ablation and that it is not nested", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  sel <- fba_select(fba_data())
  g <- glance(sel)
  expect_identical(g$ablation_blocks, 2L)
  expect_identical(g$ablation_metric, "rmse")
  expect_identical(g$ablation_nested, FALSE)
  expect_false(sel$ablation$nested)
  out <- paste(utils::capture.output(print(sel)), collapse = "\n")
  expect_match(out, "ablation", fixed = TRUE)
  expect_match(out, "not nested", fixed = TRUE)
  expect_match(out, "block_1", fixed = TRUE)
  # No refit, so no density columns are filled.
  tg <- tidy(sel, level = "group")
  expect_true(all(is.na(tg$ablation_elpd_difference)))

  plain <- fba_select(fba_data(), ablation = FALSE)
  expect_null(plain$ablation)
  gp <- glance(plain)
  expect_true(is.na(gp$ablation_blocks))
  expect_true(all(is.na(tidy(plain, level = "group")$ablation_loss_difference)))
})

test_that("plan = TRUE counts the ablation's penalised fits and engine refits", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  pl <- fba_select(fba_data(), plan = TRUE)
  ab <- pl$selection$ablation
  expect_identical(ab$blocks, c("block_1", "block_2"))
  per_path <- .fb_sel_glmnet_decision(NULL)$n_lambda + 1L
  expect_identical(ab$n_folds, 6L)
  expect_identical(ab$path_fits, 2L * 6L * as.integer(per_path))
  expect_identical(ab$protected_fits, 0L)
  expect_identical(ab$engine_refits, 0L)
  expect_identical(ab$covariance_fits, 0L)
  out <- paste(utils::capture.output(print(pl)), collapse = "\n")
  expect_match(out, "Ablation:", fixed = TRUE)

  one <- fba_select(
    fba_data(),
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
    plan = TRUE
  )$selection$ablation
  expect_identical(one$blocks, "block_1")
  expect_identical(one$path_fits, 0L)
  expect_identical(one$protected_fits, 6L)

  none <- fba_select(fba_data(), plan = TRUE, ablation = FALSE)
  expect_null(none$selection$ablation)

  # With the Bayesian refit, one engine refit per fold for the valuation and
  # one more per fold for each block.
  refit <- fb_select(
    fba_stub_base(fba_data()), fba_candidates,
    penalty = fba_two_blocks(),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 6L),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 1L, plan = TRUE,
    refit = list(method = "bayes"), ablation = TRUE
  )$selection
  expect_identical(refit$ablation$engine_refits, 12L)
  expect_identical(refit$cost$n_fits, 6L * 3L + 1L)
  out <- paste(utils::capture.output(print(structure(
    list(selection = refit), class = c("fb_select_plan", "list")
  ))), collapse = "\n")
  expect_match(out, "12 engine refits", fixed = TRUE)
})

test_that("with the Bayesian refit, the ablation's refits are held to the budget before the first", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  # The stub base records no run time, so each refit is costed at one
  # second: 6 folds x (the valuation + 2 blocks) = 18 s, over a budget of
  # 10 s that the valuation alone (6 s) would meet.
  err <- tryCatch(
    fb_select(
      fba_stub_base(fba_data()), fba_candidates,
      penalty = fba_two_blocks(),
      workflow = "penalized_optimizer",
      resamples = fb_resamples(unit = ~site, v = 6L),
      decision = fb_decision(method = "full", metric = "rmse"),
      seed = 1L, budget = 10,
      refit = list(method = "bayes"), ablation = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_penalty_runtime_budget_exceeded")
  expect_match(conditionMessage(err), "18 engine refits", fixed = TRUE)
})

# --- (vi) refusals --------------------------------------------------------- #

test_that("ablation on a route with no penalised path refuses by name before any fit", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  d <- fba_data()
  err <- tryCatch(
    fb_select(
      fba_stub_base(d), fba_candidates,
      workflow = "unpenalized_benchmark",
      resamples = fb_resamples(unit = ~site, v = 6L),
      seed = 1L, ablation = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_selection_ablation_route_undefined")
  err <- tryCatch(
    fb_select(
      fba_stub_base(d), fba_candidates,
      resamples = fb_resamples(unit = ~site, v = 6L),
      seed = 1L, plan = TRUE, ablation = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_selection_ablation_route_undefined")
  expect_error(
    fba_select(d, ablation = "yes"),
    "`ablation` must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("a block whose removal leaves a penalty the path refuses is carried with its reason", {
  skip_if_not_installed("glmnet")
  fba_quiet()
  d <- fba_data()
  d$f <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  sel <- fba_select(
    d,
    penalty = fb_penalty(term(everything()) ~ group_lasso(scale = 0.3)),
    candidates = ~ z1 + f
  )
  blocks <- sel$ablation$blocks
  expect_setequal(blocks$block, c("z1", "f"))
  # Without f, the one remaining group has one column: a grouping that
  # groups nothing, which the path refuses by name.
  expect_false(blocks$ran[blocks$block == "f"])
  expect_match(
    blocks$reason[blocks$block == "f"],
    "selection_penalty_group_structure_absent",
    fixed = TRUE
  )
  expect_true(blocks$ran[blocks$block == "z1"])
  tg <- tidy(sel, level = "group")
  expect_true(is.na(tg$ablation_loss_difference[tg$group == "f"]))
  expect_false(is.na(tg$ablation_loss_difference[tg$group == "z1"]))
  expect_identical(glance(sel)$ablation_blocks, 1L)
})

# --- (vii) the log predictive density beside the loss ---------------------- #

test_that("with the Bayesian refit, the paired elpd difference equals a hand computation from the stored per-row densities", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fba_quiet()
  d <- fba_data(n_site = 12L, per = 6L)
  base <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                     verbose = FALSE)
  sel <- fb_select(
    base, fba_candidates,
    penalty = fba_two_blocks(),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 1L,
    refit = list(method = "bayes", n_draws = 1000L),
    ablation = TRUE
  )
  expect_false(is.null(sel$value))
  full_rows <- sel$value$rows
  unit_elpd <- function(rows) {
    o <- order(rows$repeat_id, rows$row, method = "radix")
    rows <- rows[o, ]
    units <- sort(unique(rows$unit), method = "radix")
    stats::setNames(
      vapply(units, function(u) sum(rows$lpd[rows$unit == u]), numeric(1L)),
      units
    )
  }
  full <- unit_elpd(full_rows)
  tg <- tidy(sel, level = "group")
  for (block in c("block_1", "block_2")) {
    rows <- sel$ablation$rows[sel$ablation$rows$block == block, ]
    expect_false(anyNA(rows$lpd))
    ablated <- unit_elpd(rows)
    expect_identical(names(ablated), names(full))
    difference <- full - ablated
    row <- tg[tg$group == block, ]
    expect_identical(row$ablation_elpd_difference, sum(difference))
    expect_identical(
      row$ablation_elpd_se,
      sqrt(length(difference)) * stats::sd(difference)
    )
    stored <- sel$ablation$value_units[
      sel$ablation$value_units$block == block,
    ]
    expect_identical(stored$unit, names(full))
    expect_identical(stored$elpd_difference, unname(difference))
  }
  # Each ablated refit held the supports its path chose.
  for (block in c("block_1", "block_2")) {
    held <- sel$ablation$value_supports
    held <- held[held$block == block, ]
    chosen <- sel$ablation$supports
    chosen <- chosen[chosen$block == block, ]
    expect_identical(
      sort(paste(held$fold, held$column)),
      sort(paste(chosen$fold, chosen$column))
    )
  }
  # Each ablated refit ran on the value contract of the full valuation.
  expect_identical(
    sel$ablation$blocks$elpd_contract_hash,
    rep(sel$value$contract$hash, 2L)
  )
  # The loss read-out stands beside it.
  expect_false(anyNA(tg$ablation_loss_difference))
})
