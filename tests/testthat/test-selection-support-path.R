# =============================================================================
# The support history of the penalised optimiser: which columns the path holds
# at each tuning point of each fold, where each one first enters, where it last
# leaves, how often it comes back, and whether the fold could estimate it at
# all.
#
# The oracle is a hand-whitened `glmnet` path built in this file from base R
# alone: forwardsolve against the fold's Cholesky factor, `qr.resid()` for the
# Frisch-Waugh-Lovell residualisation, and `glmnet::glmnet()` at the canonical
# map. The provider is `glmnet`, so the standard is exact reproduction of the
# support at every point rather than a tolerance; the coefficients themselves
# are the provider's and are not re-derived here.
#
# What each block holds:
#
#   Reproduction. The history's membership at every benchmark point of every
#   fold equals the hand path's, for a lasso and for an elastic net.
#
#   The coordinate. Each point's canonical lambda is that fold's own
#   training-row `lambda_max` times one data-free `s` sequence, so a history
#   read on `s` compares the same fraction of each fold's null anchor. A grid
#   anchored on pooled rows breaks this block.
#
#   Entry, exit and re-entry. Read off the membership sequence, against a
#   hand-written path in which a column enters, leaves and comes back.
#
#   Availability. A column a fold's training rows cannot identify is reported
#   unavailable for that fold, which is not the same statement as a column
#   that was estimable and never entered.
#
#   The accessors. `tidy()` carries the per-column summary and the history
#   beside it; `glance()` carries the coordinate, the alpha, the chosen point
#   and the size of the path.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

# Ten sites of six rows: one protected covariate, five candidates of which two
# carry signal, and a site random effect the fold covariance whitens.
fbsp_data <- function(n_cluster = 10L, per = 6L, seed = 51L) {
  .fb_sel_with_seed(seed, {
    n <- n_cluster * per
    site <- factor(rep(seq_len(n_cluster), each = per))
    u <- stats::rnorm(n_cluster, sd = 1.4)
    d <- data.frame(site = site, x = stats::rnorm(n))
    for (j in seq_len(5L)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    d$y <- 3 + 0.8 * d$x + 1.1 * d$z1 - 0.9 * d$z2 +
      u[as.integer(site)] + stats::rnorm(n, sd = 0.5)
    d
  })
}

fbsp_declared <- paste0("z", 1:5)

fbsp_cv <- function(d = fbsp_data(), v = 4L) {
  .fb_sel_cv(
    d, "y", ~x, ~ z1 + z2 + z3 + z4 + z5, random = "site",
    resamples = fb_resamples(unit = ~site, v = v), seed = 77L
  )
}

fbsp_blocks <- function(alpha) {
  if (identical(alpha, 1)) {
    list(list(
      family = "laplace", columns = 1:5, hyperparameters = list(scale = 0.15)
    ))
  } else {
    list(list(
      family = "elastic_net", columns = 1:5,
      hyperparameters = list(lambda = 0.2, alpha = alpha)
    ))
  }
}

# Six locations of two trials of eight plots with trial-constant candidates:
# the shape in which a fold's training rows cannot identify every declared
# column, so availability and selection are different statements.
fbsp_drop_data <- function(seed = 11L) {
  .fb_sel_with_seed(seed, {
    n_loc <- 6L
    trials_per <- 2L
    plots <- 8L
    n_trial <- n_loc * trials_per
    trial_loc <- rep(sprintf("L%d", seq_len(n_loc)), each = trials_per)
    d <- data.frame(
      trial = factor(rep(sprintf("T%02d", seq_len(n_trial)), each = plots)),
      location = factor(rep(trial_loc, each = plots))
    )
    n <- nrow(d)
    tr <- as.integer(d$trial)
    for (j in 1:9) {
      d[[paste0("c", j)]] <- stats::rnorm(n_trial)[tr]
    }
    d$x <- stats::rnorm(n)
    d$w <- stats::rnorm(n)
    d$y <- 1 + 0.5 * d$x + 0.4 * d$c1 +
      stats::rnorm(n_trial, sd = 0.5)[tr] +
      stats::rnorm(n_loc, sd = 0.3)[as.integer(d$location)] +
      stats::rnorm(n, sd = 0.4)
    d
  })
}

# A base carrying the parsed model, the data and the call record the fold
# routes read, with no engine behind it.
fbsp_stub_base <- function(d) {
  fb <- fb_from_asreml(y ~ x, random = ~site, data = d)
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

fbsp_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

# The fold's whitened, protected-residualised problem, from base R alone: the
# oracle's own reduction, sharing no code with the one under test.
fbsp_hand_payload <- function(substrate, columns) {
  train_vec <- substrate$train_rows
  lower <- t(substrate$blocks$factor_tt)
  whiten <- function(m) forwardsolve(lower, m)
  candidate_vec <- colnames(substrate$x_candidate)[columns]
  kept_vec <- candidate_vec[candidate_vec %in% substrate$candidate_kept]
  x0 <- whiten(substrate$x_protected[train_vec, , drop = FALSE])
  xp <- whiten(substrate$x_candidate[train_vec, kept_vec, drop = FALSE])
  yw <- drop(whiten(matrix(substrate$response[train_vec], ncol = 1L)))
  qr0 <- qr(x0)
  xt <- qr.resid(qr0, xp)
  colnames(xt) <- kept_vec
  list(x = xt, y = qr.resid(qr0, yw), n = length(train_vec))
}

# The support of the hand path at one canonical lambda. The canonical-to-
# provider map is written out here: the L1 part is divided by the row count,
# the L2 part carries the uncentred root mean square of the response once.
fbsp_hand_support <- function(payload, lambda_canonical, alpha) {
  mixture <- alpha + (1 - alpha) * sqrt(mean(payload$y^2))
  arguments <- list(
    x = payload$x, y = payload$y, family = "gaussian",
    alpha = alpha / mixture,
    lambda = lambda_canonical * mixture / payload$n,
    standardize = FALSE, intercept = FALSE
  )
  if ("control" %in% names(formals(glmnet::glmnet))) {
    arguments$control <- list(thresh = 1e-16, maxit = 1e6)
  } else {
    arguments$thresh <- 1e-16
    arguments$maxit <- 1e6
  }
  fit <- do.call(glmnet::glmnet, arguments)
  beta <- as.numeric(stats::coef(fit))[-1L]
  sort(colnames(payload$x)[abs(beta) > 0])
}

# The membership of one fold and one column, read off the history.
fbsp_membership <- function(history, fold, column) {
  rows <- history[history$fold == fold & history$column == column, ]
  rows <- rows[order(rows$lambda_index), ]
  rows$in_support
}


# --- (i) the hand-whitened path ------------------------------------------ #

test_that("the support history reproduces a hand-whitened glmnet path at every tuning point", {
  skip_if_not_installed("glmnet")
  cv <- fbsp_cv()

  for (alpha in c(1, 0.6)) {
    path <- .fb_sel_glmnet_path(
      cv$substrates, fbsp_blocks(alpha), folds = cv$folds,
      decision = list(n_lambda = 8L), seed = 5L
    )
    history <- path$support_by_resolution$history
    expect_identical(unique(history$alpha), alpha)

    compared <- 0L
    for (label in names(cv$substrates)) {
      payload <- fbsp_hand_payload(cv$substrates[[label]], 1:5)
      rows <- history[history$fold == label, ]
      points <- sort(unique(rows$lambda_index))
      for (index in points) {
        at_point <- rows[rows$lambda_index == index, ]
        expect_identical(nrow(at_point), length(fbsp_declared))
        package_support <- sort(at_point$column[at_point$in_support])
        hand <- fbsp_hand_support(
          payload, at_point$lambda_canonical[[1L]], alpha
        )
        expect_identical(package_support, hand)
        compared <- compared + 1L
      }
    }
    # The comparison is over the whole grid of every fold, so a history that
    # held one point could not pass it.
    expect_identical(compared, 8L * length(cv$substrates))
  }

  # The spread of the check under its own conditions: the hand path computed
  # twice on identical input returns the same support, so an inexact
  # comparison would be measuring the solver rather than the history.
  payload <- fbsp_hand_payload(cv$substrates[[1L]], 1:5)
  expect_identical(
    fbsp_hand_support(payload, 4, 1),
    fbsp_hand_support(payload, 4, 1)
  )
})


# --- (ii) the tuning coordinate ------------------------------------------ #

test_that("each tuning point is the fold's own lambda_max times one data-free s sequence", {
  skip_if_not_installed("glmnet")
  cv <- fbsp_cv()
  path <- .fb_sel_glmnet_path(
    cv$substrates, fbsp_blocks(1), folds = cv$folds,
    decision = list(n_lambda = 8L), seed = 5L
  )
  history <- path$support_by_resolution$history
  lambda_max <- path$selection$lambda_max
  expect_identical(path$selection$tuning_scale, "s = lambda / lambda_max")

  sequences <- list()
  for (label in names(cv$substrates)) {
    rows <- history[history$fold == label & history$column == "z1", ]
    rows <- rows[order(rows$lambda_index), ]
    # Exact: the grid is that fold's anchor times the sequence, computed by
    # the same arithmetic on the same values, so any difference is a change
    # of the rows the anchor was read on.
    expect_identical(rows$lambda_canonical, lambda_max[[label]] * rows$s)
    sequences[[label]] <- rows$s
  }
  # One `s` on every fold, and its lambda differs between folds.
  expect_identical(length(unique(sequences)), 1L)
  expect_identical(sequences[[1L]][[1L]], 1)
  expect_gt(length(unique(as.numeric(lambda_max))), 1L)
})


# --- (iii) entry, exit and re-entry -------------------------------------- #

test_that("first entry, last exit and re-entries are read off the membership", {
  # A path written by hand: `a` enters, leaves and comes back; `b` enters and
  # stays; `c` is estimable and never enters; `d` is not estimable in this
  # fold. The fixed arm's row is on the table and must not enter the history,
  # which is the tuning grid.
  supports <- data.frame(
    arm = c("fixed", rep("benchmark", 6L)),
    lambda_index = c(1L, 1:6),
    lambda_canonical = c(5, 10, 8, 6, 4, 2, 1),
    s = c(0.5, 1, 0.8, 0.6, 0.4, 0.2, 0.1),
    fold = "f1",
    n_nonzero = c(3L, 0L, 1L, 2L, 1L, 2L, 2L),
    effective_df = 0,
    terms = c("a, b, c", "", "a", "a, b", "b", "a, b", "a, b"),
    stringsAsFactors = FALSE
  )
  built <- .fb_sel_glmnet_support_history(
    supports,
    penalised_by_fold = list(f1 = c("a", "b", "c")),
    declared_penalised = c("a", "b", "c", "d"),
    alpha = 1,
    selects = TRUE
  )

  expect_identical(nrow(built$points), 4L * 6L)
  expect_identical(sort(unique(built$points$lambda_index)), 1:6)
  expect_identical(
    fbsp_membership(built$points, "f1", "a"),
    c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
  )

  by_fold <- built$by_fold
  rownames(by_fold) <- by_fold$column
  expect_identical(by_fold["a", "first_entry_s"], 0.8)
  expect_identical(by_fold["a", "last_exit_s"], 0.4)
  expect_identical(by_fold["a", "n_reentries"], 1L)
  expect_identical(by_fold["a", "n_in_support"], 4L)
  expect_true(by_fold["a", "available"])

  expect_identical(by_fold["b", "first_entry_s"], 0.6)
  expect_identical(by_fold["b", "last_exit_s"], NA_real_)
  expect_identical(by_fold["b", "n_reentries"], 0L)

  expect_identical(by_fold["c", "n_in_support"], 0L)
  expect_identical(by_fold["c", "first_entry_s"], NA_real_)
  expect_true(by_fold["c", "available"])

  expect_false(by_fold["d", "available"])
  expect_identical(by_fold["d", "n_in_support"], 0L)
  expect_identical(by_fold["d", "n_points"], 6L)

  # The same summary, recomputed from the membership of a real path by an
  # independent route: a run-length encoding of the sequence.
  skip_if_not_installed("glmnet")
  cv <- fbsp_cv()
  path <- .fb_sel_glmnet_path(
    cv$substrates, fbsp_blocks(1), folds = cv$folds,
    decision = list(n_lambda = 8L), seed = 5L
  )
  history <- path$support_by_resolution$history
  summary_rows <- path$support_by_resolution$history_by_fold
  checked <- 0L
  for (i1 in seq_len(nrow(summary_rows))) {
    label <- summary_rows$fold[[i1]]
    column <- summary_rows$column[[i1]]
    rows <- history[history$fold == label & history$column == column, ]
    rows <- rows[order(rows$lambda_index), ]
    runs <- rle(rows$in_support)
    starts <- cumsum(c(1L, utils::head(runs$lengths, -1L)))
    entries <- starts[runs$values]
    exits <- starts[!runs$values][
      starts[!runs$values] > (if (length(entries)) entries[[1L]] else Inf)
    ]
    expect_identical(
      summary_rows$first_entry_s[[i1]],
      if (length(entries)) rows$s[[entries[[1L]]]] else NA_real_
    )
    expect_identical(
      summary_rows$last_exit_s[[i1]],
      if (length(exits)) rows$s[[exits[[length(exits)]]]] else NA_real_
    )
    expect_identical(
      summary_rows$n_reentries[[i1]], max(0L, length(entries) - 1L)
    )
    checked <- checked + 1L
  }
  expect_identical(checked, length(cv$substrates) * length(fbsp_declared))
})


# --- (iv) availability is not selection ---------------------------------- #

test_that("a column no fold can estimate is reported unavailable, not unselected", {
  skip_if_not_installed("glmnet")
  d <- fbsp_drop_data()
  cv <- .fb_sel_cv(
    d, "y", ~x, ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + w,
    random = c("trial", "location"),
    resamples = fb_resamples(unit = ~location, v = 3L), seed = 1L
  )
  path <- .fb_sel_glmnet_path(
    cv$substrates,
    list(list(
      family = "laplace", columns = 1:10, hyperparameters = list(scale = 0.4)
    )),
    folds = cv$folds, decision = list(n_lambda = 6L), seed = 1L
  )
  by_fold <- path$support_by_resolution$history_by_fold
  history <- path$support_by_resolution$history

  for (label in names(cv$substrates)) {
    dropped <- cv$substrates[[label]]$candidate_space$dropped$column
    rows <- by_fold[by_fold$fold == label, ]
    expect_identical(sort(rows$column[!rows$available]), sort(dropped))
    # An unavailable column carries the grid but no membership, and no entry.
    unavailable <- rows[!rows$available, ]
    expect_true(all(unavailable$n_in_support == 0L))
    expect_true(all(is.na(unavailable$first_entry_s)))
    expect_true(all(unavailable$n_points == 6L))
    seen <- history[history$fold == label & history$in_support, "column"]
    expect_length(intersect(seen, dropped), 0L)
  }

  # The distinction is not vacuous: a column that WAS estimable in every fold
  # and still never entered any support is on the same table, and reads
  # available with no entry.
  never <- by_fold[by_fold$available & by_fold$n_in_support == 0L, ]
  expect_gt(nrow(never), 0L)
  expect_true(all(is.na(never$first_entry_s)))
  expect_setequal(
    unique(by_fold$column[!by_fold$available]),
    path$selection$candidate_availability$column[
      path$selection$candidate_availability$n_available == 0L
    ]
  )
})


# --- (v) a family with no support ---------------------------------------- #

test_that("a ridge path reports no support history, and the assertion is not vacuous", {
  skip_if_not_installed("glmnet")
  cv <- fbsp_cv()
  ridge <- .fb_sel_glmnet_path(
    cv$substrates,
    list(list(
      family = "ridge", columns = 1:5, hyperparameters = list(scale = 0.6)
    )),
    folds = cv$folds, decision = list(n_lambda = 6L), seed = 5L
  )
  expect_false(ridge$support_by_resolution$selects)
  expect_identical(nrow(ridge$support_by_resolution$history), 0L)
  expect_identical(nrow(ridge$support_by_resolution$history_by_fold), 0L)
  expect_identical(ridge$support_by_resolution$alpha, 0)

  # The same folds under a lasso do carry a history, so the ridge statement is
  # a property of the family rather than of a path that happened to select
  # nothing.
  lasso <- .fb_sel_glmnet_path(
    cv$substrates, fbsp_blocks(1), folds = cv$folds,
    decision = list(n_lambda = 6L), seed = 5L
  )
  expect_gt(nrow(lasso$support_by_resolution$history), 0L)
  expect_identical(lasso$support_by_resolution$alpha, 1)

  # With no benchmark arm there is no grid to walk, and the history says so
  # rather than reporting the one canonical point as a path.
  fixed_only <- .fb_sel_glmnet_path(
    cv$substrates, fbsp_blocks(1), folds = cv$folds,
    decision = list(benchmark = FALSE)
  )
  expect_identical(nrow(fixed_only$support_by_resolution$history), 0L)
})


# --- (vi) the accessors --------------------------------------------------- #

test_that("tidy() carries the per-column path summary and the history beside it", {
  skip_if_not_installed("glmnet")
  fbsp_quiet()
  d <- fbsp_data()
  sel <- fb_select(
    fbsp_stub_base(d), ~ z1 + z2 + z3 + z4 + z5,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 4L),
    seed = 1L
  )
  td <- tidy(sel)
  expect_true(all(c(
    "folds_available", "folds_entered", "path_entry_s", "path_exit_s",
    "path_reentries"
  ) %in% names(td)))
  # The columns the method already returned are unchanged.
  expect_true(all(c(
    "term", "source_term", "role", "block", "family", "retained", "estimate"
  ) %in% names(td)))
  expect_identical(attr(td, "draw_source"), "none")

  # A protected column is on no penalised path, and says so rather than
  # reporting a zero.
  protected <- td[td$role == "protected", ]
  expect_gt(nrow(protected), 0L)
  expect_true(all(is.na(protected$folds_available)))
  expect_true(all(is.na(protected$path_entry_s)))

  points <- attr(td, "support_path")
  by_fold <- attr(td, "support_path_by_fold")
  expect_true(is.data.frame(points))
  expect_true(all(c(
    "alpha", "fold", "lambda_index", "s", "lambda_canonical", "column",
    "available", "in_support"
  ) %in% names(points)))
  expect_true(all(c(
    "fold", "column", "available", "n_points", "n_in_support",
    "first_entry_s", "last_exit_s", "n_reentries"
  ) %in% names(by_fold)))
  expect_identical(
    sort(unique(points$fold)),
    sort(vapply(sel$folds, function(f) f$label, character(1L)))
  )

  # Every per-column number the row reports is the stated summary over the
  # folds, recomputed here from the table beside it.
  for (column in fbsp_declared) {
    row <- td[td$term == column, ]
    rows <- by_fold[by_fold$column == column, ]
    entered <- rows[rows$n_in_support > 0L, ]
    expect_identical(row$folds_available, sum(rows$available))
    expect_identical(row$folds_entered, nrow(entered))
    expect_identical(
      row$path_entry_s,
      if (nrow(entered)) max(entered$first_entry_s) else NA_real_
    )
    expect_identical(
      row$path_exit_s,
      if (any(!is.na(entered$last_exit_s))) {
        min(entered$last_exit_s, na.rm = TRUE)
      } else {
        NA_real_
      }
    )
    expect_identical(
      row$path_reentries,
      if (nrow(entered)) max(entered$n_reentries) else NA_integer_
    )
  }
  # The chosen support and the path are two different statements, and both are
  # on the table: every retained candidate entered some fold's path.
  retained <- td[td$role == "candidate" & td$retained, ]
  expect_gt(nrow(retained), 0L)
  expect_true(all(retained$folds_entered > 0L))
})


test_that("glance() reports the tuning coordinate, the alpha and the size of the path", {
  skip_if_not_installed("glmnet")
  fbsp_quiet()
  d <- fbsp_data()
  sel <- fb_select(
    fbsp_stub_base(d), ~ z1 + z2 + z3 + z4 + z5,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 4L),
    seed = 1L
  )
  g <- glance(sel)
  expect_identical(nrow(g), 1L)
  expect_true(all(c(
    "tuning_scale", "penalty_alpha", "selected_s", "path_points",
    "path_columns_entered"
  ) %in% names(g)))
  expect_identical(g$tuning_scale, "s = lambda / lambda_max")
  expect_identical(g$penalty_alpha, 1)
  expect_identical(g$selected_s, sel$benchmark$glmnet$selection$selected_s)
  by_lambda <- as.data.frame(sel$benchmark$glmnet$support$by_lambda)
  expect_identical(
    g$path_points,
    length(unique(by_lambda$lambda_index[by_lambda$arm == "benchmark"]))
  )
  by_fold <- sel$benchmark$glmnet$support$history_by_fold
  entered <- tapply(by_fold$n_in_support, by_fold$column, function(v) any(v > 0L))
  expect_identical(g$path_columns_entered, sum(entered))
  expect_gt(g$path_columns_entered, 0L)

  # A run with no penalised optimiser behind it reports the fields as absent
  # rather than as zeros.
  bare <- structure(
    list(
      workflow = "bayes", provider = "brms", prior = "none",
      decision = "full", folds = list(),
      selection = list(chosen_terms = character(0), chosen_size = 0L),
      certificate = list(resampling = list(prediction = "new_cluster"))
    ),
    class = "fb_select"
  )
  gb <- glance(bare)
  expect_true(is.na(gb$tuning_scale))
  expect_true(is.na(gb$penalty_alpha))
  expect_true(is.na(gb$path_points))
  expect_true(is.na(gb$path_columns_entered))
})

test_that("entrance-path plot returns scored path and labelled entries", {
  skip_if_not_installed("glmnet")
  fbsp_quiet()
  d <- fbsp_data()
  sel <- fb_select(
    fbsp_stub_base(d), ~ z1 + z2 + z3 + z4 + z5,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 4L),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 1L
  )
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  plotted <- plot(sel, type = "entrance_path", max_labels = 3L)
  expect_true(is.list(plotted))
  expect_true(all(c("scores", "entries") %in% names(plotted)))
  expect_true(all(c("s", "score") %in% names(plotted$scores)))
  expect_true(all(c("variable", "entry_s", "score") %in%
                    names(plotted$entries)))
  expect_lte(nrow(plotted$entries), 3L)
  expect_true(all(plotted$entries$variable %in% fbsp_declared))
})

test_that("assessment predictions expose the selected scored rows and contract", {
  skip_if_not_installed("glmnet")
  fbsp_quiet()
  d <- fbsp_data()
  sel <- fb_select(
    fbsp_stub_base(d), ~ z1 + z2 + z3 + z4 + z5,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 4L),
    decision = fb_decision(method = "full", metric = "rmse"),
    seed = 1L
  )
  predictions <- augment(sel, type = "assessment_predictions")
  required <- c(
    "stage", "workflow", "provider", "model", "fold", "repeat_id",
    "row", "row_key", "unit", "observed_model_scale",
    "predicted_model_scale", "residual", "selected_lambda_index",
    "selected_s", "lambda_canonical", "prediction", "response",
    "contract_hash"
  )
  expect_identical(names(predictions), required)
  expect_equal(nrow(predictions), nrow(d))
  expect_equal(length(unique(predictions$row_key)), nrow(predictions))
  expect_equal(predictions$observed_model_scale, d$y[predictions$row])
  expect_equal(
    predictions$observed_model_scale - predictions$predicted_model_scale,
    predictions$residual
  )
  expect_true(all(predictions$row %in% seq_len(nrow(d))))
  expect_true(all(predictions$unit == as.character(d$site[predictions$row])))
  expect_identical(unique(predictions$model),
                   sel$benchmark$glmnet$selection$selected_model)
  expect_identical(unique(predictions$selected_lambda_index),
                   sel$benchmark$glmnet$selection$selected_lambda_index)
  expect_match(unique(predictions$contract_hash), "^[0-9a-f]{32}$")
  expect_identical(unique(predictions$prediction), "new_cluster")
  expect_identical(unique(predictions$response), "y")
})

test_that("assessment-prediction schema is stable where no optimiser ran", {
  old <- structure(list(clusters = NULL), class = "fb_select")
  predictions <- augment(old, type = "assessment_predictions")
  expect_s3_class(predictions, "data.frame")
  expect_equal(nrow(predictions), 0L)
  expect_true(all(c("row_key", "predicted_model_scale", "contract_hash") %in%
                    names(predictions)))
})

test_that("a fold's predictions do not read that fold's assessment outcomes", {
  skip_if_not_installed("glmnet")
  d <- fbsp_data(n_cluster = 8L, per = 4L)
  cv <- fbsp_cv(d, v = 4L)
  split <- cv$folds[[1L]]
  original <- cv$substrates[[split$label]]
  changed <- d
  changed$y[split$assess_rows] <- changed$y[split$assess_rows] + 1000
  perturbed <- .fb_sel_fold_substrate(
    changed,
    split,
    response = "y",
    protected = ~x,
    candidates = ~ z1 + z2 + z3 + z4 + z5,
    random = "site",
    prediction = "new_cluster",
    fold_space = original$candidate_space
  )
  fit_path <- function(substrate, data) {
    .fb_sel_glmnet_path(
      list(substrate),
      fbsp_blocks(1),
      folds = list(split),
      decision = list(rule = "min", metric = "rmse", n_lambda = 6L),
      seed = 44L,
      cluster = as.character(data$site)
    )
  }
  path_original <- fit_path(original, d)
  path_perturbed <- fit_path(perturbed, changed)
  key <- c("arm", "lambda_index", "fold", "row")
  first <- as.data.frame(path_original$fitted_values)
  second <- as.data.frame(path_perturbed$fitted_values)
  first <- first[do.call(order, first[key]), ]
  second <- second[do.call(order, second[key]), ]
  expect_identical(first[key], second[key])
  expect_equal(first$fitted, second$fitted, tolerance = 1e-12)
  expect_equal(first$lambda_canonical, second$lambda_canonical,
               tolerance = 1e-12)
  expect_true(all(second$observed - first$observed == 1000))
})


# --- (vii) the history is the grid ---------------------------------------- #

test_that("the history walks the tuning grid and not the single canonical point", {
  skip_if_not_installed("glmnet")
  cv <- fbsp_cv()
  path <- .fb_sel_glmnet_path(
    cv$substrates, fbsp_blocks(1), folds = cv$folds,
    decision = list(n_lambda = 7L), seed = 5L
  )
  history <- path$support_by_resolution$history
  expect_identical(
    nrow(history),
    7L * length(cv$substrates) * length(fbsp_declared)
  )
  expect_identical(sort(unique(history$lambda_index)), 1:7)

  # The fixed arm is one canonical point whose `s` differs between folds; it
  # is on `by_lambda` and not on the history, so an entry read on `s` is read
  # on the grid.
  by_lambda <- as.data.frame(path$support_by_resolution$by_lambda)
  fixed_rows <- by_lambda[by_lambda$arm == "fixed", ]
  expect_identical(nrow(fixed_rows), length(cv$substrates))
  expect_gt(length(unique(fixed_rows$s)), 1L)
  expect_false(any(history$lambda_canonical %in% fixed_rows$lambda_canonical))
})
