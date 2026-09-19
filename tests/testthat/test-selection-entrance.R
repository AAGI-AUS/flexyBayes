# =============================================================================
# Entrance profiles, group profiles, subsample stability and the entrance heat
# map of the penalised optimiser, read through fb_select(), tidy() and plot().
#
# The oracle is the stored per-fold coefficient table
# (`benchmark$glmnet$coefficients`) and, for stability, the stored
# per-subsample coefficient table (`benchmark$glmnet$stability$coefficients`).
# Every hand computation in this file reads those tables and the tuning grid,
# decides membership as a non-zero coefficient and availability as a
# coefficient that is not NA, and summarises with its own code. None of it
# calls the package's summarisers.
#
# What each block holds:
#
#   Column profiles. Folds available, folds entered, the earliest entry, the
#   median entry and its quartiles in `tidy()` equal the hand computation, on
#   a design where every column is estimable and on one where folds drop
#   columns.
#
#   One definition. The entrance-path labels and the heat map's markers sit
#   at `tidy()$entry_s_median`, and the columns the plot leaves unlabelled are
#   the ones whose median is NA.
#
#   Group profiles. A grouped penalty's groups and an elastic net's declared
#   blocks are profiled as a group entering when any member does.
#
#   Stability. Frequencies, the maximum over the region, the stable set, the
#   average selected count and the error bound equal the hand computation;
#   each subsample holds half of the declared units; the settings are
#   declared, printed, recorded and costed in `plan = TRUE`; and the stable
#   set separates the planted columns from noise as the signal grows.
#
#   Refusals. Stability asked for where it cannot run refuses by name.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

# Sites of rows, one protected covariate, candidates `z1`..`z<p>` of which the
# first three carry signal at `signal` times their base effect, a site random
# effect and a region column that groups the sites.
fben_data <- function(n_site = 12L, per = 6L, p = 8L, signal = 1,
                      seed = 51L) {
  .fb_sel_with_seed(seed, {
    n <- n_site * per
    site <- factor(rep(sprintf("s%02d", seq_len(n_site)), each = per))
    d <- data.frame(site = site, x = stats::rnorm(n))
    d$region <- factor(paste0("r", (as.integer(site) - 1L) %% 3L + 1L))
    for (j in seq_len(p)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    u <- stats::rnorm(n_site, sd = 1)
    d$y <- 2 + 0.5 * d$x +
      signal * (1.0 * d$z1 - 0.9 * d$z2 + 0.8 * d$z3) +
      u[as.integer(site)] + stats::rnorm(n, sd = 1)
    d
  })
}

fben_candidates <- function(p = 8L) {
  stats::as.formula(paste("~", paste0("z", seq_len(p), collapse = " + ")))
}

# Six locations of two trials of eight plots with trial-constant candidates,
# so a fold that holds a location out cannot identify every candidate.
fben_drop_data <- function(seed = 11L) {
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
    d$y <- 1 + 0.5 * d$x + 0.9 * d$c1 - 0.7 * d$c2 +
      stats::rnorm(n_trial, sd = 0.3)[tr] + stats::rnorm(n, sd = 0.4)
    d
  })
}

# A base carrying the parsed model, the data and the call record the fold
# routes read, with no engine behind it.
fben_stub_base <- function(d, random = ~site) {
  fb <- fb_from_asreml(y ~ x, random = random, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = y ~ x, random = random, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

fben_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

fben_select <- function(d, penalty, decision = fb_decision(
                          method = "full", metric = "rmse"
                        ), v = 6L, candidates = fben_candidates(),
                        seed = 1L, plan = FALSE, unit = ~site,
                        random = ~site) {
  fb_select(
    fben_stub_base(d, random = random), candidates,
    penalty = penalty,
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = unit, v = v),
    decision = decision,
    seed = seed,
    plan = plan
  )
}

# --- the hand oracle ----------------------------------------------------- #

# The stored coefficient table of the benchmark grid, candidate rows only,
# with each row's `s` read off the tuning grid by its index.
fben_grid_coefficients <- function(sel) {
  arm <- sel$benchmark$glmnet
  coef <- as.data.frame(arm$coefficients)
  coef <- coef[coef$arm == "benchmark" & coef$role == "candidate", ]
  grid <- as.data.frame(arm$selection$scores)
  coef$s <- grid$s[match(coef$lambda_index, grid$lambda_index)]
  coef
}

# One unit's (fold's or subsample's) membership sequence over the grid, in
# path order from s = 1 downwards, for a set of member columns: held where any
# member carries a non-zero coefficient, available where any member's
# coefficient is not NA.
fben_unit_summary <- function(rows, members) {
  rows <- rows[rows$term %in% members, ]
  points <- sort(unique(rows$s), decreasing = TRUE)
  held <- vapply(points, function(s) {
    e <- rows$estimate[rows$s == s]
    any(!is.na(e) & e != 0)
  }, logical(1L))
  list(
    available = any(!is.na(rows$estimate)),
    first_entry_s = if (any(held)) points[which(held)[1L]] else NA_real_,
    held = stats::setNames(held, format(points, digits = 15L))
  )
}

# The point on s by which at least a fraction `p` of the available folds had
# admitted the column: searched over the observed entries rather than read as
# an order statistic. `first` has one value per available fold, NA for a fold
# on which the column never entered.
fben_entry_at <- function(first, p) {
  if (length(first) == 0L) {
    return(NA_real_)
  }
  entries <- sort(unique(first[!is.na(first)]), decreasing = TRUE)
  admitted <- vapply(entries, function(s) {
    sum(!is.na(first) & first >= s) / length(first)
  }, numeric(1L))
  hit <- entries[admitted >= p]
  if (length(hit) == 0L) NA_real_ else max(hit)
}

# The whole column or group profile from the coefficient table and a map of
# name -> member columns.
fben_hand_profile <- function(coef, groups) {
  folds <- sort(unique(coef$fold))
  out <- lapply(names(groups), function(name) {
    per_fold <- lapply(folds, function(f) {
      fben_unit_summary(coef[coef$fold == f, ], groups[[name]])
    })
    available <- vapply(per_fold, function(u) u$available, logical(1L))
    first <- vapply(per_fold, function(u) u$first_entry_s, numeric(1L))
    first_available <- first[available]
    data.frame(
      name = name,
      folds_available = sum(available),
      folds_entered = sum(!is.na(first)),
      path_entry_s = if (any(!is.na(first))) {
        max(first, na.rm = TRUE)
      } else {
        NA_real_
      },
      entry_s_median = fben_entry_at(first_available, 0.5),
      entry_s_q25 = fben_entry_at(first_available, 0.25),
      entry_s_q75 = fben_entry_at(first_available, 0.75),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$entry_s_iqr <- out$entry_s_q25 - out$entry_s_q75
  out
}

fben_expect_profile <- function(td, hand, key = "term") {
  for (i in seq_len(nrow(hand))) {
    row <- td[td[[key]] == hand$name[[i]], ]
    expect_identical(nrow(row), 1L)
    label <- paste0(key, " ", hand$name[[i]])
    expect_identical(row$folds_available, as.integer(hand$folds_available[[i]]),
                     label = paste(label, "folds_available"))
    expect_identical(row$folds_entered, as.integer(hand$folds_entered[[i]]),
                     label = paste(label, "folds_entered"))
    for (field in c("path_entry_s", "entry_s_median", "entry_s_q25",
                    "entry_s_q75", "entry_s_iqr")) {
      expect_identical(row[[field]], hand[[field]][[i]],
                       label = paste(label, field))
    }
  }
}

# --- (i) column profiles -------------------------------------------------- #

test_that("every column's entrance profile in tidy() equals the hand computation from the coefficient table", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  sel <- fben_select(
    fben_data(),
    fb_penalty(term(everything()) ~ laplace(scale = 0.15))
  )
  td <- tidy(sel)
  coef <- fben_grid_coefficients(sel)
  columns <- paste0("z", 1:8)
  hand <- fben_hand_profile(coef, stats::setNames(as.list(columns), columns))
  fben_expect_profile(td, hand)

  # The oracle is not vacuous: on this path the earliest entry and the median
  # entry differ for some column, so a summary that confused them fails.
  expect_true(any(hand$path_entry_s != hand$entry_s_median, na.rm = TRUE))
  expect_true(all(c("entry_s_median", "entry_s_q25", "entry_s_q75",
                    "entry_s_iqr") %in% names(td)))
  # A protected column is on no penalised path.
  protected <- td[td$role == "protected", ]
  expect_true(all(is.na(protected$entry_s_median)))
  expect_true(all(is.na(protected$entry_s_iqr)))
})

test_that("a column a fold cannot estimate counts as unavailable there, not as a late entry", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  d <- fben_drop_data()
  sel <- fben_select(
    d,
    fb_penalty(term(everything()) ~ laplace(scale = 0.1)),
    v = 3L,
    candidates = stats::as.formula(paste("~", paste0("c", 1:9, collapse = " + "))),
    unit = ~location,
    random = ~trial
  )
  td <- tidy(sel)
  coef <- fben_grid_coefficients(sel)
  # A dropped column is NA at every point of its fold, never at some.
  pattern <- tapply(is.na(coef$estimate), paste(coef$fold, coef$term), mean)
  expect_true(all(pattern %in% c(0, 1)))
  expect_true(any(pattern == 1))

  columns <- paste0("c", 1:9)
  hand <- fben_hand_profile(coef, stats::setNames(as.list(columns), columns))
  expect_true(any(hand$folds_available < 3L))
  fben_expect_profile(td, hand)
})

test_that("a fold that never admitted a column counts against its median entry", {
  # Through fb_select() every estimable column enters every fold by
  # s = 0.001, so a column admitted on fewer than half of its folds is built
  # here as the per-fold record the accessors read.
  first <- list(
    zA = c(0.8, NA, NA, NA, NA),
    zB = c(0.8, 0.5, NA, NA),
    zC = c(0.5, 0.3, NA, NA),
    zD = c(0.9, 0.7, 0.5, 0.3)
  )
  available <- list(
    zA = rep(TRUE, 5L), zB = rep(TRUE, 4L), zC = c(TRUE, TRUE, FALSE, FALSE),
    zD = rep(TRUE, 4L)
  )
  history <- do.call(rbind, lapply(names(first), function(column) {
    data.frame(
      alpha = 1, fold = as.character(seq_along(first[[column]])),
      column = column, available = available[[column]], n_points = 6L,
      n_in_support = ifelse(is.na(first[[column]]), 0L, 2L),
      first_entry_s = first[[column]], last_exit_s = NA_real_,
      n_reentries = 0L, stringsAsFactors = FALSE
    )
  }))
  x <- structure(list(
    column_map = data.frame(
      column = names(first), term = names(first), role = "candidate",
      block = 1L, family = "laplace", centre = 0, scale = 1,
      stringsAsFactors = FALSE
    ),
    selection = list(chosen_terms = character(0)),
    benchmark = list(glmnet = list(
      selection = list(
        scores = data.frame(
          lambda_index = 1:6, s = c(1, 0.9, 0.8, 0.7, 0.5, 0.3),
          score = c(2, 1.9, 1.8, 1.7, 1.6, 1.65)
        ),
        metric = "rmse", selected_s = 0.5,
        tuning_scale = "s = lambda / lambda_max"
      ),
      support = list(history = data.frame(), history_by_fold = history),
      scoring_unit = "site"
    ))
  ), class = "fb_select")
  td <- tidy(x)
  for (column in names(first)) {
    row <- td[td$term == column, ]
    on_path <- first[[column]][available[[column]]]
    expect_identical(row$entry_s_median, fben_entry_at(on_path, 0.5),
                     label = paste(column, "median"))
    expect_identical(row$entry_s_q25, fben_entry_at(on_path, 0.25),
                     label = paste(column, "q25"))
    expect_identical(row$entry_s_q75, fben_entry_at(on_path, 0.75),
                     label = paste(column, "q75"))
  }
  # Written out: admitted on one fold in five has no median; on two of four
  # it has one, at the later of the two.
  expect_identical(td$entry_s_median[td$term == "zA"], NA_real_)
  expect_identical(td$path_entry_s[td$term == "zA"], 0.8)
  expect_identical(td$entry_s_median[td$term == "zB"], 0.5)
  expect_identical(td$entry_s_q75[td$term == "zB"], NA_real_)
  expect_identical(td$entry_s_q75[td$term == "zC"], 0.3)

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  plotted <- plot(x, type = "entrance_path")
  expect_setequal(plotted$entries$variable, c("zB", "zC", "zD"))
})

# --- (ii) one entry definition for tidy() and plot() ---------------------- #

test_that("the entrance-path labels and the heat-map markers sit at tidy()'s median entry", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  sel <- fben_select(
    fben_data(),
    fb_penalty(term(everything()) ~ laplace(scale = 0.15))
  )
  td <- tidy(sel)
  candidates <- td[td$role == "candidate", ]
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)

  entrance <- plot(sel, type = "entrance_path", max_labels = 100L)
  labelled <- candidates[!is.na(candidates$entry_s_median), ]
  expect_setequal(entrance$entries$variable, labelled$term)
  expect_identical(
    entrance$entries$entry_s,
    labelled$entry_s_median[match(entrance$entries$variable, labelled$term)]
  )

  heat <- plot(sel, type = "entrance_heatmap")
  expect_identical(
    heat$rows$entry_s_median,
    candidates$entry_s_median[match(heat$rows$name, candidates$term)]
  )
  expect_setequal(heat$rows$name, candidates$term)

  # Every cell of the fold heat map is the fraction of the folds that could
  # estimate the column whose coefficient is non-zero at that point.
  coef <- fben_grid_coefficients(sel)
  s_grid <- as.data.frame(sel$benchmark$glmnet$selection$scores)$s
  expect_identical(heat$s, s_grid)
  for (name in heat$rows$name) {
    for (k in seq_along(s_grid)) {
      e <- coef$estimate[coef$term == name & coef$s == s_grid[[k]]]
      expected <- if (any(!is.na(e))) {
        sum(!is.na(e) & e != 0) / sum(!is.na(e))
      } else {
        NA_real_
      }
      expect_identical(unname(heat$values[name, k]), expected,
                       label = paste0(name, " at point ", k))
    }
  }
})

test_that("labels on the entrance path no longer crowd at small s", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  sel <- fben_select(
    fben_data(),
    fb_penalty(term(everything()) ~ laplace(scale = 0.15))
  )
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  out <- plot(sel, type = "entrance_path", max_labels = 100L)
  # The grid is geometric in s, so a logarithmic axis spaces its points
  # evenly; on the linear axis half of them fell within 0.05 of zero.
  expect_true(graphics::par("xlog"))
  # Labels on adjacent grid points are drawn on opposite sides of the curve.
  s_grid <- as.data.frame(sel$benchmark$glmnet$selection$scores)$s
  index <- match(out$labels$entry_s, s_grid)
  expect_false(anyNA(index))
  adjacent <- which(abs(diff(index)) <= 1L)
  expect_gt(length(adjacent), 0L)
  for (k in adjacent) {
    expect_false(identical(out$labels$pos[[k]], out$labels$pos[[k + 1L]]))
  }
})

# --- (iii) group profiles -------------------------------------------------- #

test_that("a grouped penalty's group profile equals the hand computation over its members", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  d <- fben_data()
  d$f <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  d$g <- factor(rep(c("u", "v", "w", "t"), length.out = nrow(d)))
  sel <- fben_select(
    d,
    fb_penalty(term(everything()) ~ group_lasso(scale = 0.3)),
    v = 4L,
    candidates = ~ z1 + z2 + f + g
  )
  tg <- tidy(sel, level = "group")
  td <- tidy(sel)
  candidates <- td[td$role == "candidate", ]
  groups <- split(candidates$term, candidates$source_term)
  expect_setequal(tg$group, names(groups))
  expect_identical(
    tg$n_columns,
    as.integer(lengths(groups)[tg$group])
  )
  hand <- fben_hand_profile(fben_grid_coefficients(sel), groups)
  fben_expect_profile(tg, hand, key = "group")
  expect_true(all(tg$rule == "term"))
  expect_true(is.data.frame(attr(tg, "support_path_by_fold")))
})

test_that("an elastic net's declared blocks are profiled as groups", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  sel <- fben_select(
    fben_data(),
    fb_penalty(
      term(z1, z2, z3, z4) ~ elastic_net(lambda = 0.2, alpha = 0.5),
      term(z5, z6, z7, z8) ~ elastic_net(lambda = 0.2, alpha = 0.5)
    )
  )
  tg <- tidy(sel, level = "group")
  td <- tidy(sel)
  candidates <- td[td$role == "candidate", ]
  groups <- split(candidates$term, paste0("block_", candidates$block))
  expect_setequal(tg$group, c("block_1", "block_2"))
  expect_true(all(tg$rule == "block"))
  hand <- fben_hand_profile(fben_grid_coefficients(sel), groups)
  fben_expect_profile(tg, hand, key = "group")

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  heat <- plot(sel, type = "entrance_heatmap", level = "group")
  expect_setequal(rownames(heat$values), c("block_1", "block_2"))
  expect_identical(
    heat$rows$entry_s_median,
    tg$entry_s_median[match(heat$rows$name, tg$group)]
  )
})

# --- (iv) subsample stability --------------------------------------------- #

fben_stability_decision <- function(b = 20L, cutoff = 0.75, unit = NULL) {
  fb_decision(
    method = "full", metric = "rmse",
    stability_subsamples = b, stability_cutoff = cutoff,
    stability_unit = unit
  )
}

test_that("stability frequencies, the stable set and the bound equal the hand computation from the subsample coefficient table", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  d <- fben_data()
  sel <- fben_select(
    d,
    fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
    decision = fben_stability_decision(b = 20L, cutoff = 0.7)
  )
  st <- sel$benchmark$glmnet$stability
  expect_true(st$ran)
  expect_identical(st$subsamples, 20L)
  expect_identical(st$cutoff, 0.7)
  expect_identical(st$unit, "site")
  expect_identical(st$n_units, 12L)
  expect_identical(st$subsample_units, 6L)

  # Each subsample is half of the declared units, drawn without replacement.
  draws <- st$draws
  expect_identical(sort(unique(draws$subsample)), 1:20)
  per <- tapply(draws$unit, draws$subsample, function(u) length(unique(u)))
  expect_true(all(per == 6L))
  expect_true(all(draws$unit %in% levels(d$site)))

  coef <- as.data.frame(st$coefficients)
  s_grid <- as.data.frame(sel$benchmark$glmnet$selection$scores)$s
  # One grid for the folds and the subsamples.
  expect_identical(sort(unique(coef$s), decreasing = TRUE), s_grid)
  region <- s_grid[s_grid >= sel$benchmark$glmnet$selection$selected_s]
  expect_identical(st$region_s_min, sel$benchmark$glmnet$selection$selected_s)

  columns <- paste0("z", 1:8)
  max_frequency <- stats::setNames(numeric(length(columns)), columns)
  for (name in columns) {
    rows <- coef[coef$term == name, ]
    freq <- vapply(s_grid, function(s) {
      e <- rows$estimate[rows$s == s]
      sum(!is.na(e) & e != 0) / 20
    }, numeric(1L))
    stored <- st$path[st$path$level == "column" & st$path$name == name, ]
    stored <- stored[order(stored$lambda_index), ]
    expect_identical(stored$frequency, freq, label = paste(name, "frequency"))
    max_frequency[[name]] <- max(freq[s_grid %in% region])
  }
  summary <- st$summary[st$summary$level == "column", ]
  expect_identical(
    summary$max_frequency[match(columns, summary$name)],
    unname(max_frequency)
  )
  stable <- columns[max_frequency >= 0.7]
  expect_setequal(st$stable, stable)
  # The same draws at cutoff 1: the stable set follows the declared cutoff,
  # and on this fixture a column between one half and one tells the two apart.
  expect_true(any(max_frequency >= 0.5 & max_frequency < 1))
  strict <- fben_select(
    d,
    fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
    decision = fben_stability_decision(b = 20L, cutoff = 1)
  )$benchmark$glmnet$stability
  expect_identical(strict$coefficients, st$coefficients)
  expect_setequal(strict$stable, columns[max_frequency >= 1])

  # The average number of columns a subsample selects anywhere in the
  # region, and the Meinshausen-Buhlmann bound on the expected number of
  # false selections.
  selected_count <- vapply(1:20, function(b) {
    rows <- coef[coef$subsample == b & coef$s %in% region, ]
    length(unique(rows$term[!is.na(rows$estimate) & rows$estimate != 0]))
  }, integer(1L))
  expect_identical(st$q_hat, mean(selected_count))
  expect_identical(st$pfer_bound, mean(selected_count)^2 / ((2 * 0.7 - 1) * 8))

  # tidy() carries the same numbers per column.
  td <- tidy(sel)
  rows <- td[match(columns, td$term), ]
  expect_identical(rows$stability_max_frequency, unname(max_frequency))
  expect_identical(rows$stability_selected, unname(max_frequency >= 0.7))
  expect_identical(rows$stability_available, rep(20L, 8L))

  # The stability heat map shows the stored frequencies.
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  heat <- plot(sel, type = "entrance_heatmap", source = "stability")
  for (name in columns) {
    stored <- st$path[st$path$level == "column" & st$path$name == name, ]
    stored <- stored[order(stored$lambda_index), ]
    expect_identical(unname(heat$values[name, ]), stored$frequency)
  }
})

test_that("availability stays apart from selection in the heat map and in stability", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  sel <- fben_select(
    fben_drop_data(),
    fb_penalty(term(everything()) ~ laplace(scale = 0.1)),
    decision = fben_stability_decision(b = 10L),
    v = 3L,
    candidates = stats::as.formula(paste("~", paste0("c", 1:9, collapse = " + "))),
    unit = ~location,
    random = ~trial
  )
  columns <- paste0("c", 1:9)
  s_grid <- as.data.frame(sel$benchmark$glmnet$selection$scores)$s

  # The fold heat map: the share of the folds able to estimate the column,
  # and NA where none is.
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  heat <- plot(sel, type = "entrance_heatmap")
  coef <- fben_grid_coefficients(sel)
  for (name in columns) {
    for (k in seq_along(s_grid)) {
      e <- coef$estimate[coef$term == name & coef$s == s_grid[[k]]]
      expected <- if (any(!is.na(e))) {
        sum(!is.na(e) & e != 0) / sum(!is.na(e))
      } else {
        NA_real_
      }
      expect_identical(unname(heat$values[name, k]), expected,
                       label = paste0(name, " at point ", k))
    }
  }

  # Stability: a subsample that cannot estimate a column counts against its
  # frequency, which is over all B subsamples, and the count of subsamples
  # able to estimate it is reported beside it.
  st <- sel$benchmark$glmnet$stability
  sub <- as.data.frame(st$coefficients)
  available <- vapply(columns, function(name) {
    length(unique(sub$subsample[sub$term == name & !is.na(sub$estimate)]))
  }, integer(1L))
  expect_true(any(available < 10L))
  summary <- st$summary[st$summary$level == "column", ]
  expect_identical(summary$n_available[match(columns, summary$name)],
                   unname(available))
  for (name in columns) {
    freq <- vapply(s_grid, function(s) {
      e <- sub$estimate[sub$term == name & sub$s == s]
      sum(!is.na(e) & e != 0) / 10
    }, numeric(1L))
    stored <- st$path[st$path$level == "column" & st$path$name == name, ]
    stored <- stored[order(stored$lambda_index), ]
    expect_identical(stored$frequency, freq, label = name)
  }
  td <- tidy(sel)
  expect_identical(td$stability_available[match(columns, td$term)],
                   unname(available))
})

test_that("stability over declared blocks counts a block selected when any member is", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  sel <- fben_select(
    fben_data(),
    fb_penalty(
      term(z1, z2, z3, z4) ~ elastic_net(lambda = 0.2, alpha = 0.5),
      term(z5, z6, z7, z8) ~ elastic_net(lambda = 0.2, alpha = 0.5)
    ),
    decision = fben_stability_decision(b = 12L)
  )
  st <- sel$benchmark$glmnet$stability
  coef <- as.data.frame(st$coefficients)
  s_grid <- as.data.frame(sel$benchmark$glmnet$selection$scores)$s
  blocks <- list(block_1 = paste0("z", 1:4), block_2 = paste0("z", 5:8))
  for (name in names(blocks)) {
    freq <- vapply(s_grid, function(s) {
      held <- vapply(1:12, function(b) {
        e <- coef$estimate[coef$subsample == b & coef$s == s &
                             coef$term %in% blocks[[name]]]
        any(!is.na(e) & e != 0)
      }, logical(1L))
      sum(held) / 12
    }, numeric(1L))
    stored <- st$path[st$path$level == "group" & st$path$name == name, ]
    stored <- stored[order(stored$lambda_index), ]
    expect_identical(stored$frequency, freq, label = name)
  }
  tg <- tidy(sel, level = "group")
  expect_true(all(c("stability_max_frequency", "stability_selected") %in%
                    names(tg)))
})

test_that("the stability settings are declared, printed, recorded and costed in the plan", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  d <- fben_data()
  penalty <- fb_penalty(term(everything()) ~ laplace(scale = 0.15))
  decision <- fben_stability_decision(b = 10L, cutoff = 0.8, unit = ~region)

  printed_decision <- utils::capture.output(print(decision))
  expect_true(any(grepl("B = 10", printed_decision, fixed = TRUE)))
  expect_true(any(grepl("0.8", printed_decision, fixed = TRUE)))
  expect_true(any(grepl("region", printed_decision, fixed = TRUE)))

  # The default is no stability run, and it is recorded as such.
  expect_identical(fb_decision()$stability_subsamples, 0L)
  plain <- fben_select(d, penalty)
  expect_false(plain$benchmark$glmnet$stability$ran)
  expect_identical(plain$benchmark$glmnet$stability$subsamples, 0L)
  expect_identical(
    plain$certificate$benchmark$glmnet$stability$subsamples, 0L
  )
  expect_identical(glance(plain)$stability_subsamples, 0L)
  expect_true(all(is.na(tidy(plain)$stability_max_frequency)))

  # The plan states the cost before anything is fitted.
  plan <- fben_select(
    d, penalty,
    decision = fben_stability_decision(b = 10L, cutoff = 0.8),
    plan = TRUE
  )
  cost <- plan$selection$stability
  expect_identical(cost$subsamples, 10L)
  expect_identical(cost$cutoff, 0.8)
  expect_identical(cost$unit, "site")
  expect_identical(cost$unit_source, "resampling unit")
  expect_identical(cost$n_units, 12L)
  expect_identical(cost$subsample_units, 6L)
  expect_identical(cost$n_folds, 6L)
  expect_identical(cost$path_fits, 10L * cost$fits_per_path)
  expect_identical(cost$cv_path_fits, 6L * cost$fits_per_path)
  printed_plan <- utils::capture.output(print(plan))
  expect_true(any(grepl("Stability", printed_plan, fixed = TRUE)))
  expect_true(any(grepl("B = 10", printed_plan, fixed = TRUE)))

  # Declared on sites: it runs, and the result prints and records it.
  sel <- fben_select(
    d, penalty,
    decision = fben_stability_decision(b = 10L, cutoff = 0.8)
  )
  printed <- utils::capture.output(print(sel))
  expect_true(any(grepl("stability", printed, fixed = TRUE)))
  expect_true(any(grepl("B = 10", printed, fixed = TRUE)))
  expect_true(any(grepl("6 of 12 site", printed, fixed = TRUE)))
  record <- sel$certificate$benchmark$glmnet$stability
  expect_identical(record$subsamples, 10L)
  expect_identical(record$cutoff, 0.8)
  expect_identical(record$unit, "site")
  expect_identical(record$unit_source, "resampling unit")
  expect_identical(record$seed, 1L)
  printed_certificate <- utils::capture.output(print(sel$certificate))
  expect_true(any(grepl("B = 10", printed_certificate, fixed = TRUE)))
  g <- glance(sel)
  expect_identical(g$stability_subsamples, 10L)
  expect_identical(g$stability_cutoff, 0.8)
  expect_identical(g$stability_unit, "site")
  expect_identical(g$stability_stable, length(sel$benchmark$glmnet$stability$stable))
})

test_that("stability separates the planted columns from noise as the signal grows", {
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fben_quiet()
  # Separation is the smallest stability of a planted column minus the
  # largest of a noise column, counted in subsamples: a frequency is a
  # multiple of 1 / B, so counts compare exactly.
  b <- 30L
  signals <- c(0, 0.6, 1)
  for (seed in 1:2) {
    separation <- integer(length(signals))
    stable <- vector("list", length(signals))
    for (i in seq_along(signals)) {
      sel <- fben_select(
        fben_data(n_site = 24L, per = 5L, p = 10L, signal = signals[[i]],
                  seed = seed),
        fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
        decision = fben_stability_decision(b = b, cutoff = 0.75),
        candidates = fben_candidates(10L),
        v = 6L
      )
      st <- sel$benchmark$glmnet$stability
      summary <- st$summary[st$summary$level == "column", ]
      count <- stats::setNames(
        as.integer(round(summary$max_frequency * b)), summary$name
      )
      separation[[i]] <- min(count[paste0("z", 1:3)]) -
        max(count[paste0("z", 4:10)])
      stable[[i]] <- st$stable
    }
    label <- paste0("data seed ", seed)
    # With no signal nothing separates and nothing is stable; the separation
    # does not shrink as the signal grows, and at full signal the stable set
    # is the planted columns and nothing else.
    expect_lte(separation[[1L]], 0L, label = label)
    expect_length(stable[[1L]], 0L)
    expect_true(all(diff(separation) >= 0L), label = label)
    expect_gte(separation[[3L]] - separation[[1L]], b / 2)
    expect_setequal(stable[[3L]], paste0("z", 1:3))
  }
})

# --- (v) refusals ----------------------------------------------------------- #

fben_reason <- function(expr) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, "error")
  err$reason_code %||% NA_character_
}

test_that("stability settings refuse at construction when they cannot be read", {
  expect_identical(
    fben_reason(fb_decision(stability_subsamples = -1)),
    "selection_decision_invalid"
  )
  expect_identical(
    fben_reason(fb_decision(stability_subsamples = 2.5)),
    "selection_decision_invalid"
  )
  expect_identical(
    fben_reason(fb_decision(stability_subsamples = 10, stability_cutoff = 0.5)),
    "selection_decision_invalid"
  )
  expect_identical(
    fben_reason(fb_decision(stability_subsamples = 10, stability_cutoff = 1.2)),
    "selection_decision_invalid"
  )
  expect_s3_class(
    fb_decision(stability_subsamples = 10, stability_cutoff = 1),
    "fb_decision"
  )
})

test_that("stability refuses by name where no penalised path would read it", {
  skip_if_not_installed("glmnet")
  fben_quiet()
  d <- fben_data()
  decision <- fben_stability_decision(b = 10L)
  # No seed: the subsample draws are random.
  expect_identical(
    fben_reason(fben_select(
      d, fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
      decision = decision, seed = NULL
    )),
    "selection_seed_required"
  )
  # A ridge path has no support and no s coordinate.
  expect_identical(
    fben_reason(fben_select(
      d, fb_penalty(term(everything()) ~ ridge(scale = 1)),
      decision = decision
    )),
    "selection_decision_invalid"
  )
  # The plan refuses as the run does.
  expect_identical(
    fben_reason(fben_select(
      d, fb_penalty(term(everything()) ~ ridge(scale = 1)),
      decision = decision, plan = TRUE
    )),
    "selection_decision_invalid"
  )
  # The Bayesian route without the optimiser arm runs no penalised path, and
  # refuses before anything is fitted.
  expect_identical(
    fben_reason(fb_select(
      fben_stub_base(d), fben_candidates(),
      penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
      workflow = "bayes", benchmark = FALSE,
      resamples = fb_resamples(unit = ~site, v = 6L),
      decision = fb_decision(stability_subsamples = 10L),
      seed = 1L
    )),
    "selection_decision_invalid"
  )
  # Beside a Bayesian answer, a unit that cannot be halved refuses before the
  # reference model is sampled; the stub base has no engine to sample with.
  expect_identical(
    fben_reason(fb_select(
      fben_stub_base(d), fben_candidates(),
      penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
      workflow = "bayes", benchmark = TRUE,
      resamples = fb_resamples(unit = ~site, v = 6L),
      decision = fb_decision(stability_subsamples = 10L,
                             stability_unit = ~region),
      seed = 1L
    )),
    "selection_fold_not_usable"
  )
  # Three regions cannot be halved into subsamples of at least two units,
  # and the plan says so as the run does.
  for (plan in c(FALSE, TRUE)) {
    expect_identical(
      fben_reason(fben_select(
        d, fb_penalty(term(everything()) ~ laplace(scale = 0.15)),
        decision = fben_stability_decision(b = 10L, unit = ~region),
        plan = plan
      )),
      "selection_fold_not_usable"
    )
  }
})

test_that("the heat map says why it cannot be drawn", {
  bare <- structure(list(benchmark = list()), class = "fb_select")
  expect_message(
    out <- plot(bare, type = "entrance_heatmap"),
    "no penalised optimiser path"
  )
  expect_null(out)
  skip_if_not_installed("glmnet")
  fben_quiet()
  sel <- fben_select(
    fben_data(),
    fb_penalty(term(everything()) ~ laplace(scale = 0.15))
  )
  expect_message(
    out <- plot(sel, type = "entrance_heatmap", source = "stability"),
    "no stability subsamples"
  )
  expect_null(out)
})
