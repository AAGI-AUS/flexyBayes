# =============================================================================
# fb_select() fold compilation under a candidate-space contract: every fold
# score refers to one declared candidate space, and a fold whose training rows
# leave part of the design unidentified completes on the columns it keeps or
# refuses by name before any engine runs.
#
# The fixture reproduces the shape that killed the engine in a pilot at scale:
# candidates constant within a trial, more of them than the training trials of
# a fold that holds whole locations out can identify. Six locations of two
# trials of eight plots; folds of two locations, so each fold trains on eight
# trials. Nine trial-constant candidates are identified on all rows (twelve
# trials leave eleven) and not on any fold's training rows (eight leave
# seven), so each fold drops two. The rank depends only on the trials; eight
# plots per trial identify the residual variance well enough that an INLA
# fold refit, whose trial level is saturated once the drop is made, reaches
# the same hyperparameter mode on every run.
#
# Oracles the implementation did not write:
#   the QR rank of each fold's training design, from stats::model.matrix()
#   and qr() in this file;
#   the count of trial-constant candidates eight training trials identify,
#   from the construction;
#   stats::lm.fit() on a fold's whitened kept columns;
#   the terms each INLA fold refit was actually handed, read at the refit.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

p2_data <- function(seed = 11L) {
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

p2_candidates <- ~ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + w
p2_declared <- c(paste0("c", 1:9), "w")

p2_resamples <- function() fb_resamples(unit = ~location, v = 3L)

p2_folds <- function(d) .fb_sel_folds(p2_resamples(), d, seed = 1L)

# The rank of a fold's training design, from model.matrix() and qr() alone.
p2_rank <- function(d, rows, protected = ~x, candidates = p2_candidates) {
  X <- cbind(
    stats::model.matrix(protected, d),
    stats::model.matrix(candidates, d)[, -1L, drop = FALSE]
  )
  qr(X[rows, , drop = FALSE])$rank
}

# A base carrying the parsed model, the data and the call record the fold
# routes read, with no engine behind it.
p2_stub <- function(d, fixed, random = NULL, class = "flexybayes_inla") {
  fb <- if (is.null(random)) {
    fb_from_brms(fixed, data = d, family = "gaussian")
  } else {
    fb_from_asreml(fixed, random = random, data = d)
  }
  fb$family <- "gaussian"
  fb$link <- fb$link %||% "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = fixed, random = random, family = "gaussian")
      )
    ),
    class = c(class, "flexybayes", "list")
  )
}

p2_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

p2_cache <- new.env(parent = emptyenv())

p2_inla_base <- function() {
  if (!exists("base", envir = p2_cache, inherits = FALSE)) {
    assign(
      "base",
      flexybayes(y ~ x, random = ~ trial + location, data = p2_data(),
                 backend = "inla", verbose = FALSE),
      envir = p2_cache
    )
  }
  get("base", envir = p2_cache, inherits = FALSE)
}


# --- (i) the rule, against the analytic rank ----------------------------- #

test_that("each fold keeps the candidate columns its training rank allows, and names the rest", {
  d <- p2_data()
  folds <- p2_folds(d)
  expect_identical(p2_rank(d, seq_len(nrow(d))), 12L)

  space <- .fb_sel_candidate_space(
    ~x, p2_candidates, d, folds,
    groupings = c("trial", "location"), route = "benchmark"
  )
  expect_identical(space$declared$column, p2_declared)
  expect_identical(space$dropped_coefficient, "NA")
  expect_identical(space$tol, 1e-7)

  for (split in folds) {
    f <- space$folds[[split$label]]
    # One declared space: the kept and dropped columns partition it.
    expect_setequal(c(f$kept_columns, f$dropped$column), p2_declared)
    expect_length(intersect(f$kept_columns, f$dropped$column), 0L)
    # The analytic rank, protected rank 2.
    expect_identical(
      length(f$kept_columns),
      p2_rank(d, split$train_rows) - 2L
    )
    kept <- cbind(
      stats::model.matrix(~x, d),
      stats::model.matrix(p2_candidates, d)[, f$kept_columns, drop = FALSE]
    )[split$train_rows, , drop = FALSE]
    expect_identical(qr(kept)$rank, ncol(kept))
    # The construction: k training trials identify k - 1 trial-constant
    # candidates beside the intercept.
    n_trials <- length(unique(d$trial[split$train_rows]))
    expect_identical(sum(grepl("^c", f$kept_columns)), n_trials - 1L)
    expect_identical(f$dropped$column, c("c8", "c9"))
    expect_identical(f$dropped$fold, rep(split$label, 2L))
    expect_match(f$dropped$reason, "constant within `trial`", fixed = TRUE)
    # Two orders either side of the tolerance.
    expect_gt(f$margin$kept_min, 1e-5)
    expect_lt(f$margin$dropped_max, 1e-9)
  }
  expect_identical(
    space$availability$n_available,
    c(rep(3L, 7L), 0L, 0L, 3L)
  )
  expect_identical(unique(space$availability$n_folds), 3L)
})


# --- (ii) the benchmark route fits the kept columns ---------------------- #

test_that("the benchmark route fits each fold on its kept columns and reports NA for a dropped one", {
  d <- p2_data()
  cv <- .fb_sel_cv(
    d, "y", ~x, p2_candidates, random = c("trial", "location"),
    resamples = p2_resamples(), seed = 1L
  )
  expect_identical(cv$candidate_space$dropped$column, rep(c("c8", "c9"), 3L))
  for (s in cv$substrates) {
    beta <- s$learned$pilot_beta
    dropped <- s$candidate_space$dropped$column
    expect_identical(names(beta), c("(Intercept)", "x", p2_declared))
    expect_identical(names(beta)[is.na(beta)], dropped)
    expect_identical(s$candidate_kept, setdiff(p2_declared, dropped))

    x_kept <- cbind(
      s$x_protected,
      s$x_candidate[, s$candidate_kept, drop = FALSE]
    )[s$train_rows, , drop = FALSE]
    lower <- t(s$blocks$factor_tt)
    reference <- stats::lm.fit(
      forwardsolve(lower, x_kept),
      forwardsolve(lower, s$response[s$train_rows])
    )$coefficients
    expect_equal(
      unname(beta[!is.na(beta)]), unname(reference),
      tolerance = 1e-8
    )
  }
  expect_true(all(is.finite(cv$aggregated$score)))
})


# --- (iii) the glmnet route keeps the canonical map ---------------------- #

test_that("the glmnet path indexes every coefficient by the declared columns and reports availability beside selection", {
  skip_if_not_installed("glmnet")
  d <- p2_data()
  cv <- .fb_sel_cv(
    d, "y", ~x, p2_candidates, random = c("trial", "location"),
    resamples = p2_resamples(), seed = 1L
  )
  blocks <- list(list(
    family = "laplace",
    columns = seq_along(p2_declared),
    hyperparameters = list(scale = 0.4)
  ))
  path <- .fb_sel_glmnet_path(
    cv$substrates, blocks, folds = cv$folds,
    decision = list(n_lambda = 5L), seed = 1L
  )
  co <- as.data.frame(path$coefficients_original_scale)
  points <- unique(paste(co$arm, co$lambda_index))
  for (label in names(cv$substrates)) {
    dropped <- cv$substrates[[label]]$candidate_space$dropped$column
    for (point in points) {
      rows <- co[
        co$fold == label & paste(co$arm, co$lambda_index) == point &
          co$role == "candidate",
      ]
      expect_identical(rows$term, p2_declared)
      expect_identical(rows$term[is.na(rows$estimate)], dropped)
    }
    # The provider was handed the kept columns and nothing else.
    parts <- .fb_sel_glmnet_reduction(
      cv$substrates[[label]], seq_along(p2_declared)
    )
    expect_identical(
      colnames(parts$payload$x), setdiff(p2_declared, dropped)
    )
  }

  availability <- path$selection$candidate_availability
  expect_identical(availability$column, p2_declared)
  expect_identical(availability$n_available, c(rep(3L, 7L), 0L, 0L, 3L))
  supports <- path$support_by_resolution$by_lambda
  chosen <- supports[
    supports$arm == "benchmark" &
      supports$lambda_index == path$selection$selected_lambda_index,
  ]
  expect_identical(
    availability$n_selected,
    vapply(
      p2_declared,
      function(column) {
        sum(vapply(
          strsplit(chosen$terms, ", ", fixed = TRUE),
          function(support) column %in% support,
          logical(1L)
        ))
      },
      integer(1L),
      USE.NAMES = FALSE
    )
  )
  expect_true(all(availability$n_selected <= availability$n_available))
})


# --- (iv) the pilot's cell on INLA --------------------------------------- #

test_that("an INLA fold refit is handed the fold's kept terms, and the certificate names what it left out", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  p2_quiet()
  d <- p2_data()
  base <- p2_inla_base()
  real_refit <- .fb_sel_bayes_refit
  handed <- list()
  local_mocked_bindings(
    .fb_sel_bayes_refit = function(base, data, fixed, ...) {
      handed[[length(handed) + 1L]] <<- attr(stats::terms(fixed), "term.labels")
      real_refit(base, data, fixed, ...)
    }
  )
  sel <- fb_select(
    base, p2_candidates,
    benchmark = TRUE,
    resamples = p2_resamples(),
    seed = 1L
  )

  space <- sel$certificate$design$candidate_space
  labels <- vapply(sel$folds, function(s) s$label, character(1L))
  expect_identical(space$declared, p2_declared)
  expect_identical(space$dropped_coefficient, "NA")
  expect_identical(space$availability$n_available, c(rep(3L, 7L), 0L, 0L, 3L))

  # The cost probe and each fold's protected-only model carry no candidate;
  # the three full refits carry the kept terms of their fold, in fold order.
  full <- Filter(function(terms) any(terms %in% p2_declared), handed)
  expect_length(full, length(labels))
  for (k in seq_along(labels)) {
    left_out <- setdiff(p2_declared, full[[k]])
    expect_identical(
      left_out,
      space$dropped$column[space$dropped$fold == labels[[k]]]
    )
    expect_identical(
      length(intersect(p2_declared, full[[k]])),
      p2_rank(d, sel$folds[[k]]$train_rows) - 2L
    )
  }
  expect_true(all(is.finite(sel$scores$score)))
  # The benchmark arm dropped the same columns, so its score contract agrees
  # with the Bayesian arm's and the reference stage places the two together.
  expect_true(sel$certificate$reference_predictive$available)
  expect_setequal(
    sel$certificate$reference_predictive$arms$arm, c("bayes", "gls")
  )
  expect_match(
    paste(utils::capture.output(print(sel$certificate)), collapse = "\n"),
    "candidates : 10 declared columns", fixed = TRUE
  )
})


test_that("the benchmark and optimiser workflows carry the contract onto the object and certificate", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  p2_quiet()
  base <- p2_inla_base()
  bench <- fb_select(
    base, p2_candidates,
    workflow = "unpenalized_benchmark",
    resamples = p2_resamples(),
    seed = 1L
  )
  expect_identical(
    bench$certificate$design$candidate_space$dropped$column,
    rep(c("c8", "c9"), 3L)
  )
  opt <- fb_select(
    base, p2_candidates,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
    workflow = "penalized_optimizer",
    resamples = p2_resamples(),
    seed = 1L
  )
  availability <- opt$selection$candidate_availability
  expect_identical(availability$column, p2_declared)
  expect_identical(availability$n_available, c(rep(3L, 7L), 0L, 0L, 3L))
  expect_identical(
    opt$certificate$benchmark$glmnet$candidate_availability, availability
  )
  expect_identical(
    opt$certificate$design$candidate_space$availability$n_available,
    availability$n_available
  )
})


# --- (v) refusals before any engine runs -------------------------------- #

test_that("a protected column the training rows do not identify refuses by name before anything is fitted", {
  skip_if_not_installed("brms")
  d <- p2_data()
  folds <- p2_folds(d)
  held <- unique(as.character(d$location[folds[[1L]]$assess_rows]))
  d$irrig <- as.numeric(d$location %in% held)

  err <- tryCatch(
    .fb_sel_cv(
      d, "y", ~ x + irrig, ~ c1 + w, random = c("trial", "location"),
      resamples = p2_resamples(), seed = 1L
    ),
    error = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_selection_fold_protected_unidentified"
  )
  expect_match(conditionMessage(err), "irrig", fixed = TRUE)
  expect_match(conditionMessage(err), paste0("fold ", folds[[1L]]$label))

  base <- p2_stub(d, y ~ x + irrig, random = ~ trial + location)
  reached <- 0L
  refit <- function(...) {
    reached <<- reached + 1L
    stop("the refit was reached", call. = FALSE)
  }
  err <- tryCatch(
    .fb_sel_bayes_cv(
      base, ~ c1 + w,
      resamples = p2_resamples(), metric = "elpd", seed = 1L,
      refit_fn = refit, probe_seconds = 1
    ),
    error = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_selection_fold_protected_unidentified"
  )
  expect_identical(reached, 0L)
})


test_that("a candidate level present only in assessment rows refuses before the cost probe fits", {
  skip_if_not_installed("brms")
  d <- p2_data()
  d$soil <- factor(ifelse(
    d$location == "L1", "sand",
    ifelse(as.integer(d$trial) %% 2L == 0L, "loam", "clay")
  ))
  base <- p2_stub(d, y ~ x, random = ~ trial + location)
  reached <- 0L
  local_mocked_bindings(
    .fb_sel_bayes_refit = function(...) {
      reached <<- reached + 1L
      stop("the refit was reached", call. = FALSE)
    }
  )
  err <- tryCatch(
    fb_select(base, ~ c1 + soil, resamples = p2_resamples(), seed = 1L),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_selection_fold_level_absent")
  expect_match(conditionMessage(err), "sand", fixed = TRUE)
  expect_identical(reached, 0L)
  expect_error(
    .fb_sel_cv(
      d, "y", ~x, ~ c1 + soil, random = c("trial", "location"),
      resamples = p2_resamples(), seed = 1L
    ),
    class = "flexybayes_refusal_selection_fold_level_absent"
  )
})


test_that("a term aliased in part drops one column on the benchmark route and refuses on the Bayesian one", {
  skip_if_not_installed("brms")
  d <- p2_data()
  d$zone <- factor(c("A", "B", "C")[(as.integer(d$trial) - 1L) %% 3L + 1L])
  candidates <- ~ c1 + c2 + c3 + c4 + c5 + c6 + zone
  folds <- p2_folds(d)

  space <- .fb_sel_candidate_space(
    ~x, candidates, d, folds,
    groupings = c("trial", "location"), route = "benchmark"
  )
  for (split in folds) {
    f <- space$folds[[split$label]]
    expect_identical(f$dropped$column, "zoneC")
    expect_identical(
      length(f$kept_columns),
      p2_rank(d, split$train_rows, candidates = candidates) - 2L
    )
  }

  base <- p2_stub(d, y ~ x, random = ~ trial + location)
  reached <- 0L
  refit <- function(...) {
    reached <<- reached + 1L
    stop("the refit was reached", call. = FALSE)
  }
  err <- tryCatch(
    .fb_sel_bayes_cv(
      base, candidates,
      resamples = p2_resamples(), metric = "elpd", seed = 1L,
      refit_fn = refit, probe_seconds = 1
    ),
    error = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_selection_fold_candidate_drop_not_representable"
  )
  expect_match(conditionMessage(err), "`zone` (unidentified: zoneC)", fixed = TRUE)
  expect_identical(reached, 0L)
})


test_that("the penalised reference's out-of-fold score refuses a drop before the reference model is sampled", {
  skip_if_not_installed("brms")
  d <- p2_data()
  base <- p2_stub(
    d, y ~ x + (1 | trial) + (1 | location), class = "flexybayes_brms"
  )
  reached <- 0L
  local_mocked_bindings(
    .fb_sel_reference_fit = function(...) {
      reached <<- reached + 1L
      stop("the reference fit was reached", call. = FALSE)
    }
  )
  run <- function(candidates) {
    tryCatch(
      fb_select(
        base, candidates,
        penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
        decision = fb_decision(method = "full"),
        benchmark = TRUE,
        resamples = p2_resamples(),
        seed = 1L
      ),
      error = function(e) e
    )
  }
  err <- run(p2_candidates)
  expect_s3_class(
    err, "flexybayes_refusal_selection_fold_candidate_drop_not_representable"
  )
  expect_match(conditionMessage(err), "`c8`", fixed = TRUE)
  expect_identical(reached, 0L)
  expect_identical(
    conditionMessage(run(~ c1 + w)), "the reference fit was reached"
  )
})


test_that("a design every fold identifies drops nothing and reports full availability", {
  d <- p2_data()
  cv <- .fb_sel_cv(
    d, "y", ~x, ~ c1 + c2 + w, random = c("trial", "location"),
    resamples = p2_resamples(), seed = 1L
  )
  expect_identical(nrow(cv$candidate_space$dropped), 0L)
  expect_identical(cv$candidate_space$availability$n_available, rep(3L, 3L))
  for (s in cv$substrates) {
    expect_false(anyNA(s$learned$pilot_beta))
  }
})
