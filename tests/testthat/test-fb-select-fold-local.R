# =============================================================================
# fb_select(): the penalised reference's out-of-fold score learns its
# preprocessing on each fold's analysis rows.
#
# The defect this closes. The candidate columns were standardised once, on
# every row, and that `reference_data` was handed to every fold refit -- so an
# assessment row's covariate value set the training fold's centres and
# scales, and the compiled prior was stated against the whole sample's
# `calibration_n` and `n_train`. The score was not strictly out of fold.
#
# The falsifier puts an extreme value in one row of one assessment fold and
# requires that fold's training data, compiled prior, column map and recorded
# constants to be bitwise unchanged, while every fold that trains on the row
# does change. It runs the real `.fb_select_reference_cv()` loop; only the
# sampler and the predictive draws are replaced, by functions that record
# what they were handed.
# =============================================================================

fbfl_data <- function(seed = 5L) {
  set.seed(seed)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(12L)), each = 5L)))
  for (v in c("p1", "x1", "x2", "z1", "z2")) {
    d[[v]] <- stats::rnorm(nrow(d))
  }
  d$y <- 1 + 0.5 * d$p1 + d$x1 - 0.8 * d$x2 +
    stats::rnorm(12L, sd = 0.7)[as.integer(d$site)] +
    stats::rnorm(nrow(d), sd = 0.5)
  d
}

fbfl_base <- function(d, fixed = y ~ p1 + (1 | site)) {
  fb <- fb_from_brms(fixed, data = d, family = "gaussian")
  fb$link <- fb$link %||% "identity"
  structure(
    list(fb = fb, data = d),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
}

fbfl_resamples <- function() fb_resamples(unit = ~site, v = 3L)

fbfl_penalty <- function() {
  fb_penalty(term(everything()) ~ elastic_net(lambda = 0.1, alpha = 0.5))
}

# One run of the reference out-of-fold loop, recording what each fold's
# reference refit and prediction were handed.
fbfl_run <- function(d, penalty = fbfl_penalty()) {
  seen <- new.env(parent = emptyenv())
  seen$fits <- list()
  seen$assess <- list()
  local_mocked_bindings(
    .fb_sel_reference_fit = function(base, candidates, prior = NULL,
                                     stanvars = NULL, column_map = NULL,
                                     engine = NULL, seed = NULL, data = NULL,
                                     ...) {
      seen$fits[[length(seen$fits) + 1L]] <- list(
        data = data, prior = prior, stanvars = stanvars,
        column_map = column_map
      )
      structure(
        list(stage = "reference"),
        class = c("flexybayes_brms", "flexybayes", "list")
      )
    },
    .fb_sel_bayes_refit = function(...) {
      structure(
        list(stage = "protected"),
        class = c("flexybayes_brms", "flexybayes", "list")
      )
    },
    .fb_sel_posterior_predict_fold = function(fit, assess_rows, data, ...) {
      if (identical(fit$stage, "reference")) {
        seen$assess[[length(seen$assess) + 1L]] <- data
      }
      list(
        mu = matrix(0, nrow = 20L, ncol = length(assess_rows)),
        sigma = rep(1, 20L)
      )
    }
  )
  cv <- .fb_select_reference_cv(
    base = fbfl_base(d),
    candidates = ~ x1 + x2,
    penalty = penalty,
    protected_labels = "p1",
    resamples = fbfl_resamples(),
    prediction = "new_cluster",
    metric = "elpd",
    engine = NULL,
    seed = 7L
  )
  list(cv = cv, fits = seen$fits, assess = seen$assess)
}


test_that("an extreme value in one assessment fold leaves that fold's training prior and data unchanged", {
  skip_if_not_installed("brms")
  d <- fbfl_data()
  clean <- fbfl_run(d)
  folds <- clean$cv$folds
  expect_length(clean$fits, length(folds))

  k <- 2L
  row <- folds[[k]]$assess_rows[[1L]]
  injected <- d
  injected$x1[row] <- 1e6
  extreme <- fbfl_run(injected)

  # The partition is a deal on units and does not read covariate values.
  expect_identical(
    lapply(extreme$cv$folds, function(s) s$assess_rows),
    lapply(folds, function(s) s$assess_rows)
  )

  # The fold whose assessment rows hold the extreme value: every input to its
  # training fit, and every constant recorded for it, is bitwise unchanged.
  expect_identical(extreme$fits[[k]]$data, clean$fits[[k]]$data)
  expect_identical(extreme$fits[[k]]$prior, clean$fits[[k]]$prior)
  expect_identical(extreme$fits[[k]]$stanvars, clean$fits[[k]]$stanvars)
  expect_identical(extreme$fits[[k]]$column_map, clean$fits[[k]]$column_map)
  expect_identical(
    extreme$cv$preprocessing$by_fold[[k]],
    clean$cv$preprocessing$by_fold[[k]]
  )
  # Its other assessment rows are transformed by the same constants; only the
  # injected row moves.
  expect_identical(extreme$assess[[k]][-1L, ], clean$assess[[k]][-1L, ])
  expect_false(identical(extreme$assess[[k]]$x1[[1L]], clean$assess[[k]]$x1[[1L]]))

  # Non-vacuous: every fold that TRAINS on the row sees its constants and its
  # training data move, so the comparison above can see a leak when one is
  # there.
  for (j in setdiff(seq_along(folds), k)) {
    expect_true(row %in% folds[[j]]$train_rows)
    expect_false(identical(extreme$fits[[j]]$data, clean$fits[[j]]$data))
    expect_false(identical(
      extreme$cv$preprocessing$by_fold[[j]]$constants,
      clean$cv$preprocessing$by_fold[[j]]$constants
    ))
  }

  # The defect state, shown red on the same data: whole-sample scaling moves
  # fold k's training rows when only an assessment row changed.
  whole <- function(frame) {
    .fb_select_resolve_penalty(
      fbfl_base(frame), ~ x1 + x2, fbfl_penalty(), "p1",
      calibration = .fb_select_calibration_n(fbfl_resamples(), frame)
    )$reference_data[folds[[k]]$train_rows, , drop = FALSE]
  }
  expect_false(identical(whole(injected), whole(d)))
})


test_that("each fold's constants, calibration and n_train are its analysis rows'", {
  skip_if_not_installed("brms")
  d <- fbfl_data()
  run <- fbfl_run(d)
  folds <- run$cv$folds

  for (k in seq_along(folds)) {
    train <- folds[[k]]$train_rows
    record <- run$cv$preprocessing$by_fold[[k]]
    constants <- record$constants

    # Derived here from the training rows, not read back from the code.
    for (v in c("x1", "x2")) {
      expect_equal(
        constants$centre[constants$column == v], mean(d[[v]][train]),
        tolerance = 1e-12
      )
      expect_equal(
        constants$scale[constants$column == v], stats::sd(d[[v]][train]),
        tolerance = 1e-12
      )
      expect_equal(
        run$fits[[k]]$data[[v]],
        (d[[v]][train] - mean(d[[v]][train])) / stats::sd(d[[v]][train]),
        tolerance = 1e-12
      )
    }
    expect_true(all(constants$standardised))

    units <- length(unique(as.character(d$site[train])))
    expect_identical(record$calibration_n, units)
    expect_identical(record$n_train, units)
    expect_identical(record$n_rows, length(train))
    expect_identical(record$calibration_basis, "resampling_units")

    # The elastic-net density carries n_train as a literal, so the compiled
    # prior itself shows the fold's unit count and not the sample's twelve.
    density <- paste(run$fits[[k]]$prior$prior, collapse = " ")
    expect_match(density, paste0("target += -", units, " * "), fixed = TRUE)
    expect_false(grepl("target += -12 * ", density, fixed = TRUE))
  }
  # Sanity on the fixture: the folds' training constants really do differ
  # from the whole sample's, so agreement above is not agreement by accident.
  expect_false(isTRUE(all.equal(
    mean(d$x1[folds[[1L]]$train_rows]), mean(d$x1)
  )))
})


test_that("the horseshoe's global scale is calibrated at the fold's own units and rows", {
  skip_if_not_installed("brms")
  # Unbalanced on purpose: with the same number of rows in every site, a
  # fold's rows-per-unit equals the sample's and the fold-local ratio could
  # not be told from the whole-sample one.
  d <- fbfl_data()
  d <- d[!(d$site %in% c("S01", "S02", "S03", "S04") &
             duplicated(d$site) &
             stats::ave(seq_len(nrow(d)), d$site, FUN = seq_along) > 2L), ]
  rownames(d) <- NULL
  folds <- .fb_sel_folds(fbfl_resamples(), d, seed = 7L)
  learn <- .fb_select_fold_transform(
    ~ x1 + x2 + z1 + z2,
    fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
    character(0),
    fbfl_resamples()
  )

  # p0 = 2 of p = 4 gives the canonical ratio 1; moved to the unit it is
  # 1 * sqrt(rows / units) for THIS fold's rows and units.
  whole <- sqrt(nrow(d) / length(unique(as.character(d$site))))
  ratios <- numeric(0)
  for (k in seq_along(folds)) {
    train <- folds[[k]]$train_rows
    state <- learn(d[train, , drop = FALSE])
    h <- state$transform$blocks[[1L]]$hyperparameters
    units <- length(unique(as.character(d$site[train])))
    expect_equal(h$par_ratio_canonical, 1)
    expect_equal(h$par_ratio, sqrt(length(train) / units), tolerance = 1e-12)
    expect_match(
      paste(state$record$prior, collapse = " "),
      "horseshoe(", fixed = TRUE
    )
    ratios <- c(ratios, h$par_ratio)
  }
  # And on this design at least one fold's ratio is not the sample's.
  expect_true(any(abs(ratios - whole) > 1e-8))
})


test_that("the certificate's fold record flattens what each fold learned", {
  skip_if_not_installed("brms")
  run <- fbfl_run(fbfl_data())
  record <- .fb_select_fold_preprocessing(run$cv)
  expect_identical(record$learned_on, "analysis_rows")
  expect_identical(nrow(record$folds), length(run$cv$folds))
  expect_setequal(record$constants$fold, names(run$cv$preprocessing$by_fold))
  expect_setequal(record$constants$column, c("x1", "x2"))
  expect_true(all(grepl("target +=", record$folds$prior, fixed = TRUE)))

  # A run with no fold-local transform records nothing rather than an empty
  # claim.
  expect_null(.fb_select_fold_preprocessing(
    list(preprocessing = list(learned_on = "analysis_rows",
                              fold_local_transform = FALSE))
  ))
})
