# =============================================================================
# fb_select(protected_alias = "drop"): a protected column a fold's training
# rows do not identify is dropped as ASReml 4.2 drops it, and a held-out row
# with no estimable prediction is flagged and not scored. The default,
# "refuse", is unchanged.
#
# The fixture is an unbalanced multi-trial design (helper-alias-design.R):
# VARIETY x NAPPLIED x SEEDRATE with cells never sown and cells lost, random
# TRIAL and TRIAL:REP, folds that hold one locality out.
#
# Oracle: asreml 4.2 fitted to the same design, on all rows and on each
# fold's training rows, recorded in fixtures/asreml_alias_<shape>.rds by
# review/units/U-ALIAS_2026-09-24/asreml_alias_oracle.R: the columns it
# reports as aliased, its fixed estimates and variance components, and
# predict(classify = "VARIETY:NAPPLIED:SEEDRATE") with each cell's status.
# =============================================================================

pa_protected <- ~ VARIETY * NAPPLIED * SEEDRATE
pa_random <- c("TRIAL", "TRIAL:REP")

pa_oracle <- function(shape) {
  readRDS(test_path("fixtures", paste0("asreml_alias_", shape, ".rds")))
}

pa_resamples <- function(d) {
  fb_resamples(unit = ~LOCALITY, v = nlevels(d$LOCALITY))
}

pa_folds <- function(d) .fb_sel_folds(pa_resamples(d), d, seed = 1L)

pa_held <- function(d, split) {
  unique(as.character(d$LOCALITY[c(
    split$assess_rows,
    split$not_estimable_rows
  )]))
}

pa_cells <- function(d, rows) {
  paste(d$VARIETY, d$NAPPLIED, d$SEEDRATE, sep = "/")[rows]
}

pa_space <- function(d, folds, protected_alias = "drop") {
  .fb_sel_candidate_space(
    pa_protected,
    ~ x1 + x2,
    d,
    folds,
    groupings = c(pa_random, "LOCALITY"),
    route = "benchmark",
    protected_alias = protected_alias
  )
}

# Every fold's dropped columns and not-estimable cells against asreml's.
pa_expect_oracle <- function(shape) {
  d <- alias_design(shape)
  oracle <- pa_oracle(shape)
  folds <- pa_folds(d)
  whole <- list(list(
    fold_id = 1L,
    repeat_id = 1L,
    label = "all",
    train_rows = seq_len(nrow(d)),
    assess_rows = seq_len(nrow(d))
  ))
  all_rows <- pa_space(d, whole)$folds$all
  expect_identical(
    sort(all_rows$protected_dropped$column),
    oracle$subsets$all$aliased
  )
  expect_length(all_rows$not_estimable_rows, 0L)

  space <- pa_space(d, folds)
  for (split in folds) {
    held <- pa_held(d, split)
    f <- space$folds[[split$label]]
    expect_identical(
      sort(f$protected_dropped$column),
      oracle$subsets[[held]]$aliased
    )
    predictions <- oracle$subsets[[held]]$predictions
    aliased_cells <- predictions$cell[predictions$status != "Estimable"]
    expect_true(all(is.na(predictions$predicted[
      predictions$status != "Estimable"
    ])))
    expect_setequal(unique(pa_cells(d, f$not_estimable_rows)), aliased_cells)
    expect_identical(
      f$not_estimable_rows,
      split$assess_rows[pa_cells(d, split$assess_rows) %in% aliased_cells]
    )
  }
  space
}


# --- the rule, against asreml ------------------------------------------- #

test_that("each fold drops the protected columns asreml 4.2 reports as aliased, and flags the rows it cannot predict", {
  d <- alias_design("small")
  space <- pa_expect_oracle("small")
  # The fixture holds aliases that are not zero columns, which only the
  # column order decides between.
  folds <- pa_folds(d)
  x <- stats::model.matrix(pa_protected, d)
  nonzero <- vapply(
    folds,
    function(split) {
      dropped <- space$folds[[split$label]]$protected_dropped$column
      sum(colSums(abs(x[split$train_rows, dropped, drop = FALSE])) > 0)
    },
    integer(1L)
  )
  expect_true(all(nonzero >= 2L))
  expect_identical(
    space$not_estimable$n_not_estimable,
    vapply(
      folds,
      function(s) length(space$folds[[s$label]]$not_estimable_rows),
      integer(1L),
      USE.NAMES = FALSE
    )
  )
  expect_identical(sum(space$not_estimable$n_not_estimable), 30L)
  expect_identical(space$protected_alias, "drop")
})

test_that("on the study-shaped design each of 22 locality folds drops asreml's aliased columns and flags its unpredictable rows", {
  skip_if_not_long_tier()
  space <- pa_expect_oracle("study")
  expect_identical(sum(space$not_estimable$n_not_estimable), 27L)
})


# --- estimates and predictions, against asreml ------------------------- #

test_that("a fold's fixed estimates, variance components and held-out predictions agree with asreml", {
  d <- alias_design("small")
  oracle <- pa_oracle("small")
  cv <- .fb_sel_cv(
    d,
    "YIELD",
    pa_protected,
    ~x1,
    random = pa_random,
    resamples = pa_resamples(d),
    seed = 1L,
    protected_alias = "drop"
  )
  labels <- vapply(cv$folds, function(split) split$label, character(1L))
  for (s in cv$substrates) {
    held <- pa_held(d, cv$folds[[match(s$label, labels)]])
    a <- oracle$subsets[[held]]
    # Measured largest gaps: 5.0e-6 relative on a component, 6.4e-8 on an
    # estimate, 6.8e-8 on a prediction.
    components <- c(
      s$learned$components_terms,
      residual = s$learned$components_residual
    )[names(a$varcomp)]
    expect_equal(unname(components / a$varcomp), rep(1, 3L), tolerance = 5e-5)
    fit <- .fb_sel_gls_fit(
      s$blocks$factor_tt,
      s$x_protected[s$train_rows, , drop = FALSE],
      s$response[s$train_rows]
    )
    expect_setequal(colnames(s$x_protected), names(a$fixed))
    expect_equal(
      stats::setNames(fit$beta, colnames(s$x_protected)),
      a$fixed[colnames(s$x_protected)],
      tolerance = 1e-6
    )
    scored <- cv$rows[
      cv$rows$model == "protected_only" & cv$rows$fold == s$label,
    ]
    expect_identical(nrow(scored), length(s$assess_rows))
    predicted <- s$response[s$assess_rows] - scored$residual
    expect_equal(
      predicted,
      a$predictions$predicted[match(
        pa_cells(d, s$assess_rows),
        a$predictions$cell
      )],
      tolerance = 1e-6
    )
  }
})


# --- the default ------------------------------------------------------- #

test_that("the default still refuses a protected column the training rows do not identify", {
  d <- alias_design("small")
  err <- tryCatch(
    .fb_sel_candidate_space(
      pa_protected,
      ~ x1 + x2,
      d,
      pa_folds(d),
      groupings = c(pa_random, "LOCALITY"),
      route = "benchmark"
    ),
    error = function(e) e
  )
  expect_s3_class(
    err,
    "flexybayes_refusal_selection_fold_protected_unidentified"
  )
  expect_match(
    conditionMessage(err),
    "VARIETYV04:NAPPLIED30:SEEDRATE300",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_cv(
      d,
      "YIELD",
      pa_protected,
      ~x1,
      random = pa_random,
      resamples = pa_resamples(d),
      seed = 1L
    ),
    class = "flexybayes_refusal_selection_fold_protected_unidentified"
  )
})


test_that("the score contract leaves the flagged rows unscored and still maps each unit to its fold", {
  d <- alias_design("small")
  contract_of <- function(folds) {
    .fb_select_score_contract(
      folds = folds,
      data = d,
      response = "YIELD",
      prediction = "new_cluster",
      metric = "nlpd",
      design = matrix(1, nrow(d), 1L, dimnames = list(NULL, "(Intercept)")),
      learned_on = "analysis_rows",
      stage = "reference"
    )
  }
  folds <- .fb_sel_estimable_folds(pa_folds(d), pa_space(d, pa_folds(d)))
  contract <- contract_of(folds)
  not_estimable <- unlist(lapply(folds, function(s) s$not_estimable_rows))
  expect_length(not_estimable, 30L)
  expect_true(all(is.na(contract$row_fold[not_estimable, 1L])))
  expect_false(anyNA(contract$row_fold[-not_estimable, 1L]))
  expect_identical(contract$unit_to_fold$fold, seq_len(6L))
  # A unit whose first row is not scored keeps its fold.
  first <- folds[[1L]]$assess_rows[[1L]]
  folds[[1L]]$assess_rows <- folds[[1L]]$assess_rows[-1L]
  folds[[1L]]$not_estimable_rows <- c(first, folds[[1L]]$not_estimable_rows)
  expect_identical(contract_of(folds)$unit_to_fold$fold, seq_len(6L))
})


# --- through fb_select() ----------------------------------------------- #

pa_cache <- new.env(parent = emptyenv())

pa_base <- function() {
  if (!exists("base", envir = pa_cache, inherits = FALSE)) {
    withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
    assign(
      "base",
      flexybayes(
        YIELD ~ VARIETY * NAPPLIED * SEEDRATE,
        random = ~ TRIAL + TRIAL:REP,
        data = alias_design("small"),
        backend = "inla",
        verbose = FALSE
      ),
      envir = pa_cache
    )
  }
  get("base", envir = pa_cache, inherits = FALSE)
}

test_that("fb_select() refuses by default and, with protected_alias = \"drop\", scores every fold on the estimable rows", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- alias_design("small")
  oracle <- pa_oracle("small")
  base <- pa_base()
  expect_error(
    fb_select(
      base,
      ~ x1 + x2 + x3 + x4,
      workflow = "unpenalized_benchmark",
      resamples = pa_resamples(d),
      seed = 1L
    ),
    class = "flexybayes_refusal_selection_fold_protected_unidentified"
  )

  sel <- fb_select(
    base,
    ~ x1 + x2 + x3 + x4,
    workflow = "unpenalized_benchmark",
    resamples = pa_resamples(d),
    seed = 1L,
    protected_alias = "drop"
  )
  cs <- sel$certificate$design$candidate_space
  expect_identical(cs$protected_alias, "drop")
  for (split in sel$folds) {
    held <- pa_held(d, split)
    expect_identical(
      sort(cs$protected_dropped$column[
        cs$protected_dropped$fold == split$label
      ]),
      oracle$subsets[[held]]$aliased
    )
    predictions <- oracle$subsets[[held]]$predictions
    aliased_cells <- predictions$cell[predictions$status != "Estimable"]
    expect_setequal(
      unique(pa_cells(d, split$not_estimable_rows)),
      aliased_cells
    )
    expect_false(any(pa_cells(d, split$assess_rows) %in% aliased_cells))
    counts <- cs$not_estimable[cs$not_estimable$fold == split$label, ]
    expect_identical(counts$n_not_estimable, length(split$not_estimable_rows))
    expect_identical(
      counts$n_assess,
      length(split$assess_rows) + length(split$not_estimable_rows)
    )
    expect_identical(
      sum(sel$clusters$n[
        sel$clusters$model == "full" & sel$clusters$fold == split$label
      ]),
      length(split$assess_rows)
    )
  }
  expect_identical(sum(cs$not_estimable$n_not_estimable), 30L)
  expect_true(all(is.finite(sel$scores$score)))
  expect_match(
    paste(utils::capture.output(print(sel$certificate)), collapse = "\n"),
    "30 of 914 held-out rows not estimable and not scored",
    fixed = TRUE
  )
})

test_that("the penalised optimiser with protected_alias = \"drop\" predicts only the estimable held-out rows", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  d <- alias_design("small")
  sel <- fb_select(
    pa_base(),
    ~ x1 + x2 + x3 + x4,
    penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
    workflow = "penalized_optimizer",
    resamples = pa_resamples(d),
    seed = 1L,
    protected_alias = "drop"
  )
  scored <- sort(unlist(lapply(sel$folds, function(s) s$assess_rows)))
  expect_identical(sort(sel$assessment_predictions$row), scored)
  expect_length(scored, nrow(d) - 30L)
  expect_true("x1" %in% sel$selection$chosen_terms)
  expect_identical(
    sum(sel$certificate$design$candidate_space$not_estimable$n_not_estimable),
    30L
  )
})

test_that("protected_alias = \"drop\" refuses on the routes that write each fold as formula text", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- alias_design("small")
  expect_error(
    fb_select(
      pa_base(),
      ~ x1 + x2,
      resamples = pa_resamples(d),
      seed = 1L,
      protected_alias = "drop"
    ),
    "the Bayesian workflow writes each fold's model as formula text",
    fixed = TRUE
  )
  expect_error(
    fb_select(
      pa_base(),
      ~ x1 + x2,
      penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
      workflow = "penalized_optimizer",
      refit = list(method = "bayes"),
      resamples = pa_resamples(d),
      seed = 1L,
      protected_alias = "drop"
    ),
    "refit = list(method = \"bayes\") writes each fold's model",
    fixed = TRUE
  )
  expect_error(
    fb_select(pa_base(), ~x1, protected_alias = "keep"),
    "should be one of"
  )
})

test_that("a base-model fold refit holds out the not-estimable rows with the scored ones", {
  d <- alias_design("small")
  folds <- .fb_sel_estimable_folds(pa_folds(d), pa_space(d, pa_folds(d)))
  split <- folds[[which(vapply(
    folds,
    function(s) length(s$not_estimable_rows) > 0L,
    logical(1L)
  ))[[1L]]]]
  handed <- NULL
  local_mocked_bindings(
    .fb_sel_bayes_base_spec = function(base) {
      list(
        data = d,
        response = "YIELD",
        protected_labels = "VARIETY",
        offsets = character(0)
      )
    },
    .fb_sel_bayes_refit = function(base, data, ...) {
      handed <<- data
      stop("the refit was reached", call. = FALSE)
    }
  )
  expect_error(
    .fb_sel_cov_fold_refit(NULL, split, "inla"),
    "the refit was reached",
    fixed = TRUE
  )
  held <- setdiff(seq_len(nrow(d)), split$train_rows)
  expect_true(all(is.na(handed$YIELD[held])))
  expect_false(anyNA(handed$YIELD[split$train_rows]))
})
