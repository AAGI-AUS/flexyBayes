# =============================================================================
# fb_select()'s decision stage: the submodel each fold is scored with is
# searched on that fold's reference fit, which saw its training rows only.
#
# The defect. With decision = "projection" and benchmark = TRUE, every fold's
# reference was projected onto the terms the whole-sample search chose, and
# that search read every row, each fold's assessment rows included, so the
# decision-stage score was not out of fold.
#
# The oracle is invariance. A response built to move the ranking is placed
# only in one fold's assessment rows. That fold's searched terms stay
# bit-identical, because its reference fit and its search read its training
# rows alone; the folds that train on the changed rows are allowed to move,
# and do, so the check is not vacuous.
#
# One brms reference refit per fold per model, twice: slow, and skipped on
# CRAN.
# =============================================================================

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  options(brms.backend = "cmdstanr")
}

test_that("the decision stage searches each fold's submodel terms on that fold's training rows", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("projpred")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE
  )
  d <- fb_sel_bayes_data()
  base <- fb_sel_bayes_memo("fold_search_base", fb_sel_bayes_base(d))
  resamples <- fb_resamples(unit = ~site, v = 3L)
  seed <- 20260912L
  folds <- .fb_sel_folds(resamples, d, seed = seed)
  held <- folds[[1L]]$assess_rows
  # On the first fold's assessment rows alone, a response driven by z2. The
  # rows the whole sample holds then rank z2 ahead of x1.
  extreme <- d
  extreme$y[held] <- extreme$y[held] + 8 * extreme$z2[held]

  search <- list(
    size = 2L,
    forced = c("p1", "(1 | site)"),
    method = "forward",
    nterms_max = 5L,
    nclusters = .FB_SEL_PROJECTION_NCLUSTERS
  )
  run <- function(frame) {
    b <- base
    b$data <- frame
    .fb_select_reference_cv(
      base = b,
      candidates = ~ x1 + x2 + z1 + z2,
      penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
      protected_labels = .fb_sel_protected_labels(base$fb),
      resamples = resamples,
      prediction = "new_cluster",
      metric = "elpd",
      engine = fb_engine("brms", chains = 2L, n_samples = 1500L, warmup = 2000L),
      seed = seed,
      decision_search = search
    )
  }
  clean <- suppressWarnings(run(d))
  moved <- suppressWarnings(run(extreme))
  terms_of <- function(cv) {
    lapply(cv$derived_terms, function(x) x$projected_submodel)
  }
  label <- folds[[1L]]$label
  others <- setdiff(names(terms_of(clean)), label)

  expect_length(terms_of(clean), length(folds))
  # The forced block and one searched term, on every fold.
  for (terms in terms_of(clean)) {
    expect_true(all(search$forced %in% terms))
    expect_length(terms, length(search$forced) + 1L)
  }
  expect_identical(terms_of(moved)[[label]], terms_of(clean)[[label]])
  expect_false(identical(terms_of(moved)[others], terms_of(clean)[others]))
  expect_true("projected_submodel" %in% clean$clusters$model)
})
