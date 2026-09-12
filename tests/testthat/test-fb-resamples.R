# =============================================================================
# fb_resamples(): the resampling scheme as a public constructor.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 3.3, made
# public unchanged by SPEC_bayesian_selection_v5_2026-09-12.md section 3.3.
#
# The constructor's own refusals are exercised in test-selection-resamples.R,
# where the fold builder that reads them lives. What is tested here is the
# PUBLIC surface: that the object is exported, that it is what fb_select()
# takes by default, that its contract did not change when it became public,
# and that it prints what a reader needs in order to know which folds a
# reported number was computed on.
# =============================================================================


test_that("the constructor is exported and returns a classed object", {
  spec <- fb_resamples()
  expect_s3_class(spec, "fb_resamples")
  expect_true("fb_resamples" %in% getNamespaceExports("flexyBayes"))
  expect_identical(spec$scheme, "v_fold")
  expect_identical(spec$v, 10L)
  expect_identical(spec$repeats, 1L)
  expect_identical(spec$unit, character(0))
  expect_null(spec$fold_id)
})


test_that("it is the default fb_select() carries", {
  # The default is evaluated in fb_select()'s own frame, so a renamed or
  # unexported constructor would fail there and nowhere else.
  default <- formals(fb_select)$resamples
  expect_identical(deparse(default), "fb_resamples()")
  expect_s3_class(eval(default), "fb_resamples")
})


test_that("a unit is read from a formula or from a character vector", {
  expect_identical(fb_resamples(unit = ~site, v = 5L)$unit, "site")
  expect_identical(fb_resamples(unit = "site", v = 5L)$unit, "site")
  expect_identical(
    fb_resamples(unit = ~ site + year, v = 5L)$unit,
    c("site", "year")
  )
  expect_identical(fb_resamples(strata = ~region, v = 5L)$strata, "region")
})


test_that("a unit with no scheme is not a specification", {
  # The refusal that matters most in the public surface: a grouping unit
  # says what may not be split, not how the splitting is done, and the
  # unit-free defaults of v and repeats must not silently supply a scheme
  # to a grouped design.
  expect_error(
    fb_resamples(unit = ~site),
    class = "flexybayes_refusal_selection_resamples_underspecified"
  )
  expect_error(fb_resamples(unit = ~site), "no resampling scheme was given")
  # Negative control: naming any one of the three is enough.
  expect_s3_class(fb_resamples(unit = ~site, v = 5L), "fb_resamples")
  expect_s3_class(fb_resamples(unit = ~site, repeats = 2L), "fb_resamples")
  expect_s3_class(
    fb_resamples(unit = ~site, fold_id = rep(1:2, 10L)),
    "fb_resamples"
  )
})


test_that("an explicit partition and a scheme that would build one refuse", {
  expect_error(
    fb_resamples(v = 5L, fold_id = rep(1:2, 10L)),
    class = "flexybayes_refusal_selection_fold_id_invalid"
  )
  expect_error(
    fb_resamples(repeats = 2L, fold_id = rep(1:2, 10L)),
    "repeat the same folds"
  )
  expect_error(fb_resamples(fold_id = 1L), "one entry per row")
})


test_that("an explicit partition records the scheme it is", {
  spec <- fb_resamples(fold_id = rep(1:3, each = 4L))
  expect_identical(spec$scheme, "explicit")
  expect_true(is.na(spec$v))
  expect_length(spec$fold_id, 12L)
})


test_that("v is at least two folds", {
  expect_error(fb_resamples(v = 1L))
  expect_error(fb_resamples(v = 0L))
  expect_identical(fb_resamples(v = 2L)$v, 2L)
})


test_that("print says which folds a reported number was computed on", {
  out <- paste(
    utils::capture.output(print(fb_resamples(unit = ~site, v = 5L,
                                             repeats = 2L))),
    collapse = " "
  )
  expect_match(out, "site", fixed = TRUE)
  expect_match(out, "5 x 2 repeat(s)", fixed = TRUE)

  explicit <- paste(
    utils::capture.output(print(fb_resamples(fold_id = rep(1:3, each = 4L)))),
    collapse = " "
  )
  expect_match(explicit, "explicit", fixed = TRUE)
  expect_match(explicit, "3 explicit partition", fixed = TRUE)

  rows <- paste(
    utils::capture.output(print(fb_resamples())),
    collapse = " "
  )
  expect_match(rows, "unit    : row", fixed = TRUE)
})


test_that("the fold builder still reads the object it did before", {
  # The constructor moved file; the contract did not. This is the check
  # that the move was a move.
  set.seed(1L)
  d <- data.frame(
    site = factor(rep(sprintf("S%02d", seq_len(12L)), each = 4L)),
    y = stats::rnorm(48L)
  )
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d, seed = 7L)
  expect_length(folds, 4L)
  by_site <- tapply(
    rep(seq_along(folds), vapply(folds, function(f) length(f$assess_rows),
                                 integer(1L))),
    unlist(lapply(folds, function(f) as.character(d$site[f$assess_rows]))),
    function(f) length(unique(f))
  )
  expect_true(all(by_site == 1L))
})
