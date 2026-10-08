# Regression tests for the preflight handling of interaction random terms.
#
# A three-way interaction random term (e.g. gen:loc:yearf, IR type "combo")
# used to abort the preflight with two successive un-vectorised base-R
# errors: first "the condition has length > 1" from `.fb_dataset_levels()`
# (an `if (col %in% ...)` on the multi-column interaction), then
# "no such index at level 1" from a vector `.preflight_term_label()` used as
# a `per_term[[label]] <-` index. Both are fixed so an interaction random
# term reaches the graceful structural-refusal gate instead of crashing.

test_that("a three-way interaction random term preflights without crashing", {
  set.seed(1)
  d <- expand.grid(a = factor(1:4), b = factor(1:3), cc = factor(1:2), rep = 1:3)
  d$y <- rnorm(nrow(d))

  p <- flexybayes(y ~ 1, random = ~ a:b:cc, data = d, plan = TRUE)

  expect_s3_class(p, "fb_plan")
  # This test was written when a crossed interaction random effect refused
  # on INLA, and it asserted the refusal. The claim behind that refusal --
  # that INLA cannot represent the term class -- was false, so from
  # 0.10.0.9000 the gate accepts it and `auto` resolves to INLA. What the
  # test exists for is unchanged: a three-way interaction must preflight to
  # a typed plan rather than a base-R crash.
  expect_identical(p$gate_outcome, "accept")
  # `auto` resolves to INLA only where INLA is installed.
  skip_if_not_installed("INLA")
  expect_identical(p$backend_chosen, "inla")
  expect_true(p$will_fit)
})

test_that(".fb_dataset_levels counts distinct observed combinations for interactions", {
  d <- data.frame(
    a = factor(c(1, 1, 2, 2, 3)),
    b = factor(c(1, 2, 1, 2, 1))
  )
  ds <- .fb_dataset(d)

  # single column: the frozen-dictionary level count (unchanged behaviour)
  expect_identical(.fb_dataset_levels(ds, "a"), 3L)

  # interaction: distinct observed (a, b) combinations
  # {1:1, 1:2, 2:1, 2:2, 3:1} = 5, matching nlevels(factor(interaction(a, b)))
  expect_identical(.fb_dataset_levels(ds, c("a", "b")), 5L)

  # metadata-only dataset (no $data): NA_integer_, which the preflight
  # escalates to a graceful representation_unknown_for_preflight refusal
  meta <- .fb_dataset_metadata(ds)
  expect_true(is.na(.fb_dataset_levels(meta, c("a", "b"))))
})

test_that(".preflight_term_label collapses an interaction group to a scalar", {
  term <- list(type = "combo", vars = c("a", "b", "cc"))
  lab <- .preflight_term_label(term, "combo", kind = "random")

  expect_length(lab, 1L)
  expect_match(lab, "a:b:cc", fixed = TRUE)
})

test_that("DEF-204: the plan labels an AR1 field by its spelling", {
  d <- expand.grid(
    row = factor(1:4), col = factor(1:3), trial = factor(c("A", "B"))
  )
  d$y <- stats::rnorm(nrow(d))
  one <- d[d$trial == "A", ]
  for (case in list(
    list(random = ~ ar1(row):ar1(col), data = one,
         label = "ar1(row):ar1(col)"),
    list(random = ~ ar1(row), data = one[one$col == "1", ],
         label = "ar1(row)"),
    list(random = ~ id(trial):ar1(row):ar1(col), data = d,
         label = "id(trial):ar1(row):ar1(col)"),
    list(random = ~ at(trial):ar1(row):ar1(col), data = d,
         label = "at(trial):ar1(row):ar1(col)")
  )) {
    p <- fb_plan(y ~ 1, random = case$random, data = case$data)
    expect_true(case$label %in% names(p$representation_plan),
                label = case$label)
    expect_false(any(grepl("?", names(p$representation_plan), fixed = TRUE)))
  }
})
