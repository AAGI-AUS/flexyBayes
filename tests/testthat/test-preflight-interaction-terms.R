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
  # INLA cannot represent the crossed interaction random effect, so the
  # honest structural refusal is the correct outcome -- what matters is
  # that it is a refusal, not a base-R crash.
  expect_false(p$will_fit)
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
