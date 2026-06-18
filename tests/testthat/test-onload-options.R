# .onLoad option side effects -- the marginaleffects custom-class hook.
#
# .onLoad() appends the flexyBayes fit classes to the
# `marginaleffects_model_classes` option so marginaleffects will dispatch on
# them. That is a global option mutation, so the contract is that it must be
# (a) non-clobbering -- it never drops a user's own registered classes -- and
# (b) idempotent -- running it again adds nothing and introduces no
# duplicates. The logic is factored into .register_marginaleffects_classes()
# precisely so these two properties are unit-testable without re-triggering a
# package load.

test_that(".register_marginaleffects_classes preserves a user's own classes", {
  withr::local_options(marginaleffects_model_classes = "user_special_class")

  flexyBayes:::.register_marginaleffects_classes()
  v <- getOption("marginaleffects_model_classes")

  expect_true("user_special_class" %in% v)
  expect_true(all(c("flexybayes", "flexybayes_inla") %in% v))
})

test_that(".register_marginaleffects_classes is idempotent (no growth, no dupes)", {
  withr::local_options(marginaleffects_model_classes = NULL)

  flexyBayes:::.register_marginaleffects_classes()
  v1 <- getOption("marginaleffects_model_classes")
  flexyBayes:::.register_marginaleffects_classes()
  v2 <- getOption("marginaleffects_model_classes")

  expect_identical(sort(v2), sort(v1))
  expect_equal(sum(v2 == "flexybayes"), 1L)
  expect_equal(sum(v2 == "flexybayes_inla"), 1L)
})

test_that(".register_marginaleffects_classes seeds cleanly from an unset option", {
  withr::local_options(marginaleffects_model_classes = NULL)

  flexyBayes:::.register_marginaleffects_classes()
  v <- getOption("marginaleffects_model_classes")

  expect_setequal(v, c("flexybayes", "flexybayes_inla"))
})
