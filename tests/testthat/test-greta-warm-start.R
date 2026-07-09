# PR4: the greta emit warm-starts the intercept from the response mean
# (Francis K. C. Hui's index-approach refinement). The absent-cell half of
# that refinement is already subsumed by the emit's factor(interaction(...))
# ids, which only allocate a coefficient for an observed combination, so no
# code is needed there -- see setup_env.R.

test_that(".greta_intercept_init resolves link-scale starting values", {
  # gaussian: the raw mean
  expect_equal(
    .greta_intercept_init(c(1, 2, 3, 4, 5), list(family = "gaussian", link = "identity")),
    3
  )
  # poisson (log link): log of the mean
  expect_equal(
    .greta_intercept_init(c(1, 2, 3), list(family = "poisson", link = "log")),
    log(2)
  )
  # binomial (logit): qlogis of the mean; mean 0.5 -> 0
  expect_equal(
    .greta_intercept_init(c(0, 1, 0, 1), list(family = "binomial", link = "logit")),
    0
  )
  # unusable response -> NULL (no warm start)
  expect_null(
    .greta_intercept_init(c(NA, NaN), list(family = "gaussian", link = "identity"))
  )
})

test_that("greta emit warm-starts the intercept only when an intercept is present", {
  skip_if_no_greta()
  d <- data.frame(y = c(10, 12, 11, 13, 9, 14), x = c(1, 2, 3, 4, 5, 6))

  code <- flexybayes(
    y ~ x, data = d, backend = "greta", return_code = TRUE, verbose = FALSE
  )
  expect_match(code, "initial_values = greta::initials(mu_atg =", fixed = TRUE)

  # a no-intercept model carries no mu_atg, hence no warm start
  code0 <- flexybayes(
    y ~ 0 + x, data = d, backend = "greta", return_code = TRUE, verbose = FALSE
  )
  expect_false(grepl("initial_values", code0, fixed = TRUE))
})
