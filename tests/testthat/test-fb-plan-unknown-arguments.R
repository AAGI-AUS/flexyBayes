# DEF-196: fb_plan() refuses an argument it does not take, so a misspelt
# name cannot plan a different model from the one written.

.pua_data <- function() {
  set.seed(1)
  data.frame(
    y = stats::rpois(60, 3),
    x = stats::rnorm(60),
    g = factor(rep(1:6, 10))
  )
}

test_that("DEF-196: a misspelt argument is refused by name", {
  d <- .pua_data()
  err <- expect_error(
    fb_plan(y ~ x + (1 | g), data = d, familly = "poisson"),
    class = "flexybayes_refusal_plan_argument_unknown"
  )
  expect_match(conditionMessage(err), "`familly`", fixed = TRUE)
  expect_identical(err$arguments, "familly")
})

test_that("DEF-196: every unknown name is listed, and an unnamed extra is refused", {
  d <- .pua_data()
  err <- expect_error(
    fb_plan(y ~ x, random = ~g, data = d, fixed = y ~ x, seed = 1),
    class = "flexybayes_refusal_plan_argument_unknown"
  )
  expect_identical(err$arguments, c("fixed", "seed"))
  expect_error(
    fb_plan(
      y ~ x + (1 | g),
      d,
      NULL,
      NULL,
      "auto",
      NULL,
      list(),
      "gaussian",
      NULL,
      NULL,
      "auto",
      NULL,
      NULL,
      "auto",
      "extra"
    ),
    class = "flexybayes_refusal_plan_argument_unknown"
  )
})

test_that("DEF-196: the refusal code is registered and a valid call still plans", {
  expect_true("plan_argument_unknown" %in% fb_refusals()$reason_code)
  p <- fb_plan(y ~ x + (1 | g), data = .pua_data(), family = "poisson")
  expect_s3_class(p, "fb_plan")
})
