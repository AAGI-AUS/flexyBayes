# DEF-228, found by the 0.16.0 release gate's long tier. A stub stands in
# for the brms fit, so the test runs without the long tier.

test_that("DEF-228: proceed() fits through dispatch and keeps the call record", {
  d <- data.frame(y = c(1, 2, 3, 4), g = factor(c("a", "a", "b", "b")))
  v <- c(0.1, 0.2, 0.1, 0.2)
  dispatched <- NULL
  testthat::local_mocked_bindings(
    .dispatch_backend = function(...) {
      dispatched <<- list(...)
      structure(
        list(
          extras = list(call_info = list(backend = "brms"), fb_terms = "terms")
        ),
        class = c("flexybayes_brms", "flexybayes")
      )
    },
    .fb_fixed_alias_attach = function(fit, fb) fit
  )
  review <- flexyBayes:::.new_flexybayes_review(
    code = "model {}",
    backend = "stan_via_brms",
    ir = structure(list(), class = "fb_terms"),
    prior = NULL,
    data_name = "d",
    call = quote(flexybayes(y ~ 1)),
    seed = NULL,
    proceed_args = list(
      data = d,
      backend = "brms",
      aggregate = "auto",
      na_action = "augment"
    ),
    finish = list(
      na_meta = list(policy = "augment"),
      na_action = "augment",
      known_variance = v
    )
  )
  fit <- suppressMessages(proceed(review))
  expect_identical(dispatched$backend, "brms")
  expect_identical(dispatched$aggregate, "auto")
  expect_identical(fit$extras$fb_terms, "terms")
  expect_identical(fit$extras$call_info$known_variance, v)
  expect_identical(fit$extras$call_info$na_action, "augment")
  expect_identical(fit$extras$na_action, list(policy = "augment"))
})
