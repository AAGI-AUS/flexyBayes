# The residual-structure argument was `rcov` in ASReml 3 and `residual` in
# ASReml 4. flexyBayes adopted the ASReml 4 name in 0.9.0; `rcov` is retained
# only as a tripwire that raises a guiding error. These tests pin that contract
# across the entry surface (universal entry, ingest adapter, engine pins).

.residual_df <- function() {
  data.frame(
    y = rnorm(20),
    g = factor(rep(letters[1:4], 5)),
    env = factor(rep(c("a", "b"), 10))
  )
}

test_that("residual = parses into the residual_terms IR field", {
  df <- .residual_df()
  ir <- fb_from_asreml(y ~ env, random = ~ g, residual = ~ units, data = df)
  expect_true("residual_terms" %in% names(unclass(ir)))
  expect_identical(ir$residual_terms, list(list(type = "units")))
})

test_that("an unstated residual defaults to iid units", {
  df <- .residual_df()
  ir <- fb_from_asreml(y ~ env, data = df)
  expect_identical(ir$residual_terms, list(list(type = "units")))
})

test_that("the defunct `rcov` argument errors, pointing to `residual`", {
  df <- .residual_df()

  expect_error(
    fb_from_asreml(y ~ env, rcov = ~ units, data = df),
    "residual"
  )
  expect_error(
    flexybayes(y ~ env, rcov = ~ units, data = df, backend = "inla"),
    "residual"
  )
  # Engine pins forward to flexybayes(), so the tripwire covers them too.
  expect_error(
    fb_inla(y ~ env, rcov = ~ units, data = df),
    "residual"
  )
})
