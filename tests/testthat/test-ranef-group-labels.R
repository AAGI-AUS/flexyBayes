# ranef() names the group the user wrote, not the index INLA was given.
#
# emit_inla() builds an index column per structured random term --
# `rep_block_id` for `rep:block`, `gen_by_env_id` for an unstructured
# covariance -- and INLA names each summary.random entry after the column it
# was indexed by. The variance components were mapped back to user labels
# when that machinery was built; the random-effect table was not, so
# ranef(fit) came back with a grouping factor called `rep_block_id`. Same
# defect class as the label leak closed for summary()$varcomp, reached
# through a different accessor.

test_that("a nested term keeps its own name in ranef()", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  data(john.alpha, package = "agridat")
  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~ rep:block + gen,
    data = john.alpha, backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)
  expect_true("rep:block" %in% names(r))
  expect_true("gen" %in% names(r))
  expect_false(any(grepl("_id$", names(r))))
  # And inside the table, not only on the list name.
  expect_identical(unique(r[["rep:block"]]$group), "rep:block")
})

test_that("a plain grouping factor is left alone", {
  # The map must not rename a group whose index IS the variable the user
  # named, which is the common case and the one a blanket rewrite would
  # break.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  data(john.alpha, package = "agridat")
  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~gen, data = john.alpha,
    backend = "inla", verbose = FALSE
  ))
  expect_identical(names(ranef(fit)), "gen")
})

test_that("the label map covers every structured term type it claims", {
  # Unit-level, so it runs without INLA: build the IR for each structured
  # term and check the generated index name maps to something a reader
  # would recognise.
  d <- expand.grid(rep = factor(1:2), gen = factor(1:6), env = factor(1:3))
  d$y <- stats::rnorm(nrow(d))
  d$x <- stats::rnorm(nrow(d))

  cases <- list(
    list(random = ~ rep:gen, index = "rep_gen_id", label = "rep:gen"),
    list(random = ~ diag(env):gen, index = "gen_at_env_1",
         label = "gen (env = 1)"),
    list(random = ~ us(env):gen, index = "gen_by_env_id", label = "gen:env")
  )
  for (cs in cases) {
    fb <- fb_from_asreml(y ~ 1, random = cs$random, data = d)
    map <- .fb_inla_group_labels(fb)
    expect_true(cs$index %in% names(map), info = cs$index)
    expect_identical(unname(map[[cs$index]]), cs$label)
  }
})

test_that("the correlated slope's two index blocks share one label", {
  d <- data.frame(y = stats::rnorm(40), x = stats::rnorm(40),
                  g = factor(rep(letters[1:8], each = 5)))
  fb <- fb_from_brms(y ~ x + (x | g), data = d)
  map <- .fb_inla_group_labels(fb)
  expect_identical(unname(map[["g_cor_x_id"]]), "x | g")
  expect_identical(unname(map[["g_cor_x_slope_id"]]), "x | g")
})
