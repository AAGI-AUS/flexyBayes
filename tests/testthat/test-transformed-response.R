# A transformed response -- log(yield) -- is what every peer accepts and what
# flexyBayes refused (issue #14, fit1). These tests hold the fix and, more
# importantly, hold the two ways it can go wrong: accepting something it should
# refuse, and resolving the response outside the data.

test_that("a transformed response parses and carries its values", {
  tab <- data.frame(
    yield = c(10, 20, 40, 80, 160, 320),
    env = factor(rep(c("a", "b"), each = 3L)),
    gen = factor(rep(c("g1", "g2", "g3"), times = 2L))
  )
  fb <- fb_from_asreml(
    fixed = log(yield) ~ env,
    random = ~gen,
    data = tab,
    family = "gaussian"
  )
  expect_identical(fb$response, "log(yield)")
  expect_true(isTRUE(fb$response_is_expression))
  expect_equal(.fb_response_values(fb, tab), log(tab$yield))
})

test_that("a plain column response is unchanged and not flagged", {
  tab <- data.frame(
    yield = c(1, 2, 3, 4),
    gen = factor(c("a", "b", "a", "b"))
  )
  fb <- fb_from_asreml(
    fixed = yield ~ 1,
    random = ~gen,
    data = tab,
    family = "gaussian"
  )
  expect_identical(fb$response, "yield")
  expect_false(isTRUE(fb$response_is_expression))
  expect_equal(.fb_response_values(fb, tab), tab$yield)
})

test_that("a mistyped response still refuses by name", {
  tab <- data.frame(yield = c(1, 2, 3, 4), gen = factor(c("a", "b", "a", "b")))
  expect_error(
    fb_from_asreml(
      fixed = yeild ~ 1,
      random = ~gen,
      data = tab,
      family = "gaussian"
    ),
    "not found in data"
  )
})

test_that("a workspace object of the right length is not used as the response", {
  # The guard's whole point. The first version of this fix evaluated the
  # left-hand side with asNamespace("stats") as the enclosure, which does not
  # exclude the workspace -- every namespace's parent chain reaches globalenv
  # -- and a stray vector was silently accepted as the response. baseenv()
  # terminates at the empty environment. This test fails against that version.
  tab <- data.frame(yield = c(1, 2, 3, 4), gen = factor(c("a", "b", "a", "b")))
  notacolumn <- c(9, 9, 9, 9)
  expect_error(
    fb_from_asreml(
      fixed = notacolumn ~ 1,
      random = ~gen,
      data = tab,
      family = "gaussian"
    ),
    "not found in data"
  )
})

test_that("missing responses are found through a transformation", {
  # The consumers guarded with `resp %in% names(data)` and returned
  # integer(0) for a transformed response, which reports "no missing
  # responses" on data that has them.
  tab <- data.frame(
    yield = c(10, NA, 40, 80, NA, 320),
    env = factor(rep(c("a", "b"), each = 3L)),
    gen = factor(rep(c("g1", "g2", "g3"), times = 2L))
  )
  fb <- fb_from_asreml(
    fixed = log(yield) ~ env,
    random = ~gen,
    data = tab,
    family = "gaussian"
  )
  expect_identical(.fb_response_missing(fb, tab), c(2L, 5L))
})

test_that("streaming refuses a transformed response by name", {
  tab <- data.frame(
    yield = c(10, 20, 40, 80),
    gen = factor(c("a", "b", "a", "b"))
  )
  fb <- fb_from_asreml(
    fixed = log(yield) ~ 1,
    random = ~gen,
    data = tab,
    family = "gaussian"
  )
  expect_error(
    .fb_stream_aggregate(fb, src = tab, trials = NULL, exposure = NULL,
                         verbose = FALSE),
    "expression"
  )
})
