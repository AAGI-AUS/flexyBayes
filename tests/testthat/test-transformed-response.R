# A transformed response such as log(yield) (issue #14). The tests also cover
# the two ways it can go wrong: accepting something it should refuse, and
# resolving the response outside the data.

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
  # A namespace enclosure such as asNamespace("stats") does not exclude the
  # workspace, since every namespace's parent chain reaches globalenv.
  # baseenv() ends at the empty environment.
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
    .fb_stream_aggregate(
      fb,
      src = tab,
      trials = NULL,
      exposure = NULL,
      verbose = FALSE
    ),
    "expression"
  )
})
