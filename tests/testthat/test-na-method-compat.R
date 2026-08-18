# =============================================================================
# Accepting the missingness policy an ASReml user already writes.
#
# ASReml states two decisions through one object, `na.method(y = , x = )`:
# what to do with a missing response and what to do with a missing
# covariate. flexyBayes accepts that object, the bare list a reader
# without an asreml licence can write by hand, and its own native
# strings.
#
# The oracle for every shape assertion here is helper-asreml-shapes.R,
# recorded once from the live licensed asreml 4.2.0.392 install. asreml
# is deliberately absent from DESCRIPTION -- it is commercial and
# licensed per seat -- so every block below runs with asreml uninstalled
# and never calls it.
#
# Two recorded properties drive the normaliser, and neither is in the
# ASReml-R manual:
#
#   1. The value is a plain list with no class, so detection is by shape.
#   2. An unsupplied argument comes back as its whole default vector, so
#      the effective policy is the first element.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# 1. Shape detection                                                #
# ---------------------------------------------------------------- #

test_that("every recorded na.method() value is detected by shape", {
  cases <- c(
    "explicit", "default", "y_include", "y_omit", "y_fail",
    "x_omit", "x_include"
  )
  for (case in cases) {
    m <- .asreml_na_method_recorded(case)
    expect_true(
      flexyBayes:::.fb_is_na_method_shape(m),
      label = paste0("shape detection on the recorded case '", case, "'")
    )
    # The frozen detection rule from the Phase 0 name table holds for
    # every asreml value: a list carrying both `x` and `y`.
    expect_true(is.list(m) && all(c("x", "y") %in% names(m)), label = case)
  }
})

test_that("shape detection does not fire on things that are not policies", {
  expect_false(flexyBayes:::.fb_is_na_method_shape("augment"))
  expect_false(flexyBayes:::.fb_is_na_method_shape(list()))
  expect_false(flexyBayes:::.fb_is_na_method_shape(list("include")))
  expect_false(flexyBayes:::.fb_is_na_method_shape(list(z = "include")))
  expect_false(
    flexyBayes:::.fb_is_na_method_shape(list(y = "include", z = 1))
  )
})

test_that("detection is by shape and never by class", {
  # The recorded value carries no class of its own, so there is nothing
  # for inherits() to match. Giving it a class must not change the
  # answer either way.
  m <- .asreml_na_method_recorded("explicit")
  expect_identical(class(m), "list")
  expect_true(flexyBayes:::.fb_is_na_method_shape(m))

  disguised <- m
  class(disguised) <- c("something_else", "list")
  expect_true(flexyBayes:::.fb_is_na_method_shape(disguised))
})


# ---------------------------------------------------------------- #
# 2. The mapping, against the recorded oracle                       #
# ---------------------------------------------------------------- #

test_that("the recorded na.method() values map to the frozen policy", {
  expected <- list(
    explicit = list(y = "augment", x = "fail"),
    default = list(y = "augment", x = "fail"),
    y_include = list(y = "augment", x = "fail"),
    y_omit = list(y = "omit", x = "fail"),
    y_fail = list(y = "fail", x = "fail"),
    x_omit = list(y = "augment", x = "omit"),
    x_include = list(y = "augment", x = "include")
  )
  for (case in names(expected)) {
    got <- flexyBayes:::.fb_normalise_na_action(
      .asreml_na_method_recorded(case)
    )
    expect_identical(got, expected[[case]], label = case)
  }
})

test_that("an unreduced default vector is read as its first element", {
  # The trap: na.method(y = "include") returns x as the length-3 vector
  # c("fail", "include", "omit"). A normaliser that compared the slot to
  # a scalar would read a vector, and one that took the last element
  # would silently switch the policy to "omit".
  partial <- .asreml_na_method_recorded("y_include")
  expect_length(partial$x, 3L)
  expect_identical(
    flexyBayes:::.fb_normalise_na_action(partial)$x,
    "fail"
  )

  full <- .asreml_na_method_recorded("default")
  expect_length(full$y, 3L)
  expect_length(full$x, 3L)
  got <- flexyBayes:::.fb_normalise_na_action(full)
  expect_identical(got, list(y = "augment", x = "fail"))

  # And the first elements are exactly asreml's own documented defaults.
  defs <- .asreml_na_method_defaults()
  expect_identical(defs$y[[1L]], "include")
  expect_identical(defs$x[[1L]], "fail")
})

test_that("the fixture names the asreml build it was recorded from", {
  expect_identical(.asreml_na_method_version(), "4.2.0.392")
})


# ---------------------------------------------------------------- #
# 3. The other two accepted spellings                               #
# ---------------------------------------------------------------- #

test_that("native strings normalise to themselves", {
  for (word in c("augment", "omit", "fail")) {
    expect_identical(
      flexyBayes:::.fb_normalise_na_action(word),
      list(y = word, x = "fail"),
      label = word
    )
  }
  # The unsupplied default -- flexybayes()'s own `c("augment", "omit",
  # "fail")` formal -- reduces the way match.arg() would.
  expect_identical(
    flexyBayes:::.fb_normalise_na_action(c("augment", "omit", "fail")),
    list(y = "augment", x = "fail")
  )
})

test_that("a hand-written list needs no asreml and may omit a slot", {
  expect_identical(
    flexyBayes:::.fb_normalise_na_action(list(y = "include", x = "fail")),
    list(y = "augment", x = "fail")
  )
  expect_identical(
    flexyBayes:::.fb_normalise_na_action(list(y = "omit")),
    list(y = "omit", x = "fail")
  )
  expect_identical(
    flexyBayes:::.fb_normalise_na_action(list(x = "omit")),
    list(y = "augment", x = "omit")
  )
})

test_that("normalising is idempotent", {
  # flexybayes() normalises once for the call record and hands the same
  # value to the layer, which normalises again. The second pass must be
  # a no-op or the recorded policy and the applied policy could differ.
  for (case in c("explicit", "y_omit", "x_omit", "x_include")) {
    once <- flexyBayes:::.fb_normalise_na_action(
      .asreml_na_method_recorded(case)
    )
    expect_identical(
      flexyBayes:::.fb_normalise_na_action(once),
      once,
      label = case
    )
  }
})

test_that("an unrecognised policy word is named, not defaulted", {
  expect_error(
    flexyBayes:::.fb_normalise_na_action(list(y = "drop")),
    "not a recognised `y` policy"
  )
  expect_error(
    flexyBayes:::.fb_normalise_na_action(list(y = "include", x = "zero")),
    "not a recognised `x` policy"
  )
  expect_error(
    flexyBayes:::.fb_normalise_na_action(list(policy = "include")),
    "asreml::na.method",
    fixed = TRUE
  )
  expect_error(
    flexyBayes:::.fb_normalise_na_action("include"),
    "should be one of"
  )
})


# ---------------------------------------------------------------- #
# 4. End to end through the layer, no engine required               #
# ---------------------------------------------------------------- #

.na_method_grid <- function(n_row = 4L, n_col = 3L) {
  set.seed(11L)
  g <- expand.grid(
    row = factor(seq_len(n_row)),
    col = factor(seq_len(n_col)),
    KEEP.OUT.ATTRS = FALSE
  )
  g$x <- stats::rnorm(nrow(g))
  g$y <- stats::rnorm(nrow(g))
  g
}

.na_method_ir <- function(with_covariate = TRUE) {
  flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    link = "identity",
    fixed_terms = if (with_covariate) {
      list(list(type = "continuous", var = "x"))
    } else {
      list()
    },
    random_terms = list(),
    residual_terms = list(list(type = "units")),
    source = "asreml"
  )
}

test_that("the asreml value and the native string do the same thing", {
  g <- .na_method_grid()
  g$y[c(2L, 7L)] <- NA
  ir <- .na_method_ir()

  native <- flexyBayes:::.fb_apply_na_action(ir, g, "augment")
  asreml_shaped <- flexyBayes:::.fb_apply_na_action(
    ir, g, .asreml_na_method_recorded("explicit")
  )
  expect_identical(native$data, asreml_shaped$data)
  expect_identical(native$meta, asreml_shaped$meta)
  expect_identical(native$meta$na_action, "augment")
})

test_that("y = omit drops the missing-response rows", {
  g <- .na_method_grid()
  g$y[c(2L, 7L)] <- NA
  ir <- .na_method_ir()

  res <- flexyBayes:::.fb_apply_na_action(
    ir, g, .asreml_na_method_recorded("y_omit")
  )
  expect_equal(nrow(res$data), nrow(g) - 2L)
  expect_identical(res$meta$na_action, "omit")
  expect_equal(res$meta$n_missing_response, 2L)
  expect_equal(res$meta$n_observed, nrow(g) - 2L)
  expect_equal(res$meta$n_design, nrow(g) - 2L)
})

test_that("x = fail refuses a missing covariate by its own code", {
  g <- .na_method_grid()
  g$x[3L] <- NA
  ir <- .na_method_ir()

  err <- tryCatch(
    flexyBayes:::.fb_apply_na_action(
      ir, g, .asreml_na_method_recorded("explicit")
    ),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_missing_covariate_not_supported")
  msg <- conditionMessage(err)
  expect_match(msg, "x", fixed = TRUE)
  expect_match(msg, "ASReml", fixed = TRUE)
  # It names the remedy the user can type.
  expect_match(msg, "x = \"omit\"", fixed = TRUE)
})

test_that("x = include refuses the zero-fill by name", {
  g <- .na_method_grid()
  g$x[3L] <- NA
  ir <- .na_method_ir()

  err <- tryCatch(
    flexyBayes:::.fb_apply_na_action(
      ir, g, .asreml_na_method_recorded("x_include")
    ),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_covariate_zero_fill_not_supported")
  msg <- conditionMessage(err)
  # The message says what ASReml does, cites where it is written down,
  # and says what to do instead.
  expect_match(msg, "zero", fixed = TRUE)
  expect_match(msg, "section 3.11", fixed = TRUE)
  expect_match(msg, "x = \"omit\"", fixed = TRUE)
})

test_that("x = omit drops the covariate-NA rows and says so", {
  g <- .na_method_grid()
  g$x[c(3L, 5L)] <- NA
  ir <- .na_method_ir()

  w <- tryCatch(
    flexyBayes:::.fb_apply_na_action(
      ir, g, .asreml_na_method_recorded("x_omit")
    ),
    warning = function(w) w
  )
  expect_s3_class(w, "warning")
  msg <- conditionMessage(w)
  expect_match(msg, "dropped 2 row", fixed = TRUE)
  expect_match(msg, "x", fixed = TRUE)

  res <- suppressWarnings(flexyBayes:::.fb_apply_na_action(
    ir, g, .asreml_na_method_recorded("x_omit")
  ))
  expect_equal(nrow(res$data), nrow(g) - 2L)
  expect_false(anyNA(res$data$x))
  expect_identical(res$meta$na_covariate, "omit")
  expect_equal(res$meta$n_covariate_rows_dropped, 2L)
})

test_that("the covariate policy is a no-op when no covariate is missing", {
  g <- .na_method_grid()
  ir <- .na_method_ir()
  for (case in c("explicit", "x_omit", "x_include")) {
    res <- flexyBayes:::.fb_apply_na_action(
      ir, g, .asreml_na_method_recorded(case)
    )
    expect_equal(nrow(res$data), nrow(g), label = case)
    expect_equal(res$meta$n_covariate_rows_dropped, 0L, label = case)
  }
})

test_that("the recorded policy echoes what was asked for", {
  # Before 0.9.1 `fail` on complete data fell through into the augment
  # branch and recorded itself as `augment`, so the record was not a
  # faithful echo of the argument -- and anything reading it to repeat
  # the fit would have repeated a different one.
  g <- .na_method_grid()
  ir <- .na_method_ir()
  for (word in c("augment", "omit", "fail")) {
    res <- flexyBayes:::.fb_apply_na_action(ir, g, word)
    expect_identical(res$meta$na_action, word, label = word)
  }
})
