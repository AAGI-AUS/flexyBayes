# =============================================================================
# offset() in the fixed formula.
#
# terms() files an offset apart from the term labels, and the fixed-formula
# parser read only the labels, so `Reaction ~ Days + offset(Days)` fitted
# `Reaction ~ Days` on INLA and on brms with no message: a Days coefficient
# of 10.47 where lme4 with the offset returns 9.47. The parser now records
# every offset on the model description and evaluates it against the data;
# a route that cannot carry one refuses by name.
# =============================================================================

.off_sleepstudy <- function() {
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  e$sleepstudy
}

.off_expect_refusal <- function(expr, code, pattern = NULL) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, paste0("flexybayes_refusal_", code))
  if (inherits(err, "condition") && !is.null(pattern)) {
    for (p in pattern) {
      expect_match(conditionMessage(err), p, fixed = TRUE)
    }
  }
  invisible(err)
}

# --- ingest ---------------------------------------------------------------- #

test_that("both grammars record the offsets on the model description", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  ir_a <- fb_from_asreml(
    fixed = Reaction ~ Days + offset(Days), random = ~ Subject, data = d
  )
  ir_b <- fb_from_brms(Reaction ~ Days + offset(Days) + (1 | Subject), data = d)
  expect_identical(ir_a$offsets, "offset(Days)")
  expect_identical(ir_b$offsets, "offset(Days)")
  # The offset is not a fixed term: it carries no coefficient.
  expect_identical(
    vapply(ir_a$fixed_terms, function(t) t$label, character(1L)),
    "Days"
  )
  expect_identical(ir_a$fixed_terms, ir_b$fixed_terms)
})

test_that("several offsets are recorded in formula order and add", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  ir <- fb_from_asreml(
    fixed = Reaction ~ offset(Days) + Days + offset(log(Days + 1)),
    random = ~ Subject,
    data = d
  )
  expect_identical(ir$offsets, c("offset(Days)", "offset(log(Days + 1))"))
  # Hand computation, not the helper's own arithmetic.
  expect_equal(
    .fb_offset_values(ir, d),
    d$Days + log(d$Days + 1),
    tolerance = 0
  )
  expect_null(.fb_offset_values(fb_from_asreml(Reaction ~ Days, data = d), d))
})

test_that("a column named offset is a covariate, not an offset", {
  skip_if_not_installed("lme4")
  d <- .off_sleepstudy()
  d$offset <- d$Days
  fb <- fb_from_asreml(fixed = Reaction ~ offset, random = ~ Subject, data = d)
  expect_identical(fb$fixed[[1L]]$var, "offset")
  expect_identical(fb$fixed[[1L]]$type, "continuous")
  expect_length(fb$offsets, 0L)
})

test_that("an offset that is not a finite number on every row refuses", {
  set.seed(3L)
  d <- data.frame(
    y = stats::rnorm(20),
    x = stats::rnorm(20),
    f = factor(rep(c("a", "b"), 10)),
    z = stats::rnorm(20)
  )
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(f), data = d),
    "fixed_offset_invalid",
    c("`offset(f)`", "factor")
  )
  d_na <- d
  d_na$z[c(4, 9)] <- NA
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(z), data = d_na),
    "fixed_offset_invalid",
    c("`offset(z)`", "row(s) 4, 9")
  )
  d_inf <- d
  d_inf$z[7] <- Inf
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(z), data = d_inf),
    "fixed_offset_invalid",
    c("not finite", "row(s) 7")
  )
  # A workspace object of the right length is not a column, and the
  # baseenv() enclosure keeps it from being read as one.
  zz <- stats::rnorm(20)
  .off_expect_refusal(
    fb_from_asreml(y ~ x + offset(zz), data = d),
    "fixed_offset_invalid",
    "`offset(zz)`"
  )
  .off_expect_refusal(
    fb_from_brms(y ~ x + offset(zz), data = d),
    "fixed_offset_invalid"
  )
})

test_that("offsets in the random or residual formula keep their refusal", {
  set.seed(4L)
  d <- data.frame(
    y = stats::rnorm(20),
    g = factor(rep(1:4, 5)),
    z = stats::rnorm(20)
  )
  .off_expect_refusal(
    fb_from_asreml(y ~ 1, random = ~ g + offset(z), data = d),
    "asreml_function_not_recognised"
  )
  .off_expect_refusal(
    fb_from_asreml(y ~ 1, residual = ~ offset(z), data = d),
    "asreml_function_not_recognised"
  )
})

test_that("a design cell completed under augment needs its offset", {
  # The completed row would otherwise carry row 1's offset. `z` varies, so
  # its value at an absent cell cannot be determined and the completion
  # refuses, as it does for any other model variable.
  d <- expand.grid(row = factor(1:3), col = factor(1:3))
  d$y <- stats::rnorm(9)
  d$z <- seq_len(9) / 10
  d <- d[-5, ]
  fb <- fb_from_asreml(
    y ~ offset(z), random = ~ ar1(row):ar1(col), data = d
  )
  expect_true("z" %in% .fb_model_vars(fb))
  err <- tryCatch(
    .fb_apply_na_action(fb, d, "augment"),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "z", fixed = TRUE)
})

# --- routes ---------------------------------------------------------------- #

test_that("a fitting route that does not carry an offset refuses by name", {
  skip_if_not_installed("lme4")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .off_sleepstudy()
  for (backend in c("auto", "inla", "brms")) {
    .off_expect_refusal(
      flexybayes(
        Reaction ~ Days + offset(Days),
        random = ~ Subject,
        data = d,
        backend = backend,
        verbose = FALSE
      ),
      "fixed_offset_unsupported",
      "`offset(Days)`"
    )
  }
  cells <- data.frame(
    f = factor(rep(sprintf("f%d", 1:4), length.out = 120)),
    g = factor(rep(sprintf("g%02d", 1:12), each = 10)),
    z = stats::rnorm(120)
  )
  cells$y <- stats::rnorm(120, 10) + cells$z
  .off_expect_refusal(
    suppressMessages(flexybayes_stream(
      y ~ f + offset(z),
      random = ~ g,
      source = cells,
      backend = "inla",
      fit = FALSE,
      verbose = FALSE
    )),
    "fixed_offset_unsupported",
    "`offset(z)`"
  )
})
