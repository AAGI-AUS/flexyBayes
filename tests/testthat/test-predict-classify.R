# =============================================================================
# predict(classify = ) -- the marginal-means table.
#
# `predict(asreml_fit, classify = "Variety")` is the object a breeder takes
# to a meeting, and this is the flexyBayes equivalent. Four properties are
# asserted:
#
#   * classify = NULL leaves the existing predict() paths untouched. That
#     is the whole reason the argument exists rather than a new verb.
#   * the table has one row per level combination, in the frozen columns.
#   * the print names the estimator and says what the table is not. An
#     ASReml user reading a means table without an SED column beside it
#     should be told, not left to notice.
#   * emmeans is a suggested package, so its absence refuses by name rather
#     than through a bare requireNamespace() error.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


.pc_data <- function(seed = 20260818L, n = 60L) {
  set.seed(seed)
  d <- data.frame(
    variety = factor(rep(c("Axe", "Baudin", "Compass"), length.out = n)),
    env = factor(rep(c("E1", "E2"), each = n / 2L)),
    blk = factor(rep(letters[1:5], length.out = n))
  )
  b <- stats::rnorm(5L, sd = 0.6)[as.integer(d$blk)]
  d$yield <- 4 + c(0, 0.8, -0.5)[as.integer(d$variety)] +
    0.3 * (d$env == "E2") + b + stats::rnorm(n, sd = 0.4)
  d
}

.pc_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

.pc_cache <- new.env(parent = emptyenv())

.pc_inla_fit <- function() {
  if (is.null(.pc_cache$inla)) {
    .pc_cache$inla <- suppressMessages(fb(
      yield ~ variety + env,
      random = ~blk,
      data = .pc_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .pc_cache$inla
}

.pc_brms_fit <- function() {
  if (is.null(.pc_cache$brms)) {
    .pc_cache$brms <- suppressMessages(suppressWarnings(fb(
      yield ~ variety + env,
      random = ~blk,
      data = .pc_data(),
      backend = "brms",
      aggregate = FALSE,
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .pc_cache$brms
}

.PC_COLS <- c("estimate", "std.error", "conf.low", "conf.high")


# ---------------------------------------------------------------- #
# 1. classify = NULL is the old predict()                           #
# ---------------------------------------------------------------- #

test_that("predict() without classify is unchanged on INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  pred <- stats::predict(fit)
  expect_true(is.numeric(pred))
  expect_null(dim(pred))
  expect_identical(length(pred), nrow(.pc_data()))
  # The documented identity: fixed effects only, over the fit's own data.
  trms <- stats::delete.response(stats::terms(stats::formula(fit)))
  bhat <- stats::coef(fit)
  x_mat <- flexyBayes:::.fb_fixef_model_matrix(
    trms, fit$data, names(bhat), fit$data
  )
  expect_equal(pred, as.numeric(x_mat %*% bhat))

  # And the se.fit shape.
  out <- stats::predict(fit, se.fit = TRUE)
  expect_identical(names(out), c("fit", "se.fit"))
})

test_that("predict() on newdata is unchanged when classify is NULL", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()
  nd <- .pc_data()[1:5, , drop = FALSE]
  expect_identical(length(stats::predict(fit, newdata = nd)), 5L)
})


# ---------------------------------------------------------------- #
# 2. The table                                                      #
# ---------------------------------------------------------------- #

test_that("classify = 'variety' returns one row per level on INLA", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  tab <- stats::predict(fit, classify = "variety")
  expect_s3_class(tab, "fb_predict_classify")
  expect_s3_class(tab, "data.frame")
  expect_identical(nrow(tab), 3L)
  expect_identical(names(tab), c("variety", .PC_COLS))
  expect_identical(
    sort(as.character(tab$variety)),
    sort(levels(.pc_data()$variety))
  )
  expect_true(all(is.finite(tab$estimate)))
  expect_true(all(tab$std.error > 0))
  expect_true(all(tab$conf.low < tab$conf.high))
  # The means recover the ordering the data were generated with.
  expect_identical(
    as.character(tab$variety[which.max(tab$estimate)]),
    "Baudin"
  )
})

test_that("classify = 'variety' returns one row per level on brms", {
  skip_if_not_installed("brms")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_brms_fit()

  tab <- stats::predict(fit, classify = "variety")
  expect_s3_class(tab, "fb_predict_classify")
  expect_identical(nrow(tab), 3L)
  expect_identical(names(tab), c("variety", .PC_COLS))
  expect_true(all(tab$std.error > 0))
})

test_that("a one-sided formula and the ASReml spelling agree", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  a <- stats::predict(fit, classify = "variety")
  b <- stats::predict(fit, classify = ~variety)
  expect_identical(names(a), names(b))
  expect_identical(as.character(a$variety), as.character(b$variety))
  # The point estimate is X %*% bhat and is the same number both times.
  #
  # The standard error is not, and the band below is drawn from what its
  # noise actually measures rather than from a round number. On this
  # engine the marginal-mean standard error is sqrt(l' V l), where V is
  # vcov.flexybayes_inla()'s Monte-Carlo posterior covariance -- 2000
  # draws from INLA::inla.posterior.sample(), redrawn on every call, per
  # that method's documented default. Two calls therefore compare one
  # Monte-Carlo estimate against another and neither side is a fixed
  # reference.
  #
  # The arithmetic. A variance estimated from n independent draws carries
  # a relative standard deviation of sqrt(2 / (n - 1)) = 3.2% at n = 2000,
  # and its square root carries half of that, 1.6%. Measured on this
  # fixture over twenty independent calls the realised figure is 2.0% per
  # level -- a little above the independent-draw value because
  # inla.posterior.sample() draws across a weighted set of hyperparameter
  # configurations rather than from one fixed Gaussian. Two independent
  # calls differ with a relative standard deviation of sqrt(2) x 2.0% =
  # 2.8%, so 0.15 is a 5.3-sigma bound and 0.05 was a 1.8-sigma one --
  # which is why the 0.05 band failed roughly one comparison in ten
  # (measured: 380 ordered pairs from those twenty calls, 11.6% of them
  # over 0.05, the largest observed 0.085).
  #
  # A 5.3-sigma band still has teeth for what this test is for. A
  # spelling that resolved to a different linear combination -- a cell
  # mean where a marginal mean was asked for, say -- moves the standard
  # error by far more than 15%, while redrawing the same combination's
  # covariance does not.
  expect_equal(a$estimate, b$estimate, tolerance = 1e-10)
  expect_equal(a$std.error, b$std.error, tolerance = 0.15)
})

test_that("a two-factor classify crosses the levels", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  tab <- stats::predict(fit, classify = "variety:env")
  expect_identical(nrow(tab), 6L)
  expect_identical(names(tab), c("variety", "env", .PC_COLS))
  expect_identical(
    paste(tab$variety, tab$env),
    paste(
      rep(levels(.pc_data()$variety), times = 2L),
      rep(levels(.pc_data()$env), each = 3L)
    )
  )
  starred <- stats::predict(fit, classify = ~ variety * env)
  expect_identical(names(starred), names(tab))
  expect_equal(starred$estimate, tab$estimate, tolerance = 1e-10)
})

test_that("level widens the interval and is honoured on both engines", {
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fits <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    fits$inla <- .pc_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fits$brms <- .pc_brms_fit()
  }
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    wide <- stats::predict(fits[[nm]], classify = "variety", level = 0.99)
    narrow <- stats::predict(fits[[nm]], classify = "variety", level = 0.80)
    expect_true(
      all(wide$conf.high - wide$conf.low >
        narrow$conf.high - narrow$conf.low),
      label = nm
    )
    # The point estimate does not move with the level.
    expect_equal(wide$estimate, narrow$estimate, tolerance = 1e-8)
  }
})

test_that("an out-of-range level is refused before any fitting work", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()
  expect_error(
    stats::predict(fit, classify = "variety", level = 95),
    "strictly between 0 and 1"
  )
})


# ---------------------------------------------------------------- #
# 3. The print says what the table is, and is not                   #
# ---------------------------------------------------------------- #

test_that("the banner names the estimator and denies the SED reading", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  out <- utils::capture.output(print(stats::predict(fit, classify = "variety")))
  txt <- paste(out, collapse = "\n")
  expect_match(
    txt,
    "Estimate = posterior mean of the marginal mean. Not an SED table.",
    fixed = TRUE
  )
  # C-p2: the interval's provenance on this engine, said plainly.
  expect_match(
    txt,
    "INLA: intervals from the Gaussian approximation of the fixed effects.",
    fixed = TRUE
  )
  expect_match(txt, "95%", fixed = TRUE)
})

test_that("the brms print carries no INLA approximation line", {
  skip_if_not_installed("brms")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_brms_fit()

  txt <- paste(
    utils::capture.output(print(stats::predict(fit, classify = "variety"))),
    collapse = "\n"
  )
  expect_match(txt, "Not an SED table.", fixed = TRUE)
  expect_false(grepl("Gaussian approximation", txt, fixed = TRUE))
})

test_that("there is no sed argument and no pairwise block", {
  # D-5: the means table only. A pairwise standard error of a difference
  # is a separate object and is not smuggled in under another name.
  # (predict.flexybayes(), the bare parent method, was deleted entirely
  # at 0.9.3 -- both active engines fully override it, so the two
  # engine-specific methods are the only formals left to guard.)
  expect_false("sed" %in% names(formals(predict.flexybayes_inla)))
  expect_false("pairwise" %in% names(formals(predict.flexybayes_inla)))
  expect_false("sed" %in% names(formals(predict.flexybayes_brms)))
  expect_false("pairwise" %in% names(formals(predict.flexybayes_brms)))
})


# ---------------------------------------------------------------- #
# 4. Argument handling and the suggested-package guard              #
# ---------------------------------------------------------------- #

test_that("newdata on the classify path is ignored, not silently", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  expect_warning(
    tab <- stats::predict(
      fit,
      newdata = .pc_data()[1:3, , drop = FALSE],
      classify = "variety"
    ),
    "ignored on the classify path"
  )
  expect_identical(nrow(tab), 3L)
})

test_that("a malformed classify is refused by shape", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  expect_error(
    stats::predict(fit, classify = yield ~ variety),
    "one-sided form"
  )
  expect_error(stats::predict(fit, classify = 3L), "must be a character")
  expect_error(stats::predict(fit, classify = "  "), "names no variables")
})

test_that("emmeans absent refuses by name rather than by namespace", {
  # C-p1. Exercised without uninstalling emmeans, by mocking the named
  # availability predicate the guard consults.
  local_mocked_bindings(
    .fb_emmeans_available = function() FALSE,
    .package = "flexyBayes"
  )
  err <- tryCatch(
    flexyBayes:::.fb_predict_classify(
      structure(list(), class = "flexybayes"),
      "variety"
    ),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_classify_requires_emmeans")
  msg <- conditionMessage(err)
  expect_match(msg, "emmeans", fixed = TRUE)
  expect_match(msg, "install.packages(\"emmeans\")", fixed = TRUE)
  expect_match(msg, "coef(fit, what = \"random\")", fixed = TRUE)
})

test_that("the guard fires ahead of the classify parse", {
  # A refusal that only fires after the argument parses would let a
  # malformed classify mask the missing dependency.
  local_mocked_bindings(
    .fb_emmeans_available = function() FALSE,
    .package = "flexyBayes"
  )
  expect_error(
    flexyBayes:::.fb_predict_classify(
      structure(list(), class = "flexybayes"),
      3L
    ),
    class = "flexybayes_classify_requires_emmeans"
  )
})

test_that("the refusal code is registered in the taxonomy", {
  entries <- ls(flexyBayes:::.refusal_registry)
  expect_true("classify_requires_emmeans" %in% entries)
})


# ---------------------------------------------------------------- #
# 8. A random-effects grouping factor has no marginal mean here     #
# ---------------------------------------------------------------- #
#
# `predict(fit, classify = "Geno")` on a multi-environment trial is the
# first thing a breeder asks of a fitted model, and ASReml serves it.
# emmeans builds its reference grid from the population-level design, so a
# factor entering only as a grouping term is not in it and the call died
# inside emmeans with "No variable named Geno in the reference grid" -- a
# third-party message with no reason code, no mention of this package, and
# no route onward. The level effects are on the fit the whole time.

test_that("classify on a random-only factor refuses by name, both engines", {
  skip_on_cran()
  skip_on_ci()
  .pc_silence()

  if (requireNamespace("INLA", quietly = TRUE)) {
    expect_error(
      predict(.pc_inla_fit(), classify = "blk"),
      class = "flexybayes_refusal_classify_random_factor_not_supported"
    )
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    expect_error(
      predict(.pc_brms_fit(), classify = "blk"),
      class = "flexybayes_refusal_classify_random_factor_not_supported"
    )
  }
})

test_that("the refusal names the factor and the accessors that answer", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()

  cnd <- tryCatch(
    predict(.pc_inla_fit(), classify = "blk"),
    error = function(e) e
  )
  expect_identical(cnd$reason_code, "classify_random_factor_not_supported")
  expect_identical(cnd$factor, "blk")

  msg <- conditionMessage(cnd)
  expect_match(msg, "blk", fixed = TRUE)
  expect_match(msg, "random-effects grouping factor", fixed = TRUE)
  expect_match(msg, "reference grid", fixed = TRUE)
  expect_match(msg, "coef(fit, what = \"random\")", fixed = TRUE)
  expect_match(msg, "ranef(fit)", fixed = TRUE)
  expect_match(msg, "are planned", fixed = TRUE)
})

test_that("the named accessors do answer for that factor", {
  # The refusal is only defensible if what it points at works, so the
  # message and the accessor are asserted together.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()

  re <- stats::coef(.pc_inla_fit(), what = "random")
  expect_true("blk" %in% names(re))
  expect_gt(nrow(re$blk), 0L)
})

test_that("a fixed factor still classifies, and an unknown name is left alone", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pc_silence()
  fit <- .pc_inla_fit()

  out <- predict(fit, classify = "variety")
  expect_s3_class(out, "fb_predict_classify")
  expect_identical(nrow(out), 3L)

  # A name the model does not carry at all is not this refusal's business
  # -- emmeans's own message for that case is already about the right
  # thing, and claiming it as a random factor would be a lie.
  cnd <- tryCatch(predict(fit, classify = "nosuchvar"), error = function(e) e)
  expect_null(cnd$reason_code)
})

test_that("the random-classify refusal code is registered", {
  entries <- ls(flexyBayes:::.refusal_registry)
  expect_true("classify_random_factor_not_supported" %in% entries)
})
