# =============================================================================
# Grammar-surface refusals: the bar grammar, the smoother spellings, and the
# INLA duplicate key.
#
# The register's shared complaint across FS-6, FS-17 and FS-18 is that a
# model the package cannot fit as written raised a bare `simpleError` whose
# text sent the user in a circle -- listing bar spellings none of which
# expresses the model, or naming a function R could not find, or quoting the
# INLA subprocess's generic exit message. Where the same model has a spelling
# that DOES fit, the refusal has to name it (`DESCRIPTION`: a typed refusal
# naming the nearest implemented alternative).
#
# The grid below is the field sweep's S5 grammar-surface section plus the two
# smoother cells from S4, each row asserting the condition class and the
# sibling surface the message must name verbatim.
# =============================================================================

.gs_data <- function(seed = 20260819L) {
  set.seed(seed)
  d <- data.frame(
    g = factor(paste0("g", rep(seq_len(10L), times = 6L))),
    f = factor(rep(c("e1", "e2", "e3"), each = 20L)),
    h = factor(rep(c("b1", "b2", "b3", "b4"), times = 15L)),
    trt = factor(rep(c("a", "b"), times = 30L)),
    x = stats::rnorm(60L),
    z = stats::rnorm(60L)
  )
  d$y <- 0.6 + 0.5 * d$x + stats::rnorm(60L, 0, 0.8)
  d
}

.gs_refuse <- function(..., backend = "brms") {
  options(flexyBayes.silence_default_prior_note = TRUE)
  args <- list(..., backend = backend, verbose = FALSE)
  if (identical(backend, "brms")) {
    args <- c(
      args,
      list(chains = 1L, n_samples = 200L, warmup = 100L, seed = 1L)
    )
  }
  tryCatch(
    suppressWarnings(suppressMessages(do.call(flexybayes, args))),
    error = function(e) e
  )
}


# --- bar-grammar factor slopes (field-sweep FS-6 / finding U1) ---------- #

test_that("a correlated factor slope names `us(f):g` verbatim, on both backends", {
  d <- .gs_data()
  for (backend in c("brms", "inla")) {
    for (form in list(y ~ trt + (trt | g), y ~ trt + (0 + trt | g))) {
      err <- .gs_refuse(form, data = d, family = "gaussian", backend = backend)
      expect_s3_class(
        err,
        "flexybayes_refusal_brms_factor_random_slope_unsupported"
      )
      msg <- conditionMessage(err)
      expect_match(msg, "random = ~ us(trt):g", fixed = TRUE)
      expect_match(msg, "backend = \"brms\"", fixed = TRUE)
      expect_identical(err$asreml_surface, "random = ~ us(trt):g")
      expect_true(err$correlated)
    }
  }
})

test_that("an uncorrelated factor slope names `diag(f):g` verbatim", {
  d <- .gs_data()
  err <- .gs_refuse(y ~ trt + (trt || g), data = d, family = "gaussian")
  expect_s3_class(
    err,
    "flexybayes_refusal_brms_factor_random_slope_unsupported"
  )
  expect_match(conditionMessage(err), "random = ~ diag(trt):g", fixed = TRUE)
  expect_false(err$correlated)
  # And it points at the correlated sibling as the other option.
  expect_match(conditionMessage(err), "random = ~ us(trt):g", fixed = TRUE)
})

test_that("the named ASReml surfaces do fit -- the refusal is not a dead end", {
  # Naming an alternative that does not work is worse than naming none, so
  # both surfaces the messages above name are exercised here.
  skip_if_not_installed("brms")
  options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .gs_data()
  for (rnd in list(~ us(trt):g, ~ diag(trt):g)) {
    fit <- suppressWarnings(suppressMessages(flexybayes(
      y ~ trt,
      random = rnd,
      data = d,
      family = "gaussian",
      backend = "brms",
      chains = 1L,
      n_samples = 200L,
      warmup = 100L,
      seed = 1L,
      verbose = FALSE
    )))
    expect_s3_class(fit, "flexybayes")
  }
})

test_that("a numeric multi-slope refuses typed and claims no ASReml sibling", {
  # The honest half of the contract: this form has no surface that fits, so
  # the message must not invent one.
  d <- .gs_data()
  err <- .gs_refuse(y ~ x + z + (x + z || g), data = d, family = "gaussian")
  expect_s3_class(err, "flexybayes_refusal_brms_random_effect_form_unsupported")
  expect_false(grepl("us(", conditionMessage(err), fixed = TRUE))
  expect_false(grepl("diag(", conditionMessage(err), fixed = TRUE))
})

test_that("the supported bar forms are untouched by the factor-slope branch", {
  skip_if_not_installed("brms")
  options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .gs_data()
  for (form in list(y ~ x + (1 | g), y ~ x + (x || g), y ~ x + (1 + x || g))) {
    fit <- suppressWarnings(suppressMessages(flexybayes(
      form,
      data = d,
      family = "gaussian",
      backend = "brms",
      chains = 1L,
      n_samples = 200L,
      warmup = 100L,
      seed = 1L,
      verbose = FALSE
    )))
    expect_s3_class(fit, "flexybayes")
  }
})

test_that("a correlated NUMERIC slope keeps its own refusal", {
  d <- .gs_data()
  err <- .gs_refuse(y ~ x + (x | g), data = d, family = "gaussian")
  expect_s3_class(err, "flexybayes_correlated_slope_unsupported")
})


# --- smoother in the fixed part (field-sweep FS-18) ---------------------- #

test_that("s(x) in the fixed part refuses typed on both backends", {
  d <- .gs_data()
  for (backend in c("brms", "inla")) {
    err <- .gs_refuse(
      y ~ s(x),
      data = d,
      family = "gaussian",
      backend = backend
    )
    expect_s3_class(err, "flexybayes_refusal_fixed_smoother_not_supported")
    msg <- conditionMessage(err)
    expect_match(msg, "s(x)", fixed = TRUE)
    expect_match(msg, "random = ~ spl(x)", fixed = TRUE)
    expect_match(msg, "backend = \"inla\"", fixed = TRUE)
  }
})

test_that("every mgcv smoother spelling is refused, not just s()", {
  d <- .gs_data()
  for (head in flexyBayes:::.fb_fixed_smoother_heads()) {
    form <- stats::as.formula(sprintf("y ~ %s(x)", head))
    err <- .gs_refuse(form, data = d, family = "gaussian", backend = "inla")
    expect_s3_class(
      err,
      "flexybayes_refusal_fixed_smoother_not_supported"
    )
    expect_match(conditionMessage(err), sprintf("%s(x)", head), fixed = TRUE)
  }
})

test_that("an ordinary transformation in the fixed part still parses", {
  # The refusal is a closed list of smoother heads, not a ban on calls.
  d <- .gs_data()
  d$xpos <- abs(d$x) + 1
  for (form in list(y ~ log(xpos), y ~ I(x^2), y ~ poly(x, 2))) {
    fb <- fb_from_asreml(form, data = d)
    expect_s3_class(fb, "fb_terms")
  }
})


# --- the INLA duplicate key (field-sweep FS-17) ------------------------- #

test_that("a spline on a variable that is also fixed refuses typed on INLA", {
  # The duplicate-key guard lives in the INLA emit, so the emit has to be
  # reachable: without INLA installed the call refuses earlier, with the
  # installation refusal, and the assertions below would be graded against
  # the wrong condition. The sibling test beneath carries the same guard.
  skip_if_not_installed("INLA")
  d <- .gs_data()
  err <- .gs_refuse(
    y ~ x,
    random = ~ spl(x),
    data = d,
    family = "gaussian",
    backend = "inla"
  )
  expect_s3_class(err, "flexybayes_refusal_inla_variable_used_twice")
  msg <- conditionMessage(err)
  expect_match(msg, "`x`", fixed = TRUE)
  # Both of INLA's own remedies, plus the one specific to an rw2 smooth.
  expect_match(msg, "random = ~ spl(x)", fixed = TRUE)
  expect_match(msg, "x2 <- x", fixed = TRUE)
  expect_identical(err$variables, "x")
})

test_that("the same spline without the duplicate fixed term still fits", {
  skip_if_not_installed("INLA")
  options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .gs_data()
  fit <- suppressMessages(flexybayes(
    y ~ 1,
    random = ~ spl(x),
    data = d,
    family = "gaussian",
    backend = "inla",
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
})

test_that("the earlier fixed-and-random guard still owns the factor case", {
  # A grouping factor on both sides would also be a duplicate INLA key, but
  # it is caught upstream and backend-independently by
  # `term_in_fixed_and_random`, which explains the identifiability problem
  # rather than the engine's key restriction. The duplicate-key guard is
  # the residual case that guard does not reach: a smooth of a variable
  # that is also a fixed covariate.
  d <- .gs_data()
  err <- .gs_refuse(
    y ~ f,
    random = ~f,
    data = d,
    family = "gaussian",
    backend = "inla"
  )
  expect_s3_class(err, "flexybayes_refusal_term_in_fixed_and_random")
})

test_that("an ordinary random intercept is not caught by the duplicate-key guard", {
  skip_if_not_installed("INLA")
  options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .gs_data()
  fit <- suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    backend = "inla",
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
})


# --- the refusals are registered ---------------------------------------- #

test_that("every reason code raised here is in the registry", {
  for (code in c(
    "brms_factor_random_slope_unsupported",
    "brms_random_effect_form_unsupported",
    "brms_ingest_feature_unsupported",
    "fixed_smoother_not_supported",
    "inla_variable_used_twice"
  )) {
    reg <- fb_refusals(reason_code = code)
    expect_equal(nrow(reg), 1L, label = code)
    expect_equal(reg$since_version, "0.9.2", label = code)
  }
})
