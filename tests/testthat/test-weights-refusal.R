# =============================================================================
# Observation weights: refused by name where they are not lowered.
#
# C6 lowers `weights` for family = "gaussian" with the identity link only,
# in the ASReml / lme4 / glm(weights =) precision sense, Var(y_i) =
# sigma^2 / w_i (test-weights-lowering.R grounds the fitting side against
# lme4::lmer(weights =) and a metamorphic check on both engines). Every
# other family, and a non-identity link on Gaussian, still refuses --
# likelihood-power, frequency, trials and exposure are different models
# and cannot share this argument silently. This file covers the refusal
# side: wrong family/link (weights_requires_gaussian), the aggregated
# route (weights_not_aggregatable), and the pass-through / IR-recording
# mechanics that predate C6 and are unaffected by it.
#
# Before C6, EVERY family refused (weights_not_supported) -- this file's
# original form. That refusal is gone for the one family C6 lowers it
# for; the remaining tests below are what stayed true.
# =============================================================================

.wdat <- function(n = 8L, seed = 3L) {
  set.seed(seed)
  data.frame(y = stats::rnorm(n), x = stats::rnorm(n))
}

.wdat_count <- function(n = 40L, seed = 4L) {
  set.seed(seed)
  x <- stats::rnorm(n)
  data.frame(y = stats::rpois(n, exp(0.2 + 0.3 * x)), x = x)
}

test_that("a non-constant weight vector on Poisson is refused by name", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat_count()
  w <- as.numeric(seq_len(nrow(d)))

  for (be in c("auto", "inla", "brms")) {
    err <- tryCatch(
      suppressMessages(flexybayes(
        y ~ x, data = d, weights = w, family = "poisson", backend = be,
        return_code = TRUE
      )),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_weights_requires_gaussian")
    expect_match(conditionMessage(err), "poisson")
  }
})

test_that("a non-identity link on Gaussian is refused by name", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  w <- as.numeric(seq_len(nrow(d)))
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ x, data = d, weights = w, family = "gaussian", link = "log",
      backend = "brms", return_code = TRUE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_weights_requires_gaussian")
  expect_match(conditionMessage(err), "log")
})

test_that("Gaussian identity-link weights are NOT refused", {
  # The regression guard for C6 itself: this exact shape refused
  # (weights_not_supported) before C6. Confirms the refusal was narrowed,
  # not just renamed.
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  w <- as.numeric(seq_len(nrow(d)))
  engines <- c(
    if (requireNamespace("INLA", quietly = TRUE)) "inla",
    if (requireNamespace("brms", quietly = TRUE)) "brms"
  )
  skip_if(length(engines) == 0L, "neither engine is installed")
  for (be in engines) {
    expect_no_error(
      suppressMessages(flexybayes(
        y ~ x, data = d, weights = w, backend = be, return_code = TRUE
      ))
    )
  }
})

test_that("the refusal fires before the fit, not after", {
  # return_code = TRUE never reaches a sampler, so a refusal here proves
  # the gate sits ahead of emission rather than being a post-hoc check on
  # a fitted object.
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat_count()
  expect_error(
    suppressMessages(flexybayes(
      y ~ x, data = d, weights = as.numeric(seq_len(nrow(d))),
      family = "poisson", backend = "brms", return_code = TRUE
    )),
    class = "flexybayes_refusal_weights_requires_gaussian"
  )
})

test_that("aggregate = TRUE with weights is refused by name", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(5L)
  n <- 60L
  x <- stats::rnorm(n)
  g <- factor(rep(seq_len(6L), each = 10L))
  w <- stats::runif(n, 0.5, 3)
  y <- stats::rnorm(n, 2 + 0.3 * x, 1)
  d <- data.frame(y = y, x = x, g = g, w = w)
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ x, random = ~g, data = d, weights = w, aggregate = TRUE,
      backend = "inla", verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_weights_not_aggregatable")
})

test_that("aggregate = \"auto\" with weights falls through to per-row, unrefused", {
  # The same shape as above, but "auto" rather than an explicit TRUE:
  # the guard is documented to fall through rather than refuse.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(5L)
  n <- 60L
  x <- stats::rnorm(n)
  g <- factor(rep(seq_len(6L), each = 10L))
  w <- stats::runif(n, 0.5, 3)
  y <- stats::rnorm(n, 2 + 0.3 * x, 1)
  d <- data.frame(y = y, x = x, g = g, w = w)
  expect_no_error(
    suppressMessages(flexybayes(
      y ~ x, random = ~g, data = d, weights = w, aggregate = "auto",
      backend = "inla", verbose = FALSE
    ))
  )
})

test_that("constant weights are the unweighted model and pass through", {
  # rep(1, n) -- and any all-equal vector -- specifies the same model as
  # no weights at all, so refusing it (on any family) would be
  # gratuitous.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  expect_no_error(
    suppressMessages(flexybayes(
      y ~ x, data = d, weights = rep(1, nrow(d)),
      backend = "inla", return_code = TRUE
    ))
  )
  dc <- .wdat_count()
  expect_no_error(
    suppressMessages(flexybayes(
      y ~ x, data = dc, weights = rep(2.5, nrow(dc)), family = "poisson",
      backend = "brms", return_code = TRUE
    ))
  )
})

test_that("NULL weights are unaffected", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  expect_no_error(
    suppressMessages(flexybayes(
      y ~ x, data = d, backend = "inla", return_code = TRUE
    ))
  )
})

test_that("the ingest adapters still record the requested weights", {
  # The refusal (where it still applies) is a dispatch gate, not a
  # parser change, so the ingest adapters keep their coverage.
  d <- .wdat(5L)
  w <- c(1, 2, 3, 4, 5)
  ir <- flexyBayes:::fb_from_asreml(y ~ x, data = d, weights = w)
  types <- vapply(ir$addition_terms, function(t) t$type, character(1L))
  expect_true("weights" %in% types)
  expect_equal(ir$addition_terms[[match("weights", types)]]$values, w)
})

test_that("the guard reads the IR, not only the dispatch argument", {
  # An IR carrying weights on a non-Gaussian family must be refused even
  # when the dispatch-level `weights` argument is NULL -- otherwise a
  # caller that built the IR upstream would slip past the gate.
  d <- .wdat_count(5L)
  ir <- flexyBayes:::fb_from_asreml(
    y ~ x, data = d, weights = c(1, 2, 3, 4, 5),
    family = "poisson"
  )
  expect_error(
    flexyBayes:::.refuse_unsupported_weights(ir, weights = NULL),
    class = "flexybayes_refusal_weights_requires_gaussian"
  )
})

test_that("fb_plan() agrees with live dispatch on the family refusal", {
  # C3 (S6/FS-22) established plan == live for grammar routing; C6 adds a
  # refusal fb_plan() must also see, or "Will fit: yes" would misreport
  # a model live dispatch actually refuses.
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat_count()
  w <- as.numeric(seq_len(nrow(d)))
  expect_error(
    fb_plan(y ~ x, data = d, weights = w, family = "poisson"),
    class = "flexybayes_refusal_weights_requires_gaussian"
  )
  # And the unweighted / Gaussian-weighted cases must NOT be caught by
  # the same guard -- a regression here would silently over-refuse.
  expect_s3_class(fb_plan(y ~ x, data = d, family = "poisson"), "fb_plan")
  dg <- .wdat()
  wg <- as.numeric(seq_len(nrow(dg)))
  expect_s3_class(
    fb_plan(y ~ x, data = dg, weights = wg, family = "gaussian"), "fb_plan"
  )
})
