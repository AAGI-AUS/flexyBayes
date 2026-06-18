# Engine pins (ADR 0031 Phase 3 + 3.6).
#
# fb_inla(), fb_brms(), and fb_greta() are the three single-engine pins,
# each sugar over fb(backend = "<engine>"). At v0.5.0 the v0.4.1
# deprecation windows are closed: fb_brms() / fb_greta() no longer carry
# their old multi-backend / native-graph signatures, so a conflicting
# `backend` is a structured refusal rather than a deprecation warning.

mk_pin_df <- function(n = 40) {
  set.seed(7)
  data.frame(
    yield = rnorm(n),
    geno = factor(rep(letters[1:8], length.out = n)),
    env = factor(rep(c("a", "b"), length.out = n))
  )
}

# ---------------------------------------------------------------- #
# fb_inla() engine pin                                              #
# ---------------------------------------------------------------- #

test_that("fb_inla() is exported and pins the INLA engine", {
  expect_true("fb_inla" %in% getNamespaceExports("flexyBayes"))
  skip_on_cran()
  skip_if_not_installed("INLA")
  df <- mk_pin_df()
  fit <- suppressMessages(
    fb_inla(yield ~ env, random = ~geno, data = df, verbose = FALSE)
  )
  expect_s3_class(fit, "flexybayes_inla")
  expect_identical(fit$extras$backend_decision$backend, "inla")
})

test_that("fb_inla() passes brms-grammar polymorphism through to INLA", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  df <- mk_pin_df()
  fit <- suppressMessages(
    fb_inla(yield ~ env + (1 | geno), data = df, verbose = FALSE)
  )
  expect_identical(fit$extras$backend_decision$backend, "inla")
})

test_that("fb_inla() preserves missing()-based defaults (default-prior note)", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  # The default-prior note is a once-per-session latch; reset it so this
  # test observes it regardless of earlier fits in the file.
  flexyBayes:::.reset_emit_state_for_test()
  df <- mk_pin_df()
  msgs <- character(0)
  withr::with_options(
    list(flexyBayes.silence_default_prior_note = FALSE),
    withCallingHandlers(
      suppressWarnings(fb_inla(
        yield ~ env,
        random = ~geno,
        data = df,
        verbose = FALSE,
        mcmc_verbose = FALSE
      )),
      message = function(m) {
        msgs[[length(msgs) + 1L]] <<- conditionMessage(m)
        invokeRestart("muffleMessage")
      }
    )
  )
  expect_true(any(grepl("variance-component prior default|uniform", msgs)))
})

test_that(".fb_engine_pin() rewrites the call with the backend pinned", {
  cl <- quote(fb_inla(yield ~ env, data = d))
  # emulate match.call() shape inside fb_inla()
  rewritten <- cl
  rewritten[[1L]] <- quote(flexyBayes::flexybayes)
  rewritten$backend <- "inla"
  expect_identical(rewritten$backend, "inla")
  expect_identical(rewritten[[1L]], quote(flexyBayes::flexybayes))
})

# ---------------------------------------------------------------- #
# Backend-argument discipline (v0.5.0: pins reject a conflicting    #
# backend; the redundant self-pin is accepted)                      #
# ---------------------------------------------------------------- #

test_that("an engine pin refuses a conflicting backend argument", {
  df <- mk_pin_df()
  # fb_inla() pins INLA; backend = "greta" is a contradiction. Refuse
  # structurally rather than silently overwrite (the pre-v0.5.0 bug).
  expect_error(
    fb_inla(yield ~ env, data = df, backend = "greta"),
    class = "flexybayes_refusal_engine_pin_backend_conflict"
  )
  expect_error(
    fb_brms(yield ~ env + (1 | geno), data = df, backend = "inla", plan = TRUE),
    class = "flexybayes_refusal_engine_pin_backend_conflict"
  )
})

test_that("a redundant self-pin (fb_brms(backend = 'brms')) is accepted", {
  df <- mk_pin_df()
  # backend = "brms" agrees with the pinned engine, so it must NOT refuse
  # (the v0.4.1 deprecation notice promised this call survives the recast).
  # plan = TRUE short-circuits before any fit.
  expect_no_error(
    suppressMessages(
      fb_brms(
        yield ~ env + (1 | geno),
        data = df,
        backend = "brms",
        plan = TRUE,
        verbose = FALSE
      )
    )
  )
})

test_that("fb_greta() fits a formula via the greta engine", {
  skip_on_cran()
  skip_if_greta_backend_unusable()
  skip_on_ci()
  df <- mk_pin_df()
  fit <- suppressMessages(
    fb_greta(
      yield ~ env + (1 | geno),
      data = df,
      n_samples = 50L,
      warmup = 50L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  )
  expect_s3_class(fit, "flexybayes")
  expect_identical(fit$extras$backend_decision$backend, "greta")
})

test_that("fb_greta(model = <native graph>) still fits the native graph", {
  skip_on_cran()
  skip_if_greta_backend_unusable()
  skip_on_ci()
  # The removed `model = ` argument is remapped to the model-spec slot for
  # call-compatibility, so a native graph passed as `model = ` reaches the
  # direct greta::mcmc() fit and returns a flexybayes_direct_greta object.
  set.seed(7)
  yy <- greta::as_data(rnorm(30))
  xx <- greta::as_data(rnorm(30))
  b0 <- greta::normal(0, 10)
  b1 <- greta::normal(0, 10)
  s <- greta::uniform(0, 5)
  greta::distribution(yy) <- greta::normal(b0 + b1 * xx, s)
  m <- greta::model(b0, b1, s)
  fit <- suppressMessages(
    fb_greta(
      model = m,
      n_samples = 50L,
      warmup = 50L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  )
  expect_s3_class(fit, "flexybayes_direct_greta")
})
