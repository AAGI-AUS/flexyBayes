# Tests for triangulate()'s comparability, diagnostic and matched-prior
# gates, and for the per-parameter and overall verdicts built on them.
#
# The synthetic fits below carry a REAL fingerprint -- built by the same
# .fb_model_fingerprint() the two emits call -- over fabricated draws and
# fabricated convergence diagnostics. Fabricating the draws is what makes
# the concordance and discordance cases deterministic; fabricating the
# fingerprint would test nothing, so it is built from a real IR.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

.fixture_data <- function(n = 60L, seed = 41L) {
  set.seed(seed)
  data.frame(
    y = stats::rnorm(n),
    x = stats::rnorm(n),
    g = factor(rep(seq_len(6L), length.out = n))
  )
}

# Build the IR for `y ~ x + (1 | g)` and attach the package's own default
# prior, exactly as flexybayes() does before dispatch.
.fixture_ir <- function(data, formula = y ~ x + (1 | g)) {
  ir <- fb_from_brms(formula, data = data)
  ir$priors <- flexyBayes:::.default_uniform_prior(
    data = data,
    response = ir$response,
    family = "gaussian",
    link = NULL,
    random_groups = flexyBayes:::.fb_default_prior_targets(ir)$shared
  )
  ir
}

# A fit-shaped object: our own class first so the draws method below is
# the one that dispatches, `flexybayes_brms` behind it so the diagnostics
# gate reads the fabricated convergence slot the way it reads a real one.
mk_fp_fit <- function(draws, fingerprint, rhat = 1.001, ess = 2500,
                      n_divergent = 0L) {
  nms <- names(draws)
  psrf <- matrix(
    c(rep(rhat, length(nms)), rep(NA_real_, length(nms))),
    nrow = length(nms),
    ncol = 2L,
    dimnames = list(nms, c("Point est.", "Upper C.I."))
  )
  structure(
    list(
      draws = draws,
      extras = list(
        fingerprint = fingerprint,
        convergence = list(
          gelman = list(psrf = psrf),
          n_eff = stats::setNames(rep(ess, length(nms)), nms),
          n_divergent = n_divergent
        )
      )
    ),
    class = c("fp_fit", "flexybayes_brms", "flexybayes", "list")
  )
}

registerS3method(
  "fb_as_draws_simple",
  "fp_fit",
  function(fit, ...) fit$draws
)

# ---------------------------------------------------------------- #
# Comparability gate                                                #
# ---------------------------------------------------------------- #

test_that("two fits of the same model on the same data are comparable", {
  d <- .fixture_data()
  fp <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d), d)
  set.seed(1L)
  draws <- list("(Intercept)" = stats::rnorm(2000, 0.5, 0.2))
  tri <- triangulate(mk_fp_fit(draws, fp), mk_fp_fit(draws, fp))
  expect_true(tri$comparability$verified)
  expect_identical(tri$status, "concordant")
})

test_that("a different dataset refuses by name and names the data", {
  d1 <- .fixture_data()
  d2 <- d1
  d2$y[[3L]] <- d2$y[[3L]] + 1
  fp1 <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d1), d1)
  fp2 <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d2), d2)
  set.seed(2L)
  draws <- list("(Intercept)" = stats::rnorm(500))

  err <- tryCatch(
    triangulate(mk_fp_fit(draws, fp1), mk_fp_fit(draws, fp2)),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_triangulate_incomparable_fits")
  expect_match(conditionMessage(err), "the data contents", fixed = TRUE)
  expect_identical(err$differing_element, "data_digest")
})

test_that("a different random structure refuses and names it", {
  d <- .fixture_data()
  fp1 <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d), d)
  fp2 <- flexyBayes:::.fb_model_fingerprint(
    .fixture_ir(d, formula = y ~ x),
    d
  )
  set.seed(3L)
  draws <- list("(Intercept)" = stats::rnorm(500))

  err <- tryCatch(
    triangulate(mk_fp_fit(draws, fp1), mk_fp_fit(draws, fp2)),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_triangulate_incomparable_fits")
  expect_identical(err$differing_element, "random")
  expect_match(conditionMessage(err), "the random-effect structure",
               fixed = TRUE)
})

test_that("a different prior on a shared parameter refuses and names it", {
  d <- .fixture_data()
  ir1 <- .fixture_ir(d)
  ir2 <- .fixture_ir(d)
  # Same model, same data, a wider uniform on every SD.
  ir2$priors <- flexyBayes:::.default_uniform_prior(
    data = transform(d, y = d$y * 2),
    response = "y",
    family = "gaussian",
    link = NULL,
    random_groups = "g"
  )
  fp1 <- flexyBayes:::.fb_model_fingerprint(ir1, d)
  fp2 <- flexyBayes:::.fb_model_fingerprint(ir2, d)
  set.seed(4L)
  draws <- list(sd_g = abs(stats::rnorm(500, 1, 0.1)))

  err <- tryCatch(
    triangulate(mk_fp_fit(draws, fp1), mk_fp_fit(draws, fp2)),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_triangulate_incomparable_fits")
  expect_match(conditionMessage(err), "the prior on `sd_g`", fixed = TRUE)
  expect_false(grepl("allow_", conditionMessage(err), fixed = TRUE))
})

test_that("a fit without a fingerprint is inconclusive, not refused", {
  set.seed(5L)
  d <- list(p = stats::rnorm(500))
  tri <- triangulate(mk_fp_fit(d, NULL), mk_fp_fit(d, NULL))
  expect_false(tri$comparability$verified)
  expect_identical(tri$status, "inconclusive")
  expect_true(any(grepl("comparability unverified", tri$status_reasons)))
})

# ---------------------------------------------------------------- #
# Diagnostics gate                                                  #
# ---------------------------------------------------------------- #

test_that("a fit that failed its own diagnostics yields no verdict", {
  d <- .fixture_data()
  fp <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d), d)
  set.seed(6L)
  draws <- list("(Intercept)" = stats::rnorm(2000, 0.5, 0.2))

  tri <- triangulate(
    mk_fp_fit(draws, fp, rhat = 1.20),
    mk_fp_fit(draws, fp)
  )
  expect_identical(tri$status, "inconclusive")
  expect_true(any(grepl("failed its own diagnostics", tri$status_reasons)))
  expect_true(any(grepl("max R-hat", tri$status_reasons)))
  expect_true(all(tri$metrics$verdict == "not_compared"))
  expect_identical(tri$n_compared, 0L)
})

test_that("divergent transitions and a low ESS both fail the gate", {
  d <- .fixture_data()
  fp <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d), d)
  set.seed(7L)
  draws <- list("(Intercept)" = stats::rnorm(2000, 0.5, 0.2))

  tri_div <- triangulate(
    mk_fp_fit(draws, fp, n_divergent = 12L),
    mk_fp_fit(draws, fp)
  )
  expect_identical(tri_div$status, "inconclusive")
  expect_true(any(grepl("12 divergent transitions", tri_div$status_reasons)))

  tri_ess <- triangulate(
    mk_fp_fit(draws, fp, ess = 50),
    mk_fp_fit(draws, fp)
  )
  expect_identical(tri_ess$status, "inconclusive")
  expect_true(any(grepl("min bulk ESS", tri_ess$status_reasons)))
})

# ---------------------------------------------------------------- #
# Matched-prior gate                                                #
# ---------------------------------------------------------------- #

test_that("a variance component with no shared prior is not compared", {
  d <- data.frame(
    y = stats::rnorm(60),
    env = factor(rep(seq_len(3L), 20L)),
    gen = factor(rep(seq_len(10L), 6L))
  )
  ir <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ us(env):gen,
    data = d
  )
  ir$priors <- flexyBayes:::.default_uniform_prior(
    data = d,
    response = "y",
    family = "gaussian",
    link = NULL,
    random_groups = flexyBayes:::.fb_default_prior_targets(ir)$shared
  )
  fp <- flexyBayes:::.fb_model_fingerprint(ir, d)
  expect_true("cor_gen" %in% names(fp$engine_default_params))
  expect_true("sd_gen" %in% names(fp$priors))

  set.seed(8L)
  draws <- list(
    sd_gen = abs(stats::rnorm(2000, 1, 0.1)),
    cor_gen = stats::rnorm(2000, 0.2, 0.1)
  )
  tri <- triangulate(mk_fp_fit(draws, fp), mk_fp_fit(draws, fp))

  m <- tri$metrics
  expect_identical(m$verdict[m$param == "sd_gen"], "concordant")
  expect_identical(m$verdict[m$param == "cor_gen"], "not_compared")
  expect_match(m$reason[m$param == "cor_gen"], "LKJ")
  expect_identical(tri$n_compared, 1L)
  expect_identical(tri$status, "concordant")

  out <- capture.output(print(tri))
  expect_true(any(grepl("Not compared", out)))
  expect_true(any(grepl("status:", out)))
})

test_that("a fixed-effect coefficient needs no prior record to compare", {
  d <- .fixture_data()
  fp <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d), d)
  set.seed(9L)
  draws <- list(x = stats::rnorm(2000, 1, 0.3))
  tri <- triangulate(mk_fp_fit(draws, fp), mk_fp_fit(draws, fp))
  expect_identical(tri$metrics$verdict, "concordant")
})

# ---------------------------------------------------------------- #
# Verdicts                                                          #
# ---------------------------------------------------------------- #

test_that("a shifted posterior is discordant and says by how much", {
  d <- .fixture_data()
  fp <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d), d)
  set.seed(10L)
  a <- list(sigma = stats::rnorm(4000, 1.0, 0.10))
  b <- list(sigma = stats::rnorm(4000, 1.5, 0.10))

  tri <- triangulate(mk_fp_fit(a, fp), mk_fp_fit(b, fp))
  expect_identical(tri$status, "discordant")
  expect_identical(tri$metrics$verdict, "discordant")
  expect_gt(tri$metrics$mean_shift_sd, 1)
  expect_match(tri$status_reasons, "sigma")
})

test_that("an inflated posterior SD alone is enough to be discordant", {
  d <- .fixture_data()
  fp <- flexyBayes:::.fb_model_fingerprint(.fixture_ir(d), d)
  set.seed(11L)
  a <- list(sigma = stats::rnorm(8000, 1.0, 0.10))
  b <- list(sigma = stats::rnorm(8000, 1.0, 0.30))

  tri <- triangulate(mk_fp_fit(a, fp), mk_fp_fit(b, fp))
  expect_identical(tri$status, "discordant")
  expect_match(tri$metrics$reason, "SD ratio")
})

# ---------------------------------------------------------------- #
# Live: the collapse case the triangulation vignette teaches         #
# ---------------------------------------------------------------- #

test_that("the sleepstudy collapse is not reported as agreement", {
  skip_on_cran()
  skip_if_no_inla()
  skip_if_not_installed("brms")
  skip_if_not_installed("lme4")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE,
    # The INLA fit below lands on the boundary for sd_Subject on this
    # dataset (see the comment further down), which the 0.10.0 detector
    # correctly warns about. That is the condition under test here, not
    # an unexpected event, and it is intermittent -- expect_warning()
    # would fail on the runs where INLA finds the other mode.
    flexyBayes.silence_boundary_collapse_warning = TRUE
  )
  data(sleepstudy, package = "lme4", envir = environment())

  fit_b <- suppressMessages(fb_brms(
    Reaction ~ Days + (1 | Subject),
    data = sleepstudy,
    n_samples = 1000,
    warmup = 1000,
    chains = 4,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_i <- suppressMessages(fb_inla(
    Reaction ~ Days + (1 | Subject),
    data = sleepstudy,
    verbose = FALSE
  ))

  # Both fits carry a fingerprint and they agree: same model, same data,
  # the same injected default prior. The comparison is therefore about
  # the engines, which is the point of the vignette section it backs.
  tri <- triangulate(fit_b, fit_i, n_samples = 1000L)
  expect_true(tri$comparability$verified)
  expect_true("sd_Subject" %in% tri$common)

  # INLA's default hyperprior can leave the Subject variance at
  # numerically zero on this dataset. Either the engines disagree about
  # it, or a fit failed its own gate; what must not happen is a
  # concordant verdict on a collapsed variance component.
  expect_true(tri$status %in% c("discordant", "inconclusive"))
  row <- tri$metrics[tri$metrics$param == "sd_Subject", , drop = FALSE]
  expect_identical(nrow(row), 1L)
  expect_false(identical(row$verdict, "concordant"))
})
