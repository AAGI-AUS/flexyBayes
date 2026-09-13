# =============================================================================
# INLA containment (2026-09-13).
#
# Two refusals, each a recorded decision:
#
#   * estimated correlation matrices -- the correlated random slope (x | g)
#     and us(f):g -- are fitted on brms only
#     (correlated_covariance_requires_brms);
#   * a Gaussian fit on INLA refuses when the response standard deviation
#     exceeds .FB_INLA_RESPONSE_SD_MAX (inla_response_scale_unvalidated).
#
# A refusal protects the route it is written on and no other, so every route
# that can put a model on INLA is exercised here: the gate, an explicit
# backend = "inla" fit, the fb_inla() pin, fb_plan() under "auto", the
# aggregated route and the streamed route, the last two of which do not pass
# through the gate. The scale check reads a standard deviation recorded at
# ingest, so each ingest path is checked for recording it.
# =============================================================================

.ct_rules <- function(res) {
  if (!inherits(res, "lgm_refusal")) {
    return(character(0))
  }
  vapply(res$failures, function(f) f$rule_id, character(1L))
}

# A response with an exact standard deviation, around a level far from zero.
.ct_scaled <- function(sd_target, n = 120L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(
    g = factor(rep(sprintf("g%02d", seq_len(12L)), length.out = n)),
    x = stats::rnorm(n)
  )
  z <- stats::rnorm(n)
  d$y <- 1000 + sd_target * (z - mean(z)) / stats::sd(z)
  d
}

.ct_sleepstudy <- function(multiplier = 1) {
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  d <- e$sleepstudy
  d$Reaction <- multiplier * d$Reaction
  d
}


# --- estimated correlation matrices ------------------------------------ #

test_that("the gate refuses both correlated structures and nothing near them", {
  d <- .ct_scaled(1)
  d$env <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  refused <- list(
    fb_from_brms(y ~ x + (x | g), data = d),
    fb_from_asreml(fixed = y ~ x, random = ~ us(env):g, data = d)
  )
  for (fb in refused) {
    res <- lgm_gate(fb)
    expect_identical(res$primary_rule, "correlated_covariance_requires_brms")
    expect_identical(
      res$failures[[1L]]$reason_code,
      "correlated_covariance_requires_brms"
    )
    expect_match(res$failures[[1L]]$reason, "backend = \"brms\"", fixed = TRUE)
  }

  # The neighbours that carry no estimated correlation stay on INLA: the
  # uncorrelated slope and one variance per level.
  admitted <- list(
    fb_from_brms(y ~ x + (x || g), data = d),
    fb_from_asreml(fixed = y ~ x, random = ~ diag(env):g, data = d)
  )
  for (fb in admitted) {
    expect_false(
      "correlated_covariance_requires_brms" %in% .ct_rules(lgm_gate(fb))
    )
  }
})

test_that("an explicit INLA request for either structure raises the typed code", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ct_scaled(1)
  d$env <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  calls <- list(
    function() flexybayes(y ~ x + (x | g), data = d, backend = "inla",
                          return_code = TRUE, verbose = FALSE),
    function() flexybayes(y ~ x, random = ~ us(env):g, data = d,
                          backend = "inla", return_code = TRUE,
                          verbose = FALSE),
    function() fb_inla(y ~ x + (x | g), data = d, return_code = TRUE,
                       verbose = FALSE)
  )
  for (call in calls) {
    err <- tryCatch(suppressMessages(call()), error = function(e) e)
    expect_s3_class(
      err, "flexybayes_refusal_correlated_covariance_requires_brms"
    )
    expect_identical(err$reason_code, "correlated_covariance_requires_brms")
    expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
  }
})

test_that("auto routes both structures to brms before anything is fitted", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ct_scaled(1)
  d$env <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  p_slope <- fb_plan(y ~ x + (x | g), data = d, backend = "auto")
  expect_identical(p_slope$backend_chosen, "brms")
  expect_identical(p_slope$gate_primary_rule,
                   "correlated_covariance_requires_brms")

  p_us <- suppressMessages(flexybayes(
    y ~ x, random = ~ us(env):g, data = d, backend = "auto", plan = TRUE
  ))
  expect_identical(p_us$backend_chosen, "brms")
  expect_identical(p_us$gate_primary_rule,
                   "correlated_covariance_requires_brms")

  # An explicit INLA plan reports the same rule and does not claim a fit.
  p_inla <- fb_plan(y ~ x + (x | g), data = d, backend = "inla")
  expect_false(p_inla$will_fit)
  expect_identical(p_inla$gate_primary_rule,
                   "correlated_covariance_requires_brms")
})


# --- the response scale ------------------------------------------------ #

test_that("the scale bound sits at 150 and is read for Gaussian fits only", {
  expect_identical(.FB_INLA_RESPONSE_SD_MAX, 150)

  below <- fb_from_asreml(fixed = y ~ x, random = ~ g, data = .ct_scaled(149.9))
  above <- fb_from_asreml(fixed = y ~ x, random = ~ g, data = .ct_scaled(150.1))
  expect_false("inla_response_scale_unvalidated" %in% .ct_rules(lgm_gate(below)))
  res <- lgm_gate(above)
  expect_identical(res$primary_rule, "inla_response_scale_unvalidated")
  expect_identical(res$failures[[1L]]$reason_code,
                   "inla_response_scale_unvalidated")

  # Counts with a standard deviation in the thousands are not a Gaussian
  # fit, and the defect was measured on Gaussian fits.
  set.seed(3L)
  dp <- data.frame(g = factor(rep(letters[1:10], each = 12L)))
  dp$y <- stats::rpois(nrow(dp), lambda = rep(c(50, 5000), each = 60L))
  fp <- fb_from_asreml(fixed = y ~ 1, random = ~ g, data = dp,
                       family = "poisson")
  expect_gt(fp$data_summary$response_sd, 150)
  expect_false("inla_response_scale_unvalidated" %in% .ct_rules(lgm_gate(fp)))
})

test_that("the refusal says the bound is a proxy and names both ways forward", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  err <- tryCatch(
    suppressMessages(flexybayes(
      Reaction ~ Days + (1 | Subject), data = .ct_sleepstudy(10),
      backend = "inla", return_code = TRUE, verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_inla_response_scale_unvalidated")
  expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
  msg <- conditionMessage(err)
  expect_match(msg, "standard deviation of 563", fixed = TRUE)
  expect_match(msg, "proxy", fixed = TRUE)
  expect_match(msg, "Rescale the response", fixed = TRUE)
  expect_match(msg, "backend = \"brms\"", fixed = TRUE)

  # The engine pin reaches the same refusal.
  expect_error(
    suppressMessages(fb_inla(
      Reaction ~ Days + (1 | Subject), data = .ct_sleepstudy(10),
      return_code = TRUE, verbose = FALSE
    )),
    class = "flexybayes_refusal_inla_response_scale_unvalidated"
  )
})

test_that("sleepstudy times ten goes to brms under auto, and fits on INLA at x1", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("lme4")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  p10 <- fb_plan(Reaction ~ Days + (1 | Subject), data = .ct_sleepstudy(10),
                 backend = "auto")
  expect_identical(p10$backend_chosen, "brms")
  expect_identical(p10$gate_primary_rule, "inla_response_scale_unvalidated")

  p1 <- fb_plan(Reaction ~ Days + (1 | Subject), data = .ct_sleepstudy(1),
                backend = "auto")
  expect_identical(p1$backend_chosen, "inla")

  skip_on_cran()
  fit <- suppressWarnings(suppressMessages(flexybayes(
    Reaction ~ Days + (1 | Subject), data = .ct_sleepstudy(1),
    backend = "inla", verbose = FALSE
  )))
  expect_s3_class(fit, "flexybayes_inla")
})

test_that("the standard deviation is recorded on every ingest path", {
  skip_if_not_installed("lme4")
  d <- .ct_sleepstudy(10)
  s <- stats::sd(d$Reaction)

  expect_equal(
    fb_from_asreml(fixed = Reaction ~ Days, random = ~ Subject,
                   data = d)$data_summary$response_sd,
    s
  )
  expect_equal(
    fb_from_brms(Reaction ~ Days + (1 | Subject), data = d)$
      data_summary$response_sd,
    s
  )

  # A transformed response is measured as evaluated, not as the column it
  # names: the logarithm of a large response is on a small scale, and ten
  # times a moderate one is not.
  d1 <- .ct_sleepstudy(1)
  logged <- fb_from_asreml(fixed = log(Reaction) ~ Days, random = ~ Subject,
                           data = d1)
  expect_lt(logged$data_summary$response_sd, 1)
  tenfold <- fb_from_asreml(fixed = I(10 * Reaction) ~ Days,
                            random = ~ Subject, data = d1)
  expect_equal(tenfold$data_summary$response_sd, s)
  expect_identical(lgm_gate(tenfold)$primary_rule,
                   "inla_response_scale_unvalidated")

  # A response with gaps is still measured and still gated. The missing-
  # response policy adds or drops only rows whose response is missing, so
  # the standard deviation recorded at ingest is the one the gate reads.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d_na <- d
  d_na$Reaction[c(3L, 40L)] <- NA_real_
  p <- suppressMessages(flexybayes(
    Reaction ~ Days + (1 | Subject), data = d_na, backend = "auto",
    plan = TRUE
  ))
  expect_identical(p$gate_primary_rule, "inla_response_scale_unvalidated")
})

test_that("the aggregated route applies the bound to an explicit INLA request", {
  # An explicit backend = "inla" reaches the aggregated emit without the
  # gate, so the bound is applied there directly.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(20260913L)
  d <- expand.grid(
    rep = seq_len(20L),
    trt = factor(seq_len(5L)),
    blk = factor(seq_len(4L))
  )
  d$y <- 10 + as.numeric(d$trt) * 0.5 + stats::rnorm(nrow(d), 0, 0.3)
  d_big <- d
  d_big$y <- 1000 * d_big$y

  expect_error(
    suppressMessages(flexybayes(
      fixed = y ~ trt, random = ~ blk, data = d_big,
      backend = "inla", aggregate = TRUE, verbose = FALSE
    )),
    class = "flexybayes_refusal_inla_response_scale_unvalidated"
  )

  # The negative control: the same design at its own scale takes the
  # aggregated route.
  skip_on_cran()
  fit <- suppressMessages(flexybayes(
    fixed = y ~ trt, random = ~ blk, data = d,
    backend = "inla", aggregate = TRUE, verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
})

test_that("the streamed route reads the full data's scale, not the head's", {
  # flexybayes_stream() builds its IR from a head sample -- the first 1e5
  # rows -- and reaches the aggregated emit without the gate. The head below
  # is on a small scale and the rest of the source is not, so a check reading
  # the head's standard deviation would pass: about 1 in the head against
  # about 577 over the whole source.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(20260913L)
  n_head <- 100000L
  n <- 150000L
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
  )
  d$y <- c(stats::rnorm(n_head, 5, 1), stats::rnorm(n - n_head, 5, 1000))
  expect_lt(stats::sd(d$y[seq_len(n_head)]), 150)

  expect_error(
    suppressMessages(flexybayes_stream(
      y ~ f,
      random = ~ g,
      source = d,
      backend = "inla",
      chunk_rows = 50000L,
      verbose = FALSE
    )),
    class = "flexybayes_refusal_inla_response_scale_unvalidated"
  )

  # Without the fit, the aggregation itself is not refused: the bound is a
  # statement about INLA's Gaussian fit, not about the data.
  agg <- suppressMessages(flexybayes_stream(
    y ~ f, random = ~ g, source = d, backend = "inla",
    chunk_rows = 50000L, fit = FALSE, verbose = FALSE
  ))
  expect_s3_class(agg, "fb_aggregated")
})
