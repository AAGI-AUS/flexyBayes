# The default variance-component prior is reported with its scale. A
# Gaussian INLA fit runs on y / sd(y), and sigma ~ half-normal(a k) is
# sigma / k ~ half-normal(a). The tests convert one to the other by hand.

.ps_data <- function() {
  set.seed(20260924L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", 1:12), each = 8L)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(12L, sd = 0.4)[as.integer(d$site)] +
    stats::rnorm(nrow(d), sd = 0.3)
  d
}

# The first number after `pattern` in `text`.
.ps_number <- function(text, pattern) {
  hit <- regmatches(
    text, regexec(paste0(pattern, "([-+0-9.eE]+)"), text)
  )[[1L]]
  as.numeric(hit[[2L]])
}


test_that("the note, the banner and print() state the scale of the prior", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  local_clean_emit_state()
  withr::local_options(flexyBayes.silence_default_prior_note = FALSE)
  d <- .ps_data()
  notes <- NULL
  banner <- utils::capture.output(
    notes <- testthat::capture_messages(
      fit <- flexybayes(y ~ x, random = ~ site, data = d, backend = "inla",
                        verbose = TRUE, mcmc_verbose = FALSE)
    )
  )
  note <- grep("half_normal(scale = ", notes, fixed = TRUE, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "in the units of the response", fixed = TRUE)
  on_response <- .ps_number(note, "half_normal\\(scale = ")

  formula_line <- grep("formula:", banner, fixed = TRUE, value = TRUE)
  expect_length(formula_line, 1L)
  printed <- paste(utils::capture.output(print(fit)), collapse = "\n")
  for (txt in c(formula_line, printed)) {
    expect_match(
      txt, "on the standardised response = half-normal-SD(", fixed = TRUE
    )
    standardised <- .ps_number(txt, "<prior: half-normal-SD\\(")
    back <- .ps_number(
      txt, "on the standardised response = half-normal-SD\\("
    )
    # The same prior, converted by hand: scale on y = scale on y / sd(y)
    # times sd(y), each at the precision it is printed to.
    expect_equal(standardised * stats::sd(d$y), on_response, tolerance = 1e-3)
    expect_equal(back, on_response, tolerance = 1e-3)
  }
})


test_that("the prior INLA received is the stated prior on the response scale", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ps_data()
  fit <- flexybayes(y ~ x, random = ~ site, data = d, backend = "inla",
                    verbose = FALSE)
  stated <- prior_summary(fit)$default_scale
  sent_text <- paste(deparse(fit$extras$formula), collapse = " ")
  sent <- .ps_number(sent_text, "expression: s=")
  k <- stats::sd(d$y)
  expect_equal(stated, 5 * k, tolerance = 1e-12)
  expect_equal(sent * k, stated, tolerance = 1e-12)

  # Densities: sigma ~ half-normal(stated) on y and sigma / k ~
  # half-normal(sent) on the standardised response are one distribution.
  half_normal <- function(v, scale) 2 * stats::dnorm(v, sd = scale)
  sigma <- c(0.01, 0.2, 0.65, 1.5, 4, 10)
  expect_equal(
    half_normal(sigma, stated),
    half_normal(sigma / k, sent) / k,
    tolerance = 1e-12
  )
  expect_equal(
    stats::integrate(half_normal, 0, Inf, scale = stated)$value, 1,
    tolerance = 1e-6
  )
})


test_that("prior_summary() states the scale of the default and of each slot", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ps_data()
  fit <- flexybayes(y ~ x, random = ~ site, data = d, backend = "inla",
                    verbose = FALSE)
  out <- utils::capture.output(print(prior_summary(fit)))
  scale_line <- grep("Scale s", out, fixed = TRUE, value = TRUE)
  expect_length(scale_line, 1L)
  expect_match(scale_line, "in the units of the response", fixed = TRUE)
  slots <- grep("<prior: half-normal-SD(", out, fixed = TRUE, value = TRUE)
  expect_length(slots, 2L)
  expect_true(all(grepl(
    "on the standardised response = half-normal-SD(", slots, fixed = TRUE
  )))
})


test_that("a log-link default is stated on the log scale", {
  expect_identical(
    flexyBayes:::.fb_default_scale_units("log_link_half_normal"),
    "on the log scale of the mean"
  )
  expect_identical(
    flexyBayes:::.fb_default_scale_units("logit_half_normal"),
    "on the logit scale"
  )
  expect_identical(
    flexyBayes:::.fb_default_scale_units("identity_sd_half_normal"),
    "in the units of the response"
  )
})


test_that("a lognormal prior on the standardised response is not tagged half-normal", {
  compress <- flexyBayes:::.fb_compress_hyper_blobs
  block <- function(expr) {
    paste0(
      "f(g, model = \"iid\", hyper = list(prec = list(prior = \"",
      expr, "\")))"
    )
  }
  # The legacy lognormal(0, 1) on SD, restated on y / 2: its location moves
  # to -log(2).
  shifted <- block(
    flexyBayes:::.inla_lognormal_sd_expr(1, location = -log(2))
  )
  out <- compress(shifted)
  expect_false(grepl("half-normal", out, fixed = TRUE))
  expect_match(out, "<prior: lognormal-SD(-0.6931, 1)>", fixed = TRUE)
  expect_match(
    compress(shifted, response_scale = 2),
    "<prior: lognormal-SD(-0.6931, 1) on the standardised response = lognormal-SD(0, 1) on the response>",
    fixed = TRUE
  )
  # A uniform bound converts by the same factor.
  uniform <- block(flexyBayes:::.inla_uniform_sd_expr(upper = 2.5))
  expect_match(
    compress(uniform, response_scale = 2),
    "<prior: uniform-SD(0, 2.5) on the standardised response = uniform-SD(0, 5) on the response>",
    fixed = TRUE
  )
})
