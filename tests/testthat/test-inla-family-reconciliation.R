# =============================================================================
# The INLA family gate speaks INLA's vocabulary, and the residual hyper
# keyword is the one the likelihood declares.
#
# Two defects met on the same cells. The gate compared a flexyBayes family
# spelling against INLA's own likelihood roster, and the reconciler that
# maps between the two vocabularies runs after the gate -- so
# `negative_binomial`, `negbinom` and `binary` were refused for a naming
# reason on an engine that carries every one of those likelihoods, three of
# the eight supported spellings. And the residual hyperparameter for the
# beta likelihood was emitted as `prec` where INLA declares `phi`, so a beta
# fit died inside the engine with a raw error whose text quoted the correct
# alternatives back at us.
#
# Six of the sixteen family x backend cells the README declares as fitting
# did not fit; four of them were these.
# =============================================================================

.ifr_data <- function() {
  set.seed(20260818L)
  n <- 60L
  g <- factor(rep(seq_len(10L), each = 6L))
  x <- stats::rnorm(n)
  u <- stats::rnorm(10L, 0, 0.7)
  eta <- 0.5 + 0.4 * x + u[as.integer(g)]
  mu <- stats::plogis(eta)
  list(
    gaussian = data.frame(y = eta + stats::rnorm(n), x = x, g = g),
    poisson = data.frame(y = stats::rpois(n, exp(eta)), x = x, g = g),
    negative_binomial = data.frame(
      y = stats::rnbinom(n, mu = exp(eta), size = 2),
      x = x,
      g = g
    ),
    binomial = data.frame(y = stats::rbinom(n, 1L, mu), x = x, g = g),
    gamma = data.frame(
      y = stats::rgamma(n, shape = 2, rate = 2 / exp(eta)),
      x = x,
      g = g
    ),
    beta = data.frame(
      y = stats::rbeta(n, mu * 5, (1 - mu) * 5),
      x = x,
      g = g
    )
  )
}

# --- the gate tests the reconciled spelling ------------------------------ #

test_that("the family gate resolves the spelling before it checks the roster", {
  skip_if_not_installed("INLA")
  for (fam in c(
    "gaussian",
    "binomial",
    "binary",
    "poisson",
    "negative_binomial",
    "negbinom",
    "gamma",
    "beta"
  )) {
    check <- flexyBayes:::.lgm_check_family(list(family = fam))
    expect_true(
      isTRUE(check$pass),
      label = paste0("family_allowlist refused \"", fam, "\"")
    )
  }
})

test_that("a family INLA really does not carry is refused, naming both spellings", {
  skip_if_not_installed("INLA")
  check <- flexyBayes:::.lgm_check_family(list(family = "stdnormal"))
  expect_false(isTRUE(check$pass))
  expect_match(check$reason, "stdnormal", fixed = TRUE)
  expect_match(check$diagnostic, "resolved INLA family", fixed = TRUE)
})

# --- the residual hyper keyword is the engine's own ---------------------- #

test_that("the INLA hyper keyword table matches inla.models()", {
  skip_if_not_installed("INLA")
  roster <- INLA::inla.models()$likelihood
  keywords <- flexyBayes:::.fb_inla_residual_hyper()
  for (fam in names(keywords)) {
    entry <- roster[[fam]]
    expect_false(
      is.null(entry),
      label = paste0(fam, " absent from inla.models()")
    )
    declared <- vapply(
      entry$hyper,
      function(h) h$short.name %||% NA_character_,
      character(1)
    )
    expect_true(
      unname(keywords[[fam]]) %in% declared,
      info = paste0(
        fam,
        " declares (",
        paste(declared, collapse = ", "),
        "), table says ",
        keywords[[fam]]
      )
    )
  }
  # The one the table used to get wrong, asserted directly.
  expect_identical(unname(keywords[["beta"]]), "phi")
  expect_identical(
    as.character(roster[["beta"]]$hyper$theta$short.name),
    "phi"
  )
})

test_that("a likelihood with no residual scale takes no hyper block", {
  skip_if_not_installed("INLA")
  for (fam in c("poisson", "binomial", "nbinomial", "betabinomial")) {
    expect_null(flexyBayes:::.fb_inla_hyper_keyword(fam))
  }
  entry <- list(sigma = list(prior = "pc.prec", param = c(1, 0.01)))
  expect_length(
    flexyBayes:::.build_inla_control_family(entry, "nbinomial"),
    0L
  )
})

test_that("the control.family block carries the declared keyword", {
  entry <- list(sigma = list(prior = "pc.prec", param = c(1, 0.01)))
  gauss <- flexyBayes:::.build_inla_control_family(entry, "gaussian")
  expect_named(gauss$hyper, "prec")
  beta <- flexyBayes:::.build_inla_control_family(entry, "beta")
  expect_named(beta$hyper, "phi")
  expect_identical(beta$hyper$phi$prior, "pc.prec")
})

# --- and the cells fit ---------------------------------------------------- #

test_that("every admitted family spelling fits on INLA", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ifr_data()
  spellings <- list(
    gaussian = "gaussian",
    binomial = "binomial",
    binary = "binomial",
    poisson = "poisson",
    negative_binomial = "negative_binomial",
    negbinom = "negative_binomial",
    gamma = "gamma",
    beta = "beta"
  )
  for (fam in names(spellings)) {
    fit <- suppressMessages(flexybayes(
      y ~ x,
      random = ~g,
      data = d[[spellings[[fam]]]],
      family = fam,
      backend = "inla",
      verbose = FALSE
    ))
    expect_s3_class(fit, "flexybayes_inla")
    expect_true(
      is.finite(fit$inla$summary.fixed[["mean"]][[1L]]),
      label = paste0(fam, " intercept")
    )
  }
})

test_that("a beta fit on INLA carries a residual prior under phi", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .ifr_data()
  emitted <- suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = d$beta,
    family = "beta",
    backend = "inla",
    prior = fb_prior(sigma ~ pc(upper = 1, prob = 0.01)),
    return_code = TRUE
  ))
  expect_identical(emitted$family, "beta")
  expect_length(emitted$control_fixed, 0L)
  expect_named(emitted$control_family$hyper, "phi")
})
