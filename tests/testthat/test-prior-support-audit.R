# The prior-support audit.
#
# A uniform prior on a standard deviation reaches INLA as a density of -1e10
# outside its support, and INLA can still return a mode there
# (inst/KNOWN_ISSUES.md, issue 17). Every INLA fit is checked against INLA's
# own record of the priors it received: a fit whose mode or reported interval
# lies outside a bound is refused by name, and every fit that is returned
# carries the record of the check.

.psa_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    flexyBayes.silence_high_missingness_warning = TRUE,
    .local_envir = parent.frame()
  )
}

# A stand-in for an INLA return: a Gaussian observation precision and one
# latent term `g`, each under the same uniform bound on its standard
# deviation. `theta_g` is the log precision of `g`; the observation precision
# sits at log precision 0 (standard deviation 1).
.psa_fake_fit <- function(
  theta_g,
  q_g = c(theta_g - 0.5, theta_g + 0.5),
  upper = 3,
  lower = 0,
  prior = NULL
) {
  prior <- prior %||% .inla_uniform_sd_expr(upper = upper, lower = lower)
  one_hyper <- function() {
    list(prior = prior, param = numeric(0), short.name = "prec", fixed = FALSE)
  }
  internal <- data.frame(
    mean = c(0, theta_g),
    sd = c(0.1, 0.5),
    `0.025quant` = c(-0.2, q_g[[1L]]),
    `0.5quant` = c(0, theta_g),
    `0.975quant` = c(0.2, q_g[[2L]]),
    mode = c(0, theta_g),
    row.names = c(
      "Log precision for the Gaussian observations",
      "Log precision for g"
    ),
    check.names = FALSE
  )
  list(
    all.hyper = list(
      family = list(list(hyperid = "INLA.Data1", hyper = list(
        theta1 = one_hyper()
      ))),
      random = list(list(hyperid = "g", hyper = list(theta = one_hyper()))),
      fixed = list()
    ),
    internal.summary.hyperpar = internal,
    mode = list(theta = stats::setNames(internal$mode, rownames(internal)))
  )
}

# The N30 design: tests/testthat/test-update-carry.R .uc_count_data() before
# its `h` term was given an effect.
.psa_n30_data <- function(seed = 4242L, n = 90L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    g = factor(rep(letters[1:6], length.out = n)),
    h = factor(rep(LETTERS[1:3], length.out = n))
  )
  b_g <- stats::rnorm(6L, sd = 0.3)[as.integer(d$g)]
  eta <- 1.2 + 0.4 * (d$f == "b") - 0.3 * (d$f == "c") + b_g
  d$y <- stats::rpois(n, exp(eta))
  d
}

.psa_gaussian_data <- function(seed = 20260916L, n = 80L) {
  set.seed(seed)
  d <- data.frame(
    g = factor(rep(sprintf("g%02d", 1:8), length.out = n)),
    h = factor(rep(c("A", "B", "C", "D"), length.out = n)),
    x = stats::rnorm(n)
  )
  d$y <- 2 + 0.5 * d$x + stats::rnorm(8L, sd = 1)[as.integer(d$g)] +
    stats::rnorm(n, sd = 0.7)
  d$cnt <- stats::rpois(n, exp(1 + stats::rnorm(8L, sd = 0.5)[as.integer(d$g)]))
  d
}


# ---------------------------------------------------------------- #
# 1. Reading the bound                                               #
# ---------------------------------------------------------------- #

test_that("the bound is read from the uniform expression prior only", {
  bounds <- .fb_inla_expression_bounds(c(
    .inla_uniform_sd_expr(upper = 3),
    .inla_uniform_sd_expr(upper = 3, lower = 0.1),
    .inla_halfnormal_sd_expr(1),
    "loggamma",
    NA_character_
  ))

  expect_equal(bounds$upper, c(3, 3, NA, NA, NA))
  expect_equal(bounds$lower, c(0, 0.1, 0, 0, 0))
})


# ---------------------------------------------------------------- #
# 2. The audit on a constructed INLA return                          #
# ---------------------------------------------------------------- #

test_that("a fit inside every bound returns a record and is not refused", {
  fake <- .psa_fake_fit(theta_g = log(4.5))
  audit <- .fb_inla_support_audit(fake)

  expect_s3_class(audit, "data.frame")
  expect_equal(nrow(audit), 2L)
  expect_true(all(audit$inside))
  expect_equal(audit$mode_sd[[2L]], 1 / sqrt(4.5))
  expect_equal(audit$upper_sd, c(3, 3))
  expect_identical(.fb_inla_enforce_support(audit, "per_row"), audit)
})

test_that("a mode outside the bound is refused by name", {
  skip_if_not_installed("brms")
  # The N30 mode: precision 4.7e-07, a standard deviation near 1,460.
  fake <- .psa_fake_fit(theta_g = log(4.7e-07), q_g = c(-15, -14))
  audit <- .fb_inla_support_audit(fake)

  expect_identical(audit$inside, c(TRUE, FALSE))
  err <- expect_error(
    .fb_inla_enforce_support(audit, "per_row"),
    class = "flexybayes_refusal_prior_support_violated"
  )
  expect_match(
    conditionMessage(err),
    "g: standard deviation 1459",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "backend = \"brms\"", fixed = TRUE)
  expect_identical(err$route, "per_row")
  expect_identical(err$audit, audit)
})

test_that("an interval outside the bound is refused with the mode inside", {
  # Mode at standard deviation exp(1) = 2.72 under a bound of 3; the lower
  # quantile of the log precision implies a standard deviation of 3.67.
  fake <- .psa_fake_fit(theta_g = -2, q_g = c(-2.6, -1.5))

  expect_false(.fb_inla_support_audit(fake)$inside[[2L]])
})

test_that("the tolerance is 0.1 per cent of the bound and no wider", {
  at_bound <- -2 * log(3)
  inside <- .psa_fake_fit(
    theta_g = at_bound - 0.001,
    q_g = c(at_bound - 0.001, 1)
  )
  outside <- .psa_fake_fit(
    theta_g = at_bound - 0.01,
    q_g = c(at_bound - 0.01, 1)
  )

  expect_true(.fb_inla_support_audit(inside)$inside[[2L]])
  expect_false(.fb_inla_support_audit(outside)$inside[[2L]])
})

test_that("a two-sided bound refuses a standard deviation below it", {
  # Standard deviation exp(-3) = 0.05 under uniform(0.1, 3).
  fake <- .psa_fake_fit(theta_g = 6, q_g = c(5.5, 6.5), lower = 0.1)

  expect_false(.fb_inla_support_audit(fake)$inside[[2L]])
})

test_that("a reported value that is not finite is refused", {
  fake <- .psa_fake_fit(theta_g = log(4.5), q_g = c(NaN, 2))

  expect_false(.fb_inla_support_audit(fake)$inside[[2L]])
})

test_that("a fit with no bounded prior returns an empty record", {
  fake <- .psa_fake_fit(theta_g = log(4.7e-07), prior = "loggamma")
  audit <- .fb_inla_support_audit(fake)

  expect_equal(nrow(audit), 0L)
  expect_identical(names(audit), names(.fb_inla_support_audit_frame()))
})

test_that("hyperparameters that cannot be paired with priors are an error", {
  fake <- .psa_fake_fit(theta_g = log(4.5))
  extra <- fake$internal.summary.hyperpar[1L, , drop = FALSE]
  rownames(extra) <- "Rho for g"
  fake$internal.summary.hyperpar <- rbind(fake$internal.summary.hyperpar, extra)

  expect_error(.fb_inla_support_audit(fake), "could not be paired")
})

test_that("a bounded likelihood parameter that is not a precision is audited", {
  # The gamma likelihood's hyperparameter is INLA's "precision parameter",
  # and beta's is `phi`; the package's uniform expression prior sits on
  # their internal theta as on a Gaussian precision.
  fake <- .psa_fake_fit(theta_g = log(4.5))
  rownames(fake$internal.summary.hyperpar)[[1L]] <-
    "Intern precision-parameter for the Gamma observations"
  fake$all.hyper$family[[1L]]$hyper$theta1$short.name <- "phi"

  audit <- .fb_inla_support_audit(fake)

  expect_equal(nrow(audit), 2L)
  expect_true(all(audit$inside))
})

test_that("a likelihood row that does not name the observations is an error", {
  fake <- .psa_fake_fit(theta_g = log(4.5))
  rownames(fake$internal.summary.hyperpar)[[1L]] <- "Log precision for h"

  expect_error(.fb_inla_support_audit(fake), "could not be paired")
})

test_that("a standardised fit reports standard deviations in response units", {
  fake <- .psa_fake_fit(theta_g = log(4.5))
  spec <- list(centre = 0, scale = 2, n_obs = 10L)

  plain <- .fb_inla_support_audit(fake)
  scaled <- .fb_inla_support_audit(fake, spec)

  expect_equal(scaled$mode_sd, 2 * plain$mode_sd)
  expect_equal(scaled$upper_sd, 2 * plain$upper_sd)
  expect_identical(scaled$inside, plain$inside)
})


# ---------------------------------------------------------------- #
# 3. Live fits                                                       #
# ---------------------------------------------------------------- #

test_that("the N30 design fits inside every bound or refuses by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  .psa_silence()
  d <- .psa_n30_data()

  # The bound of 3 on each standard deviation is stated. It was the package
  # default on a log link until 0.13.0, whose half-normal default carries no
  # bound for the audit to check.
  bound_3 <- fb_prior(
    sd(group = "g") ~ uniform(lower = 0, upper = 3),
    sd(group = "h") ~ uniform(lower = 0, upper = 3)
  )
  for (i1 in seq_len(6L)) {
    res <- tryCatch(
      suppressWarnings(flexybayes(
        y ~ f,
        random = ~ g + h,
        data = d,
        family = "poisson",
        backend = "inla",
        prior = bound_3,
        verbose = FALSE
      )),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      expect_s3_class(res, "flexybayes_refusal_prior_support_violated")
    } else {
      audit <- res$extras$support_audit
      expect_true(all(audit$inside))
      # Independent of the audit: every reported precision mode implies a
      # standard deviation inside the stated bound of 3.
      hp <- res$inla$summary.hyperpar
      prec_rows <- grepl("^Precision for ", rownames(hp))
      expect_true(all(1 / sqrt(hp[prec_rows, "mode"]) <= 3 * (1 + 1e-3)))
    }
  } # ends i1, over repeats
})

test_that("every INLA fit an entry point returns has been through the audit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  .psa_silence()
  d <- .psa_gaussian_data()

  counts <- new.env(parent = emptyenv())
  counts$engine <- 0L
  counts$audit <- 0L
  real_call <- .inla_call
  real_audit <- .fb_inla_support_audit
  testthat::local_mocked_bindings(
    .inla_call = function(args) {
      counts$engine <- counts$engine + 1L
      real_call(args)
    },
    .fb_inla_support_audit = function(...) {
      counts$audit <- counts$audit + 1L
      real_audit(...)
    },
    # The restart probe calls INLA for a fit it may discard; it is switched
    # off so that every engine call is a fit an entry point returns.
    .fb_inla_restart_probe = function(...) NULL
  )

  # Each fit states a uniform bound, so every audit record has rows to show:
  # the package default is a half-normal (0.13.0), which carries none. The
  # bounds are the former defaults, 5 sd(y) and 3 on the log link.
  u_y <- 5 * stats::sd(d$y)
  bounded <- fb_prior(
    sigma ~ uniform(lower = 0, upper = u_y),
    sd(group = "g") ~ uniform(lower = 0, upper = u_y)
  )
  bounded_count <- fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 3))

  per_row <- flexybayes(
    y ~ x, random = ~g, data = d, backend = "inla", prior = bounded,
    verbose = FALSE
  )
  alias <- suppressWarnings(
    fb(y ~ x, random = ~g, data = d, backend = "inla", prior = bounded,
       verbose = FALSE),
    classes = "flexybayes_deprecated"
  )
  pinned <- suppressWarnings(
    fb_inla(y ~ x, random = ~g, data = d, prior = bounded, verbose = FALSE),
    classes = "flexybayes_deprecated"
  )
  aggregated <- flexybayes(
    cnt ~ h, random = ~g, data = d, family = "poisson", backend = "inla",
    aggregate = TRUE, prior = bounded_count, verbose = FALSE
  )
  streamed <- flexybayes_stream(
    fixed = cnt ~ h, random = ~g, source = d, family = "poisson",
    prior = bounded_count, verbose = FALSE
  )
  refit <- stats::update(per_row, random = ~ g + h)
  d_fold <- d
  d_fold$y[seq(1L, nrow(d), by = 4L)] <- NA
  fold <- .fb_sel_bayes_refit(
    per_row, d_fold, y ~ x, "inla", na_action = "augment"
  )

  fits <- list(
    flexybayes = per_row, fb = alias, fb_inla = pinned,
    aggregated = aggregated, stream = streamed, update = refit,
    fold_refit = fold
  )
  for (nm in names(fits)) {
    expect_s3_class(fits[[nm]], "flexybayes_inla")
    expect_s3_class(fits[[nm]]$extras$support_audit, "data.frame")
    expect_gt(nrow(fits[[nm]]$extras$support_audit), 0L)
  }
  expect_equal(counts$engine, length(fits))
  expect_equal(counts$audit, counts$engine)
})

test_that("an explicit INLA request surfaces the refusal", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .psa_silence()
  d <- .psa_gaussian_data()
  # The audit of the constructed N30 return, whatever the live fit found.
  violating <- .fb_inla_support_audit(
    .psa_fake_fit(theta_g = log(4.7e-07), q_g = c(-15, -14))
  )
  testthat::local_mocked_bindings(
    .fb_inla_support_audit = function(...) violating
  )

  expect_error(
    flexybayes(y ~ x, random = ~g, data = d, backend = "inla", verbose = FALSE),
    class = "flexybayes_refusal_prior_support_violated"
  )
})
