# =============================================================================
# A stated uniform bound on a standard deviation binds on INLA.
#
# A uniform-on-SD prior reaches INLA as an expression that is -1e10, flat,
# outside its support. INLA starts every precision at its own default, a
# standard deviation of 0.135 for the latent models, so a bound tighter than
# that started the optimiser outside the support with no gradient leading in:
# sd(group = "geno") ~ uniform(0, 0.05) returned a standard deviation of 0.61.
# The start now moves inside the support when, and only when, INLA's default
# lies outside it.
# =============================================================================

.bound_pois_data <- function() {
  set.seed(314L)
  n <- 3000L
  env <- sample(1:4, n, replace = TRUE)
  geno <- sample(1:20, n, replace = TRUE)
  eta <- c(0, 1.2, -0.7, 0.5)[env] + stats::rnorm(20, 0, 0.6)[geno]
  data.frame(
    env = factor(env),
    geno = factor(geno),
    y = stats::rpois(n, exp(eta - 0.5))
  )
}

.bound_gauss_data <- function() {
  set.seed(42L)
  n <- 500L
  d <- data.frame(
    f = factor(sample(letters[1:3], n, replace = TRUE)),
    g = factor(sample.int(15L, n, replace = TRUE))
  )
  d$y <- 1 + 0.5 * (d$f == "b") - 0.3 * (d$f == "c") +
    stats::rnorm(15L, 0, 0.3)[as.integer(d$g)] + stats::rnorm(n, 0, 0.15)
  d
}

.bound_varcomp <- function(fit, component) {
  tab <- summary(fit)$varcomp
  value <- tab$estimate[tab$component == component]
  stopifnot(length(value) == 1L)
  value
}

# ---------------------------------------------------------------- #
# 1. The start rule and the defaults it assumes                      #
# ---------------------------------------------------------------- #

test_that("INLA's default starting values are the ones the start rule assumes", {
  skip_if_not_installed("INLA")
  models <- INLA::inla.models()
  # INLA stores each field with an `inla.read.only` attribute, so the values
  # are compared without it
  short_name_entry <- function(hyper, short) {
    hit <- Filter(function(h) identical(as.vector(h$short.name), short), hyper)
    expect_length(hit, 1L)
    hit[[1L]]
  }

  for (model in c("iid", "rw1", "rw2", "besag", "generic0", "ar1")) {
    entry <- short_name_entry(models$latent[[model]]$hyper, "prec")
    expect_equal(
      as.vector(entry$initial),
      .FB_INLA_LATENT_PRECISION_INITIAL,
      label = model
    )
  }

  keywords <- .fb_inla_residual_hyper()
  starts <- .fb_inla_residual_initial()
  expect_setequal(names(starts), names(keywords))
  for (fam in names(keywords)) {
    entry <- short_name_entry(models$likelihood[[fam]]$hyper, keywords[[fam]])
    expect_equal(as.vector(entry$initial), starts[[fam]], label = fam)
  }
})

test_that("a start is written only when INLA's default lies outside the support", {
  tight <- list(meta = list(family = "uniform", lower = 0, upper = 0.05))
  moderate <- list(meta = list(family = "uniform", lower = 0, upper = 0.4))
  floor_above <- list(meta = list(family = "uniform", lower = 0.2, upper = 0.5))
  half_normal <- list(meta = list(family = "half_normal", scale = 0.01))

  expect_equal(.fb_inla_uniform_start(tight, 4), -2 * log(0.025))
  expect_null(.fb_inla_uniform_start(moderate, 4))
  expect_equal(.fb_inla_uniform_start(floor_above, 4), -2 * log(0.35))
  expect_null(.fb_inla_uniform_start(half_normal, 4))
})

test_that("the emitted formula carries a start for a tight bound and none for a moderate one", {
  skip_if_not_installed("INLA")
  d <- .bound_pois_data()
  emitted <- function(upper) {
    code <- suppressMessages(flexybayes(
      y ~ env,
      random = ~geno,
      data = d,
      family = "poisson",
      backend = "inla",
      # fb_prior() reads its formulas unevaluated, so the bound is substituted
      prior = eval(bquote(
        fb_prior(sd(group = "geno") ~ uniform(lower = 0, upper = .(upper)))
      )),
      return_code = TRUE
    ))
    paste(deparse(code$formula, width.cutoff = 500L), collapse = " ")
  }
  expect_match(emitted(0.05), "initial = ", fixed = TRUE)
  expect_no_match(emitted(0.4), "initial = ", fixed = TRUE)
})

# ---------------------------------------------------------------- #
# 2. Tight and moderate bounds bind in the engine                    #
# ---------------------------------------------------------------- #

test_that("tight and moderate group-SD bounds bind on a Poisson fit, per row and aggregated", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .bound_pois_data()

  fit_with <- function(upper, aggregate) {
    suppressMessages(flexybayes(
      y ~ env,
      random = ~geno,
      data = d,
      family = "poisson",
      backend = "inla",
      aggregate = aggregate,
      # fb_prior() reads its formulas unevaluated, so the bound is substituted
      prior = eval(bquote(
        fb_prior(sd(group = "geno") ~ uniform(lower = 0, upper = .(upper)))
      )),
      verbose = FALSE
    ))
  }

  for (aggregate in c(FALSE, TRUE)) {
    for (upper in c(0.05, 0.4)) {
      sd_geno <- .bound_varcomp(fit_with(upper, aggregate), "sd_geno")
      expect_lte(
        sd_geno,
        upper * 1.01,
        label = paste0("sd_geno, bound ", upper, ", aggregate = ", aggregate)
      )
    }
  }
})

test_that("tight and moderate bounds bind on a Gaussian fit, on a group and on the residual", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .bound_gauss_data()

  fit_with <- function(prior) {
    suppressMessages(flexybayes(
      y ~ f + (1 | g),
      data = d,
      backend = "inla",
      prior = prior,
      verbose = FALSE
    ))
  }

  for (upper in c(0.05, 0.2)) {
    fit <- fit_with(eval(bquote(
      fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = .(upper)))
    )))
    expect_lte(.bound_varcomp(fit, "sd_g"), upper * 1.01, label = paste0("sd_g, bound ", upper))
  }

  fit <- fit_with(fb_prior(sigma ~ uniform(lower = 0, upper = 0.1)))
  expect_lte(.bound_varcomp(fit, "sigma"), 0.1 * 1.01, label = "sigma, bound 0.1")
})

# ---------------------------------------------------------------- #
# 3. An engine failure beside a stated bound names the bound         #
# ---------------------------------------------------------------- #

test_that("an INLA program failure on a model with a stated uniform bound names the bound", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .inla_call = function(args) stop("The inla-program exited with an error.")
  )
  d <- .bound_pois_data()

  failure <- function(prior) {
    tryCatch(
      suppressMessages(flexybayes(
        y ~ env,
        random = ~geno,
        data = d,
        family = "poisson",
        backend = "inla",
        aggregate = FALSE,
        prior = prior,
        verbose = FALSE
      )),
      error = function(e) e
    )
  }

  stated <- failure(fb_prior(sd(group = "geno") ~ uniform(lower = 0, upper = 0.05)))
  expect_s3_class(stated, "flexybayes_refusal_inla_program_failed")
  expect_match(conditionMessage(stated), "states a uniform bound", fixed = TRUE)

  unstated <- failure(NULL)
  expect_s3_class(unstated, "flexybayes_refusal_inla_program_failed")
  expect_no_match(conditionMessage(unstated), "states a uniform bound", fixed = TRUE)
})
