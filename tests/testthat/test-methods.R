# Tests for S3 methods — these test the method dispatch and output structure
# without requiring a live backend fit (use mock objects where needed)

# Create a mock flexybayes object for testing methods
make_mock_flexybayes <- function() {
  n <- 50
  dat <- data.frame(
    y = rnorm(n, 50, 5),
    env = factor(rep(paste0("E", 1:5), each = 10)),
    x = rnorm(n)
  )

  # Simulate posterior draws as mcmc.list
  n_samples <- 100
  draws_mat <- matrix(rnorm(n_samples * 3), ncol = 3)
  draws_mat[, 1] <- draws_mat[, 1] + 50 # mu_atg ~ 50
  draws_mat[, 2] <- abs(draws_mat[, 2]) # sigma_geno > 0
  draws_mat[, 3] <- abs(draws_mat[, 3]) + 1 # sigma_e_atg > 0
  colnames(draws_mat) <- c("mu_atg", "sigma_geno", "sigma_e_atg")
  draws <- coda::mcmc.list(coda::mcmc(draws_mat))

  # Build GLM-compatible component
  beta <- c("(Intercept)" = 50.1)
  V <- matrix(0.5, 1, 1, dimnames = list("(Intercept)", "(Intercept)"))
  fitted_vals <- rep(50.1, n)
  resid_vals <- dat$y - fitted_vals

  glm_obj <- list(
    coefficients = beta,
    residuals = resid_vals,
    fitted.values = fitted_vals,
    family = gaussian(),
    formula = y ~ 1,
    terms = terms(y ~ 1, data = dat),
    model = model.frame(y ~ 1, data = dat),
    data = dat,
    df.residual = n - 1L,
    rank = 1L,
    y = dat$y,
    linear.predictors = fitted_vals,
    call = quote(flexybayes(fixed = y ~ 1, data = dat)),
    qr = qr(matrix(1, n, 1))
  )
  attr(glm_obj, "posterior_vcov") <- V
  class(glm_obj) <- c("flexybayes_glm", "glm", "lm")

  # Build a mock brms-shaped MCMC component. The print footer's
  # component-label table is a closed vocabulary (glm / brms / inla /
  # extras -- see .print_fit_components() in R/methods.R), so the mock
  # uses a real recognised label rather than an invented one.
  brms_out <- structure(
    list(
      model = NULL,
      draws = draws,
      env = new.env(parent = emptyenv())
    ),
    class = "flexybayes_brms_component"
  )

  # Build extras
  post_summary <- summary(draws)
  vc_table <- data.frame(
    component = c("sigma_geno", "sigma_e_atg"),
    estimate = c(1.0, 2.0),
    sd = c(0.3, 0.2),
    q2.5 = c(0.5, 1.5),
    q50 = c(1.0, 2.0),
    q97.5 = c(1.5, 2.5),
    stringsAsFactors = FALSE
  )

  extras <- structure(
    list(
      summary = post_summary,
      convergence = list(n_eff = coda::effectiveSize(draws), gelman = NULL),
      variance_comps = vc_table,
      blups = list(u_geno = rnorm(5)),
      predictions = data.frame(
        obs = 1:n,
        observed = dat$y,
        fitted = fitted_vals,
        residual = resid_vals
      ),
      code = "mu_atg <- normal(0, 10)",
      param_names = c("mu_atg", "sigma_geno", "sigma_e_atg"),
      parse_info = list(
        fixed = list(response = "y", intercept = TRUE, terms = list()),
        random = list(),
        residual = list(list(type = "units")),
        family = list(family = "gaussian", link = "identity")
      ),
      call_info = list(
        fixed = y ~ 1,
        random = NULL,
        residual = NULL,
        data_name = "dat",
        family = "gaussian",
        link = NULL,
        known_matrices = list(),
        weights = NULL,
        n_samples = 100,
        warmup = 50,
        chains = 1,
        prior_fixed_sd = 10,
        prior_vc_sd = 1
      ),
      run_time = 1.5,
      model_info = list(
        n_obs = n,
        n_fixed = 1L,
        n_random = 0L,
        n_params = 3L,
        family = "gaussian",
        link = "identity"
      )
    ),
    class = "flexybayes_extras"
  )

  structure(
    list(glm = glm_obj, brms = brms_out, extras = extras),
    class = "flexybayes"
  )
}


test_that("print.flexybayes produces output", {
  fit <- make_mock_flexybayes()
  expect_output(print(fit), "flexyBayes")
  expect_output(print(fit), "Fixed")
  expect_output(print(fit), "Sampler")
})

test_that("print.flexybayes lists only the slots the fit carries", {
  # The component footer used to name a fixed slot unconditionally, so
  # a fit shape lacking it still advertised one -- this fixture proves
  # the footer follows the object instead (within the closed glm / brms
  # / inla / extras label vocabulary).
  fit <- make_mock_flexybayes()
  brms_lines <- utils::capture.output(print(fit))
  expect_true(any(grepl("$brms", brms_lines, fixed = TRUE)))

  # Same object with the brms slot replaced by an inla slot: the
  # footer follows the object, not a hardcoded default.
  inla_shaped <- fit
  inla_shaped$brms <- NULL
  inla_shaped$inla <- list()
  inla_lines <- utils::capture.output(print(inla_shaped))
  expect_false(any(grepl("$brms", inla_lines, fixed = TRUE)))
  expect_true(any(grepl("$inla", inla_lines, fixed = TRUE)))
  expect_true(any(grepl("$extras", inla_lines, fixed = TRUE)))
})

test_that("summary.flexybayes produces output", {
  fit <- make_mock_flexybayes()
  expect_output(summary(fit), "Fixed effects")
  expect_output(summary(fit), "Variance components")
})

test_that("coef.flexybayes returns named numeric", {
  fit <- make_mock_flexybayes()
  beta <- coef(fit)
  expect_true(is.numeric(beta))
  expect_true(length(beta) > 0)
  expect_true(!is.null(names(beta)))
})

test_that("vcov.flexybayes returns matrix", {
  fit <- make_mock_flexybayes()
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), length(coef(fit)))
  expect_equal(ncol(V), length(coef(fit)))
})

test_that("confint.flexybayes refuses unconditionally on a bare-class fit", {
  # The bare method (no flexybayes_inla / flexybayes_brms override) now
  # abstains unconditionally with fit_lacks_posterior_draws -- it used
  # to compute a credible-interval matrix from the withdrawn native
  # engine's posterior draws (see NEWS.md, 0.9.3); neither active
  # engine reaches the parent, each has its own confint() method.
  fit <- make_mock_flexybayes()
  err <- tryCatch(confint(fit), error = function(e) e)
  expect_identical(err$reason_code, "fit_lacks_posterior_draws")
})

test_that("fitted.flexybayes returns numeric vector", {
  fit <- make_mock_flexybayes()
  f <- fitted(fit)
  expect_true(is.numeric(f))
  expect_equal(length(f), nobs(fit))
})

test_that("residuals.flexybayes returns numeric vector", {
  fit <- make_mock_flexybayes()
  r <- residuals(fit)
  expect_true(is.numeric(r))
  expect_equal(length(r), nobs(fit))
})

test_that("nobs.flexybayes returns integer", {
  fit <- make_mock_flexybayes()
  expect_equal(nobs(fit), 50)
})

test_that("family.flexybayes returns family object", {
  fit <- make_mock_flexybayes()
  fam <- family(fit)
  expect_true(inherits(fam, "family"))
  expect_equal(fam$family, "gaussian")
})

test_that("formula.flexybayes returns formula", {
  fit <- make_mock_flexybayes()
  f <- formula(fit)
  expect_true(inherits(f, "formula"))
})

test_that("logLik.flexybayes returns logLik object", {
  fit <- make_mock_flexybayes()
  ll <- logLik(fit)
  expect_true(inherits(ll, "logLik"))
  expect_true(is.numeric(as.numeric(ll)))
})

# predict.flexybayes() (the bare method) was deleted entirely at 0.9.3
# (see NEWS.md) -- both active engines have their own predict method
# (predict.flexybayes_inla, predict.flexybayes_brms, tested against
# real fits in test-fb-brms-stan.R and the INLA equivalent), and a
# bare-class mock like this one has no applicable predict() method any
# more. The two tests that used to live here are removed rather than
# adapted, since this fixture cannot exercise a class-specific method.

test_that("model.matrix.flexybayes returns matrix", {
  fit <- make_mock_flexybayes()
  mm <- model.matrix(fit)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), nobs(fit))
})

test_that("tidy.flexybayes returns data frame", {
  # conf.int = FALSE: the default TRUE would call confint(), which
  # refuses unconditionally on this bare-class mock (see the
  # confint.flexybayes test above) -- this test checks the tidy()
  # output shape, not the CI columns.
  fit <- make_mock_flexybayes()
  td <- tidy.flexybayes(fit, conf.int = FALSE)
  expect_true(is.data.frame(td))
  expect_true("term" %in% names(td))
  expect_true("estimate" %in% names(td))
})

test_that("tidy.flexybayes with effects='random' returns VC", {
  fit <- make_mock_flexybayes()
  td <- tidy.flexybayes(fit, effects = "random")
  expect_true(is.data.frame(td))
  expect_true(nrow(td) > 0)
})

test_that("glance.flexybayes returns one-row data frame", {
  fit <- make_mock_flexybayes()
  gl <- glance.flexybayes(fit)
  expect_true(is.data.frame(gl))
  expect_equal(nrow(gl), 1)
  expect_true("nobs" %in% names(gl))
})

test_that("augment.flexybayes returns data with .fitted and .resid", {
  fit <- make_mock_flexybayes()
  aug <- augment.flexybayes(fit)
  expect_true(is.data.frame(aug))
  expect_true(".fitted" %in% names(aug))
  expect_true(".resid" %in% names(aug))
  expect_equal(nrow(aug), nobs(fit))
})

test_that("plot.flexybayes runs without error for residuals", {
  fit <- make_mock_flexybayes()
  # Draw to a null device: without one the plot lands in an
  # Rplots.pdf in tests/testthat/, which every suite run then
  # recreates in the working tree.
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, type = "residuals"))
})

test_that("plot.flexybayes runs without error for effects", {
  fit <- make_mock_flexybayes()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, type = "effects"))
})

test_that("plot.flexybayes runs without error for variance", {
  fit <- make_mock_flexybayes()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, type = "variance"))
})

test_that("plot.flexybayes runs without error for blups", {
  fit <- make_mock_flexybayes()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit, type = "blups"))
})

test_that("plot.flexybayes refuses pp_check on a fit with no replicates", {
  # Until 0.9.1 this drew an observed-versus-fitted panel on any fit
  # carrying a response and a fitted vector, under a name that promises
  # replicated datasets. The mock carries both and no predictive draws,
  # so it now refuses by name -- catchable, and pointing at the residual
  # displays it can draw.
  fit <- make_mock_flexybayes()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot(fit, type = "pp_check"),
    class = "flexybayes_refusal_pp_check_requires_predictive_draws"
  )
})

test_that("summary.flexybayes_glm prints Bayesian summary", {
  fit <- make_mock_flexybayes()
  expect_output(summary(fit$glm), "Bayesian GLM")
  expect_output(summary(fit$glm), "posterior")
})

test_that("confint.flexybayes_glm returns credible intervals", {
  fit <- make_mock_flexybayes()
  ci <- confint(fit$glm)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)
})

test_that("anova.flexybayes compares models", {
  fit1 <- make_mock_flexybayes()
  fit2 <- make_mock_flexybayes()
  expect_output(anova(fit1, fit2), "model comparison")
})

test_that("flexybayes object has correct class", {
  fit <- make_mock_flexybayes()
  expect_true(inherits(fit, "flexybayes"))
  expect_true(inherits(fit$glm, "flexybayes_glm"))
  expect_true(inherits(fit$glm, "glm"))
  expect_true(inherits(fit$glm, "lm"))
  expect_true(inherits(fit$brms, "flexybayes_brms_component"))
  expect_true(inherits(fit$extras, "flexybayes_extras"))
})

test_that("extras contains expected components", {
  fit <- make_mock_flexybayes()
  expect_true(!is.null(fit$extras$summary))
  expect_true(!is.null(fit$extras$convergence))
  expect_true(!is.null(fit$extras$variance_comps))
  expect_true(!is.null(fit$extras$code))
  expect_true(!is.null(fit$extras$param_names))
  expect_true(!is.null(fit$extras$parse_info))
  expect_true(!is.null(fit$extras$call_info))
  expect_true(!is.null(fit$extras$run_time))
  expect_true(!is.null(fit$extras$model_info))
})
