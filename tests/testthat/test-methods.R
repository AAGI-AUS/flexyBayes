# Tests for S3 methods — these test the method dispatch and output structure
# without requiring greta (use mock objects where needed)

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

  # Build greta component (mock)
  greta_out <- structure(
    list(
      model = NULL,
      draws = draws,
      greta_arrays = list(),
      env = new.env(parent = emptyenv())
    ),
    class = "flexybayes_greta"
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
        rcov = list(list(type = "units")),
        family = list(family = "gaussian", link = "identity")
      ),
      call_info = list(
        fixed = y ~ 1,
        random = NULL,
        rcov = NULL,
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
    list(glm = glm_obj, greta = greta_out, extras = extras),
    class = "flexybayes"
  )
}


test_that("print.flexybayes produces output", {
  fit <- make_mock_flexybayes()
  expect_output(print(fit), "flexyBayes")
  expect_output(print(fit), "Fixed")
  expect_output(print(fit), "MCMC")
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

test_that("confint.flexybayes returns matrix with bounds", {
  fit <- make_mock_flexybayes()
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)
  expect_true(all(ci[, 1] < ci[, 2]))
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

test_that("predict.flexybayes returns fitted values when newdata is NULL", {
  fit <- make_mock_flexybayes()
  p <- predict(fit)
  expect_equal(length(p), nobs(fit))
})

test_that("predict.flexybayes with se.fit returns list", {
  fit <- make_mock_flexybayes()
  p <- predict(fit, se.fit = TRUE)
  expect_true(is.list(p))
  expect_true("fit" %in% names(p))
  expect_true("se.fit" %in% names(p))
})

test_that("model.matrix.flexybayes returns matrix", {
  fit <- make_mock_flexybayes()
  mm <- model.matrix(fit)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), nobs(fit))
})

test_that("tidy.flexybayes returns data frame", {
  fit <- make_mock_flexybayes()
  td <- tidy.flexybayes(fit)
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
  expect_no_error(plot(fit, type = "residuals"))
})

test_that("plot.flexybayes runs without error for effects", {
  fit <- make_mock_flexybayes()
  expect_no_error(plot(fit, type = "effects"))
})

test_that("plot.flexybayes runs without error for variance", {
  fit <- make_mock_flexybayes()
  expect_no_error(plot(fit, type = "variance"))
})

test_that("plot.flexybayes runs without error for blups", {
  fit <- make_mock_flexybayes()
  expect_no_error(plot(fit, type = "blups"))
})

test_that("plot.flexybayes runs without error for pp_check", {
  fit <- make_mock_flexybayes()
  expect_no_error(plot(fit, type = "pp_check"))
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
  expect_true(inherits(fit$greta, "flexybayes_greta"))
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
