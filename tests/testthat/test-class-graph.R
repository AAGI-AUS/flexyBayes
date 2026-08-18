# test-class-graph.R -- the shared `flexybayes` parent on INLA fits.
#
# Before 0.9.0 an INLA fit was c("flexybayes_inla", "list") while a brms fit
# was c("flexybayes_brms", "flexybayes", "list"). The split forced a
# parallel S3 method for every generic and left INLA objects reaching
# stats::*.default for the five generics with no INLA sibling. Giving INLA
# the parent closes the split; these tests pin what each of those five
# generics now does on an INLA fit, and pin that the INLA-specific methods
# still win dispatch over the parent.

skip_if_no_inla <- function() skip_if_not_installed("INLA")

mk_class_graph_data <- function(n_group = 12L, n_rep = 4L) {
  set.seed(20260815L)
  g <- factor(rep(seq_len(n_group), each = n_rep))
  x <- stats::rnorm(n_group * n_rep)
  u <- stats::rnorm(n_group, sd = 0.8)
  y <- 1.5 + 0.7 * x + u[as.integer(g)] + stats::rnorm(n_group * n_rep, sd = 0.5)
  data.frame(y = y, x = x, g = g)
}

fit_class_graph_inla <- function() {
  d <- mk_class_graph_data()
  flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
}

# ---------------------------------------------------------------- #
# The class vector itself                                          #
# ---------------------------------------------------------------- #

test_that("an INLA fit carries the shared flexybayes parent", {
  skip_if_no_inla()
  fit <- fit_class_graph_inla()

  expect_identical(class(fit), c("flexybayes_inla", "flexybayes", "list"))
  expect_true(inherits(fit, "flexybayes"))
})

test_that("INLA-specific methods still win dispatch over the parent", {
  skip_if_no_inla()
  fit <- fit_class_graph_inla()

  # logLik() must reach the INLA method (marginal-likelihood message plus
  # NA), not the parent's plug-in evaluation, which would compute a number
  # from a slot INLA does not populate the same way.
  expect_message(ll <- stats::logLik(fit), "marginal")
  expect_true(is.na(as.numeric(ll)))

  # print / summary / coef / vcov / fitted all have INLA siblings; the
  # parent must not capture any of them.
  expect_output(print(fit), "INLA")
  expect_true(is.numeric(stats::coef(fit)))
  expect_identical(
    flexyBayes:::.triangulate_source(fit),
    "inla"
  )
})

# ---------------------------------------------------------------- #
# The five parent methods with no INLA sibling                     #
# ---------------------------------------------------------------- #

test_that("nobs() answers correctly on an INLA fit through the parent", {
  skip_if_no_inla()
  fit <- fit_class_graph_inla()

  expect_identical(stats::nobs(fit), 48L)
})

test_that("model.matrix() rebuilds the INLA fit's fixed-effect design", {
  skip_if_no_inla()
  fit <- fit_class_graph_inla()

  mm <- stats::model.matrix(fit)
  expect_true(is.matrix(mm))
  expect_identical(nrow(mm), 48L)
  expect_identical(colnames(mm), c("(Intercept)", "x"))
})

test_that("confint() gives INLA marginal quantiles, not a default failure", {
  skip_if_no_inla()
  fit <- fit_class_graph_inla()

  ci <- stats::confint(fit)
  expect_true(is.matrix(ci))
  expect_identical(ncol(ci), 2L)
  expect_true(all(c("(Intercept)", "x") %in% rownames(ci)))
  expect_identical(colnames(ci), c("2.5%", "97.5%"))

  # The interval must bracket INLA's own posterior mean for each
  # coefficient, and it must widen when the credible level rises.
  cf <- stats::coef(fit)[rownames(ci)]
  expect_true(all(ci[, 1L] < cf & cf < ci[, 2L]))

  ci99 <- stats::confint(fit, level = 0.99)
  expect_true(all(ci99[, 1L] <= ci[, 1L]))
  expect_true(all(ci99[, 2L] >= ci[, 2L]))

  # Subsetting by name still works.
  one <- stats::confint(fit, parm = "x")
  expect_identical(nrow(one), 1L)
})

test_that("update() re-fits an INLA model now that the record is complete", {
  # Until 0.9.1 this block recorded the opposite: an INLA fit refused
  # `update_call_not_reconstructable`, naming `known_matrices` among the
  # arguments its record was missing. That refusal was correct about the
  # record and misleading about the cause -- the INLA emit wrote six of
  # the arguments a re-fit needs against brms's fifteen, and nothing
  # about the engine made a re-fit unsafe. The record is now complete on
  # both engines, so the refusal no longer fires here.
  skip_if_no_inla()
  fit <- fit_class_graph_inla()

  refit <- suppressMessages(stats::update(fit, random = ~g))
  expect_s3_class(refit, "flexybayes_inla")
  expect_identical(class(refit), class(fit))
})

test_that("update() still refuses by name on a genuinely short record", {
  # The refusal itself has to stay reachable: a fit read back from an
  # older object, or assembled by hand, must not be re-fitted with
  # defaults silently substituted for what it never recorded.
  skip_if_no_inla()
  fit <- fit_class_graph_inla()
  fit$extras$call_info$known_matrices <- NULL
  fit$extras$call_info$weights <- NULL

  err <- tryCatch(
    stats::update(fit, n_samples = 100L),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_update_call_not_reconstructable")
  expect_match(conditionMessage(err), "recorded call is")
  expect_match(conditionMessage(err), "known_matrices")
})

test_that("anova() refuses by name when a fit has no conditional logLik", {
  skip_if_no_inla()
  fit <- fit_class_graph_inla()

  err <- tryCatch(
    suppressMessages(stats::anova(fit)),
    condition = function(e) if (inherits(e, "error")) e else NULL
  )
  expect_s3_class(err, "flexybayes_refusal_conditional_loglik_not_available")
  expect_match(conditionMessage(err), "marginal")
})

# ---------------------------------------------------------------- #
# The parent methods refuse rather than returning an empty answer  #
# ---------------------------------------------------------------- #

test_that("confint.flexybayes refuses a fit carrying no posterior draws", {
  shell <- structure(
    list(extras = list()),
    class = c("flexybayes", "list")
  )
  err <- tryCatch(stats::confint(shell), condition = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_fit_lacks_posterior_draws")
  expect_match(conditionMessage(err), "no posterior-draw slot")
})

test_that("model.matrix.flexybayes refuses when the design cannot be built", {
  shell <- structure(
    list(extras = list()),
    class = c("flexybayes", "list")
  )
  err <- tryCatch(stats::model.matrix(shell), condition = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_model_matrix_not_recoverable")
})

test_that("logLik.flexybayes refuses instead of returning a silent NA", {
  # No response vector: the pre-0.9.0 method wrapped the failure in
  # tryCatch and returned NA with a warning, which AIC() then consumed as
  # though a log-likelihood had been computed.
  shell <- structure(
    list(glm = list(), extras = list(parse_info = list(family = list()))),
    class = c("flexybayes", "list")
  )
  err <- tryCatch(stats::logLik(shell), condition = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_conditional_loglik_not_available")

  # A family with no implemented plug-in form is named rather than
  # reported as NA.
  shell2 <- structure(
    list(
      glm = list(y = c(1, 2, 3)),
      extras = list(
        parse_info = list(family = list(family = "gamma")),
        model_info = list(n_params = 2L, n_obs = 3L)
      )
    ),
    class = c("flexybayes", "list")
  )
  err2 <- tryCatch(stats::logLik(shell2), condition = function(e) e)
  expect_s3_class(err2, "flexybayes_refusal_conditional_loglik_not_available")
  expect_match(conditionMessage(err2), "gamma")
})

test_that("update.flexybayes accepts a complete argument record", {
  # The positive control for the refusal above: with every required
  # argument recorded, the guard passes and the method proceeds to
  # rebuild the call. The rebuild is intercepted here rather than run, so
  # the test asserts the gate and not a second fit. `na_action` joined
  # the required set at 0.9.1, once both emits recorded it -- an `omit`
  # fit that re-fitted as the default would be computed on a different
  # set of rows from the one it was named after. `backend` and
  # `aggregate` joined on the same reasoning: a re-fit that took the
  # formal defaults for those came back on another engine, in another
  # representation.
  cl <- list(
    fixed = y ~ x, random = ~g, residual = NULL, data_name = "d",
    family = "gaussian", link = "identity", known_matrices = NULL,
    weights = NULL, n_samples = 100L, warmup = 100L, chains = 1L,
    prior_fixed_sd = NULL, prior_vc_sd = NULL, na_action = "augment",
    backend = "auto", aggregate = "auto", verbose = FALSE
  )
  shell <- structure(
    list(
      glm = list(data = data.frame(y = 1, x = 1, g = factor("a"))),
      extras = list(call_info = cl)
    ),
    class = c("flexybayes", "list")
  )
  err <- tryCatch(
    stats::update(shell, n_samples = 10L),
    condition = function(e) e
  )
  expect_false(
    inherits(err, "flexybayes_refusal_update_call_not_reconstructable")
  )
})
