# DEF-212 to DEF-218, DEF-221, DEF-223 and DEF-224: inputs the entry
# routes stopped on untyped, fitted as a different model, or planned
# differently from the fit.

.er_data <- function() {
  set.seed(17)
  d <- data.frame(
    g = factor(rep(seq_len(6L), each = 5L)),
    trt = factor(rep(c("a", "b", "c"), 10L)),
    x = factor(rep(c("e1", "e2", "e3"), 10L)),
    z = stats::rnorm(30L)
  )
  d$y <- 1 + stats::rnorm(30L)
  d$n <- stats::rpois(30L, 3)
  d
}

.er_quiet <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

.er_options <- function(env = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1",
    .local_envir = env
  )
}

test_that("DEF-212: a fixed or random term naming a column not in the data is refused before any engine runs", {
  skip_if_not_installed("INLA")
  .er_options()
  d <- .er_data()
  cls <- "flexybayes_refusal_model_variable_not_in_data"
  err <- expect_error(
    .er_quiet(flexybayes(
      y ~ nosuch,
      random = ~g,
      data = d,
      backend = "inla",
      verbose = FALSE
    )),
    class = cls
  )
  expect_identical(err$variables, "nosuch")
  expect_match(conditionMessage(err), "`nosuch`", fixed = TRUE)
  local_mocked_bindings(.backend_emit_fn = function(name) {
    stop("an engine was reached")
  })
  err <- expect_error(
    .er_quiet(flexybayes(
      y ~ x,
      random = ~ g + nosuch,
      data = d,
      backend = "inla",
      verbose = FALSE
    )),
    class = cls
  )
  expect_identical(err$variables, "nosuch")
  expect_error(
    .er_quiet(flexybayes(
      y ~ x,
      random = ~nosuch,
      data = d,
      backend = "brms",
      return_code = TRUE,
      verbose = FALSE
    )),
    class = cls
  )
  expect_error(
    .er_quiet(flexybayes(
      n ~ x,
      random = ~nosuch,
      data = d,
      family = "poisson",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = cls
  )
  expect_error(fb_plan(y ~ x + nosuch:z, random = ~g, data = d), class = cls)
  expect_error(fb_plan(y ~ x + (1 | nosuch), data = d), class = cls)
})

test_that("DEF-213: an aggregate value other than TRUE, FALSE or \"auto\" is refused by name", {
  d <- .er_data()
  cls <- "flexybayes_refusal_aggregate_argument_invalid"
  expect_error(
    flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      aggregate = "foo",
      verbose = FALSE
    ),
    class = cls
  )
  err <- expect_error(
    fb_plan(y ~ x, random = ~g, data = d, aggregate = NA),
    class = cls
  )
  expect_match(conditionMessage(err), "NA", fixed = TRUE)
})

test_that("DEF-214: weights of the wrong length or type are refused by name on both grammars", {
  d <- .er_data()
  cls <- "flexybayes_refusal_weights_invalid"
  err <- expect_error(
    flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      weights = c(1, 2),
      verbose = FALSE
    ),
    class = cls
  )
  expect_match(conditionMessage(err), "2 values for 30 rows", fixed = TRUE)
  expect_error(
    flexybayes(
      y ~ x + (1 | g),
      data = d,
      weights = rep("a", 30L),
      verbose = FALSE
    ),
    class = cls
  )
  expect_error(fb_plan(y ~ x, random = ~g, data = d, weights = 1), class = cls)
})

test_that("DEF-215: a prior that is not an fb_prior() object is refused by name", {
  d <- .er_data()
  cls <- "flexybayes_refusal_prior_not_fb_prior"
  err <- expect_error(
    flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      prior = "sd ~ banana(1)",
      backend = "inla",
      verbose = FALSE
    ),
    class = cls
  )
  expect_match(conditionMessage(err), "character", fixed = TRUE)
  expect_error(
    flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      prior = list(sd = 1),
      backend = "brms",
      return_code = TRUE,
      verbose = FALSE
    ),
    class = cls
  )
  expect_error(
    fb_plan(y ~ x + (1 | g), data = d, priors = "sd ~ banana(1)"),
    class = cls
  )
})

test_that("DEF-216: the aggregated routes' refusals are flexybayes refusals with a registered code", {
  skip_if_not_installed("INLA")
  .er_options()
  d <- .er_data()
  # A continuous covariate forces one cell per row.
  err <- expect_error(
    .er_quiet(flexybayes(
      n ~ z,
      random = ~g,
      data = d,
      family = "poisson",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal"
  )
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_identical(err$reason_code, "aggregation_plan_ineligible")
  expect_true("continuous_cell_key_data_dependent" %in% err$reason_codes)
  # The in-memory scope gate.
  fb <- .build_ir_polymorphic(
    fixed = n ~ x,
    random = ~g,
    residual = NULL,
    data = d,
    family = "poisson",
    link = NULL,
    weights = NULL,
    known_matrices = list(),
    prior = NULL,
    prior_fixed_sd = 100,
    prior_vc_sd = 1
  )
  err <- expect_error(
    .assert_aggregate_in_scope(fb),
    class = "flexybayes_refusal"
  )
  expect_s3_class(err, "flexybayes_aggregate_out_of_scope")
  expect_identical(err$reason_code, "aggregation_out_of_scope")
  expect_identical(err$scope_reason, "non_gaussian_family")
  # The streamed scope gate.
  err <- expect_error(
    .er_quiet(flexybayes_stream(
      y ~ x,
      random = ~ fa(trt, 2):g,
      source = d,
      fit = FALSE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal"
  )
  expect_s3_class(err, "flexybayes_aggregate_out_of_scope")
  expect_identical(err$scope_reason, "structured_random_not_aggregatable")
  codes <- c("aggregation_plan_ineligible", "aggregation_out_of_scope")
  expect_true(all(codes %in% fb_refusals()$reason_code))
})

test_that("DEF-217: a stream source that is neither a file nor a data frame is refused by name", {
  d <- .er_data()
  cls <- "flexybayes_refusal_stream_source_invalid"
  err <- expect_error(
    flexybayes_stream(
      y ~ x,
      random = ~g,
      source = "d",
      fit = FALSE,
      verbose = FALSE
    ),
    class = cls
  )
  expect_match(conditionMessage(err), "`d`", fixed = TRUE)
  expect_error(
    flexybayes_stream(
      y ~ x,
      random = ~g,
      source = 42,
      fit = FALSE,
      verbose = FALSE
    ),
    class = cls
  )
})

test_that("DEF-218: a factor or character response is refused for a family not built for one", {
  skip_if_not_installed("INLA")
  .er_options()
  d <- .er_data()
  d$trt_chr <- as.character(d$trt)
  cls <- "flexybayes_refusal_response_not_numeric"
  err <- expect_error(
    .er_quiet(flexybayes(
      trt ~ z,
      random = ~g,
      data = d,
      backend = "inla",
      verbose = FALSE
    )),
    class = cls
  )
  expect_identical(err$column, "trt")
  expect_match(conditionMessage(err), "factor", fixed = TRUE)
  expect_error(
    .er_quiet(flexybayes(
      trt ~ z,
      random = ~g,
      data = d,
      backend = "brms",
      return_code = TRUE,
      verbose = FALSE
    )),
    class = cls
  )
  expect_error(fb_plan(trt ~ z, random = ~g, data = d), class = cls)
  expect_error(
    fb_plan(trt_chr ~ z + (1 | g), data = d, family = "poisson"),
    class = cls
  )
  d$two <- factor(rep(c("no", "yes"), 15L))
  expect_error(fb_plan(two ~ z, data = d, family = "binary"), class = cls)
  # A binomial family keeps its own refusal, now in the plan as well.
  expect_error(
    fb_plan(two ~ z, data = d, family = "binomial"),
    class = "flexybayes_refusal_binomial_response_not_binary"
  )
})

test_that("DEF-221: predict() on an INLA fit reads re_formula and summary", {
  skip_if_not_installed("INLA")
  .er_options()
  set.seed(12)
  d <- data.frame(
    g = factor(rep(seq_len(6L), each = 10L)),
    f = factor(rep(c("a", "b"), 30L))
  )
  d$y <- 1 + stats::rnorm(6L)[d$g] + stats::rnorm(60L)
  fit <- .er_quiet(flexybayes(y ~ f, random = ~g, data = d, verbose = FALSE))
  # Population level: the fixed part alone, from the coefficients.
  fixed_part <- as.numeric(stats::model.matrix(~f, d) %*% coef(fit))
  expect_equal(predict(fit, re_formula = NA), fixed_part, tolerance = 1e-12)
  expect_equal(predict(fit, re_formula = NA, summary = TRUE), predict(fit))
  cls <- "flexybayes_refusal_method_argument_not_supported"
  err <- expect_error(predict(fit, re_formula = NULL), class = cls)
  expect_identical(err$argument, "re_formula")
  expect_error(predict(fit, re_formula = ~ (1 | g)), class = cls)
  err <- expect_error(predict(fit, summary = FALSE), class = cls)
  expect_identical(err$argument, "summary")
})

test_that("DEF-221: vcov() and fb_as_draws_simple() on a brms fit refuse n_samples", {
  fake <- structure(
    list(glm = structure(list(), posterior_vcov = diag(2L)), brms = NULL),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  expect_identical(vcov(fake), diag(2L))
  cls <- "flexybayes_refusal_method_argument_not_supported"
  err <- expect_error(vcov(fake, n_samples = 100L), class = cls)
  expect_identical(err$argument, "n_samples")
  err <- expect_error(fb_as_draws_simple(fake, n_samples = 100L), class = cls)
  expect_identical(err$argument, "n_samples")
})

test_that("DEF-223: fb_plan() refuses aggregate = TRUE where the fit refuses it", {
  skip_if_not_installed("INLA")
  .er_options()
  d <- .er_data()
  for (backend in c("auto", "inla")) {
    p <- fb_plan(
      n ~ z,
      random = ~g,
      data = d,
      family = "poisson",
      aggregate = TRUE,
      backend = backend
    )
    expect_false(p$will_fit)
    expect_match(p$will_fit_reason, "aggregation_plan_ineligible", fixed = TRUE)
    expect_error(
      .er_quiet(flexybayes(
        n ~ z,
        random = ~g,
        data = d,
        family = "poisson",
        aggregate = TRUE,
        backend = backend,
        verbose = FALSE
      )),
      class = "flexybayes_refusal_aggregation_plan_ineligible"
    )
  }
  # A model the aggregation compresses still plans the aggregated route.
  d4 <- d[rep(seq_len(nrow(d)), 4L), ]
  p <- fb_plan(
    n ~ x,
    random = ~g,
    data = d4,
    family = "poisson",
    aggregate = TRUE
  )
  expect_true(p$will_fit)
  expect_identical(p$path, "aggregated_inla")
})

test_that("DEF-224: the bar grammar fits a transformed response as the ASReml grammar does", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  .er_options()
  set.seed(22)
  d <- data.frame(g = factor(rep(1:6, each = 10L)), x = stats::rnorm(60L))
  d$y <- exp(
    1 + 0.3 * d$x + stats::rnorm(6L, 0, 0.3)[d$g] + stats::rnorm(60L, 0, 0.2)
  )
  bar <- .er_quiet(flexybayes(
    log(y) ~ x + (1 | g),
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
  asr <- .er_quiet(flexybayes(
    log(y) ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
  expect_identical(.fb_fit_ir(bar)$response, "log(y)")
  expect_true(isTRUE(.fb_fit_ir(bar)$response_is_expression))
  # Two identical ASReml-grammar calls differ by about 1e-8 in the
  # coefficients and 3e-6 (relative) in the precisions.
  expect_equal(coef(bar), coef(asr), tolerance = 1e-6)
  expect_equal(
    bar$inla$summary.hyperpar$mean,
    asr$inla$summary.hyperpar$mean,
    tolerance = 1e-4
  )
  # The peer: lme4's REML estimates of the same model.
  peer <- lme4::fixef(lme4::lmer(log(y) ~ x + (1 | g), data = d))
  expect_equal(unname(coef(bar)), unname(peer), tolerance = 1e-3)
})
