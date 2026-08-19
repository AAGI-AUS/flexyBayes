# =============================================================================
# Arrival: a prior the user asked for reaches the engine that fits the model.
#
# The 0.9.1 field sweep found three independent defects sharing one shape --
# an argument accepted at the door and never delivered. `half_normal(sd =
# 1.5)` fitted under `half_normal(scale = 1)`; `pc(upper = 1)` fitted under a
# tail probability the caller never wrote; `prior_fixed_sd` was documented,
# accepted, and discarded on both engines. Acceptance was tested everywhere
# in the suite. Arrival was tested nowhere.
#
# The oracle here is the engine's own account of what it received: brms's
# generated Stan program (`return_code = TRUE` runs the same
# `brms::make_stancode()` path a fit runs, so a prior row matching no model
# parameter fails here exactly as it would in `brm()`), and the INLA hyper
# specification as it is spliced into the formula INLA is handed. Neither is
# re-read from the flexyBayes object that produced it.
# =============================================================================

.arrival_data <- function(seed = 20260818L) {
  set.seed(seed)
  n <- 60L
  g <- factor(rep(seq_len(10L), each = 6L))
  x <- stats::rnorm(n)
  u <- stats::rnorm(10L, 0, 0.7)
  eta <- 0.5 + 0.4 * x + u[as.integer(g)]
  data.frame(y = eta + stats::rnorm(n), x = x, g = g)
}

# brms's own rendering of the priors it attached, as `lprior += ...` lines.
.brms_prior_lines <- function(...) {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  code <- suppressMessages(flexybayes(
    ...,
    backend = "brms",
    return_code = TRUE
  ))
  lines <- strsplit(as.character(code), "\n")[[1L]]
  trimws(grep("lprior \\+=", lines, value = TRUE))
}

# The INLA hyper specification, read back off the emitted formula rather
# than off the hyper list flexyBayes built.
.inla_emit_formula <- function(...) {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  emitted <- suppressMessages(flexybayes(
    ...,
    backend = "inla",
    return_code = TRUE
  ))
  paste(deparse(emitted$formula), collapse = " ")
}

# --- accepted specifications arrive as written --------------------------- #

test_that("an accepted fb_prior() spec reaches Stan with the values written", {
  skip_if_not_installed("brms")
  d <- .arrival_data()
  lines <- .brms_prior_lines(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    prior = fb_prior(
      sigma ~ pc(upper = 2, prob = 0.05),
      sd(group = "g") ~ half_normal(scale = 1.5),
      b("x") ~ normal(mean = 0, sd = 4)
    )
  )
  # half_normal(scale = 1.5) -> normal(0, 1.5), not the emit-time default 1.
  expect_true(any(grepl("normal_lpdf(sd_1 | 0, 1.5)", lines, fixed = TRUE)))
  # pc(upper = 2, prob = 0.05) -> exponential(-log(0.05) / 2).
  rate <- -log(0.05) / 2
  expect_true(any(grepl(
    paste0("exponential_lpdf(sigma | ", format(rate, digits = 13), ")"),
    lines,
    fixed = TRUE
  )))
  expect_true(any(grepl("normal_lpdf(b[1] | 0, 4)", lines, fixed = TRUE)))
})

test_that("a uniform-on-SD spec reaches Stan with its own bounds", {
  skip_if_not_installed("brms")
  d <- .arrival_data()
  lines <- .brms_prior_lines(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    prior = fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 4))
  )
  expect_true(any(grepl("uniform_lpdf(sd_1 | 0, 4)", lines, fixed = TRUE)))
})

test_that("an accepted fb_prior() spec reaches INLA with the values written", {
  skip_if_not_installed("INLA")
  d <- .arrival_data()
  form <- .inla_emit_formula(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    prior = fb_prior(sd(group = "g") ~ half_normal(scale = 1.5))
  )
  expect_match(form, "hyper = list(prec = list(prior =", fixed = TRUE)
  expect_match(form, "s=1.5", fixed = TRUE)

  form_pc <- .inla_emit_formula(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    prior = fb_prior(sd(group = "g") ~ pc(upper = 2, prob = 0.05))
  )
  expect_match(form_pc, "pc.prec", fixed = TRUE)
  expect_match(form_pc, "param = c(2, 0.05)", fixed = TRUE)
})

# --- the substitution the sweep caught cannot recur ---------------------- #

test_that("no accepted spec can reach an engine under a defaulted value", {
  # The three silent substitutions of FS-2, verbatim. Each is now refused at
  # construction, so none of them can reach an engine at all -- the
  # assertion is on the refusal, because there is no fit left to inspect.
  for (spec in list(
    sd(group = "g") ~ half_normal(sd = 1.5),
    sd(group = "g") ~ half_normal(scale = "1"),
    sd(group = "g") ~ half_normal(scale = 2, extra = 9),
    sigma ~ pc(upper = 1),
    sigma ~ pc(prob = 0.05)
  )) {
    expect_error(fb_prior(spec), class = "flexybayes_refusal")
  }
})

# --- the scalar prior arguments arrive, or are named as not applied ------ #

test_that("prior_fixed_sd reaches Stan when it is supplied", {
  skip_if_not_installed("brms")
  d <- .arrival_data()
  lines <- .brms_prior_lines(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    prior_fixed_sd = 7
  )
  expect_true(any(grepl("normal_lpdf(b | 0, 7)", lines, fixed = TRUE)))
  expect_true(any(grepl("normal_lpdf(Intercept | 0, 7)", lines, fixed = TRUE)))
  # It reaches the fixed effects alongside the auto-default variance
  # prior, which is the combination that used to discard it: passing the
  # scalar alone left `default_prior_active` TRUE and the legacy bridge,
  # the only consumer of the scalar, was never reached.
  expect_true(any(grepl("uniform_lpdf(sd_1 | 0,", lines, fixed = TRUE)))
})

test_that("prior_fixed_sd reaches INLA's control.fixed when it is supplied", {
  skip_if_not_installed("INLA")
  d <- .arrival_data()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  emitted <- suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    backend = "inla",
    prior_fixed_sd = 7,
    return_code = TRUE
  ))
  expect_equal(emitted$control_fixed$prec, 1 / 49)
  expect_equal(emitted$control_fixed$prec.intercept, 1 / 49)
  expect_equal(emitted$control_fixed$mean, 0)
})

test_that("an unsupplied prior_fixed_sd leaves each engine's own default", {
  skip_if_not_installed("INLA")
  d <- .arrival_data()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  emitted <- suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    backend = "inla",
    return_code = TRUE
  ))
  expect_length(emitted$control_fixed, 0L)

  skip_if_not_installed("brms")
  lines <- .brms_prior_lines(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian"
  )
  expect_false(any(grepl("normal_lpdf(b | 0, 100)", lines, fixed = TRUE)))
  expect_true(any(grepl("student_t_lpdf(Intercept |", lines, fixed = TRUE)))
})

test_that("prior_vc_sd reaches INLA as the same density brms receives", {
  skip_if_not_installed("INLA")
  d <- .arrival_data()
  form <- .inla_emit_formula(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    prior_vc_sd = 3
  )
  # lognormal(0, 3) on the SD scale, written in INLA's log-precision
  # parameterisation: log p(theta) is quadratic in theta with s = 3.
  expect_match(form, "expression: s=3;", fixed = TRUE)
  expect_match(form, "(theta*theta)/(8.0*s*s)", fixed = TRUE)

  skip_if_not_installed("brms")
  lines <- .brms_prior_lines(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    prior_vc_sd = 3
  )
  expect_true(any(grepl("lognormal_lpdf(sd_1 | 0, 3)", lines, fixed = TRUE)))
  expect_true(any(grepl("lognormal_lpdf(sigma | 0, 3)", lines, fixed = TRUE)))
})

test_that("the INLA lognormal expression prior integrates to one", {
  # The expression is a density in theta, so the check is the engine's
  # question rather than ours: does the transformed density normalise?
  expr <- flexyBayes:::.inla_lognormal_sd_expr(3)
  s <- 3
  dens <- function(theta) {
    exp(-log(2) - log(s) - 0.5 * log(2 * pi) - (theta^2) / (8 * s^2))
  }
  expect_match(expr, "expression: s=3;", fixed = TRUE)
  expect_equal(
    stats::integrate(dens, -Inf, Inf)$value,
    1,
    tolerance = 1e-6
  )
  # And it is the lognormal on the SD scale it claims to be: the density
  # of log(sigma) implied by theta = -2 log(sigma) is normal(0, s).
  theta <- -2 * log(c(0.5, 1, 2, 4))
  expect_equal(
    dens(theta) * 2,
    stats::dnorm(log(c(0.5, 1, 2, 4)), 0, s),
    tolerance = 1e-10
  )
})

test_that("an unsupplied prior_vc_sd leaves INLA's own hyperprior", {
  skip_if_not_installed("INLA")
  d <- .arrival_data()
  form <- .inla_emit_formula(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian"
  )
  expect_false(grepl("theta*theta", form, fixed = TRUE))
})

# --- prior_summary() tells one story ------------------------------------- #

test_that("prior_summary() reports the legacy scalars as applied on INLA", {
  skip_if_not_installed("INLA")
  d <- .arrival_data()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    backend = "inla",
    prior_vc_sd = 3,
    verbose = FALSE
  ))
  txt <- paste(
    utils::capture.output(print(prior_summary(fit))),
    collapse = "\n"
  )
  expect_match(txt, "prior_vc_sd    = 3", fixed = TRUE)
  # The contradiction: the header named a prior and four lines down the
  # same object said no fb_prior() was supplied and the engine used its
  # own default. Both statements cannot be true of one fit.
  expect_false(grepl("each engine used its own default", txt, fixed = TRUE))
  expect_match(txt, "prior_fixed_sd:  not supplied", fixed = TRUE)

  # The variance components are recorded as carried, not as engine
  # defaults, so triangulate() can compare them.
  expect_setequal(
    names(fit$extras$fingerprint$priors),
    c("sigma", "sd_g")
  )
  expect_length(fit$extras$fingerprint$engine_default_params, 0L)
})

test_that("a declared sigma prior on a family without one is named, not hidden", {
  not_applied <- flexyBayes:::.prior_summary_not_applied(
    priors = fb_prior(sigma ~ uniform(lower = 0, upper = 3)),
    fb_terms = list(family = "gamma"),
    backend_label = "brms"
  )
  expect_named(not_applied, "sigma")
  expect_match(not_applied[["sigma"]], "no residual scale parameter")

  # Gaussian keeps it, and the other engine is not this engine.
  expect_length(
    flexyBayes:::.prior_summary_not_applied(
      priors = fb_prior(sigma ~ uniform(lower = 0, upper = 3)),
      fb_terms = list(family = "gaussian"),
      backend_label = "brms"
    ),
    0L
  )
  expect_length(
    flexyBayes:::.prior_summary_not_applied(
      priors = fb_prior(sigma ~ uniform(lower = 0, upper = 3)),
      fb_terms = list(family = "gamma"),
      backend_label = "inla"
    ),
    0L
  )
})

# --- the named engine defaults are the engines' own ---------------------- #

test_that("the engine-default narrative matches what the engine declares", {
  skip_if_not_installed("INLA")
  defaults <- INLA::inla.set.control.fixed.default()
  expect_equal(defaults$prec, 0.001)
  expect_equal(defaults$prec.intercept, 0)
  txt <- flexyBayes:::.prior_summary_fixed_engine_default("inla")
  expect_match(txt, "prec = 0.001", fixed = TRUE)
  expect_match(txt, "prec.intercept = 0", fixed = TRUE)
})

test_that("the brms engine-default narrative matches its own prior table", {
  skip_if_not_installed("brms")
  d <- .arrival_data()
  tbl <- as.data.frame(brms::default_prior(
    y ~ x + (1 | g),
    data = d,
    family = stats::gaussian()
  ))
  b_rows <- tbl[tbl$class == "b" & !nzchar(tbl$coef), , drop = FALSE]
  int_rows <- tbl[tbl$class == "Intercept", , drop = FALSE]
  expect_true(all(!nzchar(b_rows$prior)))
  expect_true(any(grepl("student_t", int_rows$prior, fixed = TRUE)))
  txt <- flexyBayes:::.prior_summary_fixed_engine_default("brms")
  expect_match(txt, "flat on the population-level coefficients", fixed = TRUE)
  expect_match(txt, "student_t on the intercept", fixed = TRUE)
})
