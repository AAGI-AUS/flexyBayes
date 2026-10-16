# =============================================================================
# Scale-aware default priors on INLA.
#
# Every prior the package chooses for an INLA fit is written in the units of
# the quantity it governs. Every prior the caller states is left as stated
# (R/inla_default_priors.R). Expected values come from the data, from lme4, or
# from the same prior written out by hand and stated, a route that never
# passes through the code that fills or scales a default. Tolerances come from
# the spread of repeated identical fits, measured on INLA 26.8.7 and quoted
# where each is used.
# =============================================================================

# Twenty identical fits per cell, INLA 26.8.7, the numbers each tolerance
# below is built from.
#
#   15/sim block SD. Every route, with the default filled, attached or stated:
#   no error and no blow-up in 180 fits. Medians 2.946 to 2.955. The largest
#   range within one cell 0.062. Tolerance on a median: 0.2, about three
#   times that range.
#
#   sleepstudy `(Dx || Subject)`, default priors, coef(Dx) * c. Range over
#   20 calls at most 1e-13 at each c. Distance from lme4's REML 10.4672859596
#   at most 7.8e-12, where lme4 itself agrees with its own fit on rescaled
#   Days only to 1e-12. The default this replaces had a run-to-run range of
#   1.5e-10 at c = 1, and the fixed-effect multiplier was chosen against that
#   range (R/inla_default_priors.R). Tolerance: 1e-9, one order above it,
#   for engine output that differs between platforms in its last digits. A
#   multiplier of 1000 or less (7.8e-8 from lme4) fails it.
#
#   Slope SD * c. Range at most 6e-4 at each c. The means at the three scales
#   within 2.7e-3 of each other. Tolerance across scales: 0.01. Against lme4
#   the comparison is looser, because a posterior mean of a standard
#   deviation is not a REML estimate: at c = 1, where no bound binds, the two
#   stand 7.7 per cent apart. Tolerance: 10 per cent.
#
#   Poisson and Bernoulli random-intercept fits (30 groups of 8, x standard
#   normal), coefficient times multiplier, three calls per cell. Range within a
#   cell at most 1.1e-7. Between x / 1000 and x: 2.0e-7 relative (Poisson)
#   and 9.6e-6 (Bernoulli). Tolerance: 1e-4, ten times the larger. Before
#   the change the two differed by 59 and 96 per cent.
.sad_tol <- list(
  block_sd_median = 0.2,
  days_lme4 = 1e-9,
  slope_across_scales = 0.01,
  slope_reml_relative = 0.10,
  o1_across_scales = 1e-4
)

.sad_blocks <- function() {
  set.seed(3115)
  d <- expand.grid(
    rep = seq_len(20L),
    trt = factor(seq_len(5L)),
    blk = factor(seq_len(4L))
  )
  u_b <- stats::rnorm(4L, 0, 2)
  d$y <- 10 +
    0.5 * as.numeric(d$trt) +
    u_b[as.integer(d$blk)] +
    stats::rnorm(nrow(d), 0, 1)
  d
}

.sad_sleepstudy <- function(cc = 1) {
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  d <- e$sleepstudy
  d$Dx <- d$Days * cc
  d
}

.sad_component <- function(fit, component) {
  v <- summary(fit)$varcomp
  v$estimate[v$component == component][[1L]]
}

.sad_quiet <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

# A formula as one line with single spaces, so a match does not depend on
# where deparse() wraps it.
.sad_formula_text <- function(form) {
  gsub(
    "[[:space:]]+",
    " ",
    paste(deparse(form, width.cutoff = 500L), collapse = " ")
  )
}


# --- (a) an unstated random-effect SD ------------------------------------ #

test_that("an unstated block SD takes the bounded default on the per-row route", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  d <- .sad_blocks()
  sd_y <- stats::sd(d$y)
  upper <- 5 * sd_y

  # The reference: the package default written out and stated, so it reaches
  # INLA without being filled in.
  stated <- .sad_quiet(flexybayes(
    fixed = y ~ trt,
    random = ~blk,
    data = d,
    backend = "inla",
    prior = fb_prior(
      sigma ~ uniform(lower = 0, upper = upper),
      sd(group = "blk") ~ uniform(lower = 0, upper = upper)
    ),
    verbose = FALSE
  ))
  reference <- .sad_component(stated, "sd_blk")

  # The aggregated and streamed Gaussian routes now refuse by name (they
  # cannot carry the residual prior. Test-aggregation-gaussian-containment.R),
  # so the per-row route is the one a Gaussian model takes.
  sigma_only <- fb_prior(sigma ~ uniform(lower = 0, upper = upper))
  routes <- list(
    per_row = function() {
      flexybayes(
        fixed = y ~ trt,
        random = ~blk,
        data = d,
        backend = "inla",
        aggregate = FALSE,
        prior = sigma_only,
        verbose = FALSE
      )
    }
  )

  # Before the default reached these routes, twenty identical fits of this
  # design gave, per route: per-row 20 of 20 with the block SD above 1e6
  # (every one passing the numerical-confirm gate). Aggregated 10 above 1e6
  # and 5 engine errors. Streamed 10 above 1e6 and 4 errors. The fits that
  # survived sat at 1.49 against the default's 2.95. Counting only
  # blow-ups and errors, one fit fails with probability 1.0, 0.75 and 0.70,
  # so five pass with probability at most 0.25^5 = 0.001 on the aggregated
  # route and 0.30^5 = 0.0024 on the streamed one: each route below detects
  # the defect with probability above 0.99.
  for (route in names(routes)) {
    sds <- vapply(
      seq_len(5L),
      function(i) .sad_component(.sad_quiet(routes[[route]]()), "sd_blk"),
      numeric(1L)
    )
    expect_true(all(sds < 10 * sd_y), label = paste(route, "block SD bounded"))
    # Three times the largest range measured over twenty identical fits.
    expect_lt(
      abs(stats::median(sds) - reference),
      .sad_tol$block_sd_median,
      label = paste(route, "median block SD against the stated default")
    )
  }
})


# --- (b) and (c) a covariate's units --------------------------------------- #

test_that("a slope SD and its coefficient do not move with the covariate's units", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  reml <- lme4::lmer(
    Reaction ~ Days + (Days || Subject),
    data = .sad_sleepstudy(),
    REML = TRUE
  )
  reml_days <- unname(lme4::fixef(reml)[["Days"]])
  vc <- as.data.frame(lme4::VarCorr(reml))
  reml_slope_sd <- vc$sdcor[vc$var1 %in% "Days" & is.na(vc$var2)]

  cells <- lapply(c(0.001, 1, 1000), function(cc) {
    fit <- .sad_quiet(flexybayes(
      Reaction ~ Dx + (Dx || Subject),
      data = .sad_sleepstudy(cc),
      backend = "inla",
      verbose = FALSE
    ))
    # A variance-component table without the slope row is one of the ways a
    # scale-blind default shows itself, so its absence is a named failure
    # instead of an out-of-range read.
    v <- summary(fit)$varcomp
    has_slope <- any(v$component == "sd_Dx_Subject")
    expect_true(
      has_slope,
      label = paste0("c = ", cc, ": sd_Dx_Subject row in summary()$varcomp")
    )
    c(
      c = cc,
      days = stats::coef(fit)[["Dx"]] * cc,
      slope_sd = if (has_slope) {
        v$estimate[v$component == "sd_Dx_Subject"][[1L]] * cc
      } else {
        NA_real_
      }
    )
  })
  cells <- as.data.frame(do.call(rbind, cells))

  # Before: at c = 0.001 the slope SD bound 5 sd(y) = 282 sat just above a
  # fitted 279 (0.279 per day against 6.45), and the Days coefficient was
  # 8.64, which a wide stated fixed-effect prior alone restored to 10.4673.
  for (i in seq_len(nrow(cells))) {
    lab <- paste0("c = ", cells$c[[i]])
    expect_lt(
      abs(cells$days[[i]] - reml_days),
      .sad_tol$days_lme4,
      label = paste(lab, "Days against lme4")
    )
    expect_lt(
      abs(cells$slope_sd[[i]] - reml_slope_sd) / reml_slope_sd,
      .sad_tol$slope_reml_relative,
      label = paste(lab, "slope SD against lme4")
    )
  }
  # Against the fit on the covariate's own units. The posterior mean of a
  # standard deviation is not a REML estimate, at c = 1, where no bound
  # binds, the two stand 7.7 per cent apart, so the tight comparison is this
  # one.
  expect_lt(diff(range(cells$slope_sd)), .sad_tol$slope_across_scales)
  expect_lt(diff(range(cells$days)), .sad_tol$days_lme4)
})

test_that("a non-Gaussian coefficient does not move when its covariate shrinks", {
  # Poisson and Bernoulli random-intercept fits, 30 groups of 8, x standard
  # normal. Before, with x divided by 1000, the coefficient times the
  # multiplier fell to 41 per cent of its value at x on Poisson and to 4.4 per
  # cent on Bernoulli. At x * 1000 the coefficient moves by about 1 per cent
  # (Poisson) and 3 per cent (Bernoulli) whatever the prior, raw INLA does
  # the same under its own defaults and under strategy = "laplace", so that
  # scale is an engine limit, not a prior, and is not asserted here.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  count_design <- function(family) {
    n_g <- 30L
    n_per <- 8L
    set.seed(if (identical(family, "poisson")) 3201 else 3202)
    g <- factor(rep(sprintf("g%02d", seq_len(n_g)), each = n_per))
    x <- stats::rnorm(n_g * n_per)
    u <- stats::rnorm(n_g, 0, if (identical(family, "poisson")) 0.5 else 1)
    y <- if (identical(family, "poisson")) {
      stats::rpois(n_g * n_per, exp(1 + 0.5 * x + u[as.integer(g)]))
    } else {
      stats::rbinom(n_g * n_per, 1L, stats::plogis(0.5 * x + u[as.integer(g)]))
    }
    data.frame(y = y, x = x, g = g)
  }

  for (family in c("poisson", "binomial")) {
    base <- count_design(family)
    scaled <- vapply(
      c(0.001, 1),
      function(cc) {
        d <- base
        d$x <- cc * d$x
        fit <- .sad_quiet(flexybayes(
          y ~ x + (1 | g),
          data = d,
          family = family,
          backend = "inla",
          verbose = FALSE
        ))
        stats::coef(fit)[["x"]] * cc
      },
      numeric(1L)
    )
    expect_lt(
      abs(scaled[[1L]] - scaled[[2L]]) / abs(scaled[[2L]]),
      .sad_tol$o1_across_scales,
      label = family
    )
  }
})


# --- what is stated stays stated, what is filled is recorded -------------- #

test_that("a partial prior keeps what was stated and records what was filled", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  d <- .sad_sleepstudy()
  s_y <- stats::sd(d$Reaction)
  prior <- fb_prior(sd(group = "Subject") ~ half_normal(scale = 20))

  emitted <- .sad_quiet(flexybayes(
    Reaction ~ Days + (Days || Subject),
    data = d,
    backend = "inla",
    prior = prior,
    return_code = TRUE
  ))
  # Stated: mapped by the response scale and nothing else.
  expect_equal(emitted$hyper$Subject$meta$family, "half_normal")
  expect_equal(emitted$hyper$Subject$meta$scale, 20 / s_y)
  expect_null(emitted$hyper$Subject$meta$package_chosen)
  # Filled: 5 sd(y) / divisor, divided by sd(Days), then by the response scale.
  expect_identical(emitted$hyper$Days_Subject$meta$family, "half_normal")
  expect_equal(
    emitted$hyper$Days_Subject$meta$scale,
    5 / flexyBayes:::.FB_DEFAULT_SCALE_DIVISOR / stats::sd(d$Days)
  )
  expect_true(emitted$hyper$Days_Subject$meta$package_chosen)

  fit <- .sad_quiet(flexybayes(
    Reaction ~ Days + (Days || Subject),
    data = d,
    backend = "inla",
    prior = prior,
    verbose = FALSE
  ))
  # The declaration the fit carries is the caller's, unchanged.
  expect_identical(fit$extras$fb_terms$priors, prior)

  ps <- prior_summary(fit)
  expect_named(ps$filled_priors, "sd_Days_Subject")
  expect_equal(ps$covariate_sd[["Days_Subject"]], stats::sd(d$Days))
  txt <- paste(utils::capture.output(print(ps)), collapse = "\n")
  expect_match(txt, "Filled with the package default", fixed = TRUE)

  fp <- fit$extras$fingerprint
  expect_identical(fp$prior_sources[["sd_Subject"]], "explicit")
  expect_identical(fp$prior_sources[["sd_Days_Subject"]], "filled_default")
  expect_false("sd_Days_Subject" %in% names(fp$engine_default_params))

  v <- summary(fit)$varcomp
  expect_match(
    v$prior[v$component == "sd_Days_Subject"],
    "^half_normal\\(scale="
  )
})

test_that("a stated fixed-effect prior is not rescaled by its covariate", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  d <- .sad_sleepstudy(0.001)
  s_y <- stats::sd(d$Reaction)

  # A blanket scalar: sd = 7 in response units per covariate unit, on every
  # coefficient, whatever the covariate's units. On Dx, recorded in thousandths
  # of a day, it binds. That is what was asked for.
  blanket <- .sad_quiet(flexybayes(
    Reaction ~ Dx + (1 | Subject),
    data = d,
    backend = "inla",
    prior_fixed_sd = 7,
    return_code = TRUE
  ))
  expect_equal(blanket$control_fixed$prec, s_y^2 / 49)
  expect_false(is.list(blanket$control_fixed$prec))

  # A b() row keeps its own value. Only the unnamed numeric column is scaled.
  set.seed(11L)
  d$z <- stats::rnorm(nrow(d), 0, 40)
  rowwise <- .sad_quiet(flexybayes(
    Reaction ~ Dx + z + (1 | Subject),
    data = d,
    backend = "inla",
    prior = fb_prior(b("Dx") ~ normal(mean = 0, sd = 4)),
    return_code = TRUE
  ))
  expect_equal(rowwise$control_fixed$prec$Dx, s_y^2 / 16)
  # normal(0, 1e5 sd(y) / sd(z)), restated on the standardised response.
  expect_equal(rowwise$control_fixed$prec$z, stats::sd(d$z)^2 / 1e10)
})

test_that("a factor contrast keeps INLA's own fixed-effect default", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  emitted <- .sad_quiet(flexybayes(
    fixed = y ~ trt,
    random = ~blk,
    data = .sad_blocks(),
    backend = "inla",
    return_code = TRUE
  ))
  expect_length(emitted$control_fixed, 0L)
})


# --- the streamed route, the multi-way interaction, the blocks carrier ---- #

# [removed] "the streamed default bound equals the per-row bound on the same
# data", a streamed Gaussian fit now refuses by name, because the
# aggregated Gaussian route cannot carry the residual prior, and the streamed
# response-scale reading this compared went with it. The comparison returns
# with that composition.

.sad_three_way <- function() {
  set.seed(20260913L)
  d <- expand.grid(
    a = factor(1:3),
    b = factor(1:3),
    cc = factor(1:2),
    rep = 1:4
  )
  d$y <- 50 +
    stats::rnorm(18L, 0, 4)[as.integer(interaction(d$a, d$b, d$cc))] +
    stats::rnorm(nrow(d), 0, 8)
  d
}

test_that("a multi-way random interaction takes the half-normal default", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .sad_three_way()

  emitted <- .sad_quiet(flexybayes(
    fixed = y ~ 1,
    random = ~ a:b:cc,
    data = d,
    backend = "inla",
    return_code = TRUE
  ))
  # half_normal(scale = 5 sd(y) / divisor), restated on the standardised
  # response.
  expect_identical(emitted$hyper[["a:b:cc"]]$meta$family, "half_normal")
  expect_equal(
    emitted$hyper[["a:b:cc"]]$meta$scale,
    5 / flexyBayes:::.FB_DEFAULT_SCALE_DIVISOR
  )
  expect_match(
    .sad_formula_text(emitted$formula),
    "f(a_b_cc_id, model = \"iid\", hyper",
    fixed = TRUE
  )
})

test_that("a stated prior on a multi-way random interaction reaches INLA", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .sad_three_way()

  # Before, the reachability check did not list the colon-joined key of a
  # multi-way interaction, and this call was refused with
  # prior_target_not_in_model although the default fills the same key.
  emitted <- .sad_quiet(flexybayes(
    fixed = y ~ 1,
    random = ~ a:b:cc,
    data = d,
    backend = "inla",
    prior = fb_prior(sd(group = "a:b:cc") ~ uniform(lower = 0, upper = 20)),
    return_code = TRUE
  ))
  # The caller's bound restated on the standardised response, and nothing
  # else: not the default's 5 sd(y), and not marked as the package's choice.
  expected <- 20 / stats::sd(d$y)
  entry <- emitted$hyper[["a:b:cc"]]
  expect_equal(entry$meta$upper, expected, tolerance = 1e-12)
  expect_null(entry$meta$package_chosen)
  # The same bound in the f() term INLA reads.
  form <- .sad_formula_text(emitted$formula)
  pattern <- paste0(
    "^.*f\\(a_b_cc_id, model = \"iid\", hyper = list\\(prec = list\\(",
    "prior = \"expression: U=([^;]+);.*$"
  )
  expect_match(form, pattern)
  expect_equal(
    as.numeric(sub(pattern, "\\1", form)),
    expected,
    tolerance = 1e-12
  )

  # A key this model does not carry is still refused: an interaction naming a
  # factor the model does not have, and a two-way interaction of its factors.
  refused <- list(
    fb_prior(sd(group = "a:b:zz") ~ uniform(lower = 0, upper = 20)),
    fb_prior(sd(group = "a:b") ~ uniform(lower = 0, upper = 20))
  )
  for (prior in refused) {
    err <- tryCatch(
      .sad_quiet(flexybayes(
        fixed = y ~ 1,
        random = ~ a:b:cc,
        data = d,
        backend = "inla",
        prior = prior,
        return_code = TRUE
      )),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_prior_target_absent")
    expect_identical(err$reason_code, "prior_target_not_in_model")
  }
})

test_that("each block of a blocks-format vm() takes the default under a partial prior", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("Matrix")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(5L)
  d <- data.frame(
    geno = factor(rep(seq_len(5L), length.out = 60L)),
    yield = stats::rnorm(60L, 50, 5)
  )
  blocks <- list(diag(2), diag(3) + 0.1)
  # The gate refuses a blocks carrier on INLA (test-blocks-inla-varcomp.R).
  # The emitter branch it keeps for a forced override is exercised directly.
  fb <- fb_from_asreml(
    fixed = yield ~ 1,
    random = ~ vm(geno, cov = fb_cov(blocks, type = "blocks")),
    data = d,
    family = "gaussian",
    known_matrices = list(blocks = blocks),
    prior = fb_prior(sigma ~ uniform(lower = 0, upper = 30))
  )
  emitted <- .sad_quiet(emit_inla(
    fb,
    data = d,
    known_matrices = list(blocks = blocks),
    return_code = TRUE,
    verbose = FALSE
  ))
  expect_identical(emitted$hyper$geno$meta$family, "half_normal")
  expect_equal(
    emitted$hyper$geno$meta$scale,
    5 / flexyBayes:::.FB_DEFAULT_SCALE_DIVISOR
  )
  form <- .sad_formula_text(emitted$formula)
  for (k in 1:2) {
    expect_match(
      form,
      paste0(
        "f(geno_id_block_",
        k,
        ", model = \"generic0\", Cmatrix = blocks_Q_",
        k,
        ", hyper"
      ),
      fixed = TRUE
    )
  }
})

test_that("brms receives the declaration and its own defaults, unchanged", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .sad_blocks()
  code <- .sad_quiet(flexybayes(
    fixed = y ~ trt,
    random = ~blk,
    data = d,
    backend = "brms",
    prior = fb_prior(sigma ~ uniform(lower = 0, upper = 9)),
    return_code = TRUE
  ))
  lines <- strsplit(paste(code, collapse = "\n"), "\n", fixed = TRUE)[[1L]]
  sd_lines <- grep("lprior \\+= .*sd_1", lines, value = TRUE)
  expect_gt(length(sd_lines), 0L)
  expect_false(any(grepl("uniform_lpdf", sd_lines, fixed = TRUE)))
})
