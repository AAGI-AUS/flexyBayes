# Ordinal families: family = "cumulative" and "cratio" (R/ordinal.R,
# R/ordinal_methods.R, R/fb_kfold.R).
#
# The oracles are independent of the package's own code: brms's density
# functions for the category probabilities, a hand-written maximum
# likelihood fit for the coefficients and thresholds, INLA's own declared
# hyperparameter transforms, and Monte Carlo for the quadrature.

.ord_quiet <- function(expr) {
  withr::with_options(
    list(
      flexyBayes.silence_default_prior_note = TRUE,
      flexyBayes.silence_weak_identification_warning = TRUE,
      flexyBayes.silence_spatial_collapse_warning = TRUE
    ),
    expr
  )
}

# One fit per family with the field, shared by the tests below.
.ord_cache <- new.env(parent = emptyenv())
.ord_field_fit <- function(family) {
  if (is.null(.ord_cache[[family]])) {
    d <- .ord_sim_trial(if (identical(family, "cumulative")) 5L else 6L, family)
    .ord_cache[[family]] <- .ord_quiet(flexybayes(
      stage ~ depth,
      random = ~ rep + rep:depth + ar1(row):ar1(col) + plot,
      data = d,
      family = family,
      backend = "inla",
      verbose = FALSE
    ))
  }
  .ord_cache[[family]]
}


# Category probabilities ------------------------------------------------------

test_that("category probabilities equal brms's own density functions", {
  skip_if_not_installed("brms")
  set.seed(1)
  eta <- rnorm(7)
  tau <- sort(rnorm(4))
  links <- list(cumulative = c("logit", "probit"),
    cratio = c("logit", "probit", "cloglog"))
  dens <- list(
    cumulative = get("dcumulative", envir = asNamespace("brms")),
    cratio = get("dcratio", envir = asNamespace("brms"))
  )
  for (fam in names(links)) {
    for (lk in links[[fam]]) {
      ours <- .fb_ordinal_probs(eta, tau, fam, lk)
      theirs <- dens[[fam]](
        seq_len(5L),
        eta = eta,
        thres = matrix(tau, nrow = 7L, ncol = 4L, byrow = TRUE),
        link = lk
      )
      expect_equal(unname(ours), unname(theirs), tolerance = 1e-10,
        info = paste(fam, lk))
      expect_equal(rowSums(ours), rep(1, 7L), tolerance = 1e-12)
    }
  }
})

test_that("category probabilities equal a written-out reference", {
  set.seed(2)
  eta <- rnorm(9, sd = 2)
  tau <- c(-1.4, -0.2, 0.9)
  for (fam in c("cumulative", "cratio")) {
    expect_equal(
      .fb_ordinal_probs(eta, tau, fam, "logit"),
      .ord_probs_reference(eta, tau, fam),
      tolerance = 1e-12
    )
  }
})

test_that("Gauss-Hermite integration over a random effect matches Monte Carlo", {
  rule <- .fb_gauss_hermite()
  # The rule integrates x^2 exp(-x^2) exactly: sqrt(pi) / 2.
  expect_equal(sum(rule$weights * rule$nodes^2), sqrt(pi) / 2,
    tolerance = 1e-12)
  set.seed(3)
  tau <- c(-1, 0.3, 1.6)
  u <- rnorm(1e6, sd = sqrt(0.8))
  for (fam in c("cumulative", "cratio")) {
    gh <- .fb_ordinal_probs_integrated(0.4, tau, 0.8, fam, "logit")
    mc <- colMeans(.fb_ordinal_probs(0.4 + u, tau, fam, "logit"))
    expect_lt(max(abs(gh[1L, ] - mc)), 2e-3)
  }
  expect_equal(
    .fb_ordinal_probs_integrated(0.4, tau, 0, "cratio", "logit"),
    .fb_ordinal_probs(0.4, tau, "cratio", "logit")
  )
})


# Response record and coding ---------------------------------------------------

test_that("the categories are the levels observed, and an unused level is dropped", {
  d <- data.frame(
    y = factor(c("b", "d", "a", "b"), levels = c("a", "b", "c", "d"),
      ordered = TRUE),
    x = 1:4
  )
  fb <- list(family = "cumulative", response = "y")
  expect_message(rec <- .fb_ordinal_record(fb, d), "no observation at level c")
  expect_identical(rec$levels, c("a", "b", "d"))
  expect_identical(rec$K, 3L)
  fb$ordinal <- rec
  expect_identical(.fb_ordinal_codes(fb, d), c(2L, 3L, 1L, 2L))
  d$y <- c(12, 21, 12, 30)
  rec <- .fb_ordinal_record(fb, d)
  expect_identical(rec$levels, c("12", "21", "30"))
  fb$ordinal <- rec
  d$y[2L] <- NA
  expect_identical(.fb_ordinal_codes(fb, d), c(1L, NA, 1L, 3L))
  d$y[2L] <- 25
  expect_error(.fb_ordinal_codes(fb, d),
    class = "flexybayes_refusal_ordinal_response_category_unknown")
})

test_that("a response that does not state an order is refused by name", {
  fb <- list(family = "cratio", response = "y")
  expect_error(
    .fb_ordinal_record(fb, data.frame(y = factor(c("a", "b")))),
    class = "flexybayes_refusal_ordinal_response_unordered"
  )
  expect_error(
    .fb_ordinal_record(fb, data.frame(y = c("a", "b"))),
    class = "flexybayes_refusal_ordinal_response_unordered"
  )
  expect_error(
    .fb_ordinal_record(fb, data.frame(y = c(1.5, 2))),
    class = "flexybayes_refusal_ordinal_response_not_integer"
  )
  expect_error(
    .fb_ordinal_record(fb, data.frame(y = c(3L, 3L, NA))),
    class = "flexybayes_refusal_ordinal_response_too_few_categories"
  )
})

test_that("the continuation-ratio expansion writes one row per boundary reached", {
  d <- data.frame(
    y = factor(c("a", "c", "d", NA, "b"), levels = c("a", "b", "c", "d"),
      ordered = TRUE),
    x = 1:5
  )
  fb <- list(
    family = "cratio",
    response = "y",
    intercept = TRUE,
    fixed_terms = list(list(type = "continuous", var = "x", label = "x"))
  )
  setup <- .fb_ordinal_inla_setup(fb, d)
  # min(y, K - 1) rows per observed response, one for the missing one.
  expect_identical(nrow(setup$data), 1L + 3L + 3L + 1L + 2L)
  expect_identical(
    setup$data[[.FB_ORDINAL_RESPONSE]],
    c(0L, 1L, 1L, 0L, 1L, 1L, 1L, NA, 1L, 0L)
  )
  expect_identical(as.integer(setup$data[[.FB_ORDINAL_CUT]]),
    c(1L, 1L, 2L, 3L, 1L, 2L, 3L, 1L, 1L, 2L))
  expect_identical(setup$user_rows, c(1L, 2L, 5L, 8L, 9L))
  expect_false(setup$fb$intercept)
  expect_identical(setup$fb$fixed_terms[[1L]]$var, .FB_ORDINAL_CUT)
  # The expansion's Bernoulli likelihood is the category likelihood.
  tau <- c(-0.5, 0.4, 1.1)
  eta <- 0.3 * d$x
  ex <- setup$data[!is.na(setup$data[[.FB_ORDINAL_RESPONSE]]), ]
  q <- stats::plogis(0.3 * ex$x - tau[as.integer(ex[[.FB_ORDINAL_CUT]])])
  bernoulli <- sum(log(ifelse(ex[[.FB_ORDINAL_RESPONSE]] == 1L, q, 1 - q)))
  y <- as.integer(d$y)
  keep <- !is.na(y)
  categorical <- sum(log(
    .fb_ordinal_probs(eta[keep], tau, "cratio", "logit")[
      cbind(seq_len(sum(keep)), y[keep])]
  ))
  expect_equal(bernoulli, categorical, tolerance = 1e-12)
})

test_that("INLA's cumulative thresholds read through the engine's declared transforms", {
  skip_if_not_installed("INLA")
  declared <- INLA::inla.models()$likelihood$pom$hyper
  theta <- cbind(-1, matrix(c(0.2, -0.3, 0.1), nrow = 1L))
  colnames(theta) <- paste0("theta", 1:4, " for POM")
  user <- theta
  for (k in 1:4) {
    user[, k] <- declared[[paste0("theta", k)]]$from.theta(theta[, k])
  }
  expected <- cumsum(c(-1, exp(c(0.2, -0.3, 0.1))))
  expect_equal(c(.fb_ordinal_pom_alpha(theta, 5L)), expected)
  expect_equal(c(.fb_ordinal_pom_alpha(user, 5L, scale = "user")), expected)
  expect_identical(.fb_ordinal_inla_max_categories(), length(declared) + 1L)
})

test_that("the ordinal refusal codes are registered", {
  codes <- c(
    "ordinal_response_unordered", "ordinal_response_not_integer",
    "ordinal_response_too_few_categories", "ordinal_response_category_unknown",
    "ordinal_fixed_without_intercept", "ordinal_offset_not_supported",
    "ordinal_method_not_supported", "ordinal_draws_unavailable",
    "kfold_family_not_supported", "kfold_folds_invalid",
    "kfold_response_not_a_column", "kfold_category_missing_in_training"
  )
  for (code in codes) {
    expect_false(is.null(.lookup_refusal(code)), info = code)
  }
})


# Fits on INLA -----------------------------------------------------------------

test_that("both families fit the split-plot model with the AR1 x AR1 field on INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  for (fam in c("cumulative", "cratio")) {
    fit <- .ord_field_fit(fam)
    expect_s3_class(fit, c("flexybayes_ordinal", "flexybayes_inla"))
    expect_identical(nobs(fit), 900L)
    expect_identical(family(fit)$family, fam)
    expect_identical(names(coef(fit)), c("depth80", "depth120"))
    tau <- coef(fit, what = "thresholds")
    expect_identical(names(tau), paste0("threshold[Z", 1:5, "|Z", 2:6, "]"))
    expect_true(all(diff(tau) > 0), info = fam)
    vc <- summary(fit)$varcomp$component
    expect_setequal(vc, c("sd_rep", "sd_rep:depth", "sd_spatial", "rho_row",
      "rho_col", "sd_plot"))
    expect_identical(dim(vcov(fit)), c(2L, 2L))
    expect_identical(rownames(confint(fit)), c("depth80", "depth120"))
    td <- tidy(fit)
    expect_identical(td$term, c("depth80", "depth120", names(tau)))
    expect_identical(colnames(model.matrix(fit)), c("depth80", "depth120"))
    if (identical(fam, "cumulative")) {
      # The threshold prior, as INLA reports it back, is the one recorded.
      th1 <- fit$inla$all.hyper$family[[1L]]$hyper$theta1
      expect_identical(as.character(th1$prior), "dirichlet")
      expect_identical(as.numeric(th1$param), 1)
      expect_match(fit$extras$ordinal$threshold_prior, "Dirichlet(1)",
        fixed = TRUE)
      # prior_summary() lists the slot as the package's, read back as sent.
      ps <- prior_summary(fit)
      slot <- ps$slots[ps$slots$engine_slot == "family INLA.Data1 theta1", ]
      expect_identical(nrow(slot), 1L)
      expect_identical(slot$source, "package")
      expect_identical(slot$input_id, "threshold:INLA.Data1")
      expect_identical(slot$readback, "dirichlet(1)")
      expect_true(slot$readback_agrees)
      expect_false(any(ps$slots$source == "unexpected"))
      inp <- ps$inputs[ps$inputs$input_id == "threshold:INLA.Data1", ]
      expect_identical(inp$source, "structure_default")
      expect_identical(inp$n_slots, 1L)
      expect_match(ps$fixed_sd_engine_default, "intercept is pinned at zero",
        fixed = TRUE)
    } else {
      expect_match(prior_summary(fit)$fixed_sd_engine_default,
        "on each threshold row", fixed = TRUE)
    }
  }
})

test_that("threshold draws agree with the threshold table", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("posterior")
  skip_on_cran()
  for (fam in c("cumulative", "cratio")) {
    fit <- .ord_field_fit(fam)
    tab <- fit$extras$ordinal$thresholds
    dr <- posterior::as_draws_df(fit, n_draws = 2000L)
    expect_false(any(c("(Intercept)", "fb_ordinal_cut1") %in% names(dr)))
    m <- colMeans(as.data.frame(dr)[tab$term])
    # INLA's posterior sampler draws the hyperparameters from its
    # integration design, the table from their marginals: the two agree to
    # a tenth of a posterior standard deviation.
    expect_lt(max(abs(m - tab$estimate) / tab$std.error), 0.25, label = fam)
  }
})

test_that("fitted() and predict() return category probabilities", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  for (fam in c("cumulative", "cratio")) {
    fit <- .ord_field_fit(fam)
    p <- fitted(fit, n_draws = 200L, seed = 1L)
    expect_identical(dim(p), c(900L, 6L))
    expect_identical(colnames(p), paste0("Z", 1:6))
    expect_equal(rowSums(p), rep(1, 900L), tolerance = 1e-10)
    # Rows of one plot share its random effects.
    expect_equal(p[1L, ], p[10L, ])
    expect_identical(p, fitted(fit, n_draws = 200L, seed = 1L))
    nd <- data.frame(depth = factor(c(30, 80, 120)))
    zero <- predict(fit, newdata = nd, n_draws = 300L, seed = 2L)
    marg <- predict(fit, newdata = nd, marginal = TRUE, n_draws = 300L,
      seed = 2L)
    expect_equal(rowSums(marg), rep(1, 3L), tolerance = 1e-10)
    # Integrating over the random effects moves mass to the extremes.
    expect_true(all(marg[, 1L] + marg[, 6L] > zero[, 1L] + zero[, 6L]))
    # Deeper sowing raises the stage, the sign of the simulated effect.
    expect_gt(sum(zero[3L, 4:6]), sum(zero[1L, 4:6]))
    draws <- predict(fit, newdata = nd, n_draws = 50L, seed = 2L,
      summary = FALSE)
    expect_identical(dim(draws), c(50L, 3L, 6L))
  }
})

test_that("predict() at zero random effects follows the coefficients and thresholds", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .ord_field_fit("cumulative")
  dr <- .fb_ordinal_draws(fit, 400L, 3L)
  p <- predict(fit, newdata = data.frame(depth = factor(120)), n_draws = 400L,
    seed = 3L)
  by_hand <- 0
  for (s in seq_len(400L)) {
    by_hand <- by_hand +
      .ord_probs_reference(dr$beta[s, "depth120"], dr$tau[s, ], "cumulative")
  }
  expect_equal(unname(p), unname(by_hand / 400L), tolerance = 1e-10)
})

test_that("methods that read a mean response refuse by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .ord_field_fit("cratio")
  cls <- "flexybayes_refusal_ordinal_method_not_supported"
  expect_error(residuals(fit), class = cls)
  expect_error(augment(fit), class = cls)
  expect_error(plot(fit, type = "residuals"), class = cls)
  expect_error(plot(fit, type = "variogram"), class = cls)
  expect_error(recover_data.flexybayes_ordinal(fit), class = cls)
  expect_error(emm_basis.flexybayes_ordinal(fit), class = cls)
  expect_error(get_predict.flexybayes_ordinal(fit), class = cls)
  for (arg in c("type", "se.fit", "classify", "allow_new_levels")) {
    args <- stats::setNames(list(fit, TRUE), c("object", arg))
    expect_error(do.call(predict, args),
      class = "flexybayes_refusal_method_argument_not_supported", info = arg)
  }
})

test_that("an ordinal model the routes do not fit is refused by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  d <- .ord_sim_trial(9L, "cumulative", field = FALSE)
  fit_inla <- function(...) {
    .ord_quiet(flexybayes(
      ...,
      data = d,
      family = "cumulative",
      backend = "inla",
      verbose = FALSE
    ))
  }
  expect_error(fit_inla(stage ~ 0 + depth, random = ~rep),
    class = "flexybayes_refusal_ordinal_fixed_without_intercept")
  d$off <- 0.1
  expect_error(fit_inla(stage ~ depth + offset(off), random = ~rep),
    class = "flexybayes_refusal_ordinal_offset_not_supported")
  expect_error(fit_inla(stage ~ depth, random = ~rep, link = "probit"),
    class = "flexybayes_refusal_inla_gate_refused")
  expect_error(
    fit_inla(stage ~ depth, random = ~rep, residual = ~ dsum(~ units | rep)),
    class = "flexybayes_refusal_heterogeneous_residual_family_has_no_sigma"
  )
  expect_error(
    fit_inla(stage ~ depth, random = ~rep, residual = ~ ar1(row):ar1(col)),
    class = "flexybayes_refusal_ar1_residual_not_representable"
  )
  set.seed(10)
  expect_error(fit_inla(stage ~ depth, random = ~rep, weights = runif(nrow(d))),
    class = "flexybayes_refusal_weights_requires_gaussian")
  d$many <- rep(seq_len(12L), length.out = nrow(d))
  expect_error(
    .ord_quiet(flexybayes(many ~ depth, random = ~rep, data = d,
      family = "cumulative", backend = "inla", verbose = FALSE)),
    class = "flexybayes_refusal_inla_gate_refused"
  )
})

test_that("a missing response is carried, and fitted() gives its probabilities", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  d <- .ord_sim_trial(13L, "cratio")
  d$stage[c(4L, 400L)] <- NA
  fit <- .ord_quiet(flexybayes(
    stage ~ depth,
    random = ~ rep + rep:depth + ar1(row):ar1(col) + plot,
    data = d,
    family = "cratio",
    backend = "inla",
    verbose = FALSE
  ))
  p <- fitted(fit, n_draws = 100L, seed = 1L)
  expect_identical(nrow(p), 900L)
  expect_false(anyNA(p[c(4L, 400L), ]))
  # A missing plant shares its plot's random effects with the plants scored.
  expect_equal(p[4L, ], p[3L, ])
})


# Coefficients and thresholds against maximum likelihood -------------------------

test_that("with fixed effects only, INLA agrees with the maximum likelihood fit", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("MASS")
  skip_on_cran()
  set.seed(31)
  n <- 3000L
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
  eta <- 0.6 * d$x + 0.5 * (d$g == "b")
  tau <- c(-1, 0.2, 1.3)
  for (fam in c("cumulative", "cratio")) {
    p <- .ord_probs_reference(eta, tau, fam)
    d$y <- 1L + rowSums(runif(n) > t(apply(p, 1L, cumsum))[, -4L])
    fit <- .ord_quiet(flexybayes(y ~ x + g, data = d, family = fam,
      backend = "inla", verbose = FALSE))
    # The likelihood written out from the reference probabilities,
    # maximised with optim(): parameters tau_1..3, beta_x, beta_gb.
    nll <- function(par) {
      pr <- .ord_probs_reference(par[4L] * d$x + par[5L] * (d$g == "b"),
        par[1:3], fam)[cbind(seq_len(n), d$y)]
      # Thresholds out of order give the cumulative family a negative
      # probability; the optimiser is turned back from there.
      if (any(pr <= 0)) {
        return(1e10)
      }
      -sum(log(pr))
    }
    mle <- stats::optim(c(-1, 0, 1, 0, 0), nll, method = "BFGS",
      control = list(reltol = 1e-12))$par
    td <- tidy(fit)
    est <- td$estimate[match(c(paste0("threshold[", 1:3, "|", 2:4, "]"),
      "x", "gb"), td$term)]
    se <- td$std.error[match(c(paste0("threshold[", 1:3, "|", 2:4, "]"),
      "x", "gb"), td$term)]
    expect_lt(max(abs(est - mle) / se), 0.15, label = fam)
    if (identical(fam, "cumulative")) {
      polr <- MASS::polr(factor(y, ordered = TRUE) ~ x + g, data = d,
        method = "logistic")
      expect_equal(unname(c(polr$zeta, coef(polr))), mle, tolerance = 1e-4)
    }
  }
})


# Exact K-fold ---------------------------------------------------------------------

test_that("fb_kfold() scores groups jointly and refuses folds it cannot use", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  d <- .ord_sim_trial(17L, "cumulative", field = FALSE, n_plant = 4L)
  fit <- .ord_quiet(flexybayes(stage ~ depth, random = ~ rep + plot, data = d,
    family = "cumulative", backend = "inla", verbose = FALSE))
  plot_fold <- rep(1:3, length.out = 90L)
  folds <- plot_fold[as.integer(d$plot)]
  kf <- .ord_quiet(fb_kfold(fit, folds, group = "plot", n_draws = 300L,
    seed = 1L))
  expect_s3_class(kf, "fb_kfold")
  expect_identical(nrow(kf$pointwise), 90L)
  expect_true(all(kf$pointwise$n_rows == 4L))
  expect_equal(kf$estimates$elpd, sum(kf$pointwise$elpd))
  expect_true(all(is.finite(kf$pointwise$elpd)))
  expect_output(print(kf), "groups of `plot` scored jointly")
  expect_error(fb_kfold(fit, folds[-1L]),
    class = "flexybayes_refusal_kfold_folds_invalid")
  expect_error(fb_kfold(fit, rep(1L, nrow(d))),
    class = "flexybayes_refusal_kfold_folds_invalid")
  expect_error(fb_kfold(fit, rep(1:2, length.out = nrow(d)), group = "plot"),
    class = "flexybayes_refusal_kfold_folds_invalid")
  rare <- ifelse(d$stage == "Z6", 1L, 2L)
  expect_error(fb_kfold(fit, rare),
    class = "flexybayes_refusal_kfold_category_missing_in_training")
  g <- .ord_quiet(flexybayes(x ~ 1, data = data.frame(x = rnorm(20)),
    backend = "inla", verbose = FALSE))
  expect_error(fb_kfold(g, rep(1:2, 10L)),
    class = "flexybayes_refusal_kfold_family_not_supported")
})


# brms -----------------------------------------------------------------------------

test_that("the cumulative family carries INLA's threshold prior on brms", {
  skip_if_not_installed("brms")
  skip_on_cran()
  d <- .ord_sim_trial(41L, "cumulative", field = FALSE, n_plant = 4L)
  code_of <- function(formula, ...) {
    .ord_quiet(flexybayes(formula, random = ~ rep + plot, data = d,
      family = "cumulative", backend = "brms", return_code = TRUE,
      verbose = FALSE, ...))
  }
  # Dirichlet(1) on the category probabilities at a linear predictor of zero
  # is the link's density on each ordered threshold, taken where brms
  # reports b_Intercept: its centred thresholds shifted back.
  logit <- code_of(stage ~ depth)
  expect_match(logit,
    "target += logistic_lpdf(Intercept + dot_product(means_X, b) | 0, 1);",
    fixed = TRUE)
  expect_no_match(logit, "student_t_lpdf(Intercept", fixed = TRUE)
  expect_match(code_of(stage ~ depth, link = "probit"),
    "target += std_normal_lpdf(Intercept + dot_product(means_X, b));",
    fixed = TRUE)
  expect_match(code_of(stage ~ 1),
    "target += logistic_lpdf(Intercept | 0, 1);", fixed = TRUE)
  # prior_fixed_sd reaches the coefficients and leaves the thresholds alone,
  # as on INLA, where the intercept is pinned.
  scaled <- code_of(stage ~ depth, prior_fixed_sd = 1)
  expect_match(scaled, "logistic_lpdf(Intercept + dot_product", fixed = TRUE)
  expect_no_match(scaled, "normal_lpdf(Intercept", fixed = TRUE)
  expect_match(scaled, "normal_lpdf(b | 0, 1)", fixed = TRUE)
})

test_that("both families fit on brms, and agree with INLA without the field", {
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  for (fam in c("cumulative", "cratio")) {
    d <- .ord_sim_trial(41L, fam, field = FALSE)
    random <- ~ rep + rep:depth + plot
    fb <- .ord_quiet(flexybayes(stage ~ depth, random = random, data = d,
      family = fam, backend = "brms", verbose = FALSE, mcmc_verbose = FALSE,
      chains = 2L, n_samples = 1000L, warmup = 1000L, seed = 1L,
      control = list(adapt_delta = 0.95)))
    fi <- .ord_quiet(flexybayes(stage ~ depth, random = random, data = d,
      family = fam, backend = "inla", verbose = FALSE))
    expect_s3_class(fb, c("flexybayes_ordinal", "flexybayes_brms"))
    expect_identical(names(coef(fb)), c("depth80", "depth120"))
    tb <- tidy(fb)
    ti <- tidy(fi)
    gap <- abs(ti$estimate[match(tb$term, ti$term)] - tb$estimate) /
      tb$std.error
    expect_lt(max(gap), 0.3, label = fam)
    p <- fitted(fb, n_draws = 200L, seed = 1L)
    expect_identical(dim(p), c(900L, 6L))
    nd <- data.frame(depth = factor(c(30, 120)))
    expect_lt(
      max(abs(predict(fb, newdata = nd, marginal = TRUE, seed = 1L) -
        predict(fi, newdata = nd, marginal = TRUE, seed = 1L))),
      0.03,
      label = fam
    )
    expect_error(residuals(fb),
      class = "flexybayes_refusal_ordinal_method_not_supported")
    if (identical(fam, "cumulative")) {
      # The statement is read back from the program brms compiled.
      expect_match(fb$extras$ordinal$threshold_prior, "Dirichlet(1)",
        fixed = TRUE)
      ps <- prior_summary(fb)
      slot <- ps$slots[ps$slots$slot_id == "stan:threshold:Intercept", ]
      expect_identical(nrow(slot), 1L)
      expect_identical(slot$input_id, "threshold:Intercept")
      expect_true(slot$readback_agrees)
      input <- ps$inputs[ps$inputs$input_id == "threshold:Intercept", ]
      expect_identical(input$distribution, "dirichlet")
      expect_identical(input$source, "structure_default")
    }
  }
})
