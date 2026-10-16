# Heterogeneous variance over an outer factor: diag(f):g, idh(f):g, at(f):g.
#
# These assert on the emitted model: the formula that reached brms and the
# number of variance parameters it implies.
#
# The mapping diag(f):g -> (0 + f || g) was checked against ASReml and lme4
# in a separate script: all three fit exactly k variances and 0
# covariances, and the two REML fits agree to 3.5e-05.

.het_data <- function(n_site = 3L, n_gen = 8L, n_rep = 2L, seed = 1L) {
  set.seed(seed)
  d <- expand.grid(
    rep = seq_len(n_rep),
    gen = factor(seq_len(n_gen)),
    site = factor(seq_len(n_site))
  )
  d$y <- stats::rnorm(nrow(d))
  d
}

.emitted <- function(random, data = .het_data(), fixed = y ~ site) {
  code <- flexybayes(
    fixed,
    random = random,
    data = data,
    backend = "brms",
    return_code = TRUE
  )
  paste(utils::capture.output(print(code)), collapse = "\n")
}

test_that("diag(f):g emits an uncorrelated random slope", {
  skip_if_not_installed("brms")
  txt <- .emitted(~ diag(site):gen)
  # The Stan program for (0 + site || gen) carries one sd per site level and
  # no correlation block. `sd_1` with three elements, no `L_1` Cholesky.
  expect_match(txt, "sd_1", fixed = TRUE)
  expect_false(grepl("cor_1|L_1", txt))
})

test_that("at(), diag() and idh() are the same structure", {
  skip_if_not_installed("brms")
  d <- .het_data()
  a <- .emitted(~ at(site):gen, d)
  b <- .emitted(~ diag(site):gen, d)
  i <- .emitted(~ idh(site):gen, d)
  # ASReml returns the same components and standard errors for diag() and
  # idh(), so they must emit the same model.
  expect_identical(a, b)
  expect_identical(a, i)
})

test_that("us(f):g keeps the correlations that diag(f):g drops", {
  skip_if_not_installed("brms")
  d <- .het_data()
  diag_txt <- .emitted(~ diag(site):gen, d)
  us_txt <- .emitted(~ us(site):gen, d)
  # One character apart in the brms formula, k vs k(k+1)/2 parameters in the
  # model. If these two ever emit the same program, one of them is wrong.
  expect_false(identical(diag_txt, us_txt))
  expect_true(grepl("cor_1|L_1", us_txt))
})

test_that("at(f, level):g is not emitted as the diagonal structure", {
  skip_if_not_installed("brms")
  d <- .het_data()
  # at(site, 2):gen gives the effect only where site == 2. diag(site):gen
  # gives it everywhere with a per-level variance.
  at_fb <- fb_from_asreml(
    fixed = y ~ site,
    random = ~ at(site, 2):gen,
    data = d
  )
  expect_identical(
    format(flexyBayes:::.fb_to_brms_formula(at_fb)),
    "y ~ site + (0 + site2 || gen)"
  )
  expect_false(identical(
    .emitted(~ at(site, 2):gen, d),
    .emitted(~ diag(site):gen, d)
  ))
})

test_that("the refusal code is registered in the taxonomy", {
  expect_true(exists(
    "at_level_conditioning_unsupported",
    envir = flexyBayes:::.refusal_registry,
    inherits = FALSE
  ))
})

test_that("corh(f):g is refused by name, not approximated", {
  skip_if_not_installed("brms")
  d <- .het_data()
  # corh() asks for k variances plus one shared correlation. brms offers
  # only uncorrelated or fully unstructured, so it must refuse.
  err <- tryCatch(.emitted(~ corh(site):gen, d), error = function(e) e)
  expect_s3_class(
    err,
    "flexybayes_refusal_corh_no_equicorrelation_representation"
  )
  # The refusal has to be actionable: it must name both available structures.
  expect_match(conditionMessage(err), "diag\\(")
  expect_match(conditionMessage(err), "us\\(")
})

test_that("auto reaches heterogeneity without the user naming a backend", {
  skip_if_not_installed("brms")
  d <- .het_data()
  # INLA refuses these structurally, so `auto` must fall through to brms.
  # If it did not, the capability would exist but be unreachable to anyone
  # who had not read the backend documentation.
  expect_no_error(flexybayes(
    y ~ site,
    random = ~ diag(site):gen,
    data = d,
    backend = "auto",
    return_code = TRUE
  ))
})

test_that("dsum(~units|f) and at(f):units are one node with two spellings", {
  skip_if_not_installed("brms")
  d <- .het_data()
  a <- .emitted(~gen, d) # baseline, homogeneous residual
  b <- flexybayes(
    y ~ site,
    random = ~gen,
    residual = ~ dsum(~ units | site),
    data = d,
    backend = "brms",
    return_code = TRUE
  )
  cc <- flexybayes(
    y ~ site,
    random = ~gen,
    residual = ~ at(site):units,
    data = d,
    backend = "brms",
    return_code = TRUE
  )
  b_txt <- paste(utils::capture.output(print(b)), collapse = "\n")
  c_txt <- paste(utils::capture.output(print(cc)), collapse = "\n")
  expect_identical(b_txt, c_txt)
  # A heterogeneous residual becomes a linear predictor on log sigma, so the
  # Stan program gains a sigma design matrix that the homogeneous fit has not.
  expect_true(grepl("X_sigma", b_txt, fixed = TRUE))
  expect_false(grepl("X_sigma", a, fixed = TRUE))
})

test_that("a heterogeneous residual is refused for a family with no sigma", {
  skip_if_not_installed("brms")
  d <- .het_data()
  d$cnt <- stats::rpois(nrow(d), 5)
  # Poisson has no residual scale, so `sigma ~ f` names a parameter that
  # does not exist.
  err <- tryCatch(
    flexybayes(
      cnt ~ site,
      random = ~gen,
      residual = ~ dsum(~ units | site),
      data = d,
      family = "poisson",
      backend = "brms",
      return_code = TRUE
    ),
    error = function(e) e
  )
  expect_s3_class(
    err,
    "flexybayes_refusal_heterogeneous_residual_family_has_no_sigma"
  )
})

test_that("INLA refuses a heterogeneous residual for a family with no residual scale", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  # INLA split the response into one Poisson or binomial likelihood per
  # section, which is the unsectioned model under the sectioned model's
  # name. The refusal comes before the engine is called. The engine call is
  # not mocked: a mock of .inla_call here left a later brms fit in this file
  # unable to read its CmdStan output, and without the refusal the fit runs
  # and returns, which fails the class check.
  d <- .het_data()
  d$cnt <- stats::rpois(nrow(d), 5)
  d$hit <- stats::rbinom(nrow(d), 1L, 0.4)
  for (spec in list(c("cnt", "poisson"), c("hit", "binomial"))) {
    err <- tryCatch(
      flexybayes(
        stats::as.formula(paste(spec[[1L]], "~ site")),
        random = ~gen,
        residual = ~ dsum(~ units | site),
        data = d,
        family = spec[[2L]],
        backend = "inla",
        verbose = FALSE
      ),
      error = function(e) e
    )
    expect_s3_class(
      err,
      "flexybayes_refusal_heterogeneous_residual_family_has_no_sigma"
    )
  }
})

test_that("the scalar sigma prior is retargeted, not left to match nothing", {
  skip_if_not_installed("brms")
  d <- .het_data()
  # With `sigma ~ 0 + f` there is no scalar sigma. A prior aimed at it makes
  # brms refuse the fit ("The following priors do not correspond to any
  # model parameter"). The prior must land on the sigma coefficient block,
  # and the program must compile.
  txt <- paste(
    utils::capture.output(print(flexybayes(
      y ~ site,
      random = ~gen,
      residual = ~ dsum(~ units | site),
      data = d,
      backend = "brms",
      return_code = TRUE
    ))),
    collapse = "\n"
  )
  expect_true(grepl("b_sigma", txt, fixed = TRUE))
  expect_match(txt, "lprior \\+= normal_lpdf\\(b_sigma")
})

test_that("a heterogeneous fit returns one variance per level and no correlation", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  skip_if_not_long_tier()
  d <- .het_data(n_gen = 12L, n_rep = 3L, seed = 4L)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ site,
    random = ~ diag(site):gen,
    data = d,
    backend = "brms",
    chains = 1L,
    n_samples = 400L,
    warmup = 400L
  )))
  vc <- brms::VarCorr(fit$brms)
  # A diagonal over k levels is exactly k variances and zero correlations.
  # Values are not asserted: 400 draws on synthetic noise carry no
  # information.
  expect_equal(nrow(vc$gen$sd), nlevels(d$site))
  expect_null(vc$gen$cor)
})

test_that("a sectioned residual reports its per-level variances", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  skip_if_not_long_tier()
  d <- .het_data(n_gen = 10L, n_rep = 3L, seed = 6L)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ site,
    random = ~gen,
    residual = ~ dsum(~ units | site),
    data = d,
    backend = "brms",
    chains = 2L,
    n_samples = 400L,
    warmup = 400L
  )))

  tab <- flexyBayes:::.brms_residual_by_level_table(fit)
  expect_false(is.null(tab))
  # Labelled by the sectioning factor's own levels, in its own order. A
  # name such as b_sigma_site2 reads as a variance to an ASReml user.
  expect_identical(tab$level, levels(d$site))
  expect_identical(attr(tab, "factor"), "site")

  # Matches a direct computation from the draws: the median of the
  # transform, not the transform of a mean.
  draws <- as.matrix(posterior::as_draws_matrix(fit$brms))
  for (i in seq_along(tab$level)) {
    col <- paste0("b_sigma_site", tab$level[[i]])
    expect_equal(
      tab$variance[[i]],
      stats::median(exp(2 * draws[, col])),
      tolerance = 1e-12
    )
    expect_equal(
      tab$sd[[i]],
      stats::median(exp(draws[, col])),
      tolerance = 1e-12
    )
    expect_true(tab$variance_lower[[i]] < tab$variance[[i]])
    expect_true(tab$variance[[i]] < tab$variance_upper[[i]])
  }

  # There is no scalar residual variance left for the variance-component
  # table to report, which is why the block exists at all.
  expect_false("sigma" %in% colnames(draws))
  expect_false("sigma" %in% fit$extras$variance_comps$component)

  # Both user-facing surfaces carry it, and neither invents it.
  for (out in list(
    utils::capture.output(print(fit)),
    utils::capture.output(summary(fit))
  )) {
    expect_true(any(grepl("Residual by level of `site`", out, fixed = TRUE)))
    expect_true(any(grepl("level", out, fixed = TRUE)))
  }

  # A homogeneous residual must not grow the block.
  plain <- suppressWarnings(suppressMessages(flexybayes(
    y ~ site,
    random = ~gen,
    data = d,
    backend = "brms",
    chains = 1L,
    n_samples = 200L,
    warmup = 200L
  )))
  expect_null(flexyBayes:::.brms_residual_by_level_table(plain))
  expect_false(any(grepl(
    "Residual by level",
    utils::capture.output(print(plain)),
    fixed = TRUE
  )))
})

# =============================================================================
# The per-level residual reaches the object, not only the printout.
#
# The per-level residual of a dsum() fit is read from the summary object,
# without opening `$brms` or transforming `b_sigma_*` draws by hand.
# =============================================================================

test_that("the per-level residuals are rows of summary(fit)$varcomp", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  skip_if_not_long_tier()
  d <- .het_data(n_gen = 10L, n_rep = 3L, seed = 6L)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ site,
    random = ~gen,
    residual = ~ dsum(~ units | site),
    data = d,
    backend = "brms",
    chains = 2L,
    n_samples = 400L,
    warmup = 400L
  )))

  vc <- summary(fit)$varcomp
  # One row per level of the sectioning factor, plus the components the
  # table always carried.
  expect_identical(
    vc$component,
    c("sd_gen", paste0("sigma_", levels(d$site)))
  )
  expect_true(all(is.finite(vc$estimate)))
  expect_true(all(is.finite(vc$std.error)))
  expect_true(all(vc$conf.low < vc$estimate))
  expect_true(all(vc$estimate < vc$conf.high))

  # The table and the printed panel are two readings of one builder, so
  # they cannot state different numbers for one fit.
  tab <- flexyBayes:::.brms_residual_by_level_table(fit)
  resid <- vc[startsWith(vc$component, "sigma_"), , drop = FALSE]
  expect_equal(resid$estimate, tab$sd, tolerance = 1e-12)
  expect_equal(resid$conf.low, tab$sd_lower, tolerance = 1e-12)
  expect_equal(resid$conf.high, tab$sd_upper, tolerance = 1e-12)

  # The prior cell names what reached the sampler. The declared uniform on
  # the SD scale was retargeted onto the log-sigma coefficients at emit
  # time. A `sigma` spec would name a parameter this model does not have,
  # and "engine default" would claim the package set nothing.
  expect_false(any(resid$prior == "engine default"))
  expect_true(all(grepl("normal", resid$prior, fixed = TRUE)))

  # A homogeneous fit gains no such rows.
  plain <- suppressWarnings(suppressMessages(flexybayes(
    y ~ site,
    random = ~gen,
    data = d,
    backend = "brms",
    chains = 1L,
    n_samples = 200L,
    warmup = 200L
  )))
  expect_false(any(startsWith(summary(plain)$varcomp$component, "sigma_")))
})

test_that("on INLA each residual section's prior cell names the prior it received", {
  # INLA has one family hyperparameter per section, each given the residual prior.
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .het_data(n_gen = 12L, n_rep = 2L, seed = 1L)
  sections <- paste0("sigma_", levels(d$site))

  fit <- suppressWarnings(suppressMessages(fb_inla(
    y ~ site,
    random = ~gen,
    residual = ~ dsum(~ units | site),
    data = d,
    prior = fb_prior(sigma ~ half_normal(1.5), sd(gen) ~ half_normal(2))
  )))
  vc <- summary(fit)$varcomp
  expect_identical(
    vc$prior[match(sections, vc$component)],
    rep("half_normal(scale=1.5)", length(sections))
  )

  fit0 <- suppressWarnings(suppressMessages(fb_inla(
    y ~ site,
    random = ~gen,
    residual = ~ dsum(~ units | site),
    data = d
  )))
  vc0 <- summary(fit0)$varcomp
  expect_identical(
    vc0$prior[match(sections, vc0$component)],
    rep(vc0$prior[vc0$component == "sd_gen"], length(sections))
  )
})

test_that("a sectioned-residual fit answers the mean-model accessors", {
  # A bare `^b_` match also catches the distributional coefficients, which
  # are not effects on the response.
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  skip_if_not_installed("emmeans")
  d <- .het_data(n_gen = 10L, n_rep = 3L, seed = 6L)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ site,
    random = ~gen,
    residual = ~ dsum(~ units | site),
    data = d,
    backend = "brms",
    chains = 2L,
    n_samples = 400L,
    warmup = 400L
  )))

  # Mean-model coefficients only, on coef(), vcov() and confint() alike.
  expect_false(any(startsWith(names(stats::coef(fit)), "sigma_")))
  expect_identical(
    rownames(stats::confint(fit)),
    names(stats::coef(fit))
  )
  expect_identical(
    nrow(stats::vcov(fit)),
    length(stats::coef(fit))
  )

  # The response and residuals are recoverable, which they were not while
  # the formula reader was indexing a brmsformula by position.
  expect_false(all(is.na(fit$glm$y)))
  expect_false(all(is.na(fit$glm$residuals)))

  # And the marginal means the sectioning was never supposed to block.
  out <- predict(fit, classify = "site")
  expect_s3_class(out, "fb_predict_classify")
  expect_identical(nrow(out), nlevels(d$site))
  expect_true(all(is.finite(out$estimate)))
  expect_true(all(out$std.error > 0))
})
