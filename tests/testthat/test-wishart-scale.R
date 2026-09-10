# The iidKd Wishart is scaled to the data, not left on INLA's identity.
#
# iid2d and iidKd carry one joint Wishart on theta1 whose scale matrix R
# defaults to the identity -- a prior centred on unit variances. There is no
# standard-deviation slot for the prior DSL to write into, so the package's
# own default was left unspliced, which left that unit-scale prior in force
# whatever the response's scale.
#
# On lme4::sleepstudy, whose response has a standard deviation near 50, that
# made the posterior mode non-unique: eight identical calls returned intercept
# standard deviations from 0.50 to 54.8 and correlations from -1.00 to +0.75,
# none near lme4's REML answer. It was recorded as a property of INLA. It was
# the default.

test_that("the Wishart argument carries the response's scale", {
  d <- data.frame(
    y = stats::rnorm(80, 200, 50),
    x = stats::rnorm(80, 0, 4),
    g = factor(rep(letters[1:8], each = 10))
  )
  fb <- fb_from_brms(y ~ x + (x | g), data = d)
  # The prior layer's own scale basis, which this reads rather than
  # recomputing: U = 5 * sd(y).
  fb$priors <- structure(list(specs = list()), class = c("fb_prior", "list"))
  attr(fb$priors, "fb_prior_default_scale") <- 5 * 50

  sc <- .fb_inla_wishart_scale_slope(fb, fb$random_terms[[1L]])
  expect_length(sc, 2L)
  # Intercept block on var(y); slope block on var(y) / var(x) -- the slope's
  # variance in the units the slope is measured in, so it is smaller than
  # the intercept's whenever x is more spread out than one unit.
  expect_equal(sc[[1L]], 50^2, tolerance = 0.01)
  expect_equal(sc[[2L]], sc[[1L]] / stats::var(d$x), tolerance = 0.01)

  arg <- .fb_inla_wishart_arg(sc)
  expect_match(arg, "theta1", fixed = TRUE)
  # r = k + 2 = 4 at k = 2, however format() renders it.
  expect_match(arg, "param = c\\(4(\\.0+)?,")
})

test_that("a k-block Wishart takes k diagonal entries and one dof", {
  arg <- .fb_inla_wishart_arg(rep(9, 3L))
  # r = k + 2, then k diagonal scales, then k(k-1)/2 zero off-diagonals.
  expect_match(arg, "param = c(5, 9, 9, 9, 0, 0, 0)", fixed = TRUE)
})

test_that("an unusable scale falls back to INLA's own default", {
  # Better the identity than a fabricated number: the argument is dropped
  # entirely rather than guessed at.
  expect_identical(.fb_inla_wishart_arg(c(NA_real_, 1)), "")
  expect_identical(.fb_inla_wishart_arg(1), "")
})

test_that("the emitted formula carries the scaled Wishart", {
  d <- data.frame(
    y = stats::rnorm(80, 200, 50),
    x = stats::rnorm(80, 0, 4),
    g = factor(rep(letters[1:8], each = 10))
  )
  fb <- fb_from_brms(y ~ x + (x | g), data = d, family = "gaussian")
  fb$priors <- structure(list(specs = list()), class = c("fb_prior", "list"))
  attr(fb$priors, "fb_prior_default_scale") <- 5 * 50
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                collapse = " "))
  expect_match(txt, "iid2d", fixed = TRUE)
  expect_match(txt, "theta1", fixed = TRUE)
})

test_that("the correlated slope recovers a large-scale response, repeatably", {
  # The regression this exists for. On the unit-scale default these five
  # calls disagreed by more than the parameter itself.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_boundary_collapse_warning = TRUE
  )
  data(sleepstudy, package = "lme4")

  got <- t(vapply(seq_len(3L), function(i) {
    f <- suppressMessages(suppressWarnings(flexybayes(
      Reaction ~ Days + (Days | Subject), data = sleepstudy,
      backend = "inla", verbose = FALSE
    )))
    v <- summary(f)$varcomp
    g <- stats::setNames(v$estimate, v$component)
    c(g[["sd_Subject"]], g[["cor_Subject__Intercept__Days"]])
  }, numeric(2L)))

  ref <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  vc <- as.data.frame(lme4::VarCorr(ref))
  ref_sd <- vc$sdcor[vc$grp == "Subject" & vc$var1 == "(Intercept)" &
                       is.na(vc$var2)]

  # Near lme4 rather than fifty times away from it.
  expect_equal(mean(got[, 1L]), ref_sd, tolerance = 0.25)
  # And repeatable: the correlation used to span nearly the whole of
  # [-1, 1] across identical calls.
  expect_lt(diff(range(got[, 2L])), 0.5)
  expect_lt(diff(range(got[, 1L])), 5)
})
